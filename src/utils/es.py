"""Module to query ES indexes"""

from __future__ import print_function

# http://www.elasticsearch.org/guide/reference/query-dsl/custom-filters-score-query.html
# http://www.elasticsearch.org/guide/reference/query-dsl/custom-score-query.html
# http://www.elasticsearch.org/guide/reference/query-dsl/custom-boost-factor-query.html
# http://www.elasticsearch.org/guide/reference/query-dsl/boosting-query.html
import sys

import json
import re
import time
import copy
import requests

from config import (ES_INDEX_NAME_TIER1, ES_INDEX_NAME,
                    SOURCE_TRANSLATORS, GENOME_ASSEMBLY,
                    TAXONOMY, ES_HOST,  ES_INDEX_TYPE)
from biothings.utils.common import (ask, is_int, is_str,
                                    is_seq, timesofar)
from biothings.www.api.es import ESQuery, QueryError, ESQueryBuilder, \
                                 parse_facets_option
from elasticsearch import Elasticsearch
from .userfilters import UserFilters

from elasticsearch import helpers
from utils.mongo import doc_feeder


import logging
formatter = logging.Formatter("%(levelname)s:%(name)s:%(message)s")
es_logger = logging.getLogger('elasticsearch')
es_logger.setLevel(logging.WARNING)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
es_logger.addHandler(ch)

es_tracer = logging.getLogger('elasticsearch.trace')
es_tracer.setLevel(logging.WARNING)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
es_tracer.addHandler(ch)


def safe_genome_pos(s):
    '''
       safe_genome_pos(1000) = 1000
       safe_genome_pos('1000') = 1000
       safe_genome_pos('10,000') = 100000
    '''
    if isinstance(s, int):
        return s
    elif isinstance(s, str):
        return int(s.replace(',', ''))
    else:
        raise ValueError('invalid type "%s" for "save_genome_pos"' % type(s))


class ESQuery(ESQuery):
    def __init__(self):
        super(ESQuery, self).__init__()
        self._default_fields = ['name', 'symbol', 'taxid', 'entrezgene']
        self._default_species = [9606, 10090, 10116]  # human, mouse, rat
        self._tier_1_species = set(TAXONOMY.values())

    def _search(self, q, species='all', scroll_options={}, **kwargs):
        self._set_index(species)
        res = self._es.search(index=self._index, doc_type=self._doc_type,
                              body=q, **scroll_options)
        # logging.error("in search: %s" % res)
        self._index = ES_INDEX_NAME  # reset self._index
        return res

    def _msearch(self, **kwargs):
        self._set_index(kwargs.get('species', 'all'))
        # logging.debug("_msearch: %s" % kwargs['body'])
        res = super(ESQuery, self)._msearch(**kwargs)
        self._index = ES_INDEX_NAME     # reset self._index
        return res

    def _set_index(self, species):
        '''set proper index for given species parameter.'''
        if species == 'all' or len(set(species)-self._tier_1_species) > 0:
            self._index = ES_INDEX_NAME
        else:
            self._index = ES_INDEX_NAME_TIER1

    def _get_query_builder(self, **kwargs):
        return ESQueryBuilder(**kwargs)

    def _build_query(self, q, kwargs):
        # can override this function if more query types are to be added
        esqb = self._get_query_builder(**kwargs)
        return esqb.query(q)

    def _cleaned_species(self, species, default_to_none=False):
        '''return a cleaned species parameter.
           should be either "all" or a list of taxids/species_names
           or a single taxid/species_name
           returned species is always a list of taxids (even when
           only one species)
        '''
        if species is None:
            # set to default_species
            return None if default_to_none else self._default_species
        if isinstance(species, int):
            return [species]

        if is_str(species):
            if species.lower() == 'all':
                # if self.species == 'all': do not apply species filter,
                # all species is included.
                return species.lower()
            else:
                species = [s.strip().lower() for s in species.split(',')]

        if not is_seq(species):
            raise ValueError('"species" parameter must be a string, integer ' +
                             'or a list/tuple, not "{}".'.format(
                                                          type(species)))

        _species = []
        for s in species:
            if is_int(s):
                _species.append(int(s))
            elif s in TAXONOMY:
                _species.append(TAXONOMY[s])
        return _species

    def _get_options(self, options, kwargs):
        # this species parameter is added to the query, thus will
        # change facet counts.
        kwargs['species'] = self._cleaned_species(kwargs.get('species', None))
        include_tax_tree = kwargs.pop('include_tax_tree', False)
        if include_tax_tree:
            headers = {'content-type': 'application/x-www-form-urlencoded',
                       'user-agent': "Python-requests_mygene.info/%s (gzip)"
                       % requests.__version__}
            # TODO: URL as config param
            res = requests.post('http://s.biothings.io/v1/species?ids=' +
                                ','.join(['{}'.format(sid) for sid in
                                          kwargs['species']]) +
                                '&expand_species=true', headers=headers)
            if res.status_code == requests.codes.ok:
                kwargs['species'] = res.json()

        # this parameter is to add species filter without
        # changing facet counts.
        kwargs['species_facet_filter'] = self._cleaned_species(
                kwargs.get('species_facet_filter', None), default_to_none=True)
        options.kwargs = kwargs
        return options

    def metadata(self, raw=False):
        '''return metadata about the index.'''
        mapping = self._es.indices.get_mapping(self._index, self._doc_type)
        if raw:
            return mapping

        def get_fields(properties):
            for k, v in list(properties.items()):
                if 'properties' in v:
                    for f in get_fields(v['properties']):
                        yield f
                else:
                    if v.get('index', None) == 'no':
                        continue
                    f = v.get('index_name', k)
                    yield f
        mapping = list(mapping.values())[0]['mappings']
        # field_set = set(get_fields(mapping[self._doc_type]['properties']))
        metadata = {
            # 'available_fields': sorted(field_set)
            # TODO: http://mygene.info as config
            'available_fields': 'http://mygene.info/metadata/fields'
        }
        if '_meta' in mapping[self._doc_type]:
            metadata.update(mapping[self._doc_type]['_meta'])
        metadata['genome_assembly'] = GENOME_ASSEMBLY
        metadata['taxonomy'] = TAXONOMY
        if "source" not in metadata:
             # occurs when loaded from scratch, not from a change/diff file
            metadata["source"] = None
        return metadata

    def get_gene(self, geneid, **kwargs):
        '''for /gene/<geneid>'''
        options = self._get_cleaned_annotation_options(kwargs)
        qbdr = ESQueryBuilder(**options.kwargs)
        _q = qbdr.build_id_query(geneid, options.scopes)
        if options.rawquery:
            return _q
        res = self._search(_q, species=options.kwargs['species'])
        if not options.raw:
            res = self._cleaned_res(res, empty=None, single_hit=True,options=options)
        return res


class ESQueryBuilder(ESQueryBuilder):
    def __init__(self, **query_options):
        """You can pass these options:
            fields     default ['name', 'symbol', 'taxid', 'entrezgene']
            from       default 0
            size       default 10
            sort       e.g. sort='entrezgene,-symbol'
            explain    true or false
            facets     a field or a list of fields, default None

            species
            species_facet_filter
            entrezonly  default false
            ensemblonly default false
            userfilter  optional, provide the name of a saved user filter
                                  (in "userfilters" index)
            exists      optional, passing field, comma-separated fields,
                                  returned genes must have given field(s).
            missing     optional, passing field, comma-separated fields,
                                  returned genes must have NO given field(s).

        """
        super(ESQueryBuilder, self).__init__()
        self._query_options = query_options
        # species should be either 'all' or a list of taxids.
        self.species = self._query_options.pop('species', 'all')
        self.species_facet_filter = self._query_options.pop(
                    'species_facet_filter', None)
        self.entrezonly = self._query_options.pop('entrezonly', False)
        self.ensemblonly = self._query_options.pop('ensemblonly', False)
        # userfilter
        userfilter = self._query_options.pop('userfilter', None)
        self.userfilter = userfilter.split(',') if userfilter else None
        # exist filter
        existsfilter = self._query_options.pop('exists', None)
        self.existsfilter = existsfilter.split(',') if existsfilter else None
        # missing filter
        missingfilter = self._query_options.pop('missing', None)
        self.missingfilter = missingfilter.split(',') \
            if missingfilter else None

        parse_facets_option(self._query_options)
        # TODO: taken from config file ?
        self._allowed_options = ['fields', '_source', 'start', 'from', 'size',
                                 'sort', 'explain', 'version', 'aggs',
                                 'dotfield']
        for key in set(self._query_options) - set(self._allowed_options):
            del self._query_options[key]
        # convert "fields" option to "_source"
        # use "_source" instead of "fields" for ES v1.x and up
        if 'fields' in self._query_options and \
                self._query_options['fields'] is not None:
            self._query_options['_source'] = self._query_options['fields']
            del self._query_options['fields']

        # this is a fake query to make sure to return empty hits
        self._nohits_query = {
            "match": {
                'non_exist_field': ''
            }
        }

    def _translate_datasource(self, q):
        for src in SOURCE_TRANSLATORS.keys():
            q = re.sub(src, SOURCE_TRANSLATORS[src], q)
        # logging.debug(q)
        return q
        pass

    def _parse_interval_query(self, query):
        '''Check if the input query string matches interval search regex,
           if yes, return a dictionary with three key-value pairs:
              chr
              gstart
              gend
            , otherwise, return None.
        '''
        pattern = r'chr(?P<chr>\w+):(?P<gstart>[0-9,]+)-(?P<gend>[0-9,]+)'
        if query:
            mat = re.search(pattern, query)
            if mat:
                d = mat.groupdict()
                if query.startswith('hg19.'):
                    # support hg19 for human (default is hg38)
                    d['assembly'] = 'hg19'
                if query.startswith('mm9.'):
                    # support mm9 for mouse (default is mm10)
                    d['assembly'] = 'mm9'

                return d

    def dis_max_query(self, q):
        # remove '"' and '\' from q, they will break json decoder.
        q = q.replace('"', '').replace('\\', '')
        _query = {
            "dis_max": {
                "tie_breaker": 0,
                "boost": 1,
                "queries": [
                    {
                        "function_score": {
                            "query": {
                                "match": {
                                    "symbol": {
                                        "query": "%(q)s",
                                        "analyzer": "whitespace_lowercase"
                                    }
                                },
                            },
                            "weight": 5
                        }
                    },
                    {
                        "function_score": {
                            "query": {
                                # This makes phrase match of "cyclin-dependent
                                # kinase 2" appears first
                                "match_phrase": {"name": "%(q)s"},
                            },
                            "weight": 4

                        }
                    },
                    {
                        "function_score": {
                            "query": {
                                "match": {
                                    "name": {
                                        "query": "%(q)s",
                                        "operator": "and",
                                        "analyzer": "whitespace_lowercase"
                                    }
                                },
                            },
                            "weight": 3
                        }
                    },
                    {
                        "function_score": {
                            "query": {
                                "match": {
                                    "unigene": {
                                        "query": "%(q)s",
                                        "analyzer": "string_lowercase"
                                    }
                                }
                            },
                            "weight": 1.1
                        }
                    },
                    {
                        "function_score": {
                            "query": {
                                "match": {
                                    "go": {
                                        "query": "%(q)s",
                                        "analyzer": "string_lowercase"
                                    }
                                }
                            },
                            "weight": 1.1
                        }
                    },
                    # {
                    # "custom_boost_factor": {
                    #     "query" : {
                    #         "match" : { "_all" : {
                    #                     "query": "%(q)s",
                    #                     "analyzer": "whitespace_lowercase"
                    #             }
                    #         },
                    #     },
                    #     "boost_factor": 1
                    # }
                    # },
                    {
                        "function_score": {
                            "query": {
                                "query_string": {
                                    "query": "%(q)s",
                                    "default_operator": "AND",
                                    "auto_generate_phrase_queries": True
                                },
                            },
                            "weight": 1
                        }
                    },

                ]
            }
        }
        _query = json.dumps(_query)
        _query = json.loads(_query % {'q': q})

        if is_int(q):
            _query['dis_max']['queries'] = []
            _query['dis_max']['queries'].insert(
                0,
                {
                    "function_score": {
                        "query": {
                            "term": {"entrezgene": int(q)},
                        },
                        "weight": 8
                    }
                }
            )

        return _query

    def wildcard_query(self, q):
        '''q should contains either * or ?, but not the first character.'''
        _query = {
            "dis_max": {
                "tie_breaker": 0,
                "boost": 1,
                "queries": [
                    {
                        "function_score": {
                            "query": {
                                "wildcard": {
                                    "symbol": {
                                        "value": "%(q)s",
                                        # "weight": 5.0,
                                    }
                                },
                            },
                        }
                    },
                    {
                        "function_score": {
                            "query": {
                                "wildcard": {
                                    "name": {
                                        "value": "%(q)s",
                                        # "weight": 1.1,
                                    }
                                },
                            }
                        }
                    },
                    {
                        "function_score": {
                            "query": {
                                "wildcard": {
                                    "summary": {
                                        "value": "%(q)s",
                                        # "weight": 0.5,
                                    }
                                },
                            }
                        }
                    },

                ]
            }
        }
        _query = json.dumps(_query)
        try:
            _query = json.loads(_query % {'q': q.lower()})
        except ValueError:
            raise QueryError("invalid query term.")

        return _query

    def get_query_filters(self):
        '''filters added here will be applied in a filtered query,
           thus will affect the facet counts.
        '''
        filters = []
        # species filter
        if self.species and self.species != 'all':
            if len(self.species) == 1:
                filters.append({
                    "term": {"taxid": self.species[0]}
                })
            else:
                filters.append({
                    "terms": {"taxid": self.species}
                })
        if self.entrezonly:
            filters.append({
                "exists": {"field": "entrezgene"}
            })
        if self.ensemblonly:
            filters.append({
                "exists": {"field": "ensemblgene"}
            })

        if self.userfilter:
            _uf = UserFilters()
            for _fname in self.userfilter:
                _filter = _uf.get(_fname)
                if _filter:
                    filters.append(_filter['filter'])

        if self.existsfilter:
            for _filter in self.existsfilter:
                filters.append({
                    "exists": {"field": _filter}
                })
        if self.missingfilter:
            for _filter in self.missingfilter:
                filters.append({
                    "missing": {"field": _filter}
                })

        if filters:
            if len(filters) == 1:
                filters = filters[0]
            else:
                # concatenate multiple filters with "and" filter
                filters = {"and": filters}

        return filters

    def add_facet_filters(self, _query):
        """To add filters (e.g. taxid) to restrict returned hits,
            but does not change the scope for facet counts.
        """
        filters = []
        # species_facet_filter
        if self.species_facet_filter:
            if len(self.species) == 1:
                filters.append({
                    "term": {"taxid": self.species_facet_filter[0]}
                })
            else:
                filters.append({
                    "terms": {"taxid": self.species_facet_filter}
                })
        if filters:
            if len(filters) == 1:
                filters = filters[0]
            else:
                # concatenate multiple filters with "and" filter
                filters = {"and": filters}

            # this will not change facet counts
            _query["filter"] = filters

        return _query

    def add_species_custom_filters_score(self, _query):
        _query = {
            "function_score": {
                "query": _query,
                "functions": [
                    # downgrade "pseudogene" matches
                    {
                        "filter": {"term": {"name": "pseudogene"}},
                        "boost_factor": "0.5"
                    },
                    {
                        "filter": {"term": {"taxid": 9606}},
                        "boost_factor": "1.55"
                    },
                    {
                        "filter": {"term": {"taxid": 10090}},
                        "boost_factor": "1.3"
                    },
                    {
                        "filter": {"term": {"taxid": 10116}},
                        "boost_factor": "1.1"
                    },
                ],
                "score_mode": "first"
            }
        }
        return _query

    def query(self, q):
        '''mode:
                1    match query
                2    wildcard query
                3    raw_string query

               else  string_query (for test)
        '''

        # translate data source to provide back-compatibility for
        # some query fields running on ES2
        q = self._translate_datasource(q)

        # Check if special interval query pattern exists
        interval_query = self._parse_interval_query(q)
        if interval_query:
            # should also passing a "taxid" along with interval.
            if self.species != 'all':
                self.species = [self.species[0]]  # TODO: where is it used ?
                _q = self.build_genomic_pos_query(**interval_query)
                return _q
            else:
                raise QueryError('genomic interval query cannot be combined ' +
                                 'with "species=all" parameter. ' +
                                 'Specify a single species.')

        else:
            _query = self.generate_query(q)

            # TODO: this is actually not used, how useful ?
            # _query = self.string_query(q)

            _query = self.add_species_custom_filters_score(_query)
            _q = {'query': _query}
            _q = self.add_facet_filters(_q)
            if self._query_options:
                _q.update(self._query_options)
            # logging.debug("_q = %s" % json.dumps(_q))
            return _q

    # keepit (but similar)
    def build_id_query(self, id, scopes=None):
        id_is_int = is_int(id)
        if scopes is None:
            # by default search three fields:
            # ['entrezgene', 'ensemblgene', 'retired']
            if id_is_int:
                _query = {
                    "multi_match": {
                        "query": id,
                        "fields": ['entrezgene', 'retired']
                    }
                }
            else:
                _query = {
                    "match": {
                        "ensemblgene": {
                            "query": u"{}".format(id),
                            "operator": "and"
                        }
                    }
                }
        else:
            if is_str(scopes):
                _field = scopes
                if _field in ['entrezgene', 'retired']:
                    if id_is_int:
                        _query = {
                            "match": {
                                _field: id
                            }
                        }
                    else:
                        # using a fake query here to make sure
                        # return empty hits
                        _query = self._nohits_query
                else:
                    _query = {
                        "match": {
                            _field: {
                                "query": u"{}".format(id),
                                "operator": "and"
                            }
                        }
                    }
            elif is_seq(scopes):
                int_fields = []
                str_fields = copy.copy(scopes)
                if 'entrezgene' in str_fields:
                    int_fields.append('entrezgene')
                    str_fields.remove('entrezgene')
                if 'retired' in str_fields:
                    int_fields.append('retired')
                    str_fields.remove('retired')

                if id_is_int:
                    if len(int_fields) == 1:
                        _query = {
                            "match": {
                                int_fields[0]: id
                            }
                        }
                    elif len(int_fields) == 2:
                        _query = {
                            "multi_match": {
                                "query": id,
                                "fields": int_fields
                            }
                        }
                    else:
                        _query = self._nohits_query
                elif str_fields:
                    _query = {
                        "multi_match": {
                            "query": u"{}".format(id),
                            "fields": str_fields,
                            "operator": "and"
                        }
                    }
                else:
                    _query = self._nohits_query

            else:
                raise ValueError('"scopes" cannot be "%s" type' % type(scopes))

        # _query = self.add_species_filter(_query)
        _query = self.add_query_filters(_query)
        _query = self.add_species_custom_filters_score(_query)
        _q = {"query": _query}
        if self._query_options:
            _q.update(self._query_options)

        # if 'fields' in _q and _q['fields'] is not None:
        #     _q['_source'] = _q['fields']
        #     del _q['fields']
        return _q

    def build_genomic_pos_query(self, chr, gstart, gend, assembly=None):
        '''By default if assembly is None, the lastest assembly is used.
           for some species (e.g. human) we support multiple assemblies,
           exact assembly is passed as well.
        '''
        gstart = safe_genome_pos(gstart)
        gend = safe_genome_pos(gend)
        if chr.lower().startswith('chr'):
            chr = chr[3:]

        genomic_pos_field = "genomic_pos"
        if assembly:
            if assembly == 'hg19':
                genomic_pos_field = "genomic_pos_hg19"
            if assembly == 'mm9':
                genomic_pos_field = "genomic_pos_mm9"

        _query = {
            "nested": {
                "path": genomic_pos_field,
                "query": {
                    "bool": {
                        "must": [
                            {"term": {
                                genomic_pos_field + ".chr": chr.lower()}},
                            {"range": {
                                genomic_pos_field + ".start": {"lte": gend}}},
                            {"range": {
                                genomic_pos_field + ".end": {"gte": gstart}}}
                        ]
                    }
                }
            }
        }
        # _query = {
        #     'filtered': {
        #         'query': _query,
        #         'filter' : {
        #             "term" : {"taxid" : taxid}
        #         }
        #     }
        # }
        _query = self.add_query_filters(_query)
        _q = {'query': _query}
        if self._query_options:
            _q.update(self._query_options)
        return _q

# ################# #
# FROM MYGENE.HUB   #
# ################# #

def get_es(es_host=None):
    es_host = es_host or ES_HOST
    es = Elasticsearch(es_host, timeout=600, max_retries=100,retry_on_timeout=True)
    return es


def lastexception():
    exc_type, exc_value, tb = sys.exc_info()
    if exc_type is None:
        print("No exception occurs.")
        return
    print(exc_type.__name__ + ':', end='')
    try:
        excArgs = exc_value.__dict__["args"]
    except KeyError:
        excArgs = ()
    return str(exc_type)+':'+''.join([str(x) for x in excArgs])


class ESIndexer(object):
    def __init__(self, es_index_name=None, es_index_type=None, mapping=None,
                 es_host=None, step=5000):
        self.conn = get_es(es_host)
        self.ES_INDEX_NAME = es_index_name or ES_INDEX_NAME
        self.ES_INDEX_TYPE = es_index_type or ES_INDEX_TYPE
        # if self.ES_INDEX_NAME:
        #    self.conn.default_indices = [self.ES_INDEX_NAME]
        # if self.ES_INDEX_TYPE:
        #    self.conn.default_types = [self.ES_INDEX_TYPE]
        self.step = step
        # self.conn.bulk_size = self.step
        self.number_of_shards = 5      # set number_of_shards when create_index
        # optionally, can specify number of records to skip,
        # useful to continue indexing after an error.
        self.s = None
        self.use_parallel = False
        self._mapping = mapping

    def _get_es_version(self):
        info = self.conn.info()
        return info['version']['number']

    def check(self):
        '''print out ES server info for verification.'''
        # print "Servers:", self.conn.servers
        print("Servers:", '|'.join(["{host}:{port}".format(**h)
                                    for h in self.conn.transport.hosts]))
        print("Default indices:", self.ES_INDEX_NAME)
        print("ES_INDEX_TYPE:", self.ES_INDEX_TYPE)

    def create_index(self, index_name=None, mapping=None):
        index_name = index_name or self.ES_INDEX_NAME
        if not self.conn.indices.exists(index_name):
            body = {
                'settings': {
                    'number_of_shards': self.number_of_shards,
                    # set this to 0 to boost indexing
                    # after indexing, set "auto_expand_replicas": "0-all",
                    #   to make additional replicas.
                    "number_of_replicas": 0,
                }
            }
            if mapping:
                mapping = {"mappings": mapping}
                body.update(mapping)
            print(self.conn.indices.create(index=index_name, body=body))

    def exists_index(self, index):
        return self.conn.indices.exists(index)

    def delete_index_type(self, index_type, noconfirm=False):
        '''Delete all indexes for a given index_type.'''
        index_name = self.ES_INDEX_NAME
        # Check if index_type exists
        m = self.conn.indices.get_mapping(index_name, index_type)
        if not m:
            print('Error: index type "%s" does not exist in index "%s".'
                  % (index_type, index_name))
            return
        path = '/%s/%s' % (index_name, index_type)
        if noconfirm or ask(
                'Confirm to delete all data under "%s":' % path) == 'Y':
            return self.conn.indices.delete_mapping(
                    index=index_name, doc_type=index_type)

    def verify_mapping(self, update_mapping=False):
        '''Verify if index and mapping exist, update mapping if mapping does not exist,
           or "update_mapping=True" explicitly
        '''
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        doc_type = self.ES_INDEX_TYPE

        # Test if index exists
        if not self.exists_index(index_name):
            print('Error: index "%s" does not exist. Create it first.'
                  % index_name)
            return -1

        mapping = conn.indices.get_mapping(index=index_name, doc_type=doc_type)
        empty_mapping = mapping == {}

        if empty_mapping:
            # if no existing mapping available for index_type
            # force update_mapping to True
            update_mapping = True

        if update_mapping:
            print("Updating mapping...", end='')
            if not empty_mapping:
                print("\n\tRemoving existing mapping...", end='')
                print(conn.indices.delete_mapping(
                    index=index_name, doc_type=doc_type))
            self.get_field_mapping()
            print(conn.indices.put_mapping(index=index_name,
                                           doc_type=doc_type,
                                           body=self._mapping))

    def update_mapping_meta(self, meta):
        index_name = self.ES_INDEX_NAME
        doc_type = self.ES_INDEX_TYPE

        if isinstance(meta, dict) and len(set(meta) - set(
                      ['_meta', '_timestamp'])) == 0:
            body = {doc_type: meta}
            print(self.conn.indices.put_mapping(
                index=index_name,
                doc_type=doc_type,
                body=body)
            )
        else:
            raise ValueError(
                    'Input "meta" should have and only have "_meta" field.')

    def count(self, query=None, index_type=None):
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        index_type = index_type or self.ES_INDEX_TYPE
        if query:
            if isinstance(query, dict) and 'query' in query:
                return conn.count(index_name, index_type, query)
            else:
                raise ValueError("Not a valid input query")
        else:
            return conn.count(index_name, index_type)

    def get(self, id, **kwargs):
        '''get a specific doc by its id.'''
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        index_type = self.ES_INDEX_TYPE
        return conn.get(index_name, id, index_type, **kwargs)

    def get_docs(self, ids, step=None, **kwargs):
        '''return matching docs for given ids, if not found return None.
           A generator is returned and the order is perserved.
        '''
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        index_type = self.ES_INDEX_TYPE
        step = step or self.step
        for i in range(0, len(ids), step):
            _ids = ids[i:i + step]
            body = {'ids': _ids}
            res = conn.mget(
                    body=body, index=index_name, doc_type=index_type, **kwargs)
            for doc in res['docs']:
                if doc['found']:
                    yield doc['_source']
                else:
                    yield None

    def index(self, doc, id=None):
        '''add a doc to the index. If id is not None, the existing doc will be
           updated.
        '''
        return self.conn.index(
                self.ES_INDEX_NAME, self.ES_INDEX_TYPE, doc, id=id)

    def index_bulk(self, docs, step=None):
        index_name = self.ES_INDEX_NAME
        doc_type = self.ES_INDEX_TYPE
        step = step or self.step

        def _get_bulk(doc):
            doc.update({
                "_index": index_name,
                "_type": doc_type,
            })
            return doc
        actions = (_get_bulk(doc) for doc in docs)
        try:
            return helpers.bulk(self.conn, actions, chunk_size=step)
        except helpers.BulkIndexError as e:
            # try again...
            print("Bulk error, try again...")
            return self.index_bulk(docs,step)
            ##return helpers.bulk(self.conn, actions, chunk_size=step)
        except Exception as e:
            print("Err...")
            import pickle
            pickle.dump(e,open("err","wb"))

    def add_docs(self, docs, step=None):
        self.index_bulk(docs, step=step)
        self.conn.indices.flush()
        self.conn.indices.refresh()

    def delete_doc(self, index_type, id):
        '''delete a doc from the index based on passed id.'''
        return self.conn.delete(self.ES_INDEX_NAME, index_type, id)

    def delete_docs(self, ids, step=None):
        index_name = self.ES_INDEX_NAME
        doc_type = self.ES_INDEX_TYPE
        step = step or self.step

        def _get_bulk(_id):
            doc = {
                '_op_type': 'delete',
                "_index": index_name,
                "_type": doc_type,
                "_id": _id
            }
            return doc
        actions = (_get_bulk(_id) for _id in ids)
        return helpers.bulk(self.conn, actions, chunk_size=step,
                            stats_only=True, raise_on_error=False)

    def update(self, id, extra_doc, index_type=None, bulk=False):
        '''update an existing doc with extra_doc.'''
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        index_type = index_type or self.ES_INDEX_TYPE
        # old way, update locally and then push it back.
        # return self.conn.update(extra_doc, self.ES_INDEX_NAME,
        #                         index_type, id)

        if not bulk:
            body = {'doc': extra_doc}
            return conn.update(index_name, index_type, id, body)
        else:
            raise NotImplementedError
            '''
            # ES supports bulk update since v0.90.1.
            op_type = 'update'
            cmd = {op_type: {"_index": index_name,
                             "_type": index_type,
                             "_id": id}
                   }

            doc = json.dumps({"doc": extra_doc}, cls=conn.encoder)
            command = "%s\n%s" % (json.dumps(cmd, cls=conn.encoder), doc)
            conn.bulker.add(command)
            return conn.flush_bulk()
            '''

    def update_docs(self, partial_docs, **kwargs):
        index_name = self.ES_INDEX_NAME
        doc_type = self.ES_INDEX_TYPE

        def _get_bulk(doc):
            doc = {
                '_op_type': 'update',
                "_index": index_name,
                "_type": doc_type,
                "_id": doc['_id'],
                "doc": doc
            }
            return doc
        actions = (_get_bulk(doc) for doc in partial_docs)
        return helpers.bulk(self.conn, actions, chunk_size=self.step, **kwargs)

    def wait_till_all_shards_ready(self, timeout=None, interval=5):
        raise NotImplementedError

    def optimize(self):
        '''optimize the default index.'''
        return self.conn.indices.optimize(self.ES_INDEX_NAME,
                                          wait_for_merge=False,   # True,
                                          max_num_segments=5)

    def optimize_all(self):
        """optimize all indices"""
        return self.conn.indices.optimize('', wait_for_merge=False,  # True,
                                          max_num_segments=5)

    def get_field_mapping(self):
        # raise NotImplementedError
        return self._mapping

    def build_index(self, collection, update_mapping=False, verbose=False,
                    query=None, bulk=True):
        conn = self.conn
        index_name = self.ES_INDEX_NAME
        # index_type = self.ES_INDEX_TYPE

        self.verify_mapping(update_mapping=update_mapping)
        # update some settings for bulk indexing
        body = {
            "index": {
                # disable refresh temporarily
                "refresh_interval": "-1",
                "auto_expand_replicas": "0-all",
                # "number_of_replicas": 0,
                #"refresh_interval": "30s",
            }
        }
        conn.indices.put_settings(body, index_name)
        try:
            print("Building index...")
            if self.use_parallel:
                cnt = self._build_index_parallel(collection, verbose)
            else:
                cnt = self._build_index_sequential(collection, verbose,
                                                   query=query, bulk=bulk)
        finally:
            # restore some settings after bulk indexing is done.
            body = {
                "index": {
                    "refresh_interval": "1s"              # default settings
                }
            }
            conn.indices.put_settings(body, index_name)

            try:
                print("Flushing...", conn.indices.flush())
                print("Refreshing...", conn.indices.refresh())
            except:
                pass

            time.sleep(1)
            print("Validating...", end='')
            target_cnt = collection.find(query).count()
            es_cnt = self.count()['count']
            if target_cnt == es_cnt:
                print("OK [total count={}]".format(target_cnt))
            else:
                print("\nWarning: total count of gene documents does not " +
                      "match [{}, should be {}]".format(es_cnt, target_cnt))

        if cnt:
            print('Done! - {} docs indexed.'.format(cnt))
            print("Optimizing...", self.optimize())

    def _build_index_sequential(self, collection, verbose=False,
                                query=None, bulk=True):

        def rate_control(cnt, t):
            delay = 0
            if t > 90:
                delay = 30
            elif t > 60:
                delay = 10
            if delay:
                print("\tPausing for {}s...".format(delay), end='')
                time.sleep(delay)
                print("done.")

        src_docs = doc_feeder(collection, step=self.step, s=self.s,
                              batch_callback=rate_control, query=query)
        if bulk:
            res = self.index_bulk(src_docs)
            if len(res[1]) > 0:
                print("Error: {} docs failed indexing.".format(len(res[1])))
            return res[0]
        else:
            cnt = 0
            for doc in src_docs:
                self.index(doc)
                cnt += 1
                if verbose:
                    print(cnt, ':', doc['_id'])
            return cnt

    def _build_index_parallel(self, collection, verbose=False):
        raise NotImplementedError
        from utils.parallel import (run_jobs_on_ipythoncluster,
                                    collection_partition,
                                    require)
        kwargs_common = {'ES_HOST': ES_HOST,
                         'ES_INDEX_NAME': self.ES_INDEX_NAME,
                         'ES_INDEX_TYPE': self.ES_INDEX_TYPE,
                         }
        task_list = []
        for kwargs in collection_partition(collection, step=self.step):
            kwargs.update(kwargs_common)
            task_list.append(kwargs)

        @require('pymongo', 'pyes')
        def worker(kwargs):
            import pymongo
            import pyes
            server = kwargs['server']
            port = kwargs['port']
            src_db = kwargs['src_db']
            src_collection = kwargs['src_collection']
            skip = kwargs['skip']
            limit = kwargs['limit']

            mongo_conn = pymongo.MongoClient(server, port)
            src = mongo_conn[src_db]

            ES_HOST = kwargs['ES_HOST']
            ES_INDEX_NAME = kwargs['ES_INDEX_NAME']
            ES_INDEX_TYPE = kwargs['ES_INDEX_TYPE'],

            es_conn = pyes.ES(ES_HOST, default_indices=[ES_INDEX_NAME],
                              timeout=120.0, max_retries=10)

            cur = src[src_collection].find(skip=skip, limit=limit,
                                           timeout=False)
            cur.batch_size(1000)
            cnt = 0
            try:
                for doc in cur:
                    es_conn.index(doc, ES_INDEX_NAME, ES_INDEX_TYPE,
                                  doc['_id'], bulk=True)
                    cnt += 1
            finally:
                cur.close()
            es_conn.indices.flush()   # this is important to avoid missing docs
            es_conn.indices.refresh()
            return cnt

        job_results = run_jobs_on_ipythoncluster(worker, task_list)
        if job_results:
            cnt = sum(job_results)
            return cnt

    def doc_feeder(self, index_type=None, index_name=None, step=10000,
                   verbose=True, query=None, scroll='10m', **kwargs):
        conn = self.conn
        index_name = index_name or self.ES_INDEX_NAME
        doc_type = index_type or self.ES_INDEX_TYPE

        n = self.count(query=query)['count']
        cnt = 0
        t0 = time.time()
        if verbose:
            print('\ttotal docs: {}'.format(n))

        _kwargs = kwargs.copy()
        _kwargs.update(dict(size=step, index=index_name, doc_type=doc_type))
        res = helpers.scan(conn, query=query, scroll=scroll, **_kwargs)
        t1 = time.time()
        for doc in res:
            if verbose and cnt % step == 0:
                if cnt != 0:
                    print('done.[%.1f%%,%s]' % (cnt*100./n, timesofar(t1)))
                print('\t{}-{}...'.format(cnt+1, min(cnt+step, n)), end='')
                t1 = time.time()
            yield doc
            cnt += 1
        if verbose:
            print('done.[%.1f%%,%s]' % (cnt*100./n, timesofar(t1)))
            print("Finished! [{}]".format(timesofar(t0)))

    def get_id_list(self, index_type=None, index_name=None, step=100000,
                    verbose=True):
        cur = self.doc_feeder(index_type=index_type, index_name=index_name,
                              step=step, fields=[], verbose=verbose)
        id_li = [doc['_id'] for doc in cur]
        return id_li

    def get_id_list_parallel(self, taxid_li, index_type=None, index_name=None,
                             step=1000, verbose=True):
        '''return a list of all doc ids in an index_type.'''
        raise NotImplementedError
        from utils.parallel import run_jobs_on_ipythoncluster

        def _get_ids_worker(args):
            from utils.es import ESIndexer
            from pyes import MatchAllQuery
            es_kwargs, start, step = args
            q = MatchAllQuery().search()
            q.sort = [{'entrezgene': 'asc'}, {'ensembl.gene': 'asc'}]
            q.fields = []
            q.start = start
            q.size = step
            esi = ESIndexer(**es_kwargs)
            cnt = esi.count()['count']
            res = esi.conn.search_raw(q)
            assert res['hits']['total'] == cnt
            return [doc['_id'] for doc in res['hits']['hits']]

        def _get_ids_worker_by_taxid(args):
            from utils.es import ESIndexer
            from pyes import TermQuery
            es_kwargs, taxid, step = args
            q = TermQuery()
            q.add('taxid', taxid)
            q.fields = []
            q.size = step
            esi = ESIndexer(**es_kwargs)
            res = esi.conn.search(q)
            xli = [doc['_id'] for doc in res]
            assert len(xli) == res.total
            return xli

        es_kwargs = {'es_index_name': self.ES_INDEX_NAME,
                     'es_host': 'su02:9200'}
        task_li = [(es_kwargs, taxid, step) for taxid in taxid_li]
        # print task_li
        job_results = run_jobs_on_ipythoncluster(
                _get_ids_worker_by_taxid, task_li)
        return job_results

    def clone_index(self, src_index, target_index, target_es_host=None,
                    step=10000, scroll='10m', target_index_settings=None,
                    number_of_shards=None):
        '''clone src_index to target_index on the same es_host, or another one given
           by target_es_host.

           This method can now be replaced by helpers.reindex
        '''


def es_clean_indices(keep_last=2, es_host=None, verbose=True, noconfirm=False,
                     dryrun=False):
    '''clean up es indices, only keep last <keep_last> number of indices.'''
    conn = get_es(es_host)
    index_li = list(conn.indices.get_aliases().keys())
    if verbose:
        print("Found {} indices".format(len(index_li)))

    for prefix in ('genedoc_mygene', 'genedoc_mygene_allspecies'):
        pat = prefix + '_(\d{8})_\w{8}'
        _li = []
        for index in index_li:
            mat = re.match(pat, index)
            if mat:
                _li.append((mat.group(1), index))
        _li.sort()   # older collection appears first
        # keep last # of newer indices
        index_to_remove = [x[1] for x in _li[:-keep_last]]
        if len(index_to_remove) > 0:
            print("{} \"{}*\" indices will be removed.".format(
                  len(index_to_remove), prefix))
            if verbose:
                for index in index_to_remove:
                    print('\t', index)
            if noconfirm or ask("Continue?") == 'Y':
                for index in index_to_remove:
                    if dryrun:
                        print("dryrun=True, nothing is actually deleted")
                    else:
                        conn.indices.delete(index)
                print("Done.[%s indices removed]" % len(index_to_remove))
            else:
                print("Aborted.")
        else:
            print("Nothing needs to be removed.")


def get_latest_indices(es_host=None):
    conn = get_es(es_host)
    index_li = list(conn.indices.get_aliases().keys())
    print("Found {} indices".format(len(index_li)))

    latest_indices = []
    for prefix in ('genedoc_mygene', 'genedoc_mygene_allspecies'):
        pat = prefix + '_(\d{8})_\w{8}'
        _li = []
        for index in index_li:
            mat = re.match(pat, index)
            if mat:
                _li.append((mat.group(1), index))
        if not _li:
            print("Warning: no matching indices found!")
            continue
        latest_indices.append(sorted(_li)[-1])

    if len(latest_indices) == 2:
        if latest_indices[0][0] != latest_indices[1][0]:
            print("Warning: unmatched timestamp:")
            print('\n'.join([x[1] for x in latest_indices]))
        latest_indices = [x[1] for x in latest_indices]
        return latest_indices
