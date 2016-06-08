import random, os, httplib2
from biothings.tests.test_helper import BiothingTestHelperMixin, _d, TornadoRequestHelper
from nose.tools import ok_, eq_

from tornado.testing import AsyncHTTPTestCase
import www.index as index
from biothings.settings import BiothingSettings
import config


class MyGeneTest(BiothingTestHelperMixin):
    __test__ = True  # explicitly set this to be a test class

    #############################################################
    # Test functions                                            #
    #############################################################

    host = os.getenv(config.HOST_ENVAR_NAME,"")
    #if not host:
    #    raise ValueError("Missing HOST_ENVAR_NAME")
    host = host.rstrip('/')
    api = host + '/' + config.API_VERSION
    h = httplib2.Http()

    def _filter_hits(self, res, field=None):
        for hit in res.get("hits"):
            if field:
                del hit[field]
            else:
                # because can't remove elem from dict while iterate
                # need to keep track on what should be popped...
                topop = []
                for k in hit.keys():
                    if k.startswith("_"):
                        topop.append(k)
                [hit.pop(i) for i in topop]

    def json_ok(self, url, filter=False, **kwargs):
        res = super(MyGeneTest, self).json_ok(url,**kwargs)
        if filter:
            self._filter_hits(res)
        return res

    def test_main(self):
        self.get_ok(self.host)

    def test_gene_object(self):
        # test all fields are load in gene objects
        res = self.json_ok(self.get_ok(self.api + '/gene/1017'))

        attr_li = ['HGNC', 'HPRD', 'MIM', 'Vega', '_id', 'accession', 'alias',
                   'ec', 'ensembl', 'entrezgene', 'genomic_pos', 'go',
                   'homologene', 'interpro', 'ipi', 'map_location', 'name',
                   'pdb', 'pharmgkb', 'pir', 'prosite', 'reagent', 'refseq',
                   'reporter', 'summary', 'symbol', 'taxid', 'type_of_gene',
                   'unigene', 'uniprot', 'exons', 'generif']

        for attr in attr_li:
            assert res.get(attr, None) is not None, \
                'Missing field "{}" in gene "1017"'.format(attr)

        res = self.json_ok(self.get_ok(self.api + '/gene/12566'))
        attr = 'MGI'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "12566"'.format(attr)

        res = self.json_ok(self.get_ok(self.api + '/gene/245962'))
        attr = 'RGD'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "245962"'.format(attr)

        res = self.json_ok(self.get_ok(self.api + '/gene/493498'))
        attr = 'Xenbase'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "493498"'.format(attr)

        res = self.json_ok(self.get_ok(self.api + '/gene/406715'))
        attr = 'ZFIN'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "406715"'.format(attr)

        res = self.json_ok(self.get_ok(self.api + '/gene/824036'))
        attr = 'TAIR'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "824036"'.format(attr)

        res = self.json_ok(self.get_ok(self.api + '/gene/42453'))
        attr = 'FLYBASE'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "42453"'.format(attr)

        # pig
        res = self.json_ok(self.get_ok(self.api + '/gene/397593'))
        assert 'snowball' in res.get('reporter', {}), \
            'Missing field "reporter.snowball" in gene "397593"'

        # nematode
        res = self.json_ok(self.get_ok(self.api + '/gene/172677'))
        # this is not nematode, "taxid": 31234
        res = self.json_ok(self.get_ok(self.api + '/gene/9821293'))
        attr = 'WormBase'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "9821293"'.format(attr)

        # fission yeast
        res = self.json_ok(self.get_ok(self.api + '/gene/2539869'))

        # e coli.
        # res = self.json_ok(self.get_ok(self.api + '/gene/12931566'))

        # mirna
        res = self.json_ok(self.get_ok(self.api + '/gene/406881'))
        attr = 'miRBase'
        assert res.get(attr, None) is not None, \
            'Missing field "{}" in gene "406881"'.format(attr)

    def test_query(self):
        # public query api at /query via get
        # self.json_ok(self.get_ok(self.api + '/query?q=cdk2'))
        # self.json_ok(self.get_ok(self.api + '/query?q=GO:0004693'))
        # self.json_ok(self.get_ok(self.api + '/query?q=211803_at'))
        self.query_has_hits('cdk2')
        self.query_has_hits('GO:0004693')
        self.query_has_hits('211803_at')
        self.query_has_hits('IPR008351')
        self.query_has_hits('hsa-mir-503')
        self.query_has_hits('hsa-miR-503')

        # test fielded query
        self.query_has_hits('symbol:cdk2')
        # test interval query
        self.query_has_hits('chr1:151,073,054-151,383,976&species=human')

        con = self.get_ok(self.api + '/query?q=cdk2&callback=mycallback')
        ok_(con.startswith(b'mycallback('))

        # testing non-ascii character
        res = self.json_ok(self.get_ok(self.api +
                           '/query?q=54097\xef\xbf\xbd\xef\xbf\xbdmouse'))
        eq_(res['hits'], [])

        res = self.json_ok(self.get_ok(self.api + '/query'), checkerror=False)
        assert 'error' in res

        res = self.json_ok(self.get_ok(self.api + '/query?q=tRNA:Y1:85Ae'),
                           checkerror=False)
        assert 'error' in res
        # ensure returned fields by default
        res = self.json_ok(self.get_ok(self.api + '/query?q=cdk'))
        # pick one
        idx = random.randrange(0, 10)
        deffields = res["hits"][idx].keys()  # pick one...
        expected = ["_id", "_score", "taxid", "entrezgene", "name", "symbol"]
        assert sorted(list(deffields)) == sorted(expected), \
            "%s != %s" % (sorted(list(deffields)), sorted(expected))

    def test_query_post(self):
        # /query via post
        self.json_ok(self.post_ok(self.api + '/query', {'q': '1017'}))

        res = self.json_ok(self.post_ok(self.api + '/query',
                                        {'q': '1017', 'scopes': 'entrezgene'}))
        eq_(len(res), 1)
        eq_(res[0]['_id'], '1017')

        res = self.json_ok(self.post_ok(self.api + '/query',
                                        {'q': '211803_at,1018',
                                         'scopes': 'reporter,entrezgene'}))
        eq_(len(res), 2)
        eq_(res[0]['_id'], '1017')
        eq_(res[1]['_id'], '1018')

        res = self.json_ok(self.post_ok(self.api + '/query',
                                        {'q': 'CDK2',
                                         'species': 'human,10090,frog,pig',
                                         'scopes': 'symbol',
                                         'fields': 'name,symbol'}))
        assert len(res) >= 4, (res, len(res))
        res = self.json_ok(self.post_ok(self.api + '/query', {}),
                           checkerror=False)
        assert 'error' in res, res

        res = self.json_ok(self.post_ok(self.api + '/query',
                                        {'q': '[1017, "1018"]',
                                         'scopes': 'entrezgene',
                                         'jsoninput': 'true'}))
        eq_(len(res), 2)
        eq_(res[0]['_id'], '1017')
        eq_(res[1]['_id'], '1018')

    def test_query_interval(self):
        res = self.json_ok(self.get_ok(self.api +
                           '/query?q=chr1:1000-100000&species=human'))
        ok_(len(res['hits']) > 1)
        ok_('_id' in res['hits'][0])

    def test_query_size(self):
        res = self.json_ok(self.get_ok(self.api + '/query?q=cdk?'))
        eq_(len(res['hits']), 10)  # default is 10
        ok_(res['total'] > 10)

        res = self.json_ok(self.get_ok(self.api + '/query?q=cdk?&size=0'))
        eq_(len(res['hits']), 0)

        res = self.json_ok(self.get_ok(self.api + '/query?q=cdk?&limit=20'))
        eq_(len(res['hits']), 20)

        res1 = self.json_ok(self.get_ok(self.api +
                            '/query?q=cdk?&from=0&size=20'))
        res = self.json_ok(self.get_ok(self.api +
                           '/query?q=cdk?&skip=10&size=20'))
        eq_(len(res['hits']), 20)
        # print res1['hits'].index(res['hits'][0])
        # print [x['_id'] for x in res1['hits']]
        # eq_(res['hits'][0], res1['hits'][10])
        assert res['hits'][0] in res1['hits']

        # API doc says cap 1000
        res = self.json_ok(self.get_ok(self.api + '/query?q=*&size=1000'))
        eq_(len(res['hits']), 1000)
        res = self.json_ok(self.get_ok(self.api + '/query?q=*&size=1001'))
        eq_(len(res['hits']), 1000)
        res = self.json_ok(self.get_ok(self.api + '/query?q=*&size=2000'))
        eq_(len(res['hits']), 1000)

        # assert 1==0
        res = self.json_ok(self.get_ok(self.api + '/query?q=cdk?&size=1a'),
                           checkerror=False)  # invalid size parameter
        assert 'error' in res

    def test_gene(self):
        res = self.json_ok(self.get_ok(self.api + '/gene/1017'))
        eq_(res['entrezgene'], 1017)
        # testing non-ascii character
        self.get_404(self.api + '/gene/' +
                     '54097\xef\xbf\xbd\xef\xbf\xbdmouse')

        # commented out this test, as no more
        # allow dot in the geneid
        # res = self.json_ok(self.get_ok(self.api + '/gene/Y105C5B.255'))

        # testing filtering parameters
        res = self.json_ok(self.get_ok(self.api +
                           '/gene/1017?fields=symbol,name,entrezgene'))
        eq_(set(res), set(['_id', '_score', 'symbol', 'name', 'entrezgene']))
        res = self.json_ok(self.get_ok(self.api +
                           '/gene/1017?filter=symbol,go.MF'))
        eq_(set(res), set(['_id', '_score', 'symbol', 'go']))
        assert "MF" in res["go"]

        self.get_404(self.api + '/gene')
        self.get_404(self.api + '/gene/')

    def test_gene_post(self):
        res = self.json_ok(self.post_ok(self.api + '/gene', {'ids': '1017'}))
        eq_(len(res), 1)
        eq_(res[0]['entrezgene'], 1017)

        res = self.json_ok(self.post_ok(self.api + '/gene',
                                        {'ids': '1017, 1018'}))
        eq_(len(res), 2)
        eq_(res[0]['_id'], '1017')
        eq_(res[1]['_id'], '1018')

        res = self.json_ok(self.post_ok(self.api + '/gene',
                           {'ids': '1017,1018',
                            'fields': 'symbol,name,entrezgene'}))
        eq_(len(res), 2)
        for _g in res:
            eq_(set(_g), set(['_id', '_score', 'query', 'symbol',
                              'name', 'entrezgene']))

        res = self.json_ok(self.post_ok(self.api + '/gene',
                           {'ids': '1017,1018',
                            'filter': 'symbol,go.MF'}))
        eq_(len(res), 2)
        for _g in res:
            eq_(set(_g), set(['_id', '_score', 'query', 'symbol', 'go']))
            assert "MF" in _g["go"]

        # get retired gene (make sure _search ES query is run)
        res = self.json_ok(self.post_ok(self.api + '/gene',{'ids': '791256'}))
        eq_(res[0]['_id'], '50846')  # this is the corresponding _id field

    def test_status(self):
        # /status
        self.get_ok(self.host + '/status')
        self.head_ok(self.host + '/status')

    def test_metadata(self):
        root = self.json_ok(self.get_ok(self.host + '/metadata'))
        v2 = self.json_ok(self.get_ok(self.api + '/metadata'))
        eq_(root, v2)
        eq_(set(root.keys()), set(['available_fields', 'src_version',
                                   'app_revision', 'timestamp', 'taxonomy',
                                   'stats', 'genome_assembly', 'source']))
        fields = self.json_ok(self.get_ok(self.api + '/metadata/fields'))
        # test random field
        assert "refseq" in fields
        assert "accession.rna" in fields
        assert "interpro.desc" in fields
        assert "homologene" in fields
        assert "reporter.snowball" in fields

    def test_query_facets(self):
        res = self.json_ok(self.get_ok(self.api +
                                       '/query?q=cdk?&facets=taxid'))
        ok_('facets' in res)
        ok_('taxid' in res['facets'])
        eq_(res['facets']['taxid']['total'], res['total'])
        eq_(res['facets']['taxid']['other'], 0)
        eq_(res['facets']['taxid']['missing'], 0)

        u = '/query?q=cdk?&facets=taxid&species_facet_filter=human'
        res2 = self.json_ok(self.get_ok(self.api + u))
        eq_(res2['facets']['taxid']['total'], res['total'])
        eq_(res2['facets']['taxid'], res['facets']['taxid'])
        eq_([x["count"] for x in res2['facets']['taxid']['terms']
            if x["term"] == 9606][0], res2['total'])

    def test_query_userfilter(self):
        res1 = self.json_ok(self.get_ok(self.api + '/query?q=cdk'))
        res2 = self.json_ok(self.get_ok(self.api +
                            '/query?q=cdk&userfilter=bgood_cure_griffith'))
        ok_(res1['total'] > res2['total'])

        # nonexisting user filter gets ignored.
        res2 = self.json_ok(self.get_ok(self.api +
                            '/query?q=cdk&userfilter=aaaa'))
        eq_(res1['total'], res2['total'])

    def test_existsfilter(self):
        res1 = self.json_ok(self.get_ok(self.api + '/query?q=cdk'))
        res2 = self.json_ok(self.get_ok(self.api +
                            '/query?q=cdk&exists=pharmgkb'))
        ok_(res1['total'] > res2['total'])
        res3 = self.json_ok(self.get_ok(self.api +
                            '/query?q=cdk&exists=pharmgkb,pdb'))
        ok_(res2['total'] > res3['total'])

    def test_missingfilter(self):
        res1 = self.json_ok(self.get_ok(self.api + '/query?q=cdk'))
        res2 = self.json_ok(self.get_ok(self.api + '/query?q=cdk&missing=pdb'))
        ok_(res1['total'] > res2['total'])
        res3 = self.json_ok(self.get_ok(self.api +
                            '/query?q=cdk&missing=pdb,MIM'))
        ok_(res2['total'] > res3['total'])

    def test_unicode(self):
        s = u'基因'

        self.get_404(self.api + '/gene/' + s)

        res = self.json_ok(self.post_ok(self.api + '/gene', {'ids': s}))
        eq_(res[0]['notfound'], True)
        eq_(len(res), 1)
        res = self.json_ok(self.post_ok(self.api + '/gene',
                           {'ids': '1017, ' + s}))
        eq_(res[1]['notfound'], True)
        eq_(len(res), 2)

        res = self.json_ok(self.get_ok(self.api + '/query?q=' + s))
        eq_(res['hits'], [])

        res = self.json_ok(self.post_ok(self.api + '/query',
                           {"q": s, "scopes": 'symbol'}))
        eq_(res[0]['notfound'], True)
        eq_(len(res), 1)

        res = self.json_ok(self.post_ok(self.api + '/query',
                           {"q": 'cdk2+' + s}))
        eq_(res[1]['notfound'], True)
        eq_(len(res), 2)

    def test_hg19(self):
        u = '/query?q=hg19.chr12:57,795,963-57,815,592&species=human'
        res = self.json_ok(self.get_ok(self.api + u))

        ok_(len(res['hits']) == 2)
        ok_('_id' in res['hits'][0])
        u = '/query?q=chr12:57,795,963-57,815,592&species=human'
        res2 = self.json_ok(self.get_ok(self.api + u))
        ok_(res['total'] != res2['total'])

        u = '/gene/10017?fields=genomic_pos_hg19,exons_hg19'
        res = self.json_ok(self.get_ok(self.api + u))
        ok_('genomic_pos_hg19' in res)
        ok_('exons_hg19' in res)

    def test_mm9(self):
        u = '/query?q=mm9.chr12:57,795,963-57,815,592&species=mouse'
        res = self.json_ok(self.get_ok(self.api + u))
        ok_(len(res['hits']) == 2)
        ok_('_id' in res['hits'][0])
        u = '/query?q=chr12:57,795,963-57,815,592&species=mouse'
        res2 = self.json_ok(self.get_ok(self.api + u))
        ok_(res['total'] != res2['total'])
        u = '/gene/12049?fields=genomic_pos_mm9,exons_mm9'
        res = self.json_ok(self.get_ok(self.api + u))
        ok_('genomic_pos_mm9' in res)
        ok_('exons_mm9' in res)

    def test_msgpack(self):
        res = self.json_ok(self.get_ok(self.api + '/gene/1017'))
        res2 = self.msgpack_ok(self.get_ok(self.api +
                               '/gene/1017?msgpack=true'))
        ok_(res, res2)

        res = self.json_ok(self.get_ok(self.api + '/query/?q=cdk'))
        res2 = self.msgpack_ok(self.get_ok(self.api +
                               '/query/?q=cdk&msgpack=true'))
        ok_(res, res2)

        res = self.json_ok(self.get_ok(self.api + '/metadata'))
        res2 = self.msgpack_ok(self.get_ok(self.api +
                               '/metadata?msgpack=true'))
        ok_(res, res2)

    def test_taxonomy(self):
        res = self.json_ok(self.get_ok(self.api + '/species/1239'))
        ok_("lineage" in res)

        res = self.json_ok(self.get_ok(self.api +
                           '/species/46170?include_children=true'))
        ok_(len(res['children']) >= 305)

        res2 = self.json_ok(self.get_ok(self.api +
                            '/species/46170?include_children=true&has_gene=1'))
        ok_(len(res2['children']) >= 16)
        ok_(len(res2['children']) <= len(res['children']))

        u = '/query?q=lytic%20enzyme&species=1386&include_tax_tree=true'
        res = self.json_ok(self.get_ok(self.api + u))
        ok_(res['total'] >= 2)
        res2 = self.json_ok(self.get_ok(self.api +
                            '/query?q=lytic%20enzyme&species=1386'))
        eq_(res2['total'], 0)

    def test_static(self):
        self.get_ok(self.host + '/favicon.ico')
        self.get_ok(self.host + '/robots.txt')

    def test_fetch_all(self):
        res = self.json_ok(self.get_ok(self.api +
                           '/query?q=cdk2&fetch_all=true'))
        assert '_scroll_id' in res

        res2 = self.json_ok(self.get_ok(self.api +
                            '/query?scroll_id=' + res['_scroll_id']))
        assert 'hits' in res2
        ok_(len(res2['hits']) >= 2)

    def test_dotfield(self):
        # /query service
        # default dotfield=0
        rdefault = self.json_ok(self.get_ok(self.api +
                                '/query?q=ccnk&fields=refseq.rna'))
        # force no dotfield
        rfalse = self.json_ok(self.get_ok(self.api +
                              '/query?q=ccnk&fields=refseq.rna&dotfield=false'))
        # force dotfield
        rtrue = self.json_ok(self.get_ok(self.api +
                             '/query?q=ccnk&fields=refseq.rna&dotfield=true'))
        # check defaults and bool params
        # TODO: put this in json_ok as post-process filter ?
        for d in [rdefault,rfalse,rtrue]:
            for h in d["hits"]:
                del h["_score"]
        eq_(rdefault["hits"], rfalse["hits"])
        # check struct
        assert "refseq.rna" in rtrue["hits"][0].keys()
        assert "refseq" in rdefault["hits"][0].keys()
        assert "rna" in rdefault["hits"][0]["refseq"].keys()
        # TODO: no fields but dotfield => dotfield results
        # TODO: fields with dot but no dotfield => dotfield results

        # /gene service
        rdefault = self.json_ok(self.get_ok(self.api +
                                '/gene/1017?filter=symbol,go.MF'))
        rtrue = self.json_ok(self.get_ok(self.api +
                             '/gene/1017?filter=symbol,go.MF&dotfield=true'))
        rfalse = self.json_ok(self.get_ok(self.api +
                              '/gene/1017?filter=symbol,go.MF&dotfield=false'))
        eq_(rdefault, rfalse)
        assert "go.MF" in rtrue.keys()
        assert "go" in rdefault.keys()
        assert "MF" in rdefault["go"].keys()

    def test_raw(self):
        # /gene
        raw1 = self.json_ok(self.get_ok(self.api + '/gene/1017?raw=1'))
        rawtrue = self.json_ok(self.get_ok(self.api + '/gene/1017?raw=true'))
        raw0 = self.json_ok(self.get_ok(self.api + '/gene/1017?raw=0'))
        rawfalse = self.json_ok(self.get_ok(self.api + '/gene/1017?raw=false'))
        eq_(sorted(raw1), sorted(rawtrue))
        eq_(raw0, rawfalse)
        assert "_shards" in raw1
        assert "_shards" not in raw0
        assert "timed_out" in raw1
        assert "timed_out" not in raw0
        # /query
        raw1 = self.json_ok(self.get_ok(self.api + '/query?q=ccnk&raw=1'))
        rawtrue = self.json_ok(self.get_ok(self.api + '/query?q=ccnk&raw=true'))
        raw0 = self.json_ok(self.get_ok(self.api + '/query?q=ccnk&raw=0'))
        rawfalse = self.json_ok(self.get_ok(self.api + '/query?q=ccnk&raw=false'))
        # this may vary so remove in comparison
        for d in [raw1, rawtrue, raw0, rawfalse]:
            del d["took"]
        # score should be the same. approx... so remove
        for d in [raw1,rawtrue]:
            for h in d["hits"]["hits"]:
                del h["_score"]
            del d["hits"]["max_score"]
        for d in [raw0,rawfalse]:
            for h in d["hits"]:
                del h["_score"]
            del d["max_score"]
        eq_(raw1, rawtrue)
        eq_(raw0, rawfalse)
        assert "_shards" in raw1
        assert "_shards" not in raw0

    def test_species(self):
        res, con = self.h.request(self.api + "/species/9606")
        eq_(res.status, 200)
        d = _d(con.decode('utf-8'))
        eq_(set(d.keys()), set(['taxid', 'authority', 'lineage', '_id',
                                'common_name', 'genbank_common_name',
                                '_version', 'parent_taxid', 'scientific_name',
                                'has_gene', 'children', 'rank',
                                'uniprot_name']))

    def test_query_dotstar_refseq(self):
        protein = self.json_ok(self.get_ok(self.api +
                               "/query?q=refseq:NP_001670&fields=refseq"),
                               filter=True)
        u = "/query?q=refseq:NM_001679&fields=refseq"
        rna = self.json_ok(self.get_ok(self.api + u), filter=True)
        genomic = self.json_ok(self.get_ok(self.api +
                               "/query?q=refseq:NT_005612&fields=refseq"),
                               filter=True)
        u = "/query?q=refseq.protein:NP_001670&fields=refseq"
        explicit_protein = self.json_ok(self.get_ok(self.api + u), filter=True)
        self._filter_hits(explicit_protein)
        u = "/query?q=refseq.rna:NM_001679&fields=refseq"
        explicit_rna = self.json_ok(self.get_ok(self.api + u), filter=True)
        u = "/query?q=refseq.genomic:NT_005612&fields=refseq"
        explicit_genomic = self.json_ok(self.get_ok(self.api + u), filter=True)
        eq_(protein["hits"], explicit_protein["hits"])
        eq_(rna["hits"], explicit_rna["hits"])
        eq_(genomic["hits"], explicit_genomic["hits"])
        eq_(protein["hits"], rna["hits"])  # same result whatever the query
        eq_(genomic["hits"], [])  # genomic not indexed
        eq_(rna["total"], 1)
        hit = rna["hits"][0]
        eq_(hit["refseq"]["protein"], "NP_001670")
        eq_(hit["refseq"]["rna"], "NM_001679")

    def test_query_dotstar_accession(self):
        protein = self.json_ok(self.get_ok(self.api +
                               "/query?q=accession:AAH68303&fields=accession"),
                               filter=True)
        rna = self.json_ok(self.get_ok(self.api +
                           "/query?q=accession:BC068303&fields=accession"),
                           filter=True)
        genomic = self.json_ok(self.get_ok(self.api +
                               "/query?q=accession:FJ497232&fields=accession"),
                               filter=True)
        u = "/query?q=accession.protein:AAH68303&fields=accession"
        explicit_protein = self.json_ok(self.get_ok(self.api + u), filter=True)
        u = "/query?q=accession.rna:BC068303&fields=accession"
        explicit_rna = self.json_ok(self.get_ok(self.api + u), filter=True)
        u = "/query?q=accession.genomic:FJ497232&fields=accession"
        explicit_genomic = self.json_ok(self.get_ok(self.api + u), filter=True)
        eq_(protein["hits"], explicit_protein["hits"])
        eq_(rna["hits"], explicit_rna["hits"])
        eq_(genomic["hits"], explicit_genomic["hits"])
        eq_(protein["hits"], rna["hits"])  # same result whatever the query
        eq_(genomic["hits"], [])  # genomic not indexed
        eq_(rna["total"], 1)
        hit = rna["hits"][0]
        assert "AAH68303" in hit["accession"]["protein"]
        assert "BC068303" in hit["accession"]["rna"]

    def test_query_dotstar_reporter(self):
        human = self.json_ok(self.get_ok(self.api +
                             "/query?q=reporter:2842429&fields=reporter"),
                             filter=True)
        mouse = self.json_ok(self.get_ok(self.api +
                             "/query?q=reporter:1452128_a_at&fields=reporter"),
                             filter=True)
        rat = self.json_ok(self.get_ok(self.api +
                           "/query?q=reporter:1387540_at&fields=reporter"),
                           filter=True)
        # human
        eq_(human["total"], 3)
        eq_(human["hits"][0]["reporter"]["HuGene-1_1"], "8110147")
        assert "2889211" in human["hits"][0]["reporter"]["HuEx-1_0"]
        assert "TC05002114.hg.1" in human["hits"][0]["reporter"]["HTA-2_0"]
        eq_(human["hits"][0]["reporter"]["HG-U133_Plus_2"], "228805_at")
        assert "gnf1h08801_at" in human["hits"][0]["reporter"]["GNF1H"]
        eq_(human["hits"][0]["reporter"]["HuGene-1_1"], "8110147")
        eq_(human["hits"][0]["reporter"]["HuGene-2_1"], "16992761")
        # rat
        eq_(rat["total"], 1)
        eq_(rat["hits"][0]["reporter"]["RaEx-1_0"], "7082865")
        eq_(rat["hits"][0]["reporter"]["Rat230_2"], "1387540_at")
        eq_(rat["hits"][0]["reporter"]["RaGene-2_1"], "17661681")
        eq_(rat["hits"][0]["reporter"]["RaGene-1_1"], "10747640")
        assert "AF036760_at" in rat["hits"][0]["reporter"]["RG-U34A"]
        # rat
        eq_(mouse["total"], 1)
        assert "1456141_x_at" in mouse["hits"][0]["reporter"]["Mouse430_2"]
        eq_(mouse["hits"][0]["reporter"]["MTA-1_0"], "TC0X00000742.mm.1")
        assert "165150_i_at" in mouse["hits"][0]["reporter"]["MG-U74Bv2"]
        eq_(mouse["hits"][0]["reporter"]["MoEx-1_0"], "7012082")
        eq_(mouse["hits"][0]["reporter"]["GNF1M"], "gnf1m11626_at")
        eq_(mouse["hits"][0]["reporter"]["MoGene-2_1"], "17535957")
        eq_(mouse["hits"][0]["reporter"]["MoGene-1_1"], "10600512")

    def test_query_dotstar_interpro(self):
        res = self.json_ok(self.get_ok(self.api +
                           "/query?q=interpro:IPR008389&fields=interpro"))
        eq_(res["total"], 6)
        assert set([pro["id"] for hit in res["hits"]
                    for pro in hit["interpro"]]) == set(['IPR008389',
                                                         'IPR017385'])

    def test_query_dotstar_go(self):
        res = self.json_ok(self.get_ok(self.api +
                           "/query?q=GO:0016324&fields=go&sorted=true"))
        assert res["total"] > 800, \
            "Total is {}, should more than 800".format(res["total"])
        # make sure we're looking at proper
        for (i,h) in enumerate(res["hits"]):
            if h.get("go",{}).get("BP",{}).get("id") == "GO:0008150":
                break
        cc = h["go"]["CC"]
        eq_(len(cc), 4)
        assert set([c["evidence"] for c in cc]), set(["ISO", "IDA"])

    def test_query_dotstar_homologene(self):
        res = self.json_ok(self.get_ok(self.api +
                           "/query?q=homologene:44221&fields=homologene"))
        eq_(res["total"], 3)
        h = res["hits"][0]
        assert set([i[0] for i in h["homologene"]["genes"]]) == \
            set([7955, 8364, 9031, 9598, 9606, 9615, 9913, 10090, 10116])

    def test_query_dotstar_reagent(self):
        res = self.json_ok(self.get_ok(self.api +
                           "/query?q=reagent:GNF190467&fields=reagent"))
        eq_(res["total"], 1)
        hit = res["hits"][0]
        assert {"relationship": "is", "id": "GNF168655"} in \
            hit["reagent"]["GNF_Qia_hs-genome_v1_siRNA"]
        assert {"relationship": "is", "id": "GNF277345"} in \
            hit["reagent"]["GNF_mm+hs-MGC"]
        assert {"relationship": "is", "id": "GNF110093"} in \
            hit["reagent"]["NOVART_hs-genome_siRNA"]

    def test_query_dotstar_uniprot(self):
        swissid = self.json_ok(self.get_ok(self.api +
                               "/query?q=uniprot:Q8NEB7&fields=uniprot"),
                               filter=True)
        trembid = self.json_ok(self.get_ok(self.api +
                               "/query?q=uniprot:F5H2C2&fields=uniprot"),
                               filter=True)
        eq_(swissid["hits"], trembid["hits"])
        eq_(trembid["total"], 1)
        eq_(trembid["hits"][0]["uniprot"]["Swiss-Prot"], "Q8NEB7")
        assert set(trembid["hits"][0]["uniprot"]["TrEMBL"]), \
            set(["E7EP66", "F5H2C2", "F5H3P4", "F5H5S8"])

    def test_query_ensembl(self):
        u = "/query?q=ensemblprotein:ENSP00000379391&fields=ensembl"
        prot = self.json_ok(self.get_ok(self.api + u))
        u = "/query?q=ensembltranscript:ENST00000396082&fields=ensembl"
        rna = self.json_ok(self.get_ok(self.api + u))
        u = "/query?q=ensemblgene:ENSG00000100373&fields=ensembl"
        gene = self.json_ok(self.get_ok(self.api + u))
        # don' compare score, useless
        [d["hits"][0].pop("_score") for d in [prot, rna, gene]]
        eq_(prot["hits"], rna["hits"])
        eq_(rna["hits"], gene["hits"])
        eq_(rna["total"], 1)
        hit = rna["hits"][0]
        eq_(hit["ensembl"]["gene"], "ENSG00000100373")
        assert "ENSP00000216211" in hit["ensembl"]["protein"]
        assert "ENST00000216211" in hit["ensembl"]["transcript"]

    def test_sort_by_fields(self):
        res = self.json_ok(self.get_ok(self.api + "/query?q=MTFMT&sort=entrezgene"))
        hits = res["hits"]
        assert len(hits) == 3
        eq_(hits[0]["entrezgene"],69606)
        eq_(hits[1]["entrezgene"],123263)
        eq_(hits[2]["entrezgene"],315763)


# Self contained test class, used for CI tools such as Travis
# This will start a Tornado server on its own and perform tests
# against this server.
btsettings = BiothingSettings()
# force static path, as if we were in debug mode
index.settings.update({
    "static_path": btsettings.static_path
})


class MyGeneTestTornadoClient(AsyncHTTPTestCase, MyGeneTest):
    __test__ = True

    def __init__(self, methodName='runTest', **kwargs):
        super(AsyncHTTPTestCase, self).__init__(methodName, **kwargs)
        self.h = TornadoRequestHelper(self)

    def get_app(self):
        return index.get_app(index.APP_LIST)
