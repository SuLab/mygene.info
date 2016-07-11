from __future__ import print_function
import sys
import os.path
import time
import copy
from datetime import datetime
from pprint import pformat

from utils.mongo import (get_src_db, get_target_db, get_src_master,
                         get_src_build, get_src_dump, doc_feeder)
from biothings.utils.common import (timesofar, ask,
                                    dump2gridfs, get_timestamp, get_random_string)
from utils.common import safewfile, setup_logfile, loadobj
from utils.dataload import list2dict, alwayslist
from utils.es import ESIndexer
import databuild.backend
from config import LOG_FOLDER, logger as logging

'''
#Build_Config example

Build_Config = {
    "name":     "test",          #target_collection will be called "genedoc_test"
    "sources" : ['entrez_gene', 'reporter'],
    "gene_root": ['entrez_gene', 'ensembl_gene']     #either entrez_gene or ensembl_gene or both
}

#for genedoc at mygene.info
Build_Config = {
    "name":     "mygene",          #target_collection will be called "genedoc_mygene"
    "sources":  [u'ensembl_acc',
                 u'ensembl_gene',
                 u'ensembl_genomic_pos',
                 u'ensembl_interpro',
                 u'ensembl_prosite',
                 u'entrez_accession',
                 u'entrez_ec',
                 u'entrez_gene',
                 u'entrez_genesummary',
                 u'entrez_go',
                 u'entrez_homologene',
                 u'entrez_refseq',
                 u'entrez_retired',
                 u'entrez_unigene',
                 u'pharmgkb',
                 u'reagent',
                 u'reporter',
                 u'uniprot',
                 u'uniprot_ipi',
                 u'uniprot_pdb',
                 u'uniprot_pir'],
    "gene_root": ['entrez_gene', 'ensembl_gene']
}
'''


class DataBuilder():

    def __init__(self, build_config=None, backend='mongodb'):
        self.src = get_src_db()
        self.step = 10000
        self.use_parallel = False
        self.merge_logging = True     # save output into a logging file when merge is called.
        self.max_build_status = 10    # max no. of records kept in "build" field of src_build collection.

        self.using_ipython_cluster = False
        self.shutdown_ipengines_after_done = False
        self.log_folder = LOG_FOLDER

        self._build_config = build_config
        self._entrez_geneid_d = None
        self._idmapping_d_cache = {}

        self.get_src_master()

        if backend == 'mongodb':
            self.target = databuild.backend.GeneDocMongoDBBackend()
        elif backend == 'es':
            self.target = databuild.backend.GeneDocESBackend(ESIndexer())
        elif backend == 'couchdb':
            from config import COUCHDB_URL
            import couchdb
            self.target = databuild.backend.GeneDocCouchDBBackend(couchdb.Server(COUCHDB_URL))
        elif backend == 'memory':
            self.target = databuild.backend.GeneDocMemeoryBackend()
        else:
            raise ValueError('Invalid backend "%s".' % backend)

    def make_build_config_for_all(self):
        _cfg = {"sources": list(self.src_master.keys()),
                "gene_root": ['entrez_gene', 'ensembl_gene']}
        self._build_config = _cfg
        return _cfg

    def load_build_config(self, build):
        '''Load build config from src_build collection.'''
        src_build = get_src_build()
        self.src_build = src_build
        _cfg = src_build.find_one({'_id': build})
        if _cfg:
            self._build_config = _cfg
        else:
            raise ValueError('Cannot find build config named "%s"' % build)
        return _cfg

    def log_src_build(self, dict):
        '''put logging dictionary into the corresponding doc in src_build collection.
           if build_config is not loaded from src_build, nothing will be logged.
        '''
        src_build = getattr(self, 'src_build', None)
        if src_build:
            _cfg = src_build.find_one({'_id': self._build_config['_id']})
            _cfg['build'][-1].update(dict)
            src_build.update({'_id': self._build_config['_id']}, {"$set": {'build': _cfg['build']}})

    def log_building_start(self):
        if self.merge_logging:
            #setup logging
            logfile = 'databuild_{}_{}.log'.format('genedoc' + '_' + self._build_config['name'],
                                                   time.strftime('%Y%m%d'))
            logfile = os.path.join(self.log_folder, logfile)
            setup_logfile(logfile)

        src_build = getattr(self, 'src_build', None)
        if src_build:
            #src_build.update({'_id': self._build_config['_id']}, {"$unset": {"build": ""}})
            d = {'status': 'building',
                 'started_at': datetime.now(),
                 'logfile': logfile,
                 'target_backend': self.target.name}
            if self.target.name == 'mongodb':
                d['target'] = self.target.target_collection.name
            elif self.target.name == 'es':
                d['target'] = self.target.target_esidxer.ES_INDEX_NAME
            logging.info(pformat(d))
            src_build.update({'_id': self._build_config['_id']}, {"$push": {'build': d}})
            _cfg = src_build.find_one({'_id': self._build_config['_id']})
            if len(_cfg['build']) > self.max_build_status:
                #remove the first build status record
                src_build.update({'_id': self._build_config['_id']}, {"$pop": {'build': -1}})

    def _get_target_name(self):
        return 'genedoc_{}_{}_{}'.format(self._build_config['name'],
                                         get_timestamp(), get_random_string()).lower()

    def prepare_target(self, target_name=None):
        '''call self.update_backend() after validating self._build_config.'''
        if self.target.name == 'mongodb':
            _db = get_target_db()
            target_collection_name = target_name or self._get_target_name()
            self.target.target_collection = _db[target_collection_name]
            logging.info("Target: %s" % repr(target_collection_name))
        elif self.target.name == 'es':
            self.target.target_esidxer.ES_INDEX_NAME = target_name or self._get_target_name()
            self.target.target_esidxer._mapping = self.get_mapping()
        elif self.target.name == 'couchdb':
            self.target.db_name = target_name or ('genedoc' + '_' + self._build_config['name'])
        elif self.target.name == 'memory':
            self.target.target_name = target_name or ('genedoc' + '_' + self._build_config['name'])

    def get_src_master(self):
        src_master = get_src_master(self.src.client)
        self.src_master = dict([(src['_id'], src) for src in list(src_master.find())])

    def validate_src_collections(self,collection_list=None):
        if not collection_list:
            collection_list = set(self.src.collection_names())
            self.get_src_master()
            build_conf_src = self._build_config['sources']
        else:
            build_conf_src = collection_list

        logging.info("Sources: %s" % repr(build_conf_src))
        if self._build_config:
            for src in build_conf_src:
                assert src in self.src_master, '"%s" not found in "src_master"' % src
                assert src in collection_list, '"%s" not an existing collection in "%s"' % (src, self.src.name)
        else:
            raise ValueError('"build_config" cannot be empty.')

    def _load_entrez_geneid_d(self):
        self._entrez_geneid_d = loadobj(("entrez_gene__geneid_d.pyobj", self.src), mode='gridfs')

    def _load_ensembl2entrez_li(self):
        ensembl2entrez_li = loadobj(("ensembl_gene__2entrezgene_list.pyobj", self.src), mode='gridfs')
        #filter out those deprecated entrez gene ids
        logging.info(len(ensembl2entrez_li))
        ensembl2entrez_li = [(ensembl_id, self._entrez_geneid_d[int(entrez_id)]) for (ensembl_id, entrez_id) in ensembl2entrez_li
                             if int(entrez_id) in self._entrez_geneid_d]
        logging.info(len(ensembl2entrez_li))
        ensembl2entrez = list2dict(ensembl2entrez_li, 0)
        self._idmapping_d_cache['ensembl_gene'] = ensembl2entrez

    def _save_idmapping_gridfs(self):
        '''saving _idmapping_d_cache into gridfs.'''
        idmapping_gridfs_d = {}
        if self._idmapping_d_cache:
            for id_type in self._idmapping_d_cache:
                filename = 'tmp_idmapping_d_cache_' + id_type
                dump2gridfs(self._idmapping_d_cache[id_type], filename, self.src)
                idmapping_gridfs_d[id_type] = filename
        return idmapping_gridfs_d

    def make_genedoc_root(self):
        if not self._entrez_geneid_d:
            self._load_entrez_geneid_d()

        if 'ensembl_gene' in self._build_config['gene_root']:
            self._load_ensembl2entrez_li()
            ensembl2entrez = self._idmapping_d_cache['ensembl_gene']

        if "species" in self._build_config:
            _query = {'taxid': {'$in': self._build_config['species']}}
        elif "species_to_exclude" in self._build_config:
            _query = {'taxid': {'$nin': self._build_config['species_to_exclude']}}
        else:
            _query = None

        geneid_set = []
        species_set = set()
        if "entrez_gene" in self._build_config['gene_root']:
            for doc_li in doc_feeder(self.src['entrez_gene'], inbatch=True, step=self.step, query=_query):
                #target_collection.insert(doc_li, manipulate=False, check_keys=False)
                self.target.insert(doc_li)
                geneid_set.extend([doc['_id'] for doc in doc_li])
                species_set |= set([doc['taxid'] for doc in doc_li])
            cnt_total_entrez_genes = len(geneid_set)
            cnt_total_species = len(species_set)
            logging.info('# of entrez Gene IDs in total: %d' % cnt_total_entrez_genes)
            logging.info('# of species in total: %d' % cnt_total_species)

        if "ensembl_gene" in self._build_config['gene_root']:
            cnt_ensembl_only_genes = 0
            cnt_total_ensembl_genes = 0
            for doc_li in doc_feeder(self.src['ensembl_gene'], inbatch=True, step=self.step, query=_query):
                _doc_li = []
                for _doc in doc_li:
                    cnt_total_ensembl_genes += 1
                    ensembl_id = _doc['_id']
                    entrez_gene = ensembl2entrez.get(ensembl_id, None)
                    if entrez_gene is None:
                        #this is an Ensembl only gene
                        _doc_li.append(_doc)
                        cnt_ensembl_only_genes += 1
                        geneid_set.append(_doc['_id'])
                if _doc_li:
                    #target_collection.insert(_doc_li, manipulate=False, check_keys=False)
                    self.target.insert(_doc_li)
            cnt_matching_ensembl_genes = cnt_total_ensembl_genes - cnt_ensembl_only_genes
            logging.info('# of ensembl Gene IDs in total: %d' % cnt_total_ensembl_genes)
            logging.info('# of ensembl Gene IDs match entrez Gene IDs: %d' % cnt_matching_ensembl_genes)
            logging.info('# of ensembl Gene IDs DO NOT match entrez Gene IDs: %d' % cnt_ensembl_only_genes)

            geneid_set = set(geneid_set)
            logging.info('# of total Root Gene IDs: %d' % len(geneid_set))
            _stats = {'total_entrez_genes': cnt_total_entrez_genes,
                      'total_species': cnt_total_species,
                      'total_ensembl_genes': cnt_total_ensembl_genes,
                      'total_ensembl_genes_mapped_to_entrez': cnt_matching_ensembl_genes,
                      'total_ensembl_only_genes': cnt_ensembl_only_genes,
                      'total_genes': len(geneid_set)}
            self._stats = _stats
            self._src_version = self.get_src_version()
            self.log_src_build({'stats': _stats, 'src_version': self._src_version})
            return geneid_set

    def get_idmapping_d(self, src):
        if src in self._idmapping_d_cache:
            return self._idmapping_d_cache[src]
        else:
            self._load_ensembl2entrez_li()
            return self._idmapping_d_cache[src]
            #raise ValueError('cannot load "idmapping_d" for "%s"' % src)

    def merge(self, step=100000, restart_at=0,sources=None,target=None):
        t0 = time.time()
        self.validate_src_collections(sources)
        self.prepare_target(target_name=target)
        self.log_building_start()
        try:
            if self.using_ipython_cluster:
                if sources:
                    raise NotImplemented("merge speficic sources not supported when using parallel")
                self._merge_ipython_cluster(step=step)
            else:
                self._merge_local(step=step, restart_at=restart_at,src_collection_list=sources)

            if self.target.name == 'es':
                logging.info("Updating metadata...")
                self.update_mapping_meta()

            t1 = round(time.time() - t0, 0)
            t = timesofar(t0)
            self.log_src_build({'status': 'success',
                                'time': t,
                                'time_in_s': t1,
                                'timestamp': datetime.now()})

        finally:
            #do a simple validation here
            if getattr(self, '_stats', None):
                logging.info("Validating...")
                target_cnt = self.target.count()
                if target_cnt == self._stats['total_genes']:
                    logging.info("OK [total count={}]".format(target_cnt))
                else:
                    logging.info("Warning: total count of gene documents does not match [{}, should be {}]".format(target_cnt, self._stats['total_genes']))

            if self.merge_logging:
                sys.stdout.close()

    def merge_resume(self, build_config, at_collection, step=10000):
        '''resume a merging process after a failure.
             .merge_resume('mygene_allspecies', 'reporter')
        '''
        assert not self.using_ipython_cluster, "Abort. Can only resume merging in non-parallel mode."
        self.load_build_config(build_config)
        last_build = self._build_config['build'][-1]
        logging.info("Last build record:")
        logging.info(pformat(last_build))
        assert last_build['status'] == 'building', \
            "Abort. Last build does not need to be resumed."
        assert at_collection in self._build_config['sources'], \
            'Abort. Cannot resume merging from a unknown collection "{}"'.format(at_collection)
        assert last_build['target_backend'] == self.target.name, \
            'Abort. Re-initialized DataBuilder class using matching backend "{}"'.format(last_build['backend'])
        assert last_build.get('stats', None), \
            'Abort. Intital build stats are not available. You should restart the build from the scratch.'
        self._stats = last_build['stats']

        if ask('Continue to resume merging from "{}"?'.format(at_collection)) == 'Y':
            #TODO: resume logging
            target_name = last_build['target']
            self.validate_src_collections()
            self.prepare_target(target_name=target_name)
            src_cnt = 0
            for collection in self._build_config['sources']:
                if collection in ['entrez_gene', 'ensembl_gene']:
                    continue
                src_cnt += 1
                if collection == at_collection:
                    break
            self._merge_local(step=step, restart_at=src_cnt)
            if self.target.name == 'es':
                logging.info("Updating metadata...")
                self.update_mapping_meta()
            self.log_src_build({'status': 'success',
                                'timestamp': datetime.now()})

    def _merge_ipython_cluster(self, step=100000):
        '''Do the merging on ipython cluster.'''
        from ipyparallel import Client, require
        from config import CLUSTER_CLIENT_JSON

        t0 = time.time()
        src_collection_list = [collection for collection in self._build_config['sources']
                               if collection not in ['entrez_gene', 'ensembl_gene']]

        self.target.drop()
        self.target.prepare()
        geneid_set = self.make_genedoc_root()

        idmapping_gridfs_d = self._save_idmapping_gridfs()

        logging.info(timesofar(t0))

        rc = Client(CLUSTER_CLIENT_JSON)
        lview = rc.load_balanced_view()
        logging.info("\t# nodes in use: {}".format(len(lview.targets or rc.ids)))
        lview.block = False
        kwargs = {}
        target_collection = self.target.target_collection
        kwargs['server'], kwargs['port'] = target_collection.database.client.address
        kwargs['src_db'] = self.src.name
        kwargs['target_db'] = target_collection.database.name
        kwargs['target_collection_name'] = target_collection.name
        kwargs['limit'] = step

        @require('pymongo', 'time', 'types')
        def worker(kwargs):
            server = kwargs['server']
            port = kwargs['port']
            src_db = kwargs['src_db']
            target_db = kwargs['target_db']
            target_collection_name = kwargs['target_collection_name']

            src_collection = kwargs['src_collection']
            skip = kwargs['skip']
            limit = kwargs['limit']

            def load_from_gridfs(filename, db):
                import gzip
                import pickle
                import gridfs
                fs = gridfs.GridFS(db)
                fobj = fs.get(filename)
                gzfobj = gzip.GzipFile(fileobj=fobj)
                try:
                    object = pickle.load(gzfobj)
                finally:
                    gzfobj.close()
                    fobj.close()
                return object

            def alwayslist(value):
                if value is None:
                    return []
                if isinstance(value, (list, tuple)):
                    return value
                else:
                    return [value]

            conn = pymongo.MongoClient(server, port)
            src = conn[src_db]
            target_collection = conn[target_db][target_collection_name]

            idmapping_gridfs_name = kwargs.get('idmapping_gridfs_name', None)
            if idmapping_gridfs_name:
                idmapping_d = load_from_gridfs(idmapping_gridfs_name, src)
            else:
                idmapping_d = None

            cur = src[src_collection].find(skip=skip, limit=limit, timeout=False)
            cur.batch_size(1000)
            try:
                for doc in cur:
                    _id = doc['_id']
                    if idmapping_d:
                        _id = idmapping_d.get(_id, None) or _id
                    # there could be cases that idmapping returns multiple entrez_gene id.
                    for __id in alwayslist(_id): 
                        __id = str(__id)
                        doc.pop('_id', None)
                        doc.pop('taxid', None)
                        target_collection.update({'_id': __id}, doc, manipulate=False, upsert=False)
                        #target_collection.update({'_id': __id}, {'$set': doc},
            finally:
                cur.close()

        t0 = time.time()
        task_list = []
        for src_collection in src_collection_list:
            _kwargs = copy.copy(kwargs)
            _kwargs['src_collection'] = src_collection
            id_type = self.src_master[src_collection].get('id_type', None)
            if id_type:
                idmapping_gridfs_name = idmapping_gridfs_d[id_type]
                _kwargs['idmapping_gridfs_name'] = idmapping_gridfs_name
            cnt = self.src[src_collection].count()
            for s in range(0, cnt, step):
                __kwargs = copy.copy(_kwargs)
                __kwargs['skip'] = s
                task_list.append(__kwargs)

        logging.info("\t# of tasks: {}".format(len(task_list)))
        logging.info("\tsubmitting...")
        job = lview.map_async(worker, task_list)
        logging.info("done.")
        job.wait_interactive()
        logging.info("\t# of results returned: {}".format(len(job.result())))
        logging.info("\ttotal time: {}".format(timesofar(t0)))

        if self.shutdown_ipengines_after_done:
            logging.info("\tshuting down all ipengine nodes...")
            lview.shutdown()
            logging.info('Done.')

    def _merge_local(self, step=100000, restart_at=0, src_collection_list=None):
        if restart_at == 0 and src_collection_list is None:
            self.target.drop()
            self.target.prepare()
            geneid_set = self.make_genedoc_root()
        else:
            if not self._entrez_geneid_d:
                self._load_entrez_geneid_d()
            #geneid_set = set([x['_id'] for x in target_collection.find(projection=[], manipulate=False)])
            geneid_set = set(self.target.get_id_list())
            logging.info('\t', len(geneid_set))

        if not src_collection_list:
            src_collection_list = self._build_config['sources']
        src_cnt = 0
        for collection in src_collection_list:
            if collection in ['entrez_gene', 'ensembl_gene']:
                continue

            src_cnt += 1

            id_type = self.src_master[collection].get('id_type', None)
            flag_need_id_conversion = id_type is not None
            if flag_need_id_conversion:
                idmapping_d = self.get_idmapping_d(id_type)
            else:
                idmapping_d = None

            if restart_at <= src_cnt:
                if self.use_parallel:
                    self.doc_queue = []
                    self._merge_parallel_ipython(collection, geneid_set,
                                                 step=step, idmapping_d=idmapping_d)
                else:
                    self._merge_sequential(collection, geneid_set,
                                           step=step, idmapping_d=idmapping_d)
        self.target.finalize()

    def _merge_sequential(self, collection, geneid_set, step=100000, idmapping_d=None):
        for doc in doc_feeder(self.src[collection], step=step):
            _id = doc['_id']
            if idmapping_d:
                _id = idmapping_d.get(_id, None) or _id
            for __id in alwayslist(_id):    # there could be cases that idmapping returns multiple entrez_gene ids.
                __id = str(__id)
                if __id in geneid_set:
                    doc.pop('_id', None)
                    doc.pop('taxid', None)
                    # target_collection.update({'_id': __id}, {'$set': doc},
                    #                           manipulate=False,
                    #                           upsert=False) #,safe=True)
                    self.target.update(__id, doc)

    def _merge_parallel(self, collection, geneid_set, step=100000, idmapping_d=None):
        from multiprocessing import Process, Queue
        NUMBER_OF_PROCESSES = 8

        input_queue = Queue()
        input_queue.conn_pool = []

        def worker(q, target):
            while True:
                doc = q.get()
                if doc == 'STOP':
                    break
                __id = doc.pop('_id')
                doc.pop('taxid', None)
                target.update(__id, doc)
                # target_collection.update({'_id': __id}, {'$set': doc},
                #                           manipulate=False,
                #                           upsert=False) #,safe=True)

        # Start worker processes
        for i in range(NUMBER_OF_PROCESSES):
            Process(target=worker, args=(input_queue, self.target)).start()

        for doc in doc_feeder(self.src[collection], step=step):
            _id = doc['_id']
            if idmapping_d:
                _id = idmapping_d.get(_id, None) or _id
            for __id in alwayslist(_id):    # there could be cases that idmapping returns multiple entrez_gene ids.
                __id = str(__id)
                if __id in geneid_set:
                    doc['_id'] = __id
                    input_queue.put(doc)

        # Tell child processes to stop
        for i in range(NUMBER_OF_PROCESSES):
            input_queue.put('STOP')

    def _merge_parallel_ipython(self, collection, geneid_set, step=100000, idmapping_d=None):
        from IPython.parallel import Client, require

        rc = Client()
        dview = rc[:]
        #dview = rc.load_balanced_view()
        dview.block = False
        target_collection = self.target.target_collection
        dview['server'], dview['port'] = target_collection.database.client.address
        dview['database'] = target_collection.database.name
        dview['collection_name'] = target_collection.name

        def partition(lst, n):
            q, r = divmod(len(lst), n)
            indices = [q * i + min(i, r) for i in range(n + 1)]
            return [lst[indices[i]:indices[i + 1]] for i in range(n)]

        @require('pymongo', 'time')
        def worker(doc_li):
            conn = pymongo.MongoClient(server, port)
            target_collection = conn[database][collection_name]
            t0 = time.time()
            for doc in doc_li:
                __id = doc.pop('_id')
                doc.pop('taxid', None)
                target_collection.update({'_id': __id}, {'$set': doc},
                                         manipulate=False,
                                         upsert=False)  # ,safe=True)
            logging.info('Done. [%.1fs]' % (time.time() - t0))

        for doc in doc_feeder(self.src[collection], step=step):
            _id = doc['_id']
            if idmapping_d:
                _id = idmapping_d.get(_id, None) or _id
            for __id in alwayslist(_id):    # there could be cases that idmapping returns multiple entrez_gene ids.
                __id = str(__id)
                if __id in geneid_set:
                    doc['_id'] = __id
                    self.doc_queue.append(doc)

                    if len(self.doc_queue) >= step:
                        #dview.scatter('doc_li', self.doc_queue)
                        #dview.apply_async(worker)
                        dview.map_async(worker, partition(self.doc_queue, len(rc.ids)))
                        self.doc_queue = []
                        logging.info("!")

    def get_src_version(self):
        src_dump = get_src_dump(self.src.client)
        src_version = {}
        for src in src_dump.find():
            version = src.get('release', src.get('timestamp', None))
            if version:
                src_version[src['_id']] = version
        return src_version

    def get_last_src_build_stats(self):
        src_build = getattr(self, 'src_build', None)
        if src_build:
            _cfg = src_build.find_one({'_id': self._build_config['_id']})
            if _cfg['build'][-1].get('status', None) == 'success' and \
               _cfg['build'][-1].get('stats', None):
                stats = _cfg['build'][-1]['stats']
                return stats

    def get_target_collection(self):
        '''get the lastest target_collection from src_build record.'''
        src_build = getattr(self, 'src_build', None)
        if src_build:
            _cfg = src_build.find_one({'_id': self._build_config['_id']})
            if _cfg['build'][-1].get('status', None) == 'success' and \
               _cfg['build'][-1].get('target', None):
                target_collection = _cfg['build'][-1]['target']
                _db = get_target_db()
                target_collection = _db[target_collection]
                return target_collection

    def pick_target_collection(self, autoselect=True):
        '''print out a list of available target_collection, let user to pick one.'''
        target_db = get_target_db()
        target_collection_prefix = 'genedoc_' + self._build_config['name']
        target_collection_list = [target_db[name] for name in sorted(target_db.collection_names()) if name.startswith(target_collection_prefix)]
        if target_collection_list:
            logging.info("Found {} target collections:".format(len(target_collection_list)))
            logging.info('\n'.join(['\t{0:<5}{1.name:<45}\t{2}'.format(
                str(i + 1) + ':', target, target.count()) for (i, target) in enumerate(target_collection_list)]))
            logging.info()
            while 1:
                if autoselect:
                    selected_idx = input("Pick one above [{}]:".format(len(target_collection_list)))
                else:
                    selected_idx = input("Pick one above:")
                if autoselect:
                    selected_idx = selected_idx or len(target_collection_list)
                try:
                    selected_idx = int(selected_idx)
                    break
                except ValueError:
                    continue
            return target_collection_list[selected_idx - 1]
        else:
            logging.info("Found no target collections.")

    def get_mapping(self, enable_timestamp=True):
        '''collect mapping data from data sources.
           This is for GeneDocESBackend only.
        '''
        mapping = {}
        src_master = get_src_master(self.src.client)
        for collection in self._build_config['sources']:
            meta = src_master.find_one({"_id" : collection})
            if 'mapping' in meta:
                mapping.update(meta['mapping'])
            else:
                logging.info('Warning: "%s" collection has no mapping data.' % collection)
        mapping = {"properties": mapping,
                   "dynamic": False}
        if enable_timestamp:
            mapping['_timestamp'] = {
                "enabled": True,
            }
        #allow source Compression
        #Note: no need of source compression due to "Store Level Compression"
        #mapping['_source'] = {'compress': True,}
        #                      'compress_threshold': '1kb'}
        return mapping

    def update_mapping_meta(self):
        '''updating _meta field of ES mapping data, including index stats, versions.
           This is for GeneDocESBackend only.
        '''
        _meta = {}
        src_version = self.get_src_version()
        if src_version:
            _meta['src_version'] = src_version
        if getattr(self, '_stats', None):
            _meta['stats'] = self._stats

        if _meta:
            self.target.target_esidxer.update_mapping_meta({'_meta': _meta})

    def validate(self, build_config='mygene_allspecies', n=10):
        '''Validate merged genedoc, currently for ES backend only.'''
        import random
        import itertools
        import pyes

        self.load_build_config(build_config)
        last_build = self._build_config['build'][-1]
        logging.info("Last build record:")
        logging.info(pformat(last_build))
        #assert last_build['target_backend'] == 'es', '"validate" currently works for "es" backend only'

        target_name = last_build['target']
        self.validate_src_collections()
        self.prepare_target(target_name=target_name)
        logging.info("Validating...")
        target_cnt = self.target.count()
        stats_cnt = last_build['stats']['total_genes']
        if target_cnt == stats_cnt:
            logging.info("OK [total count={}]".format(target_cnt))
        else:
            logging.info("Warning: total count of gene documents does not match [{}, should be {}]".format(target_cnt, stats_cnt))

        if n > 0:
            for src in self._build_config['sources']:
                logging.info("\nSrc:", src)
                # if 'id_type' in self.src_master[src] and self.src_master[src]['id_type'] != 'entrez_gene':
                #     print "skipped."
                #     continue
                cnt = self.src[src].count()
                fdr1 = doc_feeder(self.src[src], step=10000, s=cnt - n)
                rand_s = random.randint(0, cnt - n)
                fdr2 = doc_feeder(self.src[src], step=n, s=rand_s, e=rand_s + n)
                _first_exception = True
                for doc in itertools.chain(fdr1, fdr2):
                    _id = doc['_id']
                    try:
                        es_doc = self.target.get_from_id(_id)
                    except pyes.exceptions.NotFoundException:
                        if _first_exception:
                            logging.info()
                            _first_exception = False
                        logging.info(_id, 'not found.')
                        continue
                    for k in doc:
                        if src == 'entrez_homologene' and k == 'taxid':
                            # there is occasionally known error for taxid in homologene data.
                            continue
                        assert es_doc.get(k, None) == doc[k], (_id, k, es_doc.get(k, None), doc[k])

    def build_index(self, use_parallel=True):
        target_collection = self.get_target_collection()
        if target_collection:
            es_idxer = ESIndexer(mapping=self.get_mapping())
            es_idxer.ES_INDEX_NAME = 'genedoc_' + self._build_config['name']
            es_idxer.step = 10000
            es_idxer.use_parallel = use_parallel
            #es_idxer.s = 609000
            #es_idxer.conn.indices.delete_index(es_idxer.ES_INDEX_NAME)
            es_idxer.create_index()
            es_idxer.delete_index_type(es_idxer.ES_INDEX_TYPE, noconfirm=True)
            es_idxer.build_index(target_collection, verbose=False)
            es_idxer.optimize()
        else:
            logging.info("Error: target collection is not ready yet or failed to build.")

    def build_index2(self, build_config='mygene_allspecies', last_build_idx=-1, use_parallel=False, es_host=None, es_index_name=None, noconfirm=False):
        """Build ES index from last successfully-merged mongodb collection.
            optional "es_host" argument can be used to specified another ES host, otherwise default ES_HOST.
            optional "es_index_name" argument can be used to pass an alternative index name, otherwise same as mongodb collection name
        """
        self.load_build_config(build_config)
        assert "build" in self._build_config, "Abort. No such build records for config %s" % build_config
        last_build = self._build_config['build'][last_build_idx]
        logging.info("Last build record:")
        logging.info(pformat(last_build))
        assert last_build['status'] == 'success', \
            "Abort. Last build did not success."
        assert last_build['target_backend'] == "mongodb", \
            'Abort. Last build need to be built using "mongodb" backend.'
        assert last_build.get('stats', None), \
            'Abort. Last build stats are not available.'
        self._stats = last_build['stats']
        assert last_build.get('target', None), \
            'Abort. Last build target_collection is not available.'

        # Get the source collection to build the ES index
        # IMPORTANT: the collection in last_build['target'] does not contain _timestamp field,
        #            only the "genedoc_*_current" collection does. When "timestamp" is enabled
        #            in mappings, last_build['target'] collection won't be indexed by ES correctly,
        #            therefore, we use "genedoc_*_current" collection as the source here:
        #target_collection = last_build['target']
        target_collection = "genedoc_{}_current".format(build_config)
        _db = get_target_db()
        target_collection = _db[target_collection]
        logging.info()
        logging.info('Source: ', target_collection.name)
        _mapping = self.get_mapping()
        _meta = {}
        src_version = self.get_src_version()
        if src_version:
            _meta['src_version'] = src_version
        if getattr(self, '_stats', None):
            _meta['stats'] = self._stats
        if 'timestamp' in last_build:
            _meta['timestamp'] = last_build['timestamp']
        if _meta:
            _mapping['_meta'] = _meta
        es_index_name = es_index_name or target_collection.name
        es_idxer = ESIndexer(mapping=_mapping,
                             es_index_name=es_index_name,
                             es_host=es_host,
                             step=5000)
        if build_config == 'mygene_allspecies':
            es_idxer.number_of_shards = 10   # default 5
        es_idxer.check()
        if noconfirm or ask("Continue to build ES index?") == 'Y':
            es_idxer.use_parallel = use_parallel
            #es_idxer.s = 609000
            if es_idxer.exists_index(es_idxer.ES_INDEX_NAME):
                if noconfirm or ask('Index "{}" exists. Delete?'.format(es_idxer.ES_INDEX_NAME)) == 'Y':
                    es_idxer.conn.indices.delete(es_idxer.ES_INDEX_NAME)
                else:
                    logging.info("Abort.")
                    return
            es_idxer.create_index()
            #es_idxer.delete_index_type(es_idxer.ES_INDEX_TYPE, noconfirm=True)
            es_idxer.build_index(target_collection, verbose=False)
            # time.sleep(10)    # pausing 10 second here
            # if es_idxer.wait_till_all_shards_ready():
            #     print "Optimizing...", es_idxer.optimize()

    def sync_index(self, use_parallel=True):
        from utils import diff

        sync_src = self.get_target_collection()

        es_idxer = ESIndexer(self.get_mapping())
        es_idxer.ES_INDEX_NAME = sync_src.target_collection.name
        es_idxer.step = 10000
        es_idxer.use_parallel = use_parallel
        sync_target = databuild.backend.GeneDocESBackend(es_idxer)

        changes = diff.diff_collections(sync_src, sync_target)
        return changes


def main():
    if len(sys.argv) > 1:
        config = sys.argv[1]
    else:
        config = 'mygene_allspecies'
    use_parallel = '-p' in sys.argv
    sources = None  # will build all sources
    target = None   # will generate a new collection name
    # "target_col:src_col1,src_col2" will specifically merge src_col1
    # and src_col2 into existing target_col (instead of merging everything)
    if not use_parallel and len(sys.argv) > 2:
        target,tmp = sys.argv[2].split(":")
        sources = tmp.split(",")

    t0 = time.time()
    bdr = DataBuilder(backend='mongodb')
    bdr.load_build_config(config)
    bdr.using_ipython_cluster = use_parallel
    bdr.merge(sources=sources,target=target)

    logging.info("Finished.", timesofar(t0))


if __name__ == '__main__':
    main()
