'''
Utils to compare two list of gene documents
'''
from __future__ import print_function
import time
import os.path
from biothings.utils.common import timesofar
from databuild.backend import GeneDocMongoDBBackend, GeneDocESBackend
from utils.mongo import get_target_db
from utils.es import ESIndexer


def diff_doc(doc_1, doc_2, exclude_attrs=['_timestamp']):
    diff_d = {'update': {},
              'delete': [],
              'add': {}}
    for attr in set(doc_1) | set(doc_2):
        if exclude_attrs and attr in exclude_attrs:
            continue
        if attr in doc_1 and attr in doc_2:
            _v1 = doc_1[attr]
            _v2 = doc_2[attr]
            if _v1 != _v2:
                diff_d['update'][attr] = _v2
        elif attr in doc_1 and attr not in doc_2:
            diff_d['delete'].append(attr)
        else:
            diff_d['add'][attr] = doc_2[attr]
    if diff_d['update'] or diff_d['delete'] or diff_d['add']:
        return diff_d


def full_diff_doc(doc_1, doc_2, exclude_attrs=['_timestamp']):
    diff_d = {'update': {},
              'delete': [],
              'add': {}}
    for attr in set(doc_1) | set(doc_2):
        if exclude_attrs and attr in exclude_attrs:
            continue
        if attr in doc_1 and attr in doc_2:
            _v1 = doc_1[attr]
            _v2 = doc_2[attr]
            difffound = False
            if isinstance(_v1, dict) and isinstance(_v2, dict):
                if full_diff_doc(_v1, _v2, exclude_attrs):
                    difffound = True
            elif isinstance(_v1, list) and isinstance(_v2, list):
                # there can be unhashable/unordered dict in these lists
                if len(_v1) != len(_v2):
                    difffound = True
                else:
                    for i in _v1:
                        if i not in _v2:
                            difffound = True
                            break
                    # check the other way
                    if not difffound:
                        for i in _v2:
                            if i not in _v1:
                                difffound = True
                                break
            elif _v1 != _v2:
                difffound = True

            if difffound:
                diff_d['update'][attr] = _v2

        elif attr in doc_1 and attr not in doc_2:
            diff_d['delete'].append(attr)
        else:
            diff_d['add'][attr] = doc_2[attr]
    if diff_d['update'] or diff_d['delete'] or diff_d['add']:
        return diff_d

def two_docs_iterator(b1, b2, id_list, step=10000):
    t0 = time.time()
    n = len(id_list)
    for i in range(0, n, step):
        t1 = time.time()
        print("Processing %d-%d documents..." % (i + 1, min(i + step, n)), end='')
        _ids = id_list[i:i+step]
        iter1 = b1.mget_from_ids(_ids, asiter=True)
        iter2 = b2.mget_from_ids(_ids, asiter=True)
        for doc1, doc2 in zip(iter1, iter2):
            yield doc1, doc2
        print('Done.[%.1f%%,%s]' % (i*100./n, timesofar(t1)))
    print("="*20)
    print('Finished.[total time: %s]' % timesofar(t0))


def _diff_doc_worker(args):
    #b1_target_collection, b2_es_index, ids, _path = args
    _b1, _b2, ids, _path = args
    import sys
    if _path not in sys.path:
        sys.path.append(_path)
    import utils.diff
    import importlib
    importlib.reload(utils.diff)
    from utils.diff import _diff_doc_inner_worker, get_backend

    b1 = get_backend(*_b1)
    b2 = get_backend(*_b2)

    _updates = _diff_doc_inner_worker(b1, b2, ids)
    return _updates


def _diff_doc_inner_worker(b1, b2, ids, fastdiff=False, diff_func=full_diff_doc):
    '''if fastdiff is True, only compare the whole doc,
       do not traverse into each attributes.
    '''
    _updates = []
    for doc1, doc2 in two_docs_iterator(b1, b2, ids):
        assert doc1['_id'] == doc2['_id'], repr((ids, len(ids)))
        if fastdiff:
            if doc1 != doc2:
                _updates.append({'_id': doc1['_id']})
        else:
            _diff = diff_func(doc1, doc2)
            if _diff:
                _diff['_id'] = doc1['_id']
                _updates.append(_diff)
    return _updates


# TODO: move to mongodb backend class
def get_mongodb_uri(backend):
    opt = backend.target_collection.database.client._MongoClient__options.credentials
    username = opt and opt.username or None
    password = opt and opt.password or None
    dbase = opt and opt.source or None
    uri = "mongodb://"
    if username:
        if password:
            uri += "%s:%s@" % (username,password)
        else:
            uri += "%s@" % username
    host,port = backend.target_collection.database.client.address
    uri += "%s:%s" % (host,port)
    uri += "/%s" % (dbase or backend.target_collection.database.name)
    #uri += "/%s" % backend.target_collection.name
    print("uri: %s" % uri)
    return uri


def diff_collections(b1, b2, use_parallel=True, step=10000):
    """
    b1, b2 are one of supported backend class in databuild.backend.
    e.g.,
        b1 = GeneDocMongoDBBackend(c1)
        b2 = GeneDocMongoDBBackend(c2)
    """

    id_s1 = set(b1.get_id_list())
    id_s2 = set(b2.get_id_list())
    print("Size of collection 1:\t", len(id_s1))
    print("Size of collection 2:\t", len(id_s2))

    id_in_1 = id_s1 - id_s2
    id_in_2 = id_s2 - id_s1
    id_common = id_s1 & id_s2
    print("# of docs found only in collection 1:\t", len(id_in_1))
    print("# of docs found only in collection 2:\t", len(id_in_2))
    print("# of docs found in both collections:\t", len(id_common))

    print("Comparing matching docs...")
    _updates = []
    if len(id_common) > 0:
        if not use_parallel:
            _updates = _diff_doc_inner_worker(b1, b2, list(id_common))
        else:
            from utils.parallel import run_jobs_on_ipythoncluster
            _path = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
            id_common = list(id_common)
            #b1_target_collection = b1.target_collection.name
            #b2_es_index = b2.target_esidxer.ES_INDEX_NAME
            _b1 = (get_mongodb_uri(b1), b1.target_collection.database.name, b1.target_name, b1.name)
            _b2 = (get_mongodb_uri(b2), b2.target_collection.database.name, b2.target_name, b2.name)
            #task_li = [(b1_target_collection, b2_es_index, id_common[i: i + step], _path) for i in range(0, len(id_common), step)]
            print("b1 %s" % repr(_b1))
            print("b2 %s" % repr(_b2))
            task_li = [(_b1, _b2, id_common[i: i + step], _path) for i in range(0, len(id_common), step)]
            job_results = run_jobs_on_ipythoncluster(_diff_doc_worker, task_li)
            _updates = []
            if job_results:
                for res in job_results:
                    _updates.extend(res)
            else:
                print("Parallel jobs failed or were interrupted.")
                return None

        print("Done. [{} docs changed]".format(len(_updates)))

    _deletes = []
    if len(id_in_1) > 0:
        _deletes = sorted(id_in_1)

    _adds = []
    if len(id_in_2) > 0:
        _adds = sorted(id_in_2)

    changes = {'update': _updates,
               'delete': _deletes,
               'add': _adds}
    return changes


#def get_backend(target_name, bk_type, **kwargs):
#    '''Return a backend instance for given target_name and backend type.
#        currently support MongoDB and ES backend.
#    '''
#    if bk_type == 'mongodb':
#        target_db = get_target_db()
#        target_col = target_db[target_name]
#        return GeneDocMongoDBBackend(target_col)
#    elif bk_type == 'es':
#        esi = ESIndexer(target_name, **kwargs)
#        return GeneDocESBackend(esi)

def get_backend(uri, db, col, bk_type):
    if bk_type != "mongodb":
        raise NotImplemented("Backend type '%s' not supported" % bk_type)
    from utils.mongo import MongoClient
    colobj = MongoClient(uri)[db][col]
    return GeneDocMongoDBBackend(colobj)
