import sys
import time
from utils.common import timesofar, ask
from utils.mongo import get_target_db, doc_feeder
from utils.es import ESIndexer
from utils import diff
from builder import DataBuilder
import backend


def clean_target_collection():
    bdr = DataBuilder(backend='mongodb')
    bdr.load_build_config('mygene')
    try:
        target_collection = bdr.pick_target_collection(autoselect=False)
    except KeyboardInterrupt:
        print("Aborted.")
        return

    if ask('Delete collection "{}"'.format(target_collection.name)) == 'Y':
        if ask("Double check! Are you sure?") == 'Y':
            target_collection.drop()
            print('Done, collection "{}" was dropped.'.format(target_collection.name))


def build_index(config, use_parallel=True, noconfirm=False):
    bdr = DataBuilder(backend='mongodb')
    bdr.load_build_config(config)
    target_collection = bdr.pick_target_collection()
    target_es_index = 'genedoc_' + bdr._build_config['name']

    if target_collection:
        es_idxer = ESIndexer(mapping=bdr.get_mapping())
        es_idxer.ES_INDEX_NAME = target_es_index
        es_idxer.step = 10000
        es_idxer.use_parallel = use_parallel
        es_server = es_idxer.conn.servers[0].geturl()
        print("ES target: {}/{}/{}".format(es_server,
                                           es_idxer.ES_INDEX_NAME,
                                           es_idxer.ES_INDEX_TYPE))
        if noconfirm or ask("Continue?") == 'Y':
            #es_idxer.s = 609000
            #es_idxer.conn.indices.delete_index(es_idxer.ES_INDEX_NAME)
            es_idxer.create_index()
            es_idxer.delete_index_type(es_idxer.ES_INDEX_TYPE, noconfirm=noconfirm)
            es_idxer.build_index(target_collection, verbose=False)
            es_idxer.optimize()
        else:
            print("Aborted.")
    else:
        print("Error: target collection is not ready yet or failed to build.")


def _pick_one(src_li, prompt="Pick one above: "):
    last_name = ''
    i = 0
    for src, count, name in src_li:
        i += 1
        if last_name == 'mongodb' and name == 'es':
            print
        print('\t{0:<5}{1:<45}\t{2:<9}\t{3}'.format(str(i) + ':', src, count, name))
        last_name = name

    while 1:
        selected_idx = input("\n" + prompt)
        try:
            selected_idx = int(selected_idx) - 1
            src_picked = src_li[selected_idx]
            break
        except (ValueError, IndexError):
            continue
    return src_picked


def diff2src(use_parallel=True, noconfirm=False):
    src_li = []

    target_db = get_target_db()
    src_li.extend([(name, target_db[name].count(), 'mongodb') for name in sorted(target_db.collection_names()) if name.startswith('genedoc')])

    es_idxer = ESIndexer()
    es_idxer.conn.default_indices = []
    for es_idx in es_idxer.conn.indices.get_indices():
        if es_idx.startswith('genedoc'):
            es_idxer.ES_INDEX_NAME = es_idx
            src_li.append((es_idx, es_idxer.count()['count'], 'es'))

    print("Found {} sources:".format(len(src_li)))
    src_1 = _pick_one(src_li, "Pick first source above: ")
    src_li.remove(src_1)
    print
    src_2 = _pick_one(src_li, "Pick second source above: ")

    sync_li = []
    for src in (src_1, src_2):
        if src[2] == 'mongodb':
            b = backend.GeneDocMongoDBBackend(target_db[src[0]])
        elif src[2] == 'es':
            es_idxer = ESIndexer()
            es_idxer.ES_INDEX_NAME = src[0]
            es_idxer.step = 10000
            b = backend.GeneDocESBackend(es_idxer)
        sync_li.append(b)

    sync_src, sync_target = sync_li
    print('\tsync_src:\t{:<45}{}\t{}'.format(*src_1))
    print('\tsync_target\t{:<45}{}\t{}'.format(*src_2))
    if noconfirm or ask("Continue?") == "Y":
        changes = diff.diff_collections(sync_src, sync_target, use_parallel=use_parallel)
        return changes


def test():
    target = get_target_db()
    sync_src = backend.GeneDocMongoDBBackend(target['genedoc_mygene_allspecies_20130402_uiu7bkyi'])
    idxer = ESIndexer()
    sync_target = backend.GeneDocESBackend(idxer)
    return sync_src, sync_target


def update_index(changes, sync_src, sync_target, noconfirm=False):
    # changes['_add'] = changes['delete']
    # changes['_delete'] = changes['add']
    # changes['delete'] = changes['_delete']
    # changes['add'] = changes['_add']
    # del changes['_add']
    # del changes['_delete']

    print("\t{}\trecords will be added.".format(len(changes['add'])))
    print("\t{}\trecords will be deleted.".format(len(changes['delete'])))
    print("\t{}\trecords will be updated.".format(len(changes['update'])))

    print()
    print('\tsync_src:\t{:<45}{}'.format(sync_src.target_collection.name,
                                         sync_src.name))
    print('\tsync_target\t{:<45}{}'.format(sync_target.target_esidxer.ES_INDEX_NAME,
                                           sync_target.name))

    if noconfirm or ask("Continue?") == 'Y':
        t00 = time.time()
        es_idxer = sync_target.target_esidxer

        if len(changes['add']) > 0:
            print("Adding {} new records...".format(len(changes['add'])))
            t0 = time.time()
            _q = {'_id': {'$in': changes['add']}}
            for docs in doc_feeder(sync_src.target_collection, step=1000, inbatch=True, query=_q):
                es_idxer.add_docs(docs)
            print("Done. [{}]".format(timesofar(t0)))

        if len(changes['delete']) > 0:
            print("Deleting {} old records...".format(len(changes['delete'])))
            t0 = time.time()
            es_idxer.delete_docs(changes['delete'])
            print("Done. [{}]".format(timesofar(t0)))

        if len(changes['update']) > 0:
            print("Updating {} existing records...".format(len(changes['update'])))
            t0 = time.time()
            ids = [d['_id'] for d in changes['update']]
            _q = {'_id': {'$in': ids}}
            for docs in doc_feeder(sync_src.target_collection, step=1000, inbatch=True, query=_q):
                es_idxer.add_docs(docs)
            print("Done. [{}]".format(timesofar(t0)))
        print('=' * 20)
        print('Finished. [{}]'.format(timesofar(t00)))


def sync_index(config, use_parallel=True, noconfirm=False):

    bdr = DataBuilder(backend='mongodb')
    bdr.load_build_config(config)
    target_collection = bdr.pick_target_collection()
    target_es_index = 'genedoc_' + bdr._build_config['name']

    sync_src = backend.GeneDocMongoDBBackend(target_collection)

    es_idxer = ESIndexer(bdr.get_mapping())
    es_idxer.ES_INDEX_NAME = target_es_index
    es_idxer.step = 10000
    es_idxer.use_parallel = use_parallel
    sync_target = backend.GeneDocESBackend(es_idxer)

    print('\tsync_src:\t{:<40}{}\t{}'.format(target_collection.name,
                                             sync_src.name,
                                             sync_src.count()))
    print('\tsync_target\t{:<40}{}\t{}'.format(target_es_index,
                                               sync_target.name,
                                               sync_target.count()))
    if noconfirm or ask("Continue?") == "Y":
        changes = diff.diff_collections(sync_src, sync_target)
        return changes


def main():
    if len(sys.argv) > 1:
        config = sys.argv[1]
    else:
        config = 'mygene_allspecies'
    use_parallel = '-p' in sys.argv
    noconfirm = '-b' in sys.argv
    if config == 'clean':
        clean_target_collection()
    else:
        t0 = time.time()
        build_index(config, use_parallel=use_parallel, noconfirm=noconfirm)
        print("Finished.", timesofar(t0))


if __name__ == '__main__':
    main()
