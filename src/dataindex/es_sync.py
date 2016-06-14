from __future__ import print_function
import sys
from datetime import datetime
import time
import re
import os
from pprint import pprint

from config import TARGET_ES_INDEX_SUFFIX, ES_HOST
from utils.es import ESIndexer
from utils.mongo import get_target_db, get_src_build
from utils.common import loadobj  # the one from biothings doesn't work for now
from biothings.utils.common import timesofar, ask, iter_n
from databuild.backend import GeneDocMongoDBBackend, GeneDocESBackend
from databuild.sync import get_changes_stats
from tunnel import open_tunnel, es_local_tunnel_port


class ESIndexer2(ESIndexer):
    step = 5000

    def _split_source_name(self, source):
        pat = '(\w+)_(\d{8})_\w{8}'
        mat = re.match(pat, source)
        prefix, timestamp = mat.groups()
        return prefix, timestamp

    def get_timestamp_stats(self, debug=False):
        q = {
            "aggs": {
                "timestamp": {
                    "terms": {
                        "field": "_timestamp",
                        "size": 1000,
                        "order": {"_term": "desc"}
                    }
                }
            },
            "size": 0
        }
        res = self.conn.search(index=self.ES_INDEX_NAME, doc_type=self.ES_INDEX_TYPE, body=q)
        if debug:
            return res
        return [(datetime.utcfromtimestamp(x['key'] * 1. / 1000), x['doc_count']) for x in res['aggregations']['timestamp']['buckets']]

    def get_latest_timestamp(self):
        ts_stats = self.get_timestamp_stats()
        return ts_stats[0][0]

    def pre_verify_changes(self, changes):
        src = changes['source']
        prefix = self._split_source_name(src)[0]
        _es_target = prefix + TARGET_ES_INDEX_SUFFIX    # '_current_1'
        assert self.ES_INDEX_NAME == _es_target,\
            '"{}" does not match "{}"'.format(self.ES_INDEX_NAME, _es_target)
        print('\033[34;06m{}\033[0m:'.format('[Pending changes]'))
        get_changes_stats(changes)
        print('\033[34;06m{}\033[0m:'.format('[Target ES]'))
        self.check()

    def apply_changes(self, changes, verify=True, noconfirm=False):
        if verify:
            self.pre_verify_changes(changes)

        if not (noconfirm or ask('\nContinue to apply changes?') == 'Y'):
            print("Aborted.")
            return -1
        step = self.step
        _db = get_target_db()
        source_col = _db[changes['source']]
        src = GeneDocMongoDBBackend(source_col)
        target = GeneDocESBackend(self)
        _timestamp = changes['timestamp']

        def _add_docs(ids):
            i = 0
            for _ids in iter_n(ids, step):
                t1 = time.time()
                _doc_li = src.mget_from_ids(_ids)
                for _doc in _doc_li:
                    _doc['_timestamp'] = _timestamp
                    i += 1
                target.insert(_doc_li)
                print('\t{}\t{}'.format(i, timesofar(t1)))

        t0 = time.time()
        if changes['add']:
            print("Adding {} new docs...".format(len(changes['add'])))
            t00 = time.time()
            _add_docs(changes['add'])
            print("done. [{}]".format(timesofar(t00)))
        if changes['delete']:
            print("Deleting {} discontinued docs...".format(len(changes['delete'])), end='')
            t00 = time.time()
            target.remove_from_ids(changes['delete'], step=step)
            print("done. [{}]".format(timesofar(t00)))
        if changes['update']:
            print("Updating {} existing docs...".format(len(changes['update'])))
            t00 = time.time()
            ids = [x['_id'] for x in changes['update']]
            _add_docs(ids)
            print("done. [{}]".format(timesofar(t00)))

        target.finalize()

        print("\n")
        print("Finished.", timesofar(t0))

    def get_mapping_meta(self, changes):
        src = changes['source']
        pat = 'genedoc_(\w+)_(\d{8})_\w{8}'
        mat = re.match(pat, src)
        if mat:
            build_cfg, timestamp = mat.groups()
        src_build = get_src_build()
        _cfg = src_build.find_one({'_id': build_cfg})
        _li = [x for x in _cfg['build'] if x['target'] == src]
        assert len(_li) == 1
        _build = _li[0]
        assert _build['status'] == 'success'
        meta = dict(source=src,
                    src_version=_build['src_version'],
                    stats=_build['stats'],
                    timestamp=_build['timestamp'])
        return meta

    def post_verify_changes(self, changes):
        target = GeneDocESBackend(self)
        _timestamp = changes['timestamp']
        ts_stats = self.get_timestamp_stats()

        if changes['add'] or changes['update']:
            print('Verifying "add" and "update"...', end='')
            assert ts_stats[0][0] == _timestamp, "{} != {}".format(ts_stats[0][0], _timestamp)
            _cnt = ts_stats[0][1]
            _cnt_add_update = len(changes['add']) + len(changes['update'])
            if _cnt == _cnt_add_update:
                print('...{}=={}...OK'.format(_cnt, _cnt_add_update))
            else:
                print('...{}!={}...ERROR!!!'.format(_cnt, _cnt_add_update))
        if changes['delete']:
            print('Verifying "delete"...', end='')
            _res = target.mget_from_ids(changes['delete'])
            _cnt = len([x for x in _res if x])
            if _cnt == 0:
                print('...{}==0...OK'.format(_cnt))
            else:
                print('...{}!=0...ERROR!!!'.format(_cnt))

        print("Verifying all docs have timestamp...", end='')
        _cnt = sum([x[1] for x in ts_stats])
        _cnt_all = self.count()['count']
        if _cnt == _cnt_all:
            print('{}=={}...OK'.format(_cnt, _cnt_all))
        else:
            print('ERROR!!!\n\t Should be "{}", but get "{}"'.format(_cnt_all, _cnt))

        print("Verifying all new docs have updated timestamp...")
        q = {
            'query': {
                'term': {
                    '_timestamp': _timestamp
                }
            }
        }
        cur = self.doc_feeder(query=q, fields=[], step=10000)
        _li1 = sorted(changes['add'] + [x['_id'] for x in changes['update']])
        _li2 = sorted([x['_id'] for x in cur])
        if _li1 == _li2:
            print("{}=={}...OK".format(len(_li1), len(_li2)))
        else:
            print('ERROR!!!\n\t Should be "{}", but get "{}"'.format(len(_li1), len(_li2)))


#esi = utils.es.ESIndexer2('genedoc_mygene' + TARGET_ES_INDEX_SUFFIX, es_host='su02:9500')
#esi.apply_changes(changes)
#meta = esi.get_mapping_meta(changes)
#esi.update_mapping_meta({'_meta': meta})
#esi.post_verify_changes(changes)


def _get_current_changes_fn(config):
    if config == 'genedoc_mygene_allspecies':
        pattern = 'changes_\d{8}_allspecies.pyobj'
    elif config == 'genedoc_mygene':
        pattern = 'changes_\d{8}.pyobj'

    fli = [f for f in os.listdir('.') if re.match(pattern, f)]
    if fli:
        _changes_fn = sorted(fli)[-1]
        return _changes_fn


def main():
    if len(sys.argv) > 1:
        config = sys.argv[1]
    else:
        config = 'mygene_allspecies'
    if not config.startswith('genedoc_'):
        config = 'genedoc_' + config
    assert config in ['genedoc_mygene', 'genedoc_mygene_allspecies']
    noconfirm = '-b' in sys.argv

    _changes_fn = _get_current_changes_fn(config)
    if _changes_fn:
        print("Changes file: " + _changes_fn)
    else:
        print("No changes file found. Aborted.")
        return -1
    if noconfirm or ask("Continue to load?") == 'Y':
        changes = loadobj(_changes_fn)
    else:
        print("Aborted.")
        return -2

    _es_index = config + TARGET_ES_INDEX_SUFFIX
    # ES host will be set depending on whether a tunnel is used or not
    with open_tunnel() as tunnel:
        if tunnel.ok:
            _es_host = 'localhost:' + str(es_local_tunnel_port)
        else:
            _es_host = ES_HOST

        esi = ESIndexer2(_es_index, es_host=_es_host)

        meta = esi.get_mapping_meta(changes)
        print('\033[34;06m{}\033[0m:'.format('[Metadata]'))
        pprint(meta)
        code = esi.apply_changes(changes, noconfirm=noconfirm)
        if code != -1:
            # aborted when code == -1
            _meta = {'_meta': meta}
            print(esi.conn.indices.put_mapping(esi.ES_INDEX_TYPE, _meta, [esi.ES_INDEX_NAME]))
            esi.post_verify_changes(changes)


if __name__ == '__main__':
    main()
