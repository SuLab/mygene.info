from __future__ import print_function
import sys
import os
import os.path
import re
from datetime import datetime
import time

from utils.mongo import get_target_db, doc_feeder
from .backend import GeneDocMongoDBBackend
from utils.diff import diff_collections
from utils.common import (iter_n, setup_logfile,
                          dump, send_s3_file, safewfile)
from biothings.utils.common import timesofar, ask, is_str
from config import LOG_FOLDER, logger as logging


class GeneDocSyncer:
    def __init__(self, build_config='genedoc_mygene'):
        self.build_config = build_config
        self._db = get_target_db()
        self._target_col = self._db[self.build_config + '_current']
        self.step = 10000

    def get_source_list(self):
        '''return a list of available source collections.'''
        pat = self.build_config + '_(\d{8})_\w{8}'
        _li = []
        for coll_name in self._db.collection_names():
            mat = re.match(pat, coll_name)
            if mat:
                _li.append(coll_name)
        return sorted(_li)

    def get_latest_source_col(self, n=-1):
        return self._db[self.get_source_list()[n]]

    def get_new_source_list(self):
        '''return a list of source collections have not be applied.'''
        _li = []
        target_latest_ts = self.get_target_latest_timestamp()
        for src_coll_name in self.get_source_list():
            _ts = _get_timestamp(src_coll_name)
            if _ts > target_latest_ts:
                _li.append(src_coll_name)
        return _li

    def get_changes(self, source_col, use_parallel=True):
        target_col = self._target_col
        source_col = self._db[source_col] if is_str(source_col) else source_col

        src = GeneDocMongoDBBackend(source_col)
        target = GeneDocMongoDBBackend(target_col)
        changes = diff_collections(target, src, use_parallel=use_parallel, step=self.step)
        if changes:
            changes['source'] = source_col.name
            changes['timestamp'] = _get_timestamp(source_col.name)
        return changes

    def apply_changes(self, changes):
        step = self.step
        target_col = self._target_col
        source_col = self._db[changes['source']]
        src = GeneDocMongoDBBackend(source_col)
        target = GeneDocMongoDBBackend(target_col)
        _timestamp = changes['timestamp']

        t0 = time.time()
        if changes['add']:
            logging.info("Adding {} new docs...".format(len(changes['add'])))
            t00 = time.time()
            for _ids in iter_n(changes['add'], step):
                _doc_li = src.mget_from_ids(_ids)
                for _doc in _doc_li:
                    _doc['_timestamp'] = _timestamp
                target.insert(_doc_li)
            logging.info("done. [{}]".format(timesofar(t00)))
        if changes['delete']:
            logging.info("Deleting {} discontinued docs...".format(len(changes['delete'])))
            t00 = time.time()
            target.remove_from_ids(changes['delete'], step=step)
            logging.info("done. [{}]".format(timesofar(t00)))

        if changes['update']:
            logging.info("Updating {} existing docs...".format(len(changes['update'])))
            t00 = time.time()
            i = 0
            t1 = time.time()
            for _diff in changes['update']:
                target.update_diff(_diff, extra={'_timestamp': _timestamp})
                i += 1
                if i > 1 and i % step == 0:
                    logging.info('\t{}\t{}'.format(i, timesofar(t1)))
                    t1 = time.time()
            logging.info("done. [{}]".format(timesofar(t00)))
        logging.info("\n")
        logging.info("Finished.", timesofar(t0))

    def verify_changes(self, changes):
        _timestamp = changes['timestamp']
        target = GeneDocMongoDBBackend(self._target_col)
        if changes['add']:
            logging.info('Verifying "add"...')
            # _cnt = self._target_col.find({'_id': {'$in': changes['add']}}).count()
            _cnt = target.count_from_ids(changes['add'])
            if _cnt == len(changes['add']):
                logging.info('...{}=={}...OK'.format(_cnt, len(changes['add'])))
            else:
                logging.info('...{}!={}...ERROR!!!'.format(_cnt, len(changes['add'])))
        if changes['delete']:
            logging.info('Verifying "delete"...')
            # _cnt = self._target_col.find({'_id': {'$in': changes['delete']}}).count()
            _cnt = target.count_from_ids(changes['delete'])
            if _cnt == 0:
                logging.info('...{}==0...OK'.format(_cnt))
            else:
                logging.info('...{}!=0...ERROR!!!'.format(_cnt))

        logging.info("Verifying all docs have timestamp...")
        _cnt = self._target_col.find({'_timestamp': {'$exists': True}}).count()
        _cnt_all = self._target_col.count()
        if _cnt == _cnt_all:
            logging.info('{}=={}...OK'.format(_cnt, _cnt_all))
        else:
            logging.info('ERROR!!!\n\t Should be "{}", but get "{}"'.format(_cnt_all, _cnt))

        logging.info("Verifying all new docs have updated timestamp...")
        cur = self._target_col.find({'_timestamp': {'$gte': _timestamp}}, projection={})
        _li1 = sorted(changes['add'] + [x['_id'] for x in changes['update']])
        _li2 = sorted([x['_id'] for x in cur])
        if _li1 == _li2:
            logging.info("{}=={}...OK".format(len(_li1), len(_li2)))
        else:
            logging.info('ERROR!!!\n\t Should be "{}", but get "{}"'.format(len(_li1), len(_li2)))

    def _get_cleaned_timestamp(self, timestamp):
        if is_str(timestamp):
            timestamp = datetime.strptime(timestamp, '%Y%m%d')
        assert isinstance(timestamp, datetime)
        return timestamp

    def get_change_history(self, before=None, after=None):
        _range = {}
        if before:
            before = self._get_cleaned_timestamp(before)
            _range['$lt'] = before
        if after:
            after = self._get_cleaned_timestamp(after)
            _range['$gt'] = after

        if _range:
            return self._target_col.find({'_timestamp': _range})
        else:
            raise ValueError('must provide either "before" for "after" argument.')

    def backup_timestamp(self, outfile=None, compress=True):
        '''backup "_id" and "_timestamp" fields into a output file.'''
        ts = time.strftime('%Y%m%d')
        outfile = outfile or self._target_col.name + '_tsbk_' + ts + '.txt'
        if compress:
            outfile += '.bz'
            import bz2
        logging.info('Backing up timestamps into "{}"...'.format(outfile))
        t0 = time.time()
        file_handler = bz2.BZ2File if compress else open
        with file_handler(outfile, 'w') as out_f:
            for doc in doc_feeder(self._target_col, step=100000, fields=['_timestamp']):
                out_f.write('{}\t{}\n'.format(doc['_id'], doc['_timestamp'].strftime('%Y%m%d')))
        logging.info("Done.", timesofar(t0))
        return outfile

    def get_timestamp_stats(self, returnresult=False, verbose=True):
        '''Return the count of each timestamps in _target_col.'''
        res = self._target_col.aggregate([{"$group": {"_id": "$_timestamp", "count": {"$sum": 1}}}])
        res = sorted([(x['_id'], x['count']) for x in res['result']], reverse=True)
        if verbose:
            for ts, cnt in res:
                logging.info('{}\t{}'.format(ts.strftime('%Y%m%d'), cnt))
        if returnresult:
            return res

    def get_target_latest_timestamp_0(self):
        ts_stats = self.get_timestamp_stats(returnresult=True, verbose=False)
        ts_stats.sort(reverse=True)
        latest_ts = ts_stats[0][0]
        assert ts_stats[0][1] > 0
        return latest_ts

    def get_target_latest_timestamp(self):
        cur = self._target_col.find(projection=['_timestamp']).sort([('_timestamp', -1)]).limit(1)
        try:
            doc = next(cur)
        except StopIteration:
            doc = {'_timestamp': datetime(1970, 1, 1, 0, 0, 0)}  # epoch to bootstrap if no previous
        cur.close()
        latest_ts = doc['_timestamp']
        return latest_ts


def mark_timestamp(timestamp):
    #.update({'_id': {'$in': xli1}}, {'$set': {'_timestamp': ts}}, multi=True)
    target = get_target_db()
    #genedoc_col = target.genedoc_mygene_allspecies_current
    genedoc_col = target.genedoc_mygene_xxxxx
    for doc in doc_feeder(genedoc_col):
        genedoc_col.update({'_id': doc['_id']},
                           {'$set': {'_timestamp': timestamp}},
                           manipulate=False, check_keys=False,
                           upsert=False, w=0)


def _get_timestamp(source_col, as_str=False):
    mat = re.search('_(\d{8})_\w{8}$', source_col)
    if mat:
        _timestamp = mat.group(1)
        if not as_str:
            _timestamp = datetime.strptime(_timestamp, '%Y%m%d')
        return _timestamp


def get_changes_stats(changes):
    for k in ['source', 'timestamp', 'add', 'delete', 'update']:
        if k in changes:
            v = changes[k]
            if isinstance(v, (list, dict)):
                v = len(v)
            logging.info("{}: {}".format(k, v))
    _update = changes['update']
    if _update:
        attrs = dict(add=set(), delete=set(), update=set())
        for _d in _update:
            for k in attrs.keys():
                if _d[k]:
                    attrs[k] |= set(_d[k])
        #plogging.info(attrs)
        logging.info('\n'.join(["\t{}: {} {}".format(k, len(attrs[k]), ', '.join(sorted(attrs[k]))) for k in attrs]))


def diff_two(col_1, col_2, use_parallel=True):
    target = get_target_db()
    b1 = GeneDocMongoDBBackend(target[col_1])
    b2 = GeneDocMongoDBBackend(target[col_2])
    return diff_collections(b1, b2, use_parallel=use_parallel)


def backup_timestamp_main():
    for config in ['genedoc_mygene', 'genedoc_mygene_allspecies']:
        sc = GeneDocSyncer(config)
        bkfile = sc.backup_timestamp()
        bkfile_key = 'genedoc_timestamp_bk/' + bkfile
        logging.info('Saving to S3: "{}"... '.format(bkfile_key))
        send_s3_file(bkfile, bkfile_key)
        logging.info('Done.')
        os.remove(bkfile)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == 'backup_timestamp':
        backup_timestamp_main()
        return

    if len(sys.argv) > 1:
        config = sys.argv[1]
    else:
        config = 'mygene_allspecies'
        #config = 'mygene_allspecies'
    if not config.startswith('genedoc_'):
        config = 'genedoc_' + config
    assert config in ['genedoc_mygene', 'genedoc_mygene_allspecies']
    use_parallel = '-p' in sys.argv
    no_confirm = '-b' in sys.argv

    t0 = time.time()
    sc = GeneDocSyncer(config)
    new_src_li = sc.get_new_source_list()
    if not new_src_li:
        logging.info("No new source collections need to update. Abort now.")
        return

    logging.info("Found {} new source collections need to update:".format(len(new_src_li)))
    logging.info("\n".join(['\t' + x for x in new_src_li]))

    if no_confirm or ask('Continue?') == 'Y':
        logfile = 'databuild_sync_{}_{}.log'.format(config, time.strftime('%Y%m%d'))
        logfile = os.path.join(LOG_FOLDER, logfile)
        setup_logfile(logfile)

        for src in new_src_li:
            t0 = time.time()
            logging.info("Current source collection:", src)
            ts = _get_timestamp(src, as_str=True)
            logging.info("Calculating changes... ")
            changes = sc.get_changes(src, use_parallel=use_parallel)
            logging.info("Done")
            get_changes_stats(changes)
            if no_confirm or ask("Continue to save changes...") == 'Y':
                if config == 'genedoc_mygene':
                    dumpfile = 'changes_{}.pyobj'.format(ts)
                else:
                    dumpfile = 'changes_{}_allspecies.pyobj'.format(ts)
                dump(changes, dumpfile)
                dumpfile_key = 'genedoc_changes/' + dumpfile
                logging.info('Saving to S3: "{}"... '.format(dumpfile_key))
                send_s3_file(dumpfile, dumpfile_key)
                logging.info('Done.')
                #os.remove(dumpfile)

            if no_confirm or ask("Continue to apply changes...") == 'Y':
                sc.apply_changes(changes)
                sc.verify_changes(changes)
            logging.info('=' * 20)
            logging.info("Finished.", timesofar(t0))


if __name__ == '__main__':
    main()
