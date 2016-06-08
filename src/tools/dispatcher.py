from __future__ import print_function
import time
from subprocess import Popen
from datetime import datetime
import dispatch

from biothings.utils.common import timesofar
from utils.common import src_path
from utils.mongo import src_clean_archives, target_clean_collections
from dataload.dispatch import (check_mongo, get_process_info, src_dump,
                               mark_upload_started, mark_upload_done)
from dataload.dispatch import dispatch as dispatch_src_upload

source_update_available = dispatch.Signal(providing_args=["src_to_update"])
source_upload_success = dispatch.Signal(providing_args=["src_name"])
source_upload_failed = dispatch.Signal(providing_args=["src_name"])
genedoc_merged = dispatch.Signal()
es_indexed = dispatch.Signal()

try:
    from utils.common import hipchat_msg
except:
    hipchat_msg = None


class GeneDocDispatcher:
    running_processes_upload = {}
    idle = True

    def check_src_dump(self):
        src_to_update_li = check_mongo()
        if src_to_update_li:
            print('\nDispatcher:  found pending jobs ', src_to_update_li)
            for src_to_update in src_to_update_li:
                source_update_available.send(sender=self, src_to_update=src_to_update)

    @classmethod
    def handle_src_upload(self, src_to_update, **kwargs):
        mark_upload_started(src_to_update)
        p = dispatch_src_upload(src_to_update)
        src_dump.update({'_id': src_to_update}, {"$set": {"upload.pid": p.pid}})
        p.t0 = time.time()
        self.running_processes_upload[src_to_update] = p

    def check_src_upload(self):
        running_processes = self.running_processes_upload
        jobs_finished = []
        if running_processes:
            self.idle = True
            print('Dispatcher:  {} active job(s)'.format(len(running_processes)))
            print(get_process_info(running_processes))

        for src in running_processes:
            p = running_processes[src]
            returncode = p.poll()
            if returncode is None:
                p.log_f.flush()
            else:
                t1 = round(time.time() - p.t0, 0)
                d = {
                    'upload.returncode': returncode,
                    'upload.timestamp': datetime.now(),
                    'upload.time_in_s': t1,
                    'upload.time': timesofar(p.t0),
                    'upload.logfile': p.logfile,
                    'upload.status': "success" if returncode == 0 else "failed"
                }
                mark_upload_done(src, d)
                jobs_finished.append(src)
                p.log_f.close()

                if returncode == 0:
                    msg = 'Dispatcher:  "{}" uploader finished successfully with code {} (time: {})'.format(src, returncode, timesofar(p.t0, t1=t1))
                    print(msg)
                    if hipchat_msg:
                        msg += '<a href="http://su07:8000/log/dump/{}">dump log</a>'.format(src)
                        msg += '<a href="http://su07:8000/log/upload/{}">upload log</a>'.format(src)
                        hipchat_msg(msg, message_format='html',color="green")
                    source_upload_success.send(self, src_name=src)
                else:
                    msg = 'Dispatcher:  "{}" uploader failed with code {} (time: {}s)'.format(src, returncode, t1)
                    print(msg)
                    if hipchat_msg:
                        hipchat_msg(msg,color="red")
                    source_upload_failed.send(self, src_name=src)

        for src in jobs_finished:
            del running_processes[src]

    @classmethod
    def handle_src_upload_success(self, src_name, **kwargs):
        '''when "entrez" src upload is done, trigger src_build tasks.'''
        if src_name == 'entrez':
            self.handle_src_build()

    @classmethod
    def handle_src_upload_failed(self, src_name, **kwargs):
        pass

    @classmethod
    def handle_src_build(self):

        #cleanup src and target collections
        src_clean_archives(noconfirm=True)
        target_clean_collections(noconfirm=True)

        for config in ('mygene', 'mygene_allspecies'):
            t0 = time.time()
            p = Popen(['python', '-m', 'databuild.builder', config], cwd=src_path)
            returncode = p.wait()
            t = timesofar(t0)
            if returncode == 0:
                msg = 'Dispatcher:  "{}" builder finished successfully with code {} (time: {})'.format(config, returncode, t)
                color = "green"
            else:
                msg = 'Dispatcher:  "{}" builder failed successfully with code {} (time: {})'.format(config, returncode, t)
                color = "red"
            print(msg)
            if hipchat_msg:
                msg += '<a href="http://su07:8000/log/build/{}">build log</a>'.format(config)
                hipchat_msg(msg, message_format='html',color=color)

            assert returncode == 0, "Subprocess failed. Check error above."
        genedoc_merged.send(self)

    @classmethod
    def handle_genedoc_merged(self, **kwargs):
        for config in ('mygene', 'mygene_allspecies'):
            t0 = time.time()
            p = Popen(['python', '-m', 'databuild.sync', config, '-p', '-b'], cwd=src_path)
            returncode = p.wait()
            t = timesofar(t0)
            if returncode == 0:
                msg = 'Dispatcher:  "{}" syncer finished successfully with code {} (time: {})'.format(config, returncode, t)
                color = "green"
            else:
                msg = 'Dispatcher:  "{}" syncer failed successfully with code {} (time: {})'.format(config, returncode, t)
                color = "red"
            print(msg)
            if hipchat_msg:
                msg += '<a href="http://su07:8000/log/sync/{}">sync log</a>'.format(config)
                hipchat_msg(msg, message_format='html',color=color)

            assert returncode == 0, "Subprocess failed. Check error above."

    def check_src_build(self):
        pass

    def check_src_index(self):
        pass

    def main(self):
        import sys
        _flag = len(sys.argv) == 2
        while 1:
            self.idle = True
            self.check_src_dump()
            self.check_src_upload()
            if _flag:
                if sys.argv[1] == 'flag1':
                    # resend source_upload_success signal
                    source_upload_success.send(self, src_name='entrez')
                elif sys.argv[1] == 'flag2':
                    # resend genedoc_merged signal
                    genedoc_merged.send(self)
                _flag = False
            self.check_src_build()
            self.check_src_index()

            if self.idle:
                print('\b' * 50, end='')
                for i in range(100):
                    print('\b' * 2 + [chr(8212), '\\', '|', '/'][i % 4], end='')
                    time.sleep(0.1)
            else:
                time.sleep(10)


source_update_available.connect(GeneDocDispatcher.handle_src_upload)
source_upload_success.connect(GeneDocDispatcher.handle_src_upload_success)
source_upload_failed.connect(GeneDocDispatcher.handle_src_upload_failed)
genedoc_merged.connect(GeneDocDispatcher.handle_genedoc_merged)

if __name__ == '__main__':
    GeneDocDispatcher().main()
