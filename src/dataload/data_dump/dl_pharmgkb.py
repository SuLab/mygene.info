# Copyright [2010-2013] [Chunlei Wu]
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import sys
import os
import os.path
import time
from datetime import datetime

import requests
from biothings.utils.common import ask, timesofar, safewfile

src_path = os.path.split(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])[0]
sys.path.append(src_path)
from utils.common import LogPrint, rmdashfr
from utils.mongo import get_src_dump
from config import DATA_ARCHIVE_ROOT, ARCHIVE_DATA


TIMESTAMP = time.strftime('%Y%m%d')
if ARCHIVE_DATA:
    DATA_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/pharmgkb', TIMESTAMP)
else:
    DATA_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/pharmgkb/latest')

#GENES_URL = 'http://www.pharmgkb.org/commonFileDownload.action?filename=genes.zip'
GENES_URL = 'http://www.pharmgkb.org/download.do?objId=genes.zip&dlCls=common'


def download(no_confirm=False):
    orig_path = os.getcwd()
    try:
        os.chdir(DATA_FOLDER)
        filename = 'genes.zip'
        url = GENES_URL
        if os.path.exists(filename):
            if no_confirm or ask('Remove existing file "%s"?' % filename) == 'Y':
                os.remove(filename)
            else:
                print("Skipped!")
                return
        print('Downloading "%s"...' % filename)
        cmdline = 'wget "%s" -O %s' % (url, filename)
        #cmdline = 'axel -a -n 5 %s' % url   #faster than wget using 5 connections
        return_code = os.system(cmdline)
        if return_code == 0:
            print("Success.")
        else:
            print("Failed with return code (%s)." % return_code)
        print("=" * 50)
    finally:
        os.chdir(orig_path)


def check_header():
    req = requests.Request('HEAD', GENES_URL)
    res = requests.session().send(req.prepare())
    assert res.status_code == 200, "Error: fail to access download url."
    lastmodified = res.headers.get('last-modified', '')
    if lastmodified:
        # an example: 'last-modified': 'Thu, 06 Dec 2012 11:01:50 GMT'
        lastmodified = datetime.strptime(lastmodified, "%a, %d %b %Y %H:%M:%S %Z")

    return lastmodified


if __name__ == '__main__':
    no_confirm = True   # set it to True for running this script automatically without intervention.

    src_dump = get_src_dump()
    lastmodified = check_header()
    doc = src_dump.find_one({'_id': 'pharmgkb'})
    if doc and 'lastmodified' in doc and lastmodified <= doc['lastmodified']:
        data_file = os.path.join(doc['data_folder'], 'genes.zip')
        if os.path.exists(data_file):
            print("No newer file found. Abort now.")
            sys.exit(0)

    if not ARCHIVE_DATA:
        rmdashfr(DATA_FOLDER)

    if not os.path.exists(DATA_FOLDER):
        os.makedirs(DATA_FOLDER)
    else:
        if not (no_confirm or len(os.listdir(DATA_FOLDER)) == 0 or ask('DATA_FOLDER (%s) is not empty. Continue?' % DATA_FOLDER) == 'Y'):
            sys.exit(0)

    log_f, logfile = safewfile(os.path.join(DATA_FOLDER, 'pharmgkb_dump.log'), prompt=(not no_confirm), default='O')
    sys.stdout = LogPrint(log_f, timestamp=True)

    #mark the download starts
    doc = {'_id': 'pharmgkb',
           'timestamp': TIMESTAMP,
           'data_folder': DATA_FOLDER,
           'lastmodified': lastmodified,
           'logfile': logfile,
           'status': 'downloading'}
    src_dump.save(doc)
    t0 = time.time()
    try:
        download(no_confirm)
    finally:
        sys.stdout.close()
    #mark the download finished successfully
    _updates = {
        'status': 'success',
        'time': timesofar(t0),
        'pending_to_upload': True    # a flag to trigger data uploading
    }
    src_dump.update({'_id': 'pharmgkb'}, {'$set': _updates})
