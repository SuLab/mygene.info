# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/

# Copyright [2014-] [Chunlei Wu]
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
from ftplib import FTP
from io import StringIO
from dataload.data_dump.dl_entrez import _get_ascp_cmdline, _expand_wildchar_urls
from utils.common import safewfile, setup_logfile
from biothings.utils.common import timesofar, ask
from utils.mongo import get_src_dump
from config import DATA_ARCHIVE_ROOT, logger as logging


src_path = os.path.split(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])[0]
sys.path.append(src_path)

timestamp = time.strftime('%Y%m%d')
# DATA_FOLDER=os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/ucsc', timestamp)
REFSEQ_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/refseq')

FTP_SERVER = 'ftp.ncbi.nih.gov'
BASE_PATH = '/refseq/release/complete/'
RELEASE_FILE = '/refseq/release/RELEASE_NUMBER'
DATA_FILE = 'complete.*.rna.gbff.gz'


def get_refseq_release():
    ftp = FTP(FTP_SERVER)
    ftp.login()
    release = StringIO()
    ftp.retrlines('RETR ' + RELEASE_FILE, release.write)
    ftp.quit()
    return release.getvalue()


def check_refseq_release():
    refseq_release = get_refseq_release()
    src_dump = get_src_dump()
    doc = src_dump.find_one({'_id': 'refseq'})
    if doc and 'release' in doc and refseq_release <= doc['release']:
        data_file = os.path.join(doc['data_folder'], 'complete.109.rna.gbff.gz')
        if os.path.exists(data_file):
            logging.info("No newer release found. Abort now.")
            sys.exit(0)


def download(path, release, no_confirm=False):
    out = []
    orig_path = os.getcwd()
    try:
        data_folder = os.path.join(path, release)
        if not os.path.exists(data_folder):
            os.mkdir(data_folder)

        _url = 'ftp://' + FTP_SERVER + BASE_PATH + DATA_FILE
        url_li = _expand_wildchar_urls(_url)
        logging.info('Found {} "{}" files to download.'.format(len(url_li), DATA_FILE))

        for url in url_li:
            os.chdir(data_folder)
            filename = os.path.split(url)[1]
            if os.path.exists(filename):
                if no_confirm or ask('Remove existing file "%s"?' % filename) == 'Y':
                    os.remove(filename)
                else:
                    logging.info("Skipped!")
                    continue
            logging.info('Downloading "%s"...' % filename)
            #cmdline = 'wget %s' % url
            #cmdline = 'axel -a -n 5 %s' % url   #faster than wget using 5 connections
            cmdline = _get_ascp_cmdline(url)
            return_code = os.system(cmdline)
            #return_code = 0;print cmdline    #for testing
            if return_code == 0:
                logging.info("Success.")
            else:
                logging.info("Failed with return code (%s)." % return_code)
                out.append((url, return_code))
            logging.info("=" * 50)
    finally:
        os.chdir(orig_path)

    return out


def main_cron():
    no_confirm = True   # set it to True for running this script automatically without intervention.

    logging.info("Checking latest refseq release:\t", end='')
    refseq_release = get_refseq_release()
    logging.info(refseq_release)

    src_dump = get_src_dump()
    doc = src_dump.find_one({'_id': 'refseq'})
    if doc and 'release' in doc and refseq_release <= doc['release']:
        data_file = os.path.join(doc['data_folder'], 'complete.109.rna.gbff.gz')
        if os.path.exists(data_file):
            logging.info("No newer release found. Abort now.")
            sys.exit(0)

    DATA_FOLDER = os.path.join(REFSEQ_FOLDER, str(refseq_release))
    if not os.path.exists(DATA_FOLDER):
        os.makedirs(DATA_FOLDER)
    else:
        if not (no_confirm or len(os.listdir(DATA_FOLDER)) == 0 or ask('DATA_FOLDER (%s) is not empty. Continue?' % DATA_FOLDER) == 'Y'):
            sys.exit(0)

    logfile = os.path.join(DATA_FOLDER, 'refseq_dump.log')
    setup_logfile(logfile)

    #mark the download starts
    doc = {'_id': 'refseq',
           'release': refseq_release,
           'timestamp': time.strftime('%Y%m%d'),
           'data_folder': DATA_FOLDER,
           'logfile': logfile,
           'status': 'downloading'}
    src_dump.save(doc)
    t0 = time.time()

    try:
        download(DATA_FOLDER, refseq_release, no_confirm=no_confirm)
    finally:
        sys.stdout.close()

    #mark the download finished successfully
    _updates = {
        'status': 'success',
        'time': timesofar(t0),
        'pending_to_upload': True    # a flag to trigger data uploading
    }
    src_dump.update({'_id': 'refseq'}, {'$set': _updates})


if __name__ == '__main__':
    main_cron()
