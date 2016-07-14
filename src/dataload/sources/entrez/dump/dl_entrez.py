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
import re
from ftplib import FTP, error_temp

from biothings.utils.common import ask, timesofar, safewfile

src_path = os.path.split(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])[0]
sys.path.append(src_path)

from utils.common import setup_logfile, rmdashfr, hipchat_msg
from utils.mongo import get_src_dump
from config import DATA_ARCHIVE_ROOT, ASCP_ROOT, ARCHIVE_DATA, logger as logging

import glob
from parse_refseq_gbff import main as parse_refseq_gbff

TIMESTAMP = time.strftime('%Y%m%d')
if ARCHIVE_DATA:
    DATA_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/entrez', TIMESTAMP)
else:
    DATA_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/entrez/latest')


FILE_LIST = {
    'gene': {
        'url': 'ftp://ftp.ncbi.nih.gov/gene/DATA/',
        'files': [
            'gene_info.gz',
            'gene2accession.gz',
            'gene2refseq.gz',
            'gene2unigene',
            'gene2go.gz',
            'gene_history.gz',
            'gene2ensembl.gz'
        ]
    },

    'refseq': {
        'url': 'ftp://ftp.ncbi.nih.gov/refseq/',
        'files': [
            'H_sapiens/mRNA_Prot/human.*.rna.gbff.gz',
            'M_musculus/mRNA_Prot/mouse.*.rna.gbff.gz',
            'R_norvegicus/mRNA_Prot/rat.*.rna.gbff.gz',
            'D_rerio/mRNA_Prot/zebrafish.*.rna.gbff.gz',
            'X_tropicalis/mRNA_Prot/frog.*.rna.gbff.gz',
            'B_taurus/mRNA_Prot/cow.*.rna.gbff.gz',
            'S_scrofa/mRNA_Prot/pig.*.rna.gbff.gz'
        ]
    },

    'Homologene': {
        'url': 'ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/',
        'files': ['homologene.data']
    },

    'generif': {
        'url': 'ftp://ftp.ncbi.nih.gov/gene/GeneRIF/',
        'files': ['generifs_basic.gz']
    },
}


def _get_ascp_cmdline(url):
    '''
    ~/opt/aspera_connect/bin/ascp -QT -l640M -i  \
      ~/opt/aspera_connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nih.gov:/refseq/H_sapiens/mRNA_Prot/human.rna.gbff.gz .
    '''
    execpath = ASCP_ROOT + '/bin/ascp'
    keypath = ASCP_ROOT + '/etc/asperaweb_id_dsa.openssh'
    cmd = execpath + ' -QT -l640M -i ' + keypath + ' anonftp@'
    _url = url[6:]   # remove 'ftp://'
    _url = _url.replace('.gov/', '.gov:/')
    cmd = cmd + _url + ' .'
    return cmd


def _expand_wildchar_urls(url):
    '''
    do globbing on ftp url with wildchar to get the list of matching urls.
    if url contains no wildchar, return url as it is.
    '''
    if url.find('*') != -1 or url.find('?') != -1:
        mat = re.match('^ftp://([\w\.]+)(/.+$)', url)
        host, path = mat.groups()
        ftp = FTP(host)
        ftp.login()
        try:
            url_list = ftp.nlst(path)
            ftp.quit()
        except error_temp:
            # not found
            url_list = []
        url_list = ['ftp://' + host + fn for fn in url_list]
        url_list.sort()
    else:
        url_list = [url]
    return url_list


def _expand_refseq_files():
    '''expand refseq url list with wildchar.'''
    base = FILE_LIST['refseq']['url']
    _files = []
    for f in FILE_LIST['refseq']['files']:
        _files.extend(_expand_wildchar_urls(base + f))
    _i = len(base)
    _files = [fn[_i:] for fn in _files]    # remove leading "ftp://ftp.ncbi.nih.gov/refseq/"
    FILE_LIST['refseq']['files'] = _files
    return _files


def download(path, no_confirm=False):
    out = []
    orig_path = os.getcwd()
    try:
        _expand_refseq_files()
        for subfolder in FILE_LIST:
            filedata = FILE_LIST[subfolder]
            baseurl = filedata['url']
            data_folder = os.path.join(path, subfolder)
            if not os.path.exists(data_folder):
                os.mkdir(data_folder)

            for f in filedata['files']:
                url = baseurl + f
                os.chdir(data_folder)
                filename = os.path.split(f)[1]
                if os.path.exists(filename):
                    if no_confirm or ask('Remove existing file "%s"?' % filename) == 'Y':
                        os.remove(filename)
                    else:
                        logging.info("Skipped!")
                        continue
                logging.info('Downloading "%s"...' % f)
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


def parse_gbff(path):
    refseq_folder = os.path.join(path, 'refseq')
    gbff_files = glob.glob(os.path.join(refseq_folder, '*.rna.gbff.gz'))
    assert len(gbff_files) >= 15, 'Missing "*.gbff.gz" files? Found %d in %s:\n%s' % (len(gbff_files), refseq_folder, '\n'.join(gbff_files))
    parse_refseq_gbff(refseq_folder)


def redo_parse_gbff(path):
    '''call this function manually to re-start the parsing step and set src_dump.
       This is used when main() is broken at parsing step, then parsing need to be re-started
       after the fix.
    '''
    #mark the download starts
    src_dump = get_src_dump()

    t0 = time.time()
    t_download = timesofar(t0)
    t1 = time.time()
    #mark parsing starts
    src_dump.update({'_id': 'entrez'}, {'$set': {'status': 'parsing'}})
    parse_gbff(path)
    t_parsing = timesofar(t1)
    t_total = timesofar(t0)

    #mark the download finished successfully
    _updates = {
        'status': 'success',
        'time': {
            'download': t_download,
            'parsing': t_parsing,
            'total': t_total
        },
        'pending_to_upload': True    # a flag to trigger data uploading
    }

    src_dump.update({'_id': 'entrez'}, {'$set': _updates})


def main():
    no_confirm = True   # set it to True for running this script automatically without intervention.

    if not ARCHIVE_DATA:
        rmdashfr(DATA_FOLDER)

    if not os.path.exists(DATA_FOLDER):
        os.makedirs(DATA_FOLDER)
    else:
        if not (no_confirm or len(os.listdir(DATA_FOLDER)) == 0 or ask('DATA_FOLDER (%s) is not empty. Continue?' % DATA_FOLDER) == 'Y'):
            sys.exit()

    logfile = os.path.join(DATA_FOLDER, 'entrez_dump.log')
    setup_logfile(logfile)

    #mark the download starts
    src_dump = get_src_dump()
    doc = {'_id': 'entrez',
           'timestamp': TIMESTAMP,
           'data_folder': DATA_FOLDER,
           'logfile': logfile,
           'status': 'downloading'}
    src_dump.save(doc)
    t0 = time.time()
    download(DATA_FOLDER, no_confirm=no_confirm)
    t_download = timesofar(t0)
    t1 = time.time()
    #mark parsing starts
    src_dump.update({'_id': 'entrez'}, {'$set': {'status': 'parsing'}})
    parse_gbff(DATA_FOLDER)
    t_parsing = timesofar(t1)
    t_total = timesofar(t0)

    #mark the download finished successfully
    _updates = {
        'status': 'success',
        'time': {
            'download': t_download,
            'parsing': t_parsing,
            'total': t_total
        },
        'pending_to_upload': True    # a flag to trigger data uploading
    }

    src_dump.update({'_id': 'entrez'}, {'$set': _updates})


if __name__ == '__main__':
    try:
        main()
        hipchat_msg('"entrez" downloader finished successfully',color='green')
    except Exception as e:
        import traceback
        logging.error("Error while downloading: %s" % traceback.format_exc())
        hipchat_msg('"entrez" downloader failed: %s' % e,color='red')
        sys.exit(255)
