# Copyright [2010-2015] [Chunlei Wu]
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


'''
Usage:
     python dl_ensembl.py check      # Check the lastest Ensembl/BioMart version
     python dl_ensembl.py <ensembl_ver>   # perform the actual download
'''
from __future__ import print_function
import sys
import os
import time
from ftplib import FTP

import requests
from biothings.utils.common import ask, timesofar, safewfile

src_path = os.path.split(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])[0]
sys.path.append(src_path)
from utils.common import LogPrint, rmdashfr
from utils.mongo import get_src_dump
from utils.dataload import tab2list
from config import DATA_ARCHIVE_ROOT, ARCHIVE_DATA



ENSEMBL_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/ensembl')

MART_URL = "http://www.biomart.org/biomart/martservice"
MART_URL = "http://uswest.ensembl.org/biomart/martservice"

XML_QUERY_TEMPLATE_EXAMPLE = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >

    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "ensembl_transcript_id" />
        <Attribute name = "ensembl_peptide_id" />
    </Dataset>
</Query>
'''

XML_QUERY_TEMPLATE = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "%(virtual_schema)s" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "%(dataset)s" interface = "default" >
        %(filters)s
        %(attributes)s
    </Dataset>
</Query>
'''

DataSet_D = {
    'human': 'hsapiens_gene_ensembl',
    'mouse': 'mmusculus_gene_ensembl',
    'rat': 'rnorvegicus_gene_ensembl',
    'fruitfly': 'dmelanogaster_gene_ensembl',
    'nematode': 'celegans_gene_ensembl',
    'zebrafish': 'drerio_gene_ensembl',
    'thale-cress': '',   # not available
    'frog': 'xtropicalis_gene_ensembl',
    'pig': 'sscrofa_gene_ensembl',
}


def chk_latest_mart_version():
    ftp = FTP('ftp.ensembl.org')
    ftp.login()
    #release_li = ftp.nlst('/pub/release-*')
    release_li = [x for x in ftp.nlst('/pub') if x.startswith('/pub/release-')]
    return sorted([int(fn.split('-')[-1]) for fn in release_li])[-1]


def _to_int(taxid):
    try:
        return int(taxid)
    except:
        return None


def get_all_species(release):
    import tempfile
    outfile = tempfile.mktemp() + '.txt.gz'
    try:
        print('Downloading "species.txt.gz"...', end=' ')
        out_f = open(outfile, 'wb')
        ftp = FTP('ftp.ensembl.org')
        ftp.login()
        species_file = '/pub/release-%s/mysql/ensembl_production_%s/species.txt.gz' % (release, release)
        ftp.retrbinary("RETR " + species_file, out_f.write)
        out_f.close()
        print('Done.')

        #load saved file
        print('Parsing "species.txt.gz"...', end=' ')
        species_li = tab2list(outfile, (1, 2, 7), header=0)   # db_name,common_name,taxid
        species_li = [x[:-1] + [_to_int(x[-1])] for x in species_li]
        print('Done.')
    finally:
        os.remove(outfile)
        pass

    return species_li


class MartException(Exception):
    pass


class BioMart(object):
    url = MART_URL
    template = XML_QUERY_TEMPLATE

    def __init__(self, release=None, species_li=None):
        self.species_li = species_li
        self.release = release
        self.no_confirm = False
        if not self.species_li and self.release:
            self.get_species_list()

    def _query(self, *args, **kwargs):
        req = requests.Request(*args, **kwargs)
        res = requests.session().send(req.prepare())
        if res.status_code != 200:
            raise MartException(res)
        if res.text.startswith('Query ERROR:'):
            raise MartException(res.text)
        return res.text

    def _make_query_xml(self, dataset, attributes, filters=None):
        attrib_xml = '\n'.join(['<Attribute name = "%s" />' % attrib for attrib in attributes])
        if filters:
            filter_xml = '\n'.join(['<Filter name = "%s" excluded = "0"/>' % filter for filter in filters])
        else:
            filter_xml = ''
        return XML_QUERY_TEMPLATE % {'virtual_schema': self.get_virtual_schema(),
                                     'dataset': dataset,
                                     'attributes': attrib_xml,
                                     'filters': filter_xml}

    def get_species_list(self):
        self.species_li = get_all_species(self.release)
        return self.species_li

    def get_virtual_schema(self):
        return 'default'

    def _get_species_table_prefix(self, species):
        x = species.split('_')
        return x[0][0] + x[1]

    def get_dataset_name(self, species):
        return '%s_gene_ensembl' % self._get_species_table_prefix(species[0])

    def chk_latest_mart_version(self):
        pass

    def query_mart(self, xml):
        return self._query('POST', self.url, data='query=%s\n' % xml)

    def get_registry(self):
        return self._query('GET', self.url + '?type=registry')

    def get_datasets(self):
        con = self._query('GET', self.url + '?type=datasets&mart=ensembl')
        out = []
        for line in con.split('\n'):
            line = line.strip()
            if line != '':
                ld = line.split('\t')
                out.append((ld[1], ld[2], ld[4]))
        return out

    def _fetch_data(self, outfile, attributes, filters='', header=None, debug=False):
        cnt_all = 0
        out_f, outfile = safewfile(outfile, prompt=(not self.no_confirm), default='O')
        if header:
            out_f.write('\t'.join(header) + '\n')
        print('Dumping "%s"...' % os.path.split(outfile)[1])
        for species in self.species_li:
            dataset = self.get_dataset_name(species)
            taxid = species[2]
            if not dataset:
                continue
            xml = self._make_query_xml(dataset, attributes=attributes, filters=filters)
            if debug:
                print(xml)
            try:
                con = self.query_mart(xml)
            except MartException:
                import traceback
                err_msg = traceback.format_exc()
                print(species[0], err_msg)
                continue
            cnt = 0
            for line in con.split('\n'):
                if line.strip() != '':
                    out_f.write(str(taxid) + '\t' + line + '\n')
                    cnt += 1
                    cnt_all += 1
            print(species[0], cnt)
        out_f.close()
        print("Total: %d" % cnt_all)

    def get_gene__main(self, outfile, debug=False):
        header = ['taxonomy_id',
                  'ensembl_gene_id',
                  'symbol',
                  'gene_chrom_start', 'gene_chrom_end', 'chr_name', 'chrom_strand',
                  'description']
        attributes = ["ensembl_gene_id",
                      "external_gene_name",   # symbols, called "external_gene_id" before release 76
                      "start_position", "end_position", "chromosome_name", "strand",
                      "description"]
        self._fetch_data(outfile, attributes, header=header, debug=debug)

    def get_translation__main(self, outfile, debug=False):
        header = ['taxonomy_id',
                  'gene_stable_id',
                  'transcript_stable_id',
                  'translation_stable_id']
        attributes = ["ensembl_gene_id",
                      "ensembl_transcript_id",
                      "ensembl_peptide_id"]
        self._fetch_data(outfile, attributes, header=header, debug=debug)

    def get_xref_entrezgene(self, outfile, debug=False):
        header = ['taxonomy_id',
                  'gene_stable_id',
                  'dbprimary_id']
        attributes = ["ensembl_gene_id",
                      "entrezgene"]
        filters = ["with_entrezgene"]
        self._fetch_data(outfile, attributes, filters, header=header, debug=debug)

    def get_profile(self, outfile, debug=False):
        header = ['taxonomy_id',
                  'gene_stable_id',
                  'transcript_stable_id',
                  'translation_stable_id',
                  'profile_id']
        attributes = ["ensembl_gene_id",
                      "ensembl_transcript_id",
                      "ensembl_peptide_id",
                      "profile"]
        filters = ["with_profile"]
        self._fetch_data(outfile, attributes, filters, header=header, debug=debug)

    def get_interpro(self, outfile, debug=False):
        header = ['taxonomy_id',
                  'gene_stable_id',
                  'transcript_stable_id',
                  'translation_stable_id',
                  'interpro_id', 'short_description', 'description']
        attributes = ["ensembl_gene_id",
                      "ensembl_transcript_id",
                      "ensembl_peptide_id",
                      "interpro", "interpro_short_description", "interpro_description"]
        filters = ["with_interpro"]
        self._fetch_data(outfile, attributes, filters, header=header, debug=debug)

    def get_pfam(self, outfile, debug=False):
        header = ['taxonomy_id',
                  'gene_stable_id',
                  'transcript_stable_id',
                  'translation_stable_id',
                  'pfam']
        attributes = ["ensembl_gene_id",
                      "ensembl_transcript_id",
                      "ensembl_peptide_id",
                      "pfam"]
        filters = ["with_protein_feature_pfam"]
        self._fetch_data(outfile, attributes, filters, header=header, debug=debug)


def main():
    if len(sys.argv) == 2 and sys.argv[1] == 'check':
        print("Checking latest mart_version:\t", end=' ')
        mart_version = chk_latest_mart_version()
        print(mart_version)
        return

    if len(sys.argv) > 1:
        mart_version = sys.argv[1]
    else:
        print("Checking latest mart_version:\t", end=' ')
        mart_version = chk_latest_mart_version()
        print(mart_version)

    BM = BioMart()
    BM.species_li = get_all_species(mart_version)
    if ARCHIVE_DATA:
        DATA_FOLDER = os.path.join(ENSEMBL_FOLDER, mart_version)
    else:
        DATA_FOLDER = os.path.join(ENSEMBL_FOLDER,"latest")
        rmdashfr(DATA_FOLDER)

    if not os.path.exists(DATA_FOLDER):
        os.makedirs(DATA_FOLDER)
    else:
        if not (len(os.listdir(DATA_FOLDER)) == 0 or ask('DATA_FOLDER (%s) is not empty. Continue?' % DATA_FOLDER) == 'Y'):
            return
    log_f, logfile = safewfile(os.path.join(DATA_FOLDER, 'ensembl_mart_%s.log' % mart_version))
    sys.stdout = LogPrint(log_f, timestamp=True)

    BM.get_gene__main(os.path.join(DATA_FOLDER, 'gene_ensembl__gene__main.txt'))
    BM.get_translation__main(os.path.join(DATA_FOLDER, 'gene_ensembl__translation__main.txt'))
    BM.get_xref_entrezgene(os.path.join(DATA_FOLDER, 'gene_ensembl__xref_entrezgene__dm.txt'))

    BM.get_profile(os.path.join(DATA_FOLDER, 'gene_ensembl__prot_profile__dm.txt'))
    BM.get_interpro(os.path.join(DATA_FOLDER, 'gene_ensembl__prot_interpro__dm.txt'))
    sys.stdout.close()


def main_cron():
    no_confirm = True   # set it to True for running this script automatically without intervention.

    src_dump = get_src_dump()
    print("Checking latest mart_version:\t", end=' ')
    mart_version = chk_latest_mart_version()
    print(mart_version)

    doc = src_dump.find_one({'_id': 'ensembl'})
    if doc and 'release' in doc and mart_version <= doc['release']:
        data_file = os.path.join(doc['data_folder'], 'gene_ensembl__gene__main.txt')
        if os.path.exists(data_file):
            print("No newer release found. Abort now.")
            sys.exit(0)

    DATA_FOLDER = os.path.join(ENSEMBL_FOLDER, str(mart_version))
    if not os.path.exists(DATA_FOLDER):
        os.makedirs(DATA_FOLDER)
    else:
        if not (no_confirm or len(os.listdir(DATA_FOLDER)) == 0 or ask('DATA_FOLDER (%s) is not empty. Continue?' % DATA_FOLDER) == 'Y'):
            sys.exit(0)

    log_f, logfile = safewfile(os.path.join(DATA_FOLDER, 'ensembl_mart_%s.log' % mart_version), prompt=(not no_confirm), default='O')
    sys.stdout = LogPrint(log_f, timestamp=True)

    #mark the download starts
    doc = {'_id': 'ensembl',
           'release': mart_version,
           'timestamp': time.strftime('%Y%m%d'),
           'data_folder': DATA_FOLDER,
           'logfile': logfile,
           'status': 'downloading'}
    src_dump.save(doc)
    t0 = time.time()

    try:
        BM = BioMart()
        BM.species_li = get_all_species(mart_version)
        BM.get_gene__main(os.path.join(DATA_FOLDER, 'gene_ensembl__gene__main.txt'))
        BM.get_translation__main(os.path.join(DATA_FOLDER, 'gene_ensembl__translation__main.txt'))
        BM.get_xref_entrezgene(os.path.join(DATA_FOLDER, 'gene_ensembl__xref_entrezgene__dm.txt'))

        BM.get_profile(os.path.join(DATA_FOLDER, 'gene_ensembl__prot_profile__dm.txt'))
        BM.get_interpro(os.path.join(DATA_FOLDER, 'gene_ensembl__prot_interpro__dm.txt'))
        BM.get_pfam(os.path.join(DATA_FOLDER, 'gene_ensembl__prot_pfam__dm.txt'))
    finally:
        sys.stdout.close()

    #mark the download finished successfully
    _updates = {
        'status': 'success',
        'time': timesofar(t0),
        'pending_to_upload': True    # a flag to trigger data uploading
    }
    src_dump.update({'_id': 'ensembl'}, {'$set': _updates})


if __name__ == '__main__':
    main_cron()
