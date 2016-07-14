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
import os.path
import glob
import time
from Bio import SeqIO

from utils.common import SubStr
from utils.dataload import anyfile
from config import DATA_ARCHIVE_ROOT

src_path = os.path.split(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0])[0]
sys.path.append(src_path)

# TODO: what if process is done day after ?
timestamp = time.strftime('%Y%m%d')
DATA_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/entrez', timestamp, 'refseq')


class GBFFParser():
    def __init__(self, infile):
        self.infile = infile
        self.in_f = anyfile(self.infile)

    def parse(self):
        out_li = []
        for rec in SeqIO.parse(self.in_f, 'genbank'):
            geneid = self.get_geneid(rec)
            if geneid:
                summary = self.get_summary(rec)
                ec_list = self.get_ec_numbers(rec)
                if summary or ec_list:
                    out_li.append((geneid, summary, ec_list))
        return out_li

    def get_geneid(self, rec):
        '''Return geneid as integer, None if not found.'''
        geneid = None
        gene_feature = [x for x in rec.features if x.type == 'gene']
        # NCBI has now fixed this issue (https://twitter.com/kdpru/status/474673626730741761)
        # if len(gene_feature) == 0 and rec.id == 'NR_001526.1':
        #     print "Known error for NR_001526.1. Fixed."
        #     return '252949'         # a temp fix for this wrong rec from NCBI
        assert len(gene_feature) == 1, '#: {}, id: {}'.format(len(gene_feature), rec.id)
        gene_feature = gene_feature[0]
        db_xref = gene_feature.qualifiers.get('db_xref', None)
        if db_xref:
            x = [x for x in db_xref if x.startswith('GeneID:')]
            if len(x) == 1:
                geneid = int(SubStr(x[0], 'GeneID:'))
        return geneid

    def get_summary(self, rec):
        '''Return summary string if available, return '' otherwise.'''
        summary = ''
        comment = rec.annotations.get('comment', None)
        if comment:
            if comment.find('Summary:') != -1:
                summary = SubStr(comment, 'Summary: ',).replace('\n', ' ')
                for end_str in [# '[provided by RefSeq].',
                                '[provided by ',
                                '[supplied by ',
                                '##',
                                # '[RGD',
                                'COMPLETENESS:',
                                'Sequence Note:',
                                'Transcript Variant:',
                                'Publication Note:']:
                    if summary.find(end_str) != -1:
                        summary = SubStr(summary, end_string=end_str)
                summary = summary.strip()
        return summary

    def get_ec_numbers(self, rec):
        '''Return a list of EC numbers if available, return [] if not found.'''
        ec_list = []
        cds_feature = [x for x in rec.features if x.type == 'CDS']
        if cds_feature:
            assert len(cds_feature) == 1, self.get_geneid(rec)
            cds_feature = cds_feature[0]
            ec_list = cds_feature.qualifiers.get('EC_number', [])
            assert isinstance(ec_list, list)
#            ec_list = [SubStr(x, 'EC_number="', '"') for x in ec_qualifiers]
        return ec_list


def dump_object(obj, outfile):
    '''Dump a python object to a output file
    '''
    import pickle
    import bz2

    protocol = 2
    print("Dumping into \"%s\"..." % outfile, end='')
    out_f = open(outfile, 'wb')
    out_f.write(bz2.compress(pickle.dumps(obj, protocol)))
    out_f.close()
    print("Done!")


def output_gene2summary(out_d, outfile):
    '''Output tab delimited file for gene summary, with two column
       gene summary
       (no header line)
    '''
    out_f = open(outfile, 'w')
    for species in out_d:
        out_li = []
        for rec in out_d[species]:
            geneid = rec[0]
            summary = rec[1]
            if summary:
                out_li.append((geneid, summary))

        for geneid, summary in sorted(set(out_li)):
            out_f.write('%s\t%s\n' % (geneid, summary))
    out_f.close()


def output_gene2ec(out_d, outfile):
    '''Output tab delimited file for gene EC numbers, with two column
       gene ec_number
       (multiple ec_numbers are comma-seperated)
       (no header line)
    '''
    out_f = open(outfile, 'w')
    for species in out_d:
        dd = {}
        for rec in out_d[species]:
            geneid = rec[0]
            ec_list = rec[2]
            if ec_list:
                if geneid in dd:
                    dd[geneid].extend(ec_list)
                else:
                    dd[geneid] = ec_list
        for geneid in sorted(dd.keys()):
            out_f.write('%s\t%s\n' % (geneid,
                                      ','.join(sorted(set(dd[geneid])))))
    out_f.close()


def main(data_folder):
    # if '-b' in sys.argv:
    #     batch_mode = True
    #     sys.argv.remove('-b')
    # else:
    #     batch_mode = False

    out_d = {}
    gbff_files = glob.glob(os.path.join(data_folder, '*.rna.gbff.gz'))
    if len(gbff_files) > 0:
        for infile in gbff_files:
            filename = os.path.split(infile)[1]
            species = filename.split('.')[0]
            print('Parsing "%s" (%s)...' % (filename, species), end='')
            gp = GBFFParser(infile)
            out_li = gp.parse()
            print('Done! [%s]' % len(out_li))
            if species in out_d:
                out_d[species].extend(out_li)
            else:
                out_d[species] = out_li

        #Now dump out_d to disk
        outfile = os.path.join(data_folder, 'rna.gbff.parsed.pyobj')
        dump_object(out_d, outfile)

        #output gene2summary text file
        print('Output "gene2summary_all.txt" file...', end='')
        output_gene2summary(out_d, os.path.join(data_folder,
                                                'gene2summary_all.txt'))
        print("Done!")

        #output gene2ec text file
        print('Output "gene2ec_all.txt" file...', end='')
        output_gene2ec(out_d, os.path.join(data_folder, 'gene2ec_all.txt'))
        print("Done!")
    else:
        print("Error: cannot found any '*.rna.gbff.gz' files in '%s'" % data_folder)

if __name__ == '__main__':
    data_folder = sys.argv[1] if len(sys.argv) > 1 else DATA_FOLDER
    main(data_folder)
