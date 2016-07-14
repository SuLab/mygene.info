'''
Populates MICROBE gene entries with genomic position data
Currently updates the 120 microbial taxids that are NCBI Reference Sequences

run get_ref_microbe_taxids function to get an updated file for TAXIDS_FILE
when it's necessary.
'''
import os.path
from biothings.utils.common import (dump, loadobj, get_timestamp)
from utils.dataload import (tab2list, load_start, load_done)
from dataload import get_data_folder


DATA_FOLDER = get_data_folder('entrez')
print('DATA_FOLDER: ' + DATA_FOLDER)

__metadata__ = {
    '__collection__': 'entrez_genomic_pos',
}

TAXIDS_FILE = os.path.join(DATA_FOLDER, "../ref_microbe_taxids_20151014.pyobj")
DATAFILE = os.path.join(DATA_FOLDER, 'gene/gene2refseq.gz')


def load_genedoc(self):
    """
    Loads gene data from NCBI's refseq2gene.gz file.
    Parses it based on genomic position data and refseq status provided by the
    list of taxids from get_ref_microbe_taxids() as lookup table
    :return:
    """
    taxids = loadobj(TAXIDS_FILE)
    taxid_set = set(taxids)
    load_start(DATAFILE)

    def _includefn(ld):
        return ld[0] in taxid_set  # match taxid from taxid_set

    cols_included = [0, 1, 7, 9, 10, 11]  # 0-based col idx
    gene2genomic_pos_li = tab2list(DATAFILE, cols_included, header=1,
                                   includefn=_includefn)
    count = 0
    last_id = None
    for gene in gene2genomic_pos_li:
        count += 1
        strand = 1 if gene[5] == '+' else -1
        _id = gene[1]

        mgi_dict = {
            '_id': _id,
            'genomic_pos': {
                'start': int(gene[3]),
                'end': int(gene[4]),
                'chr': gene[2],
                'strand': strand
            }
        }
        if _id != last_id:
            # rows with dup _id will be skipped
            yield mgi_dict
        last_id = _id

    load_done('[%d]' % count)


def get_mapping(self):
    mapping = {
        "genomic_pos": {
            "dynamic": False,
            "type": "nested",
            "properties": {
                "chr": {"type": "string"},
                "start": {"type": "long"},
                "end": {"type": "long"},
                "strand": {
                    "type": "byte",
                    "index": "no"
                },
            },
        },
    }

    return mapping


def get_ref_microbe_taxids():
    """
    Downloads the latest bacterial genome assembly summary from the NCBI genome
    ftp site and generate a list of taxids of the bacterial reference genomes.

    :return:
    """
    import urllib
    import csv

    urlbase = 'ftp://ftp.ncbi.nlm.nih.gov'
    urlextension = '/genomes/refseq/bacteria/assembly_summary.txt'
    assembly = urllib.urlopen(urlbase + urlextension)
    datareader = csv.reader(assembly.read().splitlines(), delimiter="\t")
    taxid = []

    for row in datareader:
        if row[4] == 'reference genome':
            taxid.append(row[5])

    ts = get_timestamp()
    dump(taxid, "ref_microbe_taxids_{}.pyobj".format(ts))

    return taxid
