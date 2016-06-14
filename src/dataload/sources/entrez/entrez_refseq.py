from .entrez_base import Gene2RefseqParser


__metadata__ = {
    '__collection__': 'entrez_refseq',
}


def load_genedoc(self):
    parser = Gene2RefseqParser()
    parser.set_all_species()
    gene2refseq = parser.load()
    return gene2refseq


def get_mapping(self):
    mapping = {
        "refseq": {
            "dynamic": False,
            #"path": "just_name",      #make both fields, refseq, refseq.rna, work
            "properties": {
                "genomic": {
                    "type": "string",
                    "index": "no",
                    "include_in_all": False,
                },
                "rna": {
                    "type": "string",
                    "analyzer": "refseq_analyzer",
                    "copy_to": "refseq",
                },
                'protein': {
                    "type": "string",
                    "analyzer": "refseq_analyzer",
                    "copy_to": "refseq",
                }
            }
        }
    }
    return mapping
