from .entrez_base import Gene2AccessionParser

__metadata__ = {
    '__collection__': 'entrez_accession',
}


def load_genedoc(self):
    parser = Gene2AccessionParser()
    parser.set_all_species()
    gene2acc = parser.load()
    return gene2acc


def get_mapping(self):
    mapping = {
            "accession": {
                "dynamic": False,
                #"path": "just_name",      #make both fields, accession.rna and rna, work
                "properties": {
                    "genomic": {
                        "type": "string",
			"index": "no",
                        "include_in_all": False,
                        },
                    "rna": {
                        "type": "string",
                        "analyzer": "refseq_analyzer",
                        },
                    'protein': {
                        "type": "string",
                        "analyzer": "refseq_analyzer",
                        #"index_name": "accession",
                        },
                    'translation': {
                        "type": "object",
			"enabled": False,
                        "include_in_all": False,
                        },
                    }
                }
            }
    return mapping
