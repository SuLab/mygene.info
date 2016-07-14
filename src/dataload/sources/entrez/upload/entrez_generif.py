from .entrez_base import Gene2GeneRifParser

__metadata__ = {
    '__collection__': 'entrez_generif',
}


def load_genedoc(self=None):
    gene2generif = Gene2GeneRifParser().load()
    return gene2generif


def get_mapping(self=None):
    mapping = {
        # do not index generif
        "generif": {
            "properties": {
                "pubmed": {
                    "type": "long",
                    "index": "no",
                },
                "text": {
                    "type": "string",
                    "index": "no",
                }
            }
        }
    }
    return mapping
