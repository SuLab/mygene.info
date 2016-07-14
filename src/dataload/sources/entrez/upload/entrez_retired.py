from .entrez_base import Gene2RetiredParser

__metadata__ = {
    '__collection__': 'entrez_retired',
}


def load_genedoc(self):
    parser = Gene2RetiredParser()
    parser.set_all_species()
    gene2retired = parser.load()
    return gene2retired


def get_mapping(self):
    mapping = {
        "retired": {"type": "long"},
    }
    return mapping
