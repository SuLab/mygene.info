from .entrez_base import Gene2UnigeneParser


__metadata__ = {
    '__collection__': 'entrez_unigene',
}


def load_genedoc(self):
    self.parser = Gene2UnigeneParser()
    self.parser.set_all_species()
    gene2unigene = self.parser.load()
    return gene2unigene


def get_mapping(self):
    mapping = {
        "unigene":  {
            "type": "string",
            "analyzer": "string_lowercase"
        }
    }
    return mapping
