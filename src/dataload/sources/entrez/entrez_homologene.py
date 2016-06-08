from .entrez_base import HomologeneParser

__metadata__ = {
    '__collection__': 'entrez_homologene',
}


def load_genedoc(self):
    parser = HomologeneParser()
    parser.set_all_species()
    gene2homologene = parser.load()
    return gene2homologene


def get_mapping(self):
    mapping = {
        "homologene": {
            "dynamic": False,
            #"path": "just_name",
            "properties": {
                "genes": {
                    "type": "long",
                    "index": "no",
                    "include_in_all": False,
                },
                "id": {
                    "type": "long",
                    "include_in_all": False,
                    "copy_to": "homologene"
                }
            }
        }
    }
    return mapping
