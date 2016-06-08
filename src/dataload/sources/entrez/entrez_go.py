from .entrez_base import Gene2GOParser

__metadata__ = {
    '__collection__': 'entrez_go',
}


def load_genedoc(self=None):
    parser = Gene2GOParser()
    parser.set_all_species()
    gene2go = parser.load()
    return gene2go


def get_mapping(self=None):
    mapping = {
        "go": {
            "dynamic": False,
            "properties": {
                "MF": {
                    "dynamic": False,
                    "properties": {
                        "term": {
                            "type": "string",
                            "include_in_all": False
                        },
                        "id": {
                            "type": "string",
                            "analyzer": "string_lowercase",
                        },
                        "evidence": {
                            "type": "string",
                            "index": "no"
                        },
                        "pubmed": {
                            "type": "long",
                            "index": "no"
                        }
                    }
                },
                "CC": {
                    "dynamic": False,
                    "properties": {
                        "term": {
                            "type": "string",
                            "include_in_all": False
                        },
                        "id": {
                            "type": "string",
                            "analyzer": "string_lowercase",
                        },
                        "evidence": {
                            "type": "string",
                            "index": "no"
                        },
                        "pubmed": {
                            "type": "long",
                            "index": "no"
                        }
                    }
                },
                "BP": {
                    "dynamic": False,
                    "properties": {
                        "term": {
                            "type": "string",
                            "include_in_all": False
                        },
                        "id": {
                            "type": "string",
                            "analyzer": "string_lowercase",
                        },
                        "evidence": {
                            "type": "string",
                            "index": "no"
                        },
                        "pubmed": {
                            "type": "long",
                            "index": "no"
                        }
                    }
                }
            }
        }
    }

    return mapping
