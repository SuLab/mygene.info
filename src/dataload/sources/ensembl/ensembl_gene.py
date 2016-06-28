from .ensembl_base import EnsemblParser


__metadata__ = {
    '__collection__': 'ensembl_gene',
    'required_fields': ['taxid'],
    'ENSEMBL_GENEDOC_ROOT': True,
    'id_type': 'ensembl_gene',

}


def load_genedoc(self=None):
    ep = EnsemblParser(load_ensembl2entrez=False)
    ensembl_genes = ep.load_ensembl_main()
    return ensembl_genes


def get_mapping(self=None):
    mapping = {
        "taxid":  {"type": "integer",
                   "include_in_all": False},
        "symbol": {"type": "string",
                   "analyzer": "string_lowercase",
                   "boost": 5.0},
        "name":   {"type": "string",
                   "boost": 0.8},    # downgrade name field a little bit

    }
    return mapping


def get_mapping_to_entrez(self=None):
    ep = EnsemblParser()
    ep._load_ensembl2entrez_li()
    return ep.ensembl2entrez_li
