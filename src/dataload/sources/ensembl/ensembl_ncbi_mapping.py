from __future__ import print_function
import os.path
from collections import defaultdict

import mygene

from dataload import get_data_folder
from utils.dataload import anyfile, tabfile_feeder
from utils.common import safewfile


ENSEMBL_DATA_FOLDER = get_data_folder('ensembl')
print('Ensembl DATA_FOLDER: ' + ENSEMBL_DATA_FOLDER)
Entrez_DATA_FOLDER = get_data_folder('entrez')
print('Ensembl DATA_FOLDER: ' + Entrez_DATA_FOLDER)


gene_ensembl_1_xref_dm_file = os.path.join(ENSEMBL_DATA_FOLDER, "gene_ensembl__xref_entrezgene__dm.txt")
gene_ensembl_2_main_file = os.path.join(ENSEMBL_DATA_FOLDER, "gene_ensembl__gene__main.txt")
gene2ensembl_file = os.path.join(Entrez_DATA_FOLDER, "gene/gene2ensembl.gz")
gene_main_file = os.path.join(Entrez_DATA_FOLDER, "gene/gene_info.gz")

outfile = os.path.join(ENSEMBL_DATA_FOLDER, "gene_ensembl__gene__extra.txt")


def find_multiple_mappings_from_entrezgene_file(gene_ensembl_entrezgene_dm_file):
    """Input gene_ensembl_entrezgene_dm_file, and identify how many NCBI gene IDs there are for
    each ensembl gene ID. Lines in input file are:

    'gene_ensembl__xref_entrezgene__dm.txt' (useful columns in input_file):

    col1: Ensembl gene ID
    col2: NCBI gene ID

    If there is > 1 NCBI gene ID, we need to process further.
    """
    print("ensembl_ncbi_mapping.py script is running now...")
    print("step 1 start: find where multiple NCBI IDs map to one Ensembl ID")
    ensembl_dict_with_entrez = defaultdict(list)

    with open(gene_ensembl_entrezgene_dm_file) as file_in:
        next(file_in)
        for line in file_in:
            split_line = line.split("\t")
            ensembl_gene_id = split_line[1].strip()
            ncbi_gene_id = split_line[2].strip()
            ensembl_dict_with_entrez[ensembl_gene_id].append(ncbi_gene_id)

    print("number of Ensembl gene IDs mapped to NCBI gene ID: ", len(ensembl_dict_with_entrez))
    multi_mapping_dict = {}
    for key in ensembl_dict_with_entrez:
        if len(ensembl_dict_with_entrez[key]) > 1:
            multi_mapping_dict[key] = ensembl_dict_with_entrez[key]
    print("number of Ensembl IDs with > 1 NCBI gene ID: ", len(multi_mapping_dict))
    print("step 1 end")
    return multi_mapping_dict, len(ensembl_dict_with_entrez)


def create_ensembl_gene_id_dict(gene_ensembl_main_file, multi_mapping_dict):
    """Using gene_ensembl_main_file, identify correct ensembl symbol for each
    ensembl gene ID. Add this information to a new dictionary.

    'gene_ensembl__gene__main.txt' (useful columns in input file):

    col1: Ensembl gene ID
    col2: Ensembl symbol
    """
    print("step 2 start: get Ensembl symbol from Ensembl main file")
    ensembl_dict = defaultdict(list)
    with open(gene_ensembl_main_file) as file_in:
        next(file_in)
        for line in file_in:
            split_line = line.split("\t")
            ensembl_id_from_main = str(split_line[1].strip())
            ensembl_symbol_from_main = split_line[2].strip()

            if ensembl_id_from_main in multi_mapping_dict:
                ensembl_id_dict = {}
                ensembl_id_dict['data'] = {'ncbi_list': multi_mapping_dict[ensembl_id_from_main],
                                           'symbol': ensembl_symbol_from_main,
                                           'gene2ensembl': []}
                ensembl_dict[ensembl_id_from_main] = ensembl_id_dict
    print("step 2 end")
    return ensembl_dict


def find_ncbi_ids_from_gene2ensembl(ensembl_dict, gene2ensembl_file):
    """Input is gene2ensembl_file; maps NCBI gene ID to one Ensembl gene ID.

    'gene2ensembl' (useful columns in input file):

    col1: NCBI gene ID
    col2: Ensembl gene ID
    """
    print("step 3 start: find NCBI IDs from gene2ensembl file")
    file_in = anyfile(gene2ensembl_file)
    #with open(gene2ensembl_file) as file_in:
    try:
        next(file_in)
        for line in file_in:
            split_line = line.split("\t")
            ensembl_gene_id_from_gene2ensembl = split_line[2].strip()
            ncbi_gene_id_from_gene2ensembl = split_line[1].strip()

            if ensembl_gene_id_from_gene2ensembl in ensembl_dict:
                ensembl_dict[ensembl_gene_id_from_gene2ensembl]['data']['gene2ensembl'].append(ncbi_gene_id_from_gene2ensembl)
    finally:
        file_in.close()

    count = 0
    for key in ensembl_dict:
        if len(ensembl_dict[key]['data']['gene2ensembl']) == 1:
            count += 1

    print("Total number of Ensembl gene IDs mapping uniquely with gene2ensembl: ", count)
    print("step 3 end")
    return ensembl_dict, count


def find_ncbi_symbols(gene_info_file,ensembl_dict):
    print("step 4 start: read NCBI gene symbol")
    ncbi_list_to_find = {}
    for key in ensembl_dict:
        ncbi_list = ensembl_dict[key]['data']['ncbi_list']
        for e in ncbi_list:
            ncbi_list_to_find[e] = True
        #gene2ensembl_ncbi_gene_id_match_list = ensembl_dict[key]['data']['gene2ensembl']
        #if len(gene2ensembl_ncbi_gene_id_match_list) != 1:
        #    ncbi_list_to_find.append(ncbi_list)

    for e in list(set([item for sublist in ncbi_list_to_find for item in sublist])):
        ncbi_list_to_find[e] = True

    ncbi_id_symbols = {}
    print("todiufd %s" % repr(ncbi_list_to_find))
    for ld in tabfile_feeder(gene_info_file):
        if ld[1] in ncbi_list_to_find:
            ncbi_id_symbols[ld[1]] = ld[2]

    print("number of unique NCBI gene IDs to be queried using mygene.info: ", len(ncbi_list_to_find))
    print("number symbols found from querying mygene.info: ", len(ncbi_id_symbols))
    print("step 4 end")
    return ncbi_id_symbols


def merge_mapping(ensembl_dict, mygene_website_dict, add_source=False):
    """First use gene2ensembl as single match NCBI gene ID (if == 1 match).
    Next, if no gene2ensembl match, then look at mygene.info to find which NCBI
    ID from the NCBI multi mapping list returns the same ensembl symbol as the
    ensembl main file, and use corresponding NCBI gene ID as single match.

    OUTPUT generator:
    ---------------------
    Tuple with ensembl gene ID and NCBI gene ID
    """
    print("step 5 start: Generator-decide whether to use gene2ensembl or mygene.info for mapping")
    for key in ensembl_dict:
        ncbi_list = ensembl_dict[key]['data']['ncbi_list']
        ensembl_symbol = ensembl_dict[key]['data']['symbol'].upper()
        gene2ensembl_ncbi_gene_id_match_list = ensembl_dict[key]['data']['gene2ensembl']

        if len(gene2ensembl_ncbi_gene_id_match_list) == 1:
            if add_source is False:
                yield (key, int(gene2ensembl_ncbi_gene_id_match_list[0]))
            else:
                yield (key, int(gene2ensembl_ncbi_gene_id_match_list[0]), '1')

        else:
            ensembl_symbol_list_from_mygene = []
            for ncbi_id in ncbi_list:
                try:
                    ensembl_symbol_list_from_mygene.append(mygene_website_dict[ncbi_id].upper())
                except KeyError:
                    # need this here; keeps list size/order will never match with ensembl_symbol)
                    ensembl_symbol_list_from_mygene.append('symbol_not_found')

            if ensembl_symbol in ensembl_symbol_list_from_mygene:
                if ensembl_symbol_list_from_mygene.count(ensembl_symbol) == 1:
                    ncbi_idx = ensembl_symbol_list_from_mygene.index(ensembl_symbol)
                    if add_source is False:
                        yield (key, int(ncbi_list[ncbi_idx]))
                    else:
                        yield (key, int(ncbi_list[ncbi_idx]), '2')

    print("step 5 end")


def write_mapping_file(mapping_generator, confirm=True):
    """OUTPUT is mapping file:
    -------------------------
    Note: you will not know the source of the mapping unless you use
    the optional parameter "add_source=True" to merge_mapping() function
    col0: Ensembl gene ID
    col2 "add_source" == 1: NCBI ID gene ID from gene2ensembl
    col2 "add_source" == 2: NCBI ID gene ID from ncbi_list if mygene.info symbol == ensembl symbol
        (i.e. iterate through ncbi list (for each Ensembl ID) on mygene.info
        (ex: http://mygene.info/v2/gene/100894237?fields=symbol )
        and when the symbol found matches the ensembl symbol use this
        NCBI ID if symbols match only once)
    """
    print("step 6 start: write file from mapping generator of tuples")
    mapping_file, mapping_filename = safewfile(outfile, prompt=confirm,default='O')

    count = 0
    for item in mapping_generator:
        count += 1
        split_item = list(item)
        split_item = '\t'.join([str(i) for i in split_item])
        mapping_file.write(split_item + "\n")

    print("total Ensembl IDs uniquely mapped to NCBI gene ID:", count)
    mapping_file.close()
    print("Output file: \"{}\"".format(mapping_filename))
    print("step 6 end\n")
    return count


def run_stats(total_ensembl_IDs, ensembl_dict, ensembl_map_count, total_mapped):
    print("Final Summary:")
    print("--------------")
    print("Total Ensembl gene IDs", total_ensembl_IDs)
    print("Total Ensembl gene IDs with multiple NCBI gene IDs: ", len(ensembl_dict))
    print("Percent of Ensembl gene IDs with multiple NCBI gene IDs: ", round((len(ensembl_dict) * 1. / (total_ensembl_IDs)) * 100, 1))
    print("Total Ensembl gene IDs successfully and uniquely mapped to 1 NCBI gene ID: ", total_mapped)
    print("Total mapped using gene2ensembl: ", ensembl_map_count)
    print("Total mapped from mygene.info: ", total_mapped - ensembl_map_count)
    print("Percent of Ensembl IDs uniquely mapped out of Ensembl IDs with > 1 NCBI gene ID: ", round((total_mapped * 1. / (len(ensembl_dict))) * 100, 1))


# def main(gene_ensembl_1, gene_ensembl_2, gene2ensembl):
def main(confirm=True):
    multi_mapping_dict, total_ensembl_IDs = find_multiple_mappings_from_entrezgene_file(gene_ensembl_1_xref_dm_file)
    ensembl_dict = create_ensembl_gene_id_dict(gene_ensembl_2_main_file, multi_mapping_dict)
    ensembl_dict, ensembl_match_count = find_ncbi_ids_from_gene2ensembl(ensembl_dict, gene2ensembl_file)
    ncbi_id_symbols = find_ncbi_symbols(gene_main_file, ensembl_dict)
    mapping_generator = merge_mapping(ensembl_dict, ncbi_id_symbols, add_source=False)
    total_mapped = write_mapping_file(mapping_generator, confirm=confirm)
    run_stats(total_ensembl_IDs, ensembl_dict, ensembl_match_count, total_mapped)

# use this for testing purposes:
main(False)#gene_ensembl_1_xref_dm_file, gene_ensembl_2_main_file, gene2ensembl_file)
