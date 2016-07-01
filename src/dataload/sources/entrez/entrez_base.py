from __future__ import print_function
import os.path
from config import SPECIES_LI, TAXONOMY
from utils.common import file_newer, loadobj, dump
from dataload import get_data_folder
from utils.dataload import (load_start, load_done,
                            tab2dict, tab2list, value_convert,
                            normalized_value, dict_convert, dict_to_list,
                            tab2dict_iter
                            )

DATA_FOLDER = get_data_folder('entrez')
print('DATA_FOLDER: ' + DATA_FOLDER)


class EntrezParserBase(object):
    def __init__(self):
        # if species_li is None, include all species
        self.set_species_li(SPECIES_LI)
        self.DATA_FOLDER = DATA_FOLDER
        self.datafile = os.path.join(self.DATA_FOLDER, self.DATAFILE)

    def set_all_species(self):
        '''To load all species.'''
        self.species_li = None
        self.taxid_set = None
        self.species_filter = None

    def set_species_li(self, species_li):
        '''To load only specied species if species_li is not None.'''
        if species_li:
            self.species_li = species_li
            self.taxid_set = set([TAXONOMY[species] for species in species_li])
            self.species_filter = lambda ld: int(ld[0]) in self.taxid_set
        else:
            self.set_all_species()

    def load(self, aslist=False):
        raise NotImplementedError


class GeneInfoParser(EntrezParserBase):
    '''Parser for NCBI gene_info.gz file.'''
    DATAFILE = 'gene/gene_info.gz'

    def format(self,doc):
        gid, info = list(doc.items())[0]
        info['entrezgene'] = int(gid)
        info["_id"] = gid
        return info

    def load(self, aslist=False):
        '''
        loading ncbi "gene_info" file
        This must be called first to create basic gene documents
        with all basic fields, e.g., name, symbol, synonyms, etc.

        format of gene_info file:
        #Format: tax_id GeneID Symbol LocusTag Synonyms dbXrefs
                 map_location description type_of_gene Symbol_from
                 nomenclature_authority Full_name_from_nomenclature_authority
        Nomenclature_status Other_designations Modification_da
        te (tab is used as a separator, pound sign - start of a comment)

        '''
        load_start(self.datafile)
        gene_d = tab2dict(self.datafile, (0, 1, 2, 3, 4, 5, 7, 8, 9), key=1,
                          alwayslist=0, includefn=self.species_filter)

        def _ff(d):
            (
                taxid, symbol, locus_tag, synonyms,
                dbxrefs, map_location,
                description, type_of_gene
            ) = d
            out = dict(taxid=int(taxid),
                       symbol=symbol,
                       name=description)
            if map_location != '-':
                out['map_location'] = map_location
            if type_of_gene != '-':
                out['type_of_gene'] = type_of_gene
            if synonyms != '-':
                out['alias'] = normalized_value(synonyms.split('|'))
            if locus_tag != '-':
                out['locus_tag'] = locus_tag

            for x in dbxrefs.split('|'):
                if x == '-':
                    continue
                xd = x.split(':')
                if len(xd) == 3 and xd[0] == xd[1] and \
                        xd[0] in ['HGNC', 'MGI']:
                    # a fix for NCBI bug for dup xref prefix, 'HGNC:HGNC:36328'
                    xd = xd[1:]
                try:
                    _db, _id = xd
                except:
                    print(x)
                    raise
                # we don't need ensembl xref from here, we will get it from
                # Ensembl directly
                if _db.lower() in ['ensembl', 'imgt/gene-db']:
                    # we don't need 'IMGT/GENE-DB" xref either, because they
                    # are mostly the same as gene symbol
                    continue
                # add "MGI:" prefix for MGI ids.
                if _db.lower() == 'mgi':
                    _id = "MGI:" + _id
                out[_db] = _id
            return out


        # add entrezgene field
        cnt = 0
        for d in gene_d:
            d = value_convert(d, _ff)
            yield self.format(d)
            cnt += 1

        load_done('[%d]' % cnt)


def get_geneid_d(species_li=None, load_cache=True, save_cache=True,only_for={}):
    '''return a dictionary of current/retired geneid to current geneid mapping.
       This is useful, when other annotations were mapped to geneids may
       contain retired gene ids.

       if species_li is None, genes from all species are loaded.

       Note that all ids are int type.
    '''
    if species_li:
        taxid_set = set([TAXONOMY[species] for species in species_li])
    else:
        taxid_set = None

    orig_cwd = os.getcwd()
    os.chdir(DATA_FOLDER)

    # check cache file
    _cache_file = 'gene/geneid_d.pyobj'
    if load_cache and os.path.exists(_cache_file) and \
       file_newer(_cache_file, 'gene/gene_info.gz') and \
       file_newer(_cache_file, 'gene/gene_history.gz'):

        print('Loading "geneid_d" from cache file...', end='')
        _taxid_set, out_d = loadobj(_cache_file)
        assert _taxid_set == taxid_set
        print('Done.')
        os.chdir(orig_cwd)
        return out_d

    DATAFILE = os.path.join(DATA_FOLDER, 'gene/gene_info.gz')
    load_start(DATAFILE)
    if species_li:
        species_filter = lambda ld: int(ld[0]) in taxid_set and (only_for and ld[1] in only_for)
    else:
        species_filter = lambda ld: only_for and ld[1] in only_for
    geneid_li = set(tab2list(DATAFILE, 1, includefn=species_filter))
    load_done('[%d]' % len(geneid_li))

    DATAFILE = os.path.join(DATA_FOLDER, 'gene/gene_history.gz')
    load_start(DATAFILE)

    if species_li:
        _includefn = lambda ld: int(ld[0]) in taxid_set and ld[1] in geneid_li
    else:
        _includefn = lambda ld: ld[1] in geneid_li    # include all species
    retired2gene = tab2dict(DATAFILE, (1, 2), 1, alwayslist=0, includefn=_includefn)
    # includefn above makes sure taxid is for species_li and filters out those
    # mapped_to geneid exists in gene_info list

    load_done('[%d]' % len(retired2gene))
    # convert key/value to int
    out_d = dict_convert(retired2gene, keyfn=int, valuefn=int)
    for g in geneid_li:
        _g = int(g)
        out_d[_g] = _g

    if save_cache:
        if species_li:
            dump((taxid_set, out_d), _cache_file)
        else:
            dump((None, out_d), _cache_file)

    os.chdir(orig_cwd)
    return out_d


class HomologeneParser(EntrezParserBase):
    '''Parser for NCBI homologenes.data file.'''
    DATAFILE = 'Homologene/homologene.data'

    def _sorted_homologenes(self, homologenes):
        '''sort list of homologenes [(taxid, geneid),...] based on the order
            defined in species_li.
        '''
        d = {}
        for i, species in enumerate(SPECIES_LI):
            d[TAXONOMY[species]] = i
        gene_li = [(d.get(taxid, taxid), taxid, geneid)
                   for taxid, geneid in homologenes]
        return [g[1:] for g in sorted(gene_li)]

    def load(self, aslist=False):
        '''
        loading ncbi "homologene.data" file
        adding "homologene" field in gene doc
        '''

        load_start(self.datafile)
        homo_d = tab2dict(self.datafile,(2,1),0,header=0)
        DATAFILE = os.path.join(DATA_FOLDER, 'gene/gene_history.gz')
        retired2gene = tab2dict(DATAFILE, (1, 2), 1, alwayslist=0,includefn=lambda ld: ld[1] != '-')
        for id in list(homo_d.keys()):
            homo_d[retired2gene.get(id,id)] = homo_d[id]

        with open(self.datafile) as df:
            homologene_d = {}
            doc_li = []
            print()
            geneid_d = get_geneid_d(self.species_li,load_cache=False,save_cache=False,only_for=homo_d)

            for line in df:
                ld = line.strip().split('\t')
                hm_id, tax_id, geneid = [int(x) for x in ld[:3]]
                if (self.taxid_set is None or tax_id in self.taxid_set) and geneid in geneid_d:
                    # for selected species only
                    # and also ignore those geneid does not match any
                    # existing gene doc
                    # in case of orignal geneid is retired, replaced with the
                    # new one, if available.
                    geneid = geneid_d[geneid]
                    genes = homologene_d.get(hm_id, [])
                    genes.append((tax_id, geneid))
                    homologene_d[hm_id] = genes

                    doc_li.append(dict(_id=str(geneid), taxid=tax_id,
                                       homologene={'id': hm_id}))

            for i, gdoc in enumerate(doc_li):
                gdoc['homologene']['genes'] = self._sorted_homologenes(
                    set(homologene_d[gdoc['homologene']['id']]))
                doc_li[i] = gdoc

            load_done('[%d]' % len(doc_li))

        if aslist:
            return doc_li
        else:
            gene_d = dict([(d['_id'], d) for d in doc_li])
            return gene_d


class GeneSummaryParser(EntrezParserBase):
    '''Parser for gene2summary_all.txt, adding "summary" field in gene doc'''
    DATAFILE = 'refseq/gene2summary_all.txt'

    def load(self, aslist=False):
        load_start(self.datafile)
        with open(self.datafile) as df:
            geneid_set = set()
            doc_li = []
            for line in df:
                geneid, summary = line.strip().split('\t')
                if geneid not in geneid_set:
                    doc_li.append(dict(_id=geneid, summary=str(summary)))
                    geneid_set.add(geneid)
        load_done('[%d]' % len(doc_li))

        if aslist:
            return doc_li
        else:
            gene_d = dict([(d['_id'], d) for d in doc_li])
            return gene_d


class Gene2AccessionParserBase(EntrezParserBase):
    DATAFILE = 'to_be_specified'
    fieldname = 'to_be_specified'

    def format(self,doc):
        gid, acc = list(doc.items())[0]
        d = {"_id" : gid}
        d.update(acc)
        return d

    def load(self, aslist=False):
        load_start(self.datafile)
        gene2acc = tab2dict_iter(self.datafile, (1, 3, 5, 7), 0, alwayslist=1,
                             includefn=self.species_filter)

        def _ff(d):
            out = {
                'rna': [],
                'protein': [],
                'genomic': [],
                'mapping': []
            }
            for rna, prot, dna in d:
                if rna == '-': rna = None
                if prot == '-': prot = None
                if dna == '-': dna = None
                if rna is not None:
                    out['rna'].append(rna)
                if prot is not None:
                    out['protein'].append(prot)
                if dna is not None:
                    out['genomic'].append(dna)
                if rna is not None or prot is not None:
                    out['mapping'].append((rna,prot))
            # remove dup
            for k in out:
                out[k] = normalized_value(out[k])
            # remove empty rna/protein/genomic field
            _out = {}
            for k, v in out.items():
                if v:
                    _out[k] = v
            if _out:
                _out = {self.fieldname: _out}
            return _out

        #gene/2acc = dict_convert(gene2acc, valuefn=_ff)
        cnt = 0
        for gd in gene2acc:
            convd = self.format(dict_convert(gd, valuefn=_ff))
            yield convd
            cnt += 1
        load_done('[%d]' % cnt)

        if aslist:
            return dict_to_list(gene2acc)
        else:
            return gene2acc


class Gene2AccessionParser(Gene2AccessionParserBase):
    DATAFILE = 'gene/gene2accession.gz'
    fieldname = 'accession'


class Gene2RefseqParser(Gene2AccessionParserBase):
    DATAFILE = 'gene/gene2refseq.gz'
    fieldname = 'refseq'


class Gene2UnigeneParser(EntrezParserBase):
    DATAFILE = 'gene/gene2unigene'

    def format(self,doc):
        gid, unigene = list(doc.items())[0]
        return {"_id" : gid, "unigene" : unigene} 

    def load(self, aslist=False):
        load_start(self.datafile)
        print()
        uni_d = tab2dict(self.datafile, (0, 1), 0, alwayslist=0)
        DATAFILE = os.path.join(DATA_FOLDER, 'gene/gene_history.gz')
        retired2gene = tab2dict(DATAFILE, (1, 2), 1, alwayslist=0,includefn=lambda ld: ld[1] != '-')
        for id in list(uni_d.keys()):
            uni_d[retired2gene.get(id,id)] = uni_d[id]
        geneid_d = get_geneid_d(self.species_li,load_cache=False,save_cache=False,only_for=uni_d)
        #geneid_d = get_geneid_d(self.species_li)
        gene2unigene = tab2dict_iter(self.datafile, (0, 1), 0, alwayslist=0,
                                 includefn=lambda ld: int(ld[0]) in geneid_d)
        cnt = 0
        for doc in gene2unigene:
            yield self.format(doc)
            cnt += 1
        load_done('[%d]' % cnt)



class Gene2GOParser(EntrezParserBase):
    DATAFILE = 'gene/gene2go.gz'

    def load(self, aslist=False):
        load_start(self.datafile)
        gene2go = tab2dict_iter(self.datafile, (1, 2, 3, 4, 5, 6, 7), 0, alwayslist=1,
                           includefn=self.species_filter)
        category_d = {'Function': 'MF',
                      'Process': 'BP',
                      'Component': 'CC'}

        def _ff(d):
            out = {}
            for goid, evidence, qualifier, goterm, pubmed, gocategory in d:
                _gocategory = category_d[gocategory]
                _d = out.get(_gocategory, [])
                _rec = dict(id=goid, term=goterm)
                if evidence != '-':
                    _rec['evidence'] = evidence
                if qualifier != '-':
                    # here I also fixing some inconsistency issues in NCBI data
                    # Colocalizes_with -> colocalizes_with
                    # Contributes_with -> contributes_with
                    # Not -> NOT
                    _rec['qualifier'] = qualifier.replace('Co', 'co').replace('Not', 'NOT')
                if pubmed != '-':
                    if pubmed.find('|') != -1:
                        pubmed = [int(pid) for pid in pubmed.split('|')]
                    else:
                        pubmed = int(pubmed)
                    _rec['pubmed'] = pubmed
                _d.append(_rec)
                out[_gocategory] = _d
            for k in out:
                if len(out[k]) == 1:
                    out[k] = out[k][0]
            return out

        for gd in gene2go:
            convd = dict_convert(gd, valuefn=_ff)
            assert len(list(convd.items())) == 1, "nope: %s" % list(convd.items())
            gid, go = list(convd.items())[0]
            gene_d = {}
            gene_d[gid] = {'go': go}
            yield gene_d


class Gene2RetiredParser(EntrezParserBase):
    '''
    loading ncbi gene_history file, adding "retired" field in gene doc
    '''

    DATAFILE = 'gene/gene_history.gz'

    def load(self, aslist=False):
        load_start(self.datafile)
        if self.species_li:
            _includefn = lambda ld: int(ld[0]) in self.taxid_set and ld[1] != '-'
        else:
            _includefn = lambda ld: ld[1] != '-'
        gene2retired = tab2dict(self.datafile, (1, 2), 0, alwayslist=1,
                                includefn=_includefn)
        gene2retired = dict_convert(gene2retired, valuefn=lambda x: normalized_value([int(xx) for xx in x]))

        gene_d = {}
        for gid, retired in gene2retired.items():
            gene_d[gid] = {'retired': retired}
        load_done('[%d]' % len(gene_d))

        if aslist:
            return dict_to_list(gene_d)
        else:
            return gene_d


class Gene2ECParser(EntrezParserBase):
    '''
    loading gene2ec data, adding "ec" field in gene doc

    Sample lines for input file:
        24159   2.3.3.8
        24161   3.1.3.2,3.1.3.48
    '''
    DATAFILE = 'refseq/gene2ec_all.txt'

    def load(self, aslist=False):
        load_start(self.datafile)
        with open(self.datafile) as df:
            geneid_set = set()
            doc_li = []
            for line in df:
                geneid, ec = line.strip().split('\t')
                if ec.find(',') != -1:
                    # there are multiple EC numbers
                    ec = [str(x) for x in ec.split(',')]
                else:
                    ec = str(ec)
                if geneid not in geneid_set:
                    doc_li.append(dict(_id=geneid, ec=ec))
                    geneid_set.add(geneid)
        load_done('[%d]' % len(doc_li))

        if aslist:
            return doc_li
        else:
            gene_d = dict([(d['_id'], d) for d in doc_li])
            return gene_d


class Gene2GeneRifParser(EntrezParserBase):
    '''
    '''
    DATAFILE = 'generif/generifs_basic.gz'

    def _cvt_pubmed(self, pubmed_str):
        """input is a string of one or multiple pubmed ids, sep by comma"""
        _li = [int(x) for x in pubmed_str.split(',')]
        if len(_li) == 1:
            return _li[0]
        else:
            return _li

    def load(self):
        load_start(self.datafile)
        cnt = 0
        for datadict in tab2dict_iter(self.datafile, (1, 2, 4), 0, alwayslist=1):
            datadict = dict_convert(datadict, valuefn=lambda v: {
                            'generif': [dict(pubmed=self._cvt_pubmed(x[0]), text=x[1]) for x in v]})

            for id,doc in datadict.items():
                cnt += 1
                doc['_id'] = id
                yield doc
        #gene2generif = tab2dict(self.datafile, (1, 2, 4), 0, alwayslist=1)
        #gene2generif = dict_convert(gene2generif, valuefn=lambda v: {
        #    'generif': [dict(pubmed=self._cvt_pubmed(x[0]), text=x[1]) for x in v]})
        load_done('[%d]' % cnt)
