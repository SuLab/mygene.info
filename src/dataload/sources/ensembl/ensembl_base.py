import os.path
import copy
#from config import DATA_ARCHIVE_ROOT
from dataload import get_data_folder
from utils.common import SubStr
from utils.dataload import (load_start, load_done,
                            tab2dict, tab2list, value_convert, normalized_value,
                            list2dict, dict_nodup, dict_attrmerge, tab2dict_iter
                            )

#DATA_FOLDER = os.path.join(DATA_ARCHIVE_ROOT, 'by_resources/ensembl/69')
DATA_FOLDER = get_data_folder('ensembl')
print('DATA_FOLDER: ' + DATA_FOLDER)


#fn to skip lines with LRG records.'''
def _not_LRG(ld):
    return not ld[1].startswith("LRG_")

def map_id(hdocs,mapdict):
    res = []
    for k,v in hdocs.items():
        entrez_ids = mapdict.get(k)
        if entrez_ids:
            for eid in entrez_ids:
                d = {"_id" : eid}
                d.update(v)
                res.append(d)
        else:
            d = {"_id" : k}
            d.update(v)
            res.append(d)

    return res


class EnsemblParser:
    def __init__(self,load_ensembl2entrez=True):
        self.ensembl2entrez_li = None
        self.ensembl_main = None
        if load_ensembl2entrez:
            self._load_ensembl2entrez_li()
            self.ensembl2entrez = list2dict(self.ensembl2entrez_li, 0,alwayslist=True)


    #TODO: not used
    def _load_ensembl_2taxid(self):
        """ensembl2taxid"""
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__translation__main.txt')
        load_start(DATAFILE)
        ensembl2taxid = dict_nodup(tab2dict(DATAFILE, (0, 1), 1, includefn=_not_LRG))
        # need to convert taxid to integer here
        ensembl2taxid = value_convert(ensembl2taxid, lambda x: int(x))
        load_done('[%d]' % len(ensembl2taxid))
        return ensembl2taxid

    #TODO: not used
    def _load_ensembl2name(self):
        """loading ensembl gene to symbol+name mapping"""
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__gene__main.txt')
        load_start(DATAFILE)
        ensembl2name = tab2dict(DATAFILE, (1, 2, 7), 0, includefn=_not_LRG)

        def _fn(x):
            out = {}
            if x[0].strip() not in ['', '\\N']:
                out['symbol'] = x[0].strip()
            if x[1].strip() not in ['', '\\N']:
                _name = SubStr(x[1].strip(), '', ' [Source:').strip()
                if _name:
                    out['name'] = _name
            return out
        ensembl2name = value_convert(ensembl2name, _fn)
        load_done('[%d]' % len(ensembl2name))
        return ensembl2name

    def _load_ensembl2entrez_li(self):
        """gene_ensembl__xref_entrezgene__dm"""
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__xref_entrezgene__dm.txt')
        load_start(DATAFILE)
        ensembl2entrez_li = tab2list(DATAFILE, (1, 2), includefn=_not_LRG)   # [(ensembl_gid, entrez_gid),...]
        load_done('[%d]' % len(ensembl2entrez_li))
        self.ensembl2entrez_li = ensembl2entrez_li


    def load_ensembl_main(self):

        """loading ensembl gene to symbol+name mapping"""
        def _fn(x):
            out = {'taxid' : int(x[0])}
            if x[1].strip() not in ['', '\\N']:
                out['symbol'] = x[1].strip()
            if x[2].strip() not in ['', '\\N']:
                _name = SubStr(x[2].strip(), '', ' [Source:').strip()
                if _name:
                    out['name'] = _name
            return out

        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__gene__main.txt')
        load_start(DATAFILE)
        for datadict in tab2dict_iter(DATAFILE, (0, 1, 2, 7), 1, includefn=_not_LRG):
            datadict = value_convert(datadict, _fn)
            for id,doc in datadict.items():
                doc['_id'] = id
                yield doc

        #ensembl2name = value_convert(ensembl2name, _fn)
        #load_done('[%d]' % len(ensembl2name))
        #return ensembl2name

        #em2name = self._load_ensembl2name()
        #em2taxid = self._load_ensembl_2taxid()
        #assert set(em2name) == set(em2taxid)   # should have the same ensembl ids

        ##merge them together
        #ensembl_main = em2name
        #for k in ensembl_main:
        #    ensembl_main[k].update({'taxid': em2taxid[k]})
        #return ensembl_main

    def load_ensembl2acc(self):
        """
        loading ensembl to transcripts/proteins data
        """
        #Loading all ensembl GeneIDs, TranscriptIDs and ProteinIDs
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__translation__main.txt')
        load_start(DATAFILE)

        def _fn(x, eid):
            out = {'gene': eid, 'mapping' : []}
            def mapping(transcript_id, protein_id):
                trid = transcript_id and transcript_id != '\\N' and transcript_id or None
                pid = protein_id and protein_id != '\\N' and protein_id or None
                out['mapping'].append((trid,pid))

            if isinstance(x, list):
                transcript_li = []
                protein_li = []
                for _x in x:
                    if _x[0] and _x[0] != '\\N':
                        transcript_li.append(_x[0])
                    if _x[1] and _x[1] != '\\N':
                        protein_li.append(_x[1])
                    mapping(_x[0],_x[1])

                if transcript_li:
                    out['transcript'] = normalized_value(transcript_li)
                if protein_li:
                    out['protein'] = normalized_value(protein_li)
            else:
                if x[0] and x[0] != '\\N':
                    out['transcript'] = x[0]
                if x[1] and x[1] != '\\N':
                    out['protein'] = x[1]
                mapping(x[0],x[1])

            return out

        #ensembl2acc = tab2dict(DATAFILE, (1, 2, 3), 0, includefn=_not_LRG)
        for datadict in tab2dict_iter(DATAFILE, (1, 2, 3), 0, includefn=_not_LRG):
            for k in datadict:
                datadict[k] = {'ensembl': _fn(datadict[k], k), '__aslistofdict__' : 'ensembl'}
            for doc in map_id(datadict,self.ensembl2entrez):
                yield doc

        #for k in ensembl2acc:
        #    ensembl2acc[k] = {'ensembl': _fn(ensembl2acc[k], k)}

        #return self.convert2entrez(ensembl2acc)

    def load_ensembl2pos(self):
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__gene__main.txt')
        load_start(DATAFILE)
        for datadict in tab2dict_iter(DATAFILE, (1, 3, 4, 5, 6), 0, includefn=_not_LRG):
            datadict = dict_nodup(datadict)
            datadict = value_convert(datadict, lambda x: {'chr': x[2], 'start': int(x[0]), 'end': int(x[1]), 'strand': int(x[3])})
            datadict = value_convert(datadict, lambda x: {'genomic_pos': x, '__aslistofdict__' : 'genomic_pos'}, traverse_list=False) 
            for doc in map_id(datadict,self.ensembl2entrez):
                yield doc

    def load_ensembl2prosite(self):
        #Prosite
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__prot_profile__dm.txt')
        load_start(DATAFILE)
        for datadict in tab2dict_iter(DATAFILE, (1, 4), 0):
            datadict = dict_nodup(datadict)
            datadict = value_convert(datadict, lambda x: {'prosite': x}, traverse_list=False)
            for doc in map_id(datadict,self.ensembl2entrez):
                yield doc

    def load_ensembl2interpro(self):
        #Interpro
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__prot_interpro__dm.txt')
        load_start(DATAFILE)
        for datadict in tab2dict_iter(DATAFILE, (1, 4, 5, 6), 0):
            datadict = dict_nodup(datadict)
            # optimize with on call/convert
            datadict = value_convert(datadict, lambda x: {'id': x[0], 'short_desc': x[1], 'desc': x[2]})
            # __aslistofdict__ : merge to 'interpro' as list of dict, not merging keys as list
            # (these are merging instructions for later called merge_struct)
            # 'interpro' : {'a': 1, 'b': 2} and 'interpro' : {'a': 3, 'b': 4} should result in:
            # => 'interpro' : [{'a': 1, 'b': 2},{'a': 3, 'b': 4}]
            # or not:
            # => 'interpro' : {'a': [1,3], 'b': [2,4]}
            datadict = value_convert(datadict, lambda x: {'interpro': x, '__aslistofdict__' : 'interpro'}, traverse_list=False)
            for doc in map_id(datadict,self.ensembl2entrez):
                yield doc

    def load_ensembl2pfam(self):
        #Prosite
        DATAFILE = os.path.join(DATA_FOLDER, 'gene_ensembl__prot_pfam__dm.txt')
        load_start(DATAFILE)
        for datadict in tab2dict_iter(DATAFILE, (1, 4), 0):
            datadict = dict_nodup(datadict)
            datadict = value_convert(datadict, lambda x: {'pfam': x}, traverse_list=False)
            for doc in map_id(datadict,self.ensembl2entrez):
                yield doc

    #TODO: not used
    def convert2entrez(self, ensembl2x):
        '''convert a dict with ensembl gene ids as the keys to matching entrezgene ids as the keys.'''
        if not self.ensembl2entrez_li:
            self._load_ensembl2entrez_li()

        if not self.ensembl_main:
            self.ensembl_main = self.load_ensembl_main()

        ensembl2entrez = list2dict(self.ensembl2entrez_li, 0)
        entrez2ensembl = list2dict(self.ensembl2entrez_li, 1)

        #Now make a dictionary indexed by entrez gene id
        print('# of ensembl IDs in total: %d' % len(set(ensembl2x) | set(ensembl2entrez)))
        print('# of ensembl IDs match entrez Gene IDs: %d' % len(set(ensembl2x) & set(ensembl2entrez)))
        print('# of ensembl IDs DO NOT match entrez Gene IDs: %d' % len(set(ensembl2x) - set(ensembl2entrez)))

        #all genes with matched entrez
        def _fn(eid, taxid=None):
            d = copy.copy(ensembl2x.get(eid, {}))   # need to make a copy of the value here.
            return d                                # otherwise, it will cause issue when multiple entrezgene ids
                                                    # match the same ensembl gene, for example,
                                                    #      ENSMUSG00000027104 --> (11909, 100047997)

        data = value_convert(entrez2ensembl, _fn)

        #add those has no matched entrez geneid, using ensembl id as the key
        for eid in (set(ensembl2x) - set(ensembl2entrez)):
            _g = ensembl2x[eid]
            #_g.update(self.ensembl_main.get(eid, {}))
            data[eid] = _g

        for id in data:
            if isinstance(data[id], dict):
                _doc = dict_nodup(data[id], sort=True)
            else:
                #if one entrez gene matches multiple ensembl genes
                _doc = dict_attrmerge(data[id], removedup=True, sort=True)
            data[id] = _doc

        return data
