# Copyright [2010-2012] [Chunlei Wu]
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import print_function
import os
import os.path
import io
import itertools
import csv
import json
from utils.common import safewfile
from biothings.utils.common import ask

csv.field_size_limit(10000000)   # default is 131072, too small for some big files


#===============================================================================
# Misc. Utility functions
#===============================================================================
def load_start(datafile):
    print('Loading "%s"...' % os.path.split(datafile)[1], end='')


def load_done(msg=''):
    print("Done." + msg)


#===============================================================================
# List Utility functions
#===============================================================================
def llist(list, sep='\t'):
    '''Nicely output the list with each item a line.'''
    for x in list:
        if isinstance(x, (list, tuple)):
            xx = sep.join([str(i) for i in x])
        else:
            xx = str(x)
        print(xx)


def listitems(a_list, *idx):
    '''Return multiple items from list by given indexes.'''
    if isinstance(a_list, tuple):
        return tuple([a_list[i] for i in idx])
    else:
        return [a_list[i] for i in idx]


def list2dict(a_list, keyitem, alwayslist=False):
    '''Return a dictionary with specified keyitem as key, others as values.
       keyitem can be an index or a sequence of indexes.
       For example: li=[['A','a',1],
                        ['B','a',2],
                        ['A','b',3]]
                    list2dict(li,0)---> {'A':[('a',1),('b',3)],
                                         'B':('a',2)}
       if alwayslist is True, values are always a list even there is only one item in it.
                    list2dict(li,0,True)---> {'A':[('a',1),('b',3)],
                                              'B':[('a',2),]}
    '''
    _dict = {}
    for x in a_list:
        if isinstance(keyitem, int):      # single item as key
            key = x[keyitem]
            value = tuple(x[:keyitem] + x[keyitem + 1:])
        else:
            key = tuple([x[i] for i in keyitem])
            value = tuple([x[i] for i in range(len(a_list)) if i not in keyitem])
        if len(value) == 1:      # single value
            value = value[0]
        if key not in _dict:
            if alwayslist:
                _dict[key] = [value, ]
            else:
                _dict[key] = value
        else:
            current_value = _dict[key]
            if not isinstance(current_value, list):
                current_value = [current_value, ]
            current_value.append(value)
            _dict[key] = current_value
    return _dict


def list_nondup(a_list):
    x = {}
    for item in a_list:
        x[item] = None
    return x.keys()


def listsort(a_list, by, reverse=False, cmp=None, key=None):
    '''Given list is a list of sub(list/tuple.)
       Return a new list sorted by the ith(given from "by" item)
       item of each sublist.'''
    new_li = [(x[by], x) for x in a_list]
    new_li.sort(cmp=cmp, key=key, reverse=reverse)
    return [x[1] for x in new_li]


def list_itemcnt(a_list):
    '''Return number of occurrence for each type of item in the list.'''
    x = {}
    for item in a_list:
        if item in x:
            x[item] += 1
        else:
            x[item] = 1
    return [(i, x[i]) for i in x]


def alwayslist(value):
    """If input value if not a list/tuple type, return it as a single value list."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return value
    else:
        return [value]


def merge_struct(v1, v2):

    if isinstance(v1, list):
        #print("v1 is list ", end="")
        if isinstance(v2, list):
            #print("v2 is list -> extend")
            v1.extend(v2)
        else:
            #print("v2 not list -> append")
            v1.append(v2)

    elif isinstance(v1, dict):
        assert isinstance(v2, dict)
        #print("v1 & v2 is dict")
        for k in list(v1.keys()):
            #print("v1[%s]" % k)
            if k in v2:
                #print("%s is both dict -> merge/update v1[%s],v2[%s]" % (k, k, k))
                #v1.update(merge(v1[k], v2[k]))
                v1[k] = merge_struct(v1[k], v2[k])
            else:
                #print("%s not in v2 -> update v1 with v2" % k)
                v1.update(v2)
        for k in v2:
            #print("v2[%s]" % k)
            if k in v1:
                pass  # already done
            else:
                v1[k] = v2[k]

    elif isinstance(v1, str) or isinstance(v1, int) or isinstance(v1, float):
        if isinstance(v2, str) or isinstance(v2, int) or isinstance(v2, float):
            if v1 != v2:
                #print("v1 & v2 not iterable -> list of 2")
                v1 = [v1, v2]
            else:
                pass
                #print("v1 == v2, skip")
        else:
            #print("v2 iterable, reverse merge")
            return merge_struct(v2, v1)
    else:
        raise TypeError("dunno how to merge type %s" % type(v1))

    #print("return %s" % v1)
    return v1

#===============================================================================
# File Utility functions
#===============================================================================
def anyfile(infile, mode='r'):
    '''
    return a file handler with the support for gzip/zip comppressed files
    if infile is a two value tuple, then first one is the compressed file;
      the second one is the actual filename in the compressed file.
      e.g., ('a.zip', 'aa.txt')

    '''
    if isinstance(infile, tuple):
        infile, rawfile = infile[:2]
    else:
        rawfile = os.path.splitext(infile)[0]
    filetype = os.path.splitext(infile)[1].lower()
    if filetype == '.gz':
        import gzip
        in_f = io.TextIOWrapper(gzip.GzipFile(infile, 'r'))
    elif filetype == '.zip':
        import zipfile
        in_f = io.TextIOWrapper(zipfile.ZipFile(infile, 'r').open(rawfile, 'r'))
    else:
        in_f = open(infile, mode)
    return in_f


def tabfile_tester(datafile, header=1, sep='\t'):
    reader = csv.reader(anyfile(datafile), delimiter=sep)
    lineno = 0
    try:
        for i in range(header):
            next(reader)
            lineno += 1

        for ld in reader:
            lineno += 1
    except:
        print("Error at line number:", lineno)
        raise


def dupline_seperator(dupline, dup_sep, dup_idx=None, strip=False):
    '''
    for a line like this:
        a   b1,b2  c1,c2

    return a generator of this list (breaking out of the duplicates in each field):
        [(a,b1,c1),
         (a,b2,c1),
         (a,b1,c2),
         (a,b2,c2)]
    example:
         dupline_seperator(dupline=['a', 'b1,b2', 'c1,c2'],
                           dup_idx=[1,2],
                           dup_sep=',')
    if dup_idx is None, try to split on every field.
    if strip is True, also tripe out of extra spaces.
    '''
    value_li = list(dupline)
    for idx, value in enumerate(value_li):
        if dup_idx:
            if idx in dup_idx:
                value = value.split(dup_sep)
                if strip:
                    value = [x.strip() for x in value]
            else:
                value = [value]
        else:
            value = value.split(dup_sep)
            if strip:
                value = [x.strip() for x in value]
        value_li[idx] = value
    return itertools.product(*value_li)    # itertools.product fits exactly the purpose here


def tabfile_feeder(datafile, header=1, sep='\t',
                   includefn=None,
                   coerce_unicode=True,
                   assert_column_no=None):
    '''a generator for each row in the file.'''

    in_f = anyfile(datafile)
    reader = csv.reader(in_f, delimiter=sep)
    lineno = 0
    try:
        for i in range(header):
            next(reader)
            lineno += 1

        for ld in reader:
            if assert_column_no:
                if len(ld) != assert_column_no:
                    err = "Unexpected column number:" \
                          " got {}, should be {}".format(len(ld), assert_column_no)
                    raise ValueError(err)
            if not includefn or includefn(ld):
                lineno += 1
                if coerce_unicode:
                    yield [str(x) for x in ld]
                else:
                    yield ld
    except ValueError:
        print("Error at line number:", lineno)
        raise


def tab2list(datafile, cols, **kwargs):
    if os.path.exists(datafile):
        if isinstance(cols, int):
            return [ld[cols] for ld in tabfile_feeder(datafile, **kwargs)]
        else:
            return [listitems(ld, *cols) for ld in tabfile_feeder(datafile, **kwargs)]
    else:
        print('Error: missing "%s". Skipped!' % os.path.split(datafile)[1])
        return {}


def tab2dict(datafile, cols, key, alwayslist=False, **kwargs):
    if isinstance(datafile, tuple):
        _datafile = datafile[0]
    else:
        _datafile = datafile
    if os.path.exists(_datafile):
        return list2dict([listitems(ld, *cols) for ld in tabfile_feeder(datafile, **kwargs)], key, alwayslist=alwayslist)
    else:
        print('Error: missing "%s". Skipped!' % os.path.split(_datafile)[1])
        return {}


def tab2dict_iter(datafile, cols, key, alwayslist=False, **kwargs):
    if isinstance(datafile, tuple):
        _datafile = datafile[0]
    else:
        _datafile = datafile
    if os.path.exists(_datafile):
        bulk = []
        prev_id = None
        for ld in tabfile_feeder(datafile, **kwargs):
            #print(ld)
            li = listitems(ld, *cols)
            #print("key %s len bulk %s prev %s" % (li[key],len(bulk),prev_id))
            if prev_id == None or (li[key] == prev_id):
                #print("\t\tfound same")
                bulk.append(li)
                prev_id = li[key]
            else:
                #print("bulk size: %s" % bulk)
                di = list2dict(bulk, key, alwayslist=alwayslist)
                bulk = []
                # changed key, init next bulk
                bulk.append(li)
                prev_id = li[key]
                #print("on yield: %s" % di)
                yield di
        # flush remaining bulk
        if bulk:
            di = list2dict(bulk, key, alwayslist=alwayslist)
            yield di
    else:
        print('Error: missing "%s". Skipped!' % os.path.split(_datafile)[1])
        return {}


def file_merge(infiles, outfile=None, header=1, verbose=1):
    '''merge a list of input files with the same format.
       if header will be removed from the 2nd files in the list.
    '''
    outfile = outfile or '_merged'.join(os.path.splitext(infiles[0]))
    out_f, outfile = safewfile(outfile)
    if verbose:
        print("Merging...")
    cnt = 0
    for i, fn in enumerate(infiles):
        print(os.path.split(fn)[1], '...', end='')
        line_no = 0
        in_f = anyfile(fn)
        if i > 0:
            for k in range(header):
                in_f.readline()
        for line in in_f:
            out_f.write(line)
            line_no += 1
        in_f.close()
        cnt += line_no
        print(line_no)
    out_f.close()
    print("=" * 20)
    print("Done![total %d lines output]" % cnt)


#===============================================================================
# Dictionary Utility functions
#===============================================================================
def value_convert(_dict, fn, traverse_list=True):
    '''For each value in _dict, apply fn and then update
       _dict with return the value.
       if traverse_list is True and a value is a list,
       apply fn to each item of the list.
    '''
    for k in _dict:
        if traverse_list and isinstance(_dict[k], list):
            _dict[k] = [fn(x) for x in _dict[k]]
        else:
            _dict[k] = fn(_dict[k])
    return _dict


def dict_convert(_dict, keyfn=None, valuefn=None):
    '''Return a new dict with each key converted by keyfn (if not None),
       and each value converted by valuefn (if not None).
    '''
    if keyfn is None and valuefn is not None:
        for k in _dict:
            _dict[k] = valuefn(_dict[k])
        return _dict

    elif keyfn is not None:
        out_dict = {}
        for k in _dict:
            out_dict[keyfn(k)] = valuefn(_dict[k]) if valuefn else _dict[k]
        return out_dict
    else:
        return _dict


def updated_dict(_dict, attrs):
    '''Same as dict.update, but return the updated dictionary.'''
    out = _dict.copy()
    out.update(attrs)
    return out


def merge_dict(dict_li, attr_li, missingvalue=None):
    '''
    Merging multiple dictionaries into a new one.
    Example:
    In [136]: d1 = {'id1': 100, 'id2': 200}
    In [137]: d2 = {'id1': 'aaa', 'id2': 'bbb', 'id3': 'ccc'}
    In [138]: merge_dict([d1,d2], ['number', 'string'])
    Out[138]:
    {'id1': {'number': 100, 'string': 'aaa'},
     'id2': {'number': 200, 'string': 'bbb'},
     'id3': {'string': 'ccc'}}
    In [139]: merge_dict([d1,d2], ['number', 'string'], missingvalue='NA')
    Out[139]:
    {'id1': {'number': 100, 'string': 'aaa'},
     'id2': {'number': 200, 'string': 'bbb'},
     'id3': {'number': 'NA', 'string': 'ccc'}}
    '''
    dd = dict(zip(attr_li, dict_li))
    key_set = set()
    for attr in dd:
        key_set = key_set | set(dd[attr])

    out_dict = {}
    for k in key_set:
        value = {}
        for attr in dd:
            if k in dd[attr]:
                value[attr] = dd[attr][k]
            elif missingvalue is not None:
                value[attr] = missingvalue
        out_dict[k] = value
    return out_dict


def normalized_value(value, sort=True):
    '''Return a "normalized" value:
           1. if a list, remove duplicate and sort it
           2. if a list with one item, convert to that single item only
           3. if a list, remove empty values
           4. otherwise, return value as it is.
    '''
    if isinstance(value, list):
        value = [x for x in value if x]   # remove empty values
        try:
            _v = list(set(value))
        except TypeError:
            #use alternative way
            _v = [json.loads(x) for x in set([json.dumps(x) for x in value])]
        if len(_v) and sort:
            # py3 won't sort dict anymore...
            if isinstance(_v[0],dict):
                _v = sorted(_v,key=lambda x: sorted(x.keys()))
            else:
                _v = sorted(_v)
        if len(_v) == 1:
            _v = _v[0]
    else:
        _v = value

    return _v


def dict_nodup(_dict, sort=True):
    for k in _dict:
        _dict[k] = normalized_value(_dict[k], sort=sort)
    return _dict


def dict_attrmerge(dict_li, removedup=True, sort=True, special_fns={}):
    '''
        dict_attrmerge([{'a': 1, 'b':[2,3]},
                        {'a': [1,2], 'b':[3,5], 'c'=4}])
        sould return
             {'a': [1,2], 'b':[2,3,5], 'c'=4}

        special_fns is a dictionary of {attr:  merge_fn}
         used for some special attr, which need special merge_fn
         e.g.,   {'uniprot': _merge_uniprot}
    '''
    out_dict = {}
    keys = []
    for d in dict_li:
        keys.extend(d.keys())
    keys = set(keys)
    for k in keys:
        _value = []
        for d in dict_li:
            if d.get(k, None):
                if isinstance(d[k], list):
                    _value.extend(d[k])
                else:
                    _value.append(d[k])
        if len(_value) == 1:
            out_dict[k] = _value[0]
        else:
            out_dict[k] = _value

        if k in special_fns:
            out_dict[k] = special_fns[k](out_dict[k])

    if removedup:
        out_dict = dict_nodup(out_dict, sort=sort)
    return out_dict


def dict_apply(dict, key, value, sort=True):
    '''

    '''
    if key in dict:
        _value = dict[key]
        if not isinstance(_value, list):
            _value = [_value]
        if isinstance(value, list):
            _value.extend(value)
        else:
            _value.append(value)
    else:
        _value = value

    dict[key] = normalized_value(_value, sort=sort)


def dict_to_list(gene_d):
    '''return a list of genedoc from genedoc dictionary and
       make sure the "_id" field exists.
    '''
    doc_li = [updated_dict(gene_d[k], {'_id': str(k)}) for k in sorted(gene_d.keys())]
    return doc_li


#===============================================================================
# Network Utility functions
#===============================================================================
def download(url, output_folder, output_file, no_confirm=False, use_axel=False):
    orig_path = os.getcwd()
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  # create output_folder if doesn not exist
    try:
        os.chdir(output_folder)
        if os.path.exists(output_file):
            if no_confirm or ask('Remove existing file "%s"?' % output_file) == 'Y':
                os.remove(output_file)
            else:
                print("Skipped!")
                return
        print('Downloading "%s"...' % output_file)
        if use_axel:
            #faster than wget using 5 connections
            cmdline = 'axel -a -n 5 "{}" -o "{}"'.format(url, output_file)
        else:
            cmdline = 'wget "{}" -O "{}"'.format(url, output_file)
        return_code = os.system(cmdline)
        if return_code == 0:
            print("Success.")
        else:
            print("Failed with return code (%s)." % return_code)
        print("="*50)
    finally:
        os.chdir(orig_path)
