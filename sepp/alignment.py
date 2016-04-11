#!/usr/bin/env python
# -*- coding: utf-8 -*-
###########################################################################
##    Copyright 2012 Siavash Mirarab, Nam Nguyen, and Tandy Warnow.
##    This file is part of SEPP, and some parts of the code are
##    adopted from SATe software package 
##   (https://github.com/sate-dev/sate-core/blob/master/sate/alignment.py).
##
##    SEPP is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    SEPP is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with SEPP.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

import re
#from sepp import #_LOG
from sepp.filemgr import open_with_intermediates

from collections import Mapping
from abc import ABCMeta
import copy
from sepp import get_logger
import pdb
_INDEL = re.compile(r"[-]")
_DANGEROUS_NAME_CHARS = re.compile(r"[^a-zA-Z0-9]")

DATASET_TAXA_ATTR = "taxon_sets"
DATASET_CHAR_ATTR = "char_matrices"

_LOG = get_logger(__name__)

def get_pdistance(distances, leaves, stat = 'mean'):
    """Returns the mean p-distance given a distance matrix and a set of names"""
    counts = 0
    pdistance = 0
    for i in xrange(0,len(leaves)-1):
        for j in xrange(i+1,len(leaves)):
            name = "".join([leaves[i],leaves[j]])
            if stat == 'mean':
                pdistance = pdistance+distances[name]
                counts = counts + 1
            elif stat == 'max':
                if (distance[name] > pdistance):
                    pdistance = distances[name]

    if (stat == 'mean'):
        pdistance = pdistance / counts
    return pdistance
          
def hamming_distance(seq1,seq2):
    """Returns hamming distance between two sequences"""
    #xors = [ord(a) ^ ord(b) for a,b in zip(seq1,seq2)]
    #ors = [ord(a) | ord(b) for a,b in zip(seq1,seq2)]
    #num_mismatch = 0.0
    #non_gap = 0.0
    #for (a,b) in zip(xors,ors):
        #if (a != 0 & a <= 127):
            #num_mismatch+=1
        #if (b != 255):
            #non_gap+=1

    #if non_gap == 0:
        #return -1
    #return num_mismatch/non_gap;
    length = len(seq1)
    num_mismatch = 0.0
    non_gap = 0.0
    for i in xrange(0,length):
        if seq1[i] == "-" or seq2[i] == "-":
            continue
        non_gap+=1
        if seq1[i].lower() !=  seq2[i].lower():
            num_mismatch+=1
    if non_gap == 0:
        return 0
    return num_mismatch/non_gap;


def is_sequence_legal(seq):
    """Check for illegal characters -- TODO:, currently returns True"""
    return True

def _read_fasta(src):
    """generator that returns (name, sequence) tuples from either a FASTA
    formatted file or file object.
    Use MutableAlignment class instead of directly using this method
    """
    file_obj = None
    if isinstance(src, str):
        try:
            file_obj = open(src, "rU")
        except IOError:
            print("The file `%s` does not exist, exiting gracefully" % src)
    elif isinstance(src, file):
            file_obj = src
    else:
        raise TypeError('FASTA reader cannot recognize the source of %s' % src)
    name = None
    seq_list = list()
    for line_number, i in enumerate(file_obj):
        if i.startswith('>'):
            if name:
                yield name, ''.join(seq_list)
                seq_list = list()
            name = i[1:].strip()
        else:
            seq = ''.join(i.strip().upper().split())
            if not is_sequence_legal(seq):
                raise Exception("Error: illegal characeters in sequence at line %d" % line_number)
            seq_list.append(seq)
    yield name, ''.join(seq_list)
    if isinstance(src, str):
        file_obj.close()        

def _write_fasta(alignment, dest):
    """Writes the `alignment` in FASTA format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    for name, seq in alignment.iteritems():
        file_obj.write('>%s\n%s\n' % (name, seq))
    if isinstance(dest, str):
        file_obj.close()

#def is_taxon_a_full_length_seq(x):
#    if isinstance(x, str):
#        return False if "_" in x else True
#    return False if " " in x.label else True

#def get_taxon_label(taxon):
#    return taxon.label.replace(" ", "_")

class ReadOnlyAlignment(Mapping, object):
    '''
    An abstract class that provide all the read-only operations for an alignment
    '''
    __metaclass__ = ABCMeta    
    def get_num_taxa(self):
        """returns the number sequences"""
        return len(self.get_sequence_names())
    
    def get_p_distance(self):
        """returns the average and max p-distance of alignment"""
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        max = 0
        average = 0.0
        count = 0.0        
        copy = MutableAlignment()
        copy.set_alignment(self)
        names = copy.get_sequence_names()
        for name in copy.keys():
            copy[name] = copy[name].replace('-',chr(255)).lower()            
        
        for x in xrange(0,self.get_num_taxa()):
            for y in xrange(x+1,self.get_num_taxa()):
                distance = hamming_distance(copy[names[x]],copy[names[y]])
                if distance >= 0:
                    if distance > max:
                        max = distance
                    count+=1
                    average+=distance
        average = average/count
        return (average,max)
        
    def get_length(self):
        """returns the number of columns in the alignment"""
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            return len(self.values()[0])

    def get_sequence_names(self):
        """returns a list of sequence names"""
        return self.keys()

    def write_to_path(self, filename, schema='FASTA'):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        _LOG.debug("Writing alignment of type %s with length %d to file %s" %(str(self.__class__),len(self),filename))
        file_obj = open_with_intermediates(filename, 'w')
        self.write(file_obj, file_format=schema)
        file_obj.close()

    def write(self, file_obj, file_format):
        """Writes the sequence data in the specified `file_format` to `file_obj`"""
        if (file_format.upper() == 'FASTA'):
            write_func = _write_fasta
#        elif (file_format.upper() == 'NEXUS'):
#            write_func = write_nexus
#        elif (file_format.upper() == 'PHYLIP'):
#            write_func = write_phylip
        else:
            write_func = _write_fasta
        write_func(self, file_obj)

    def write_unaligned_fasta(self, filename):
        """Writes the sequence data without gaps as FASTA, but note that the
        lines may bet "ragged".
        """
        file_obj = open_with_intermediates(filename, 'w')
        for name, seq in self.items():
            file_obj.write('>%s\n%s\n' % (name, re.sub(_INDEL, '', seq)))
        file_obj.close()

    def is_empty(self):
        return self.__len__() < 1

    def is_aligned(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            l0 = self.get_length()
            return all([len(i) == l0 for i in self.values()])                

    def is_all_gap(self, pos):
        '''Checks to see if a column is all gap column at position x'''
        for seq in self.itervalues():
            if (seq[pos] != '-'):
                return False
        return True    
        
    def __str__(self):
        return '\n'.join([">%s\n%s" %(k, self[k]) for k in sorted(self.keys())])
    
    def divide_to_equal_chunks(self,chunks):
        names = self.get_sequence_names()        
        ret = []
        for i in xrange(0,chunks):
            subset = names[i:len(names):chunks]
            if subset:
                subset_alg = self.get_soft_sub_alignment(subset)
            else:
                subset_alg = None
            ret.append(subset_alg)            
        return ret     

    def get_soft_sub_alignment(self, sub_key):
        '''
        Returns a read-only sub-alignment, which won't consume extra memory, 
        since it will not hold a separate copy of the alignment.
        '''
        return ReadonlySubalignment(sub_key, self)       
    
class MutableAlignment(dict, ReadOnlyAlignment, object):
    """ An alignment object, that can be modified. This is the class that should
    be used mainly for holding alignments.  
    """
    def __init__(self):
        "creates an empty matrix"
        dict.__init__(self)
        self.datatype = None
    
    def set_alignment(self, alignment):
        for name, seq in alignment.iteritems():
            self[name] = seq.upper()

    def read_filepath(self, filename, file_format='FASTA'):
        """Augments the matrix by reading the filepath.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        file_obj = open(filename, 'rU')
        return self.read_file_object(file_obj, file_format=file_format)

    def read_file_object(self, file_obj, file_format='FASTA'):
        """Augments the matrix by reading the file object.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        if (file_format.upper() == 'FASTA'):
            read_func = _read_fasta
#        elif (file_format.upper() == 'NEXUS'):
#            read_func = read_nexus
#        elif (file_format.upper() == 'PHYLIP'):
#            read_func = read_phylip
        else:
            raise NotImplementedError("Unknown file format (%s) is not supported" % file_format)
        for name, seq in read_func(file_obj):
            self[name] = seq.upper()
        return self
        
    def add_column(self, pos, char='-'):
        #_LOG.debug("Added a column to reference alignment at position %d" %pos)
        for name, seq in self.iteritems():
            if hasattr(char, "get"):
                c = char.get(name,'-')
            else:
                c = char
            seq = seq[:pos] + c + seq[pos:]
            self[name] = seq
    
    def remove_column(self, pos):
        for name, seq in self.iteritems():            
            seq = seq[:pos] + seq[pos + 1:]
            self[name] = seq

    def remove_columns(self, indexes):
        for name, seq in self.iteritems():
            self[name] = ''.join((char for idx, char in enumerate(seq) if idx not in indexes))
    
    def delete_all_gap(self):
        '''
        Delete all sites that consists of nothing but gaps
        '''
        #pdb.set_trace()
        pos = 0
        i = 0
        subset = []        
        while (pos < self.get_length()):            
            if (self.is_all_gap(pos)):
                self.remove_column(pos)
            else:
                subset.append(i)
                pos = pos + 1
            i+=1
            
        _LOG.debug("Alignment length reduced to %d" %len(subset))
        return subset

    def get_hard_sub_alignment(self, sub_keys):
        "Creates a new alignment with a subset of the taxa."
        new_alignment = MutableAlignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            new_alignment[key] = self[key]
        new_alignment.delete_all_gap()
        return new_alignment         

class ReadonlySubalignment(ReadOnlyAlignment):
    '''
    A class that emulates a subalignment of a given alignment. This is a 
    readonly alignment and does not actually hold sequences in memory. It
    simply keeps a set of sequences that belong to the sub alignment, 
    and implements methods of alignment class (and dictionaries) such that 
    only subalignment keys are returned. 
    '''
    def __init__(self, keys, parent_alignment):
        self.seq_names = set(keys)
        self.parent_alignment = parent_alignment

    def __getitem__(self, key):
        if key in self.seq_names:
            return self.parent_alignment[key]
        else:
            raise KeyError(key)

    def __len__(self):
        return len(self.seq_names)


    def __iter__(self):
        for key in self.seq_names:
            yield key

    def get_mutable_alignment(self):
        ret = MutableAlignment()
        ret.set_alignment(self)
        return ret

class _AlignmentLookupHelper(object):
    '''
    Internal helper calss
    '''
    def __init__(self,pos,ref):
        self.pos = pos
        self.ref = ref        
        
    def get(self, key, default=None):
        if self.ref.has_key(key):
            return self.ref[key][self.pos]
        else:
            return default

class ExtendedAlignment(MutableAlignment):
    '''
    This is used to keep an extended alignment. An extended alignment 
    has a bunch of sequences that are labeles as fragments. More importantly,
    columns of extended alignments are labeled with numbers and also it is 
    known whether a column is an "insertion" column or a normal column. 
    1) these can be read from .sto files. 
    2) these alignments can be merged together.
    '''
    def __init__(self, fragment_names):
        MutableAlignment.__init__(self)
        self.fragments = set(fragment_names)
        self._col_labels = []

    def set_alignment(self, alignment):
        MutableAlignment.set_alignment(self, alignment)
        self._reset_col_names()        

    def read_file_object(self, file_obj, file_format='FASTA'):
        ''' currently supports only fasta'''
        ret = MutableAlignment.read_file_object(self, file_obj, file_format)
        self._reset_col_names()
        return ret

    def add_column(self, pos, char='-', new_label=None):
        '''
        A new column is added to the alignment. The new label can be value, 
        or one of three directives (see below). 
        '''
        MutableAlignment.add_column(self, pos, char)
        if new_label == "MAX":
            self._col_labels.insert(pos, max(self._col_labels) + 1)
        elif new_label == "INC_LAST":
            self._col_labels.add(max(self._col_labels) + 1)
        elif new_label == "RESET":
            self._reset_col_names()
        else:
            self._col_labels.insert(pos, new_label)
        
    def remove_column(self, pos, labels="REMOVE"):
        '''
        Remove a column and potentially adjust column names.
        '''
        MutableAlignment.remove_column(self, pos)
        if labels == "RESET":
            self._reset_col_names()
        elif labels == "REMOVE":
            self._col_labels = self._col_labels[:pos] + self._col_labels[pos + 1:]
        
    def _get_col_labels(self):
        return self._col_labels
    col_labels = property(_get_col_labels)
    
    def _reset_col_names(self):
        ''' sequentially label columns'''
        self._col_labels = range(0,self.get_length())
    
    def get_fragments_readonly_alignment(self):
        '''
        Return a readonly alignment that contains only the fragments. 
        '''        
        return ReadonlySubalignment(self.fragments, self)
    
    def get_base_readonly_alignment(self):
        '''
        Returns a readonly subalignment that does not contain fragments. 
        '''
        return ReadonlySubalignment(self.get_base_seq_names(), self)
    
    def get_fragment_names(self):
        return self.fragments
    
    def get_base_seq_names(self):
        return list(set(self.keys()) - self.fragments)

    def _read_sto(self, handle):
        '''
        Reads a sto file and populates this current object. Figures out
        insertion columns by finding lower case letters, asterisks, and dots. 
        '''
        p = re.compile(r'[a-z*.]')
        insertions = []
        lastStartInd = -1
        for line in handle:
            line = line.strip() 
            if line == "//":
                # End of the alignment. ignore meta-data
                break
            elif line == "" or line[0] == "#":                
                pass
            else: # not a comment
                key, seq = line.split()
                current = self.get(key,"")
                startind = len(current)
                if startind != lastStartInd: # finding insertion columns in limited to first sequence
                    #s = sum(len(x) for x in current)
                    insertions.extend([m.start()+startind for m in p.finditer(seq)])
                    lastStartInd = startind                                            
                self[key] = current + seq.replace(".","-")
#        for k,v in self.items():
#            self[k] = "".join(v)
        self._reset_col_names() 
        return set(insertions)

    def build_extended_alignment(self, base_alignment, path_to_sto_extension,
                                 convert_to_string=True):
        '''
        Given a base alignment and a path to an .sto file (or a list of paths), 
        this methods populates self with an extended alignment by first reading
        the base alignment, and then merging in all the .sto extension alignments.         
        Note that the .sto alignments should contain only fragments, and also
        note that there should be no taxon overlap among extensions, or between
        extensions and the base.
        '''
        if isinstance(base_alignment, ReadOnlyAlignment):
            self.set_alignment(copy.deepcopy(base_alignment))            
        elif isinstance(base_alignment, str):
            _LOG.info("Reading base alignment: %s." %(base_alignment))
            self.read_filepath(base_alignment, "FASTA")
        
        if isinstance(path_to_sto_extension, str):
            paths = [path_to_sto_extension]
        else:
            paths = path_to_sto_extension
        
        self.from_string_to_bytearray()
        for path in paths:
            _LOG.info("Reading sto extended alignment: %s." %(path))
            ext = ExtendedAlignment(self.fragments)            
            ext.read_extended_alignment(path)  
            _LOG.info("Merging extension sto file (%s) into base alignment (%s)." %(path,base_alignment))             
            self.merge_in(ext,False)
            _LOG.debug("Finished merging extension sto file (%s) into base alignment (%s)." %(path,base_alignment))
            del ext
        if convert_to_string:
            self.from_bytearray_to_string() 
            
    
    def read_extended_alignment(self, path, aformat = "stockholm", assertion = False):
        ''' Reads alignment from given path and figures out "insertion"
        columns. Labels insertion columns with special labels and labels the 
        rest of columns (i.e. original columns) sequentially. 
        '''         
        handle = open(path,'rU')
        insertions = None
        if aformat.lower() == "stockholm":
            insertions = self._read_sto(handle)
            _LOG.debug("%s insertions: %d" %(path,len(insertions)))
        else:
            raise ValueError("format %s is not supported yet" %aformat)
        assert insertions is not None        
        
        '''Assert that insertion columns have only gaps in original seqs and
        give them appropriate labels'''
        insertion = -1
        for c in insertions:            
            if assertion:
                assert not any ([self[k][c] != "-" for k in self.get_base_seq_names()]), (
                            "Insertion column has sequence among original "
                            "sequences. An error? column: %d k= %s" %(c,k))
            self.col_labels[c] = insertion
            insertion -= 1
        ''' Adjust labels of other columns '''
        i = 0
        for c in xrange(0,self.get_length()):
            if self.col_labels[c] >= 0:
                self.col_labels[c] = i
                i += 1
            
    def _is_insertion_label(self, i):
        ''' This differentiates between an insertion and an original column        
        '''
        return i < 0
    
    def is_insertion_column(self, col):
        return self._is_insertion_label(self.col_labels[col])
    
    def relabel_original_columns(self, original_labels):
        '''
        This methods relabels non-insertion columns in self based on the 
        input lables. Insertion column labels will not be affected.  
        '''
        j = 0
        _LOG.debug("Relabeling %d (%d) with %d labels." %(self.get_length(),len(self._col_labels),len(original_labels)))
        for i in xrange(0,self.get_length()):
            if not self._is_insertion_label(self.col_labels[i]):
                assert j < len(original_labels), ("Not enough labels"
                           " %d %d.\n %s" %(j,len(original_labels),str(self._col_labels)))
                self.col_labels[i] = original_labels[j]
                j += 1
        assert j == len(original_labels), ("Some of original labels are unused."
                           " Some columns from original alignment went missing? %d %d" %(j,len(original_labels)))    

    def get_insertion_columns(self):	
	return [i for (i,x) in enumerate(self.col_labels) if self._is_insertion_label(x)]

    def get_insertion_column_ranges(self):
        pos=[]
        strt = None
        prev = None
        for p in self.get_insertion_columns():
            if prev is None or prev+1 != p:
                if strt is not None:
                    pos.append((strt,prev))
                strt = p
            prev = p
        if strt is not None:
            pos.append((strt,prev))
        return pos
    
    def write_insertion_column_indexes(self,path):
        file_obj = open(path,'w')
        r = ','.join("%d-%d" %(pair[0],pair[1]) for pair in self.get_insertion_column_ranges())
        file_obj.write(r);
        file_obj.write("\n")
        file_obj.close() 

    def remove_insertion_columns(self):
        '''
        Outputs a new alignment with insertion columns masked out.
        '''
        cols = self.get_insertion_columns()
        s=[]
        a=0
        for b in cols: 
            if b > a:
                s.append((a,b)); 
            a=b+1;
        s.append((a,len(self.col_labels)))
        for name, seq in self.items():
            news = []
            for c in s:
                news.append(seq[c[0]:c[1]])
            self[name] = "".join(news)
        
    def write_insertion_maked_to_file(self,path):
        cols = self.get_insertion_columns()
        s=[]
        a=0
        for b in cols: 
            if b > a:
                s.append((a,b)); 
            a=b+1;
        s.append((a,len(self.col_labels)))
        file_obj = open(path,'w')
        for name, seq in self.items():
            file_obj.write('>%s\n' % name)
            for c in s:
                file_obj.write(seq[c[0]:c[1]])
            file_obj.write("\n")
        file_obj.close()
    
    def from_bytearray_to_string(self):
        for k,v in self.iteritems():
            self[k] = str(v)

    def from_string_to_bytearray(self):
        for k,v in self.iteritems():
            self[k] = bytearray(v)   
                                
    def merge_in(self, other, convert_to_string = True):
        '''
        Merges another alignment in with the current alignment. The other
        alignment needs to be an ExtendedAlignment as well. Since both 
        alignments are extended alignments, we know column labels, and we also
        know which columns are insertions. Columns with the same labels are
        merged together. In other cases gaps are introduced as necessary to
        merge the two alignments. Insertion columns are considered different
        even when they have the same label. 
        
        convert_to_string is by default true, meaning that in the beginning
        sequences of self are assumed to be string objects. These are converted
        to bytearray for better speed, but at the end of the merge, are converted
        back to string. When convert_to_string is False, no conversion back and
        from string is performed, *but* everything is assumed to be in bytearrays.
        So, self is assumed to be in bytearray before calling merge, and will
        remain in bytearray after calling merge. This is useful in cases where
        multiple alignments are merged in with self. Initially, everything
        should be turned into bytearray, and after all merging is finished, 
        everything can be converted back to string (using from_bytearray_to_string
        and from_string_to_bytearray).        
        '''        
        assert isinstance(other, ExtendedAlignment)
        _LOG.debug("Merging started ...")
        if other.is_empty():
            return
        me = 0
        she = 0 # Assumption: alignments are female!
        me_len = self.get_length() if not self.is_empty() else 0 
        she_len = other.get_length()   
        insertion = -1
        
        merged_insertion_columns = 0
        
        ''' Add sequences from her to my alignment '''
        for f in other.fragments:
            self.fragments.add(f)
        if convert_to_string:
            self.from_string_to_bytearray()
        
        selfother={}
        for k,v in other.iteritems():
            #assert k not in self, "Merging overlapping alignments not implemented"
            if (k not in self):
              selfother[k] = bytearray(v)
        while True:
            ''' Check exit conditions'''
            if me == me_len and she == she_len:
                break
            
            ''' Check the 5 possible statuses between she and I '''            
            if she != she_len and other.is_insertion_column(she):
                if me != me_len and self.is_insertion_column(me):
                    ''' We both have a series of insertion columns''' 
                    start = me
                    while me != me_len and self.is_insertion_column(me) and she != she_len and other.is_insertion_column(she):
                        me += 1
                        she += 1
                        merged_insertion_columns += 1
                    run = me - start
                    self.col_labels[start:me] = range(insertion,insertion-run,-1)
                else:
                    ''' Hers is a series of insertion columns'''                
                    start = she
                    while she != she_len and other.is_insertion_column(she):
                        she += 1
                    run = she - start                    
                    ins = bytearray(b"-") * run 
                    for seq in self.itervalues():               
                        seq[me:me] = ins
                    self._col_labels[me:me] = range(insertion,insertion-run,-1)
                    insertion -= run   
                    me += run             
                    me_len += run            
            elif me != me_len and self.is_insertion_column(me):
                ''' Mine is a series of insertion column'''
                start = me
                while me != me_len and self.is_insertion_column(me):
                    me += 1
                run = me - start
                ins = bytearray(b"-") * run
                for v in selfother.itervalues():
                    v[start:start] =ins
                self.col_labels[start:me] = range(insertion,insertion-run,-1)
                insertion -= run
            elif she == she_len or (me != me_len and self.col_labels[me] < other.col_labels[she]):
                ''' My column is not present (i.e. was allgap) in the "other"'''
                start = me
                while me < me_len and (she == she_len or me != me_len and self.col_labels[me] < other.col_labels[she]):
                    me += 1
                run = me - start
                ins = bytearray(b"-") * run
                for v in selfother.itervalues():
                    v[start:start] = ins                         
            elif me == me_len or (she != she_len and self.col_labels[me] > other.col_labels[she]):
                ''' Her column is not present (i.e. was allgap) in "me"'''
                start = she
                while she < she_len and (me == me_len  or she != she_len and self.col_labels[me] > other.col_labels[she]):             
                    she += 1
                run = she - start
                ins = bytearray(b"-") * run 
                for seq in self.itervalues():               
                    seq[me:me] = ins
                self._col_labels[me:me] = other.col_labels[start:she]
                me += run
                me_len += run                                            
            elif self.col_labels[me] == other.col_labels[she]:
                ''' A shared column'''
                while me < me_len and she < she_len and self.col_labels[me] == other.col_labels[she]:
                    she += 1
                    me += 1
            else:
                raise "hmmm, we thought this should be impossible? %d %d" %(me, she)
        
        self.update(selfother)
        
        if convert_to_string:
            self.from_bytearray_to_string()
        _LOG.debug("Merging finished ...")
        
        return merged_insertion_columns
