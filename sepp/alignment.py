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

from collections import Mapping, Container
from abc import ABCMeta
import copy
from sepp import get_logger

_INDEL = re.compile(r"[-]")
_DANGEROUS_NAME_CHARS = re.compile(r"[^a-zA-Z0-9]")

DATASET_TAXA_ATTR = "taxon_sets"
DATASET_CHAR_ATTR = "char_matrices"

_LOG = get_logger(__name__)


def is_sequence_legal(seq):
    """Check for illegal characters -- TODO, currently returns True"""
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
    for name, seq in alignment.items():
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
        for name in self.keys():
            if (self[name][pos] != '-'):
                return False
        return True    
        
    def __str__(self):
        return '\n'.join([">%s\n%s" %(k, self[k]) for k in sorted(self.keys())])     

class MutableAlignment(dict, ReadOnlyAlignment, object):
    """ An alignment object, that can be modified. This is the class that should
    be used mainly for holding alignments.  
    """
    def __init__(self):
        "creates an empty matrix"
        dict.__init__(self)
        self.datatype = None
    
    def set_alignment(self, alignment):
        for name in alignment.keys():
            self[name] = alignment[name].upper()

    def read_filepath(self, filename, file_format='FASTA'):
        """Augments the matrix by reading the filepath.
        If duplicate sequence names are encountered then the old name will be replaced.
        """
        file_obj = open(filename, 'r')
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
        for name in self.keys():
            seq = self.get(name)
            if hasattr(char, "get"):
                c = char.get(name,'-')
            else:
                c = char
            seq = seq[:pos] + c + seq[pos:]
            self[name] = seq
    
    def remove_column(self, pos):
        for name in self.keys():
            seq = self.get(name)
            seq = seq[:pos] + seq[pos + 1:]
            self[name] = seq
    
    def delete_all_gap(self):
        '''
        Delete all sites that consists of nothing but gaps
        '''
        pos = 0
        i = 0
        subset = []
        name = self.keys()[0]
        while (pos < len(self[name])):            
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

    def get_soft_sub_alignment(self, sub_key):
        '''
        Returns a read-only sub-alignment, which won't consume extra memory, 
        since it will not hold a separate copy of the alignment.
        '''
        return ReadonlySubalignment(sub_key, self)                

class ReadonlySubalignment(ReadOnlyAlignment):
    '''
    A class that emulates a subalignment of a given alignment. This is a 
    readonly alignment and does not actually hold sequences in memory. It
    simply keeps the list of sequences that belong to the sub alignment, 
    and implements methods of alignment class (and dictionaries) such that 
    only subalignment keys are returned. 
    '''
    def __init__(self, keys, parent_alignment):
        self.seq_names = keys
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
        insertion columns by finding lower case letters and dots. 
        '''
        insertions = set()
        for line in handle:
            line = line.strip() 
            if line == '# STOCKHOLM 1.0': # If multiple alignments, just treat them as one. 
                pass
            elif line == "//":
                # End of the alignment. ignore meta-data
                break
            elif line == "":                
                pass
            elif line[0] != "#": # not a comment
                parts = [x.strip() for x in line.split(" ",1)]
                if len(parts) != 2:
                    raise ValueError("Could not split line into identifier " \
                                      + "and sequence:\n" + line)
                key, seq = parts
                if key not in self.keys():
                    self[key] = ""
                startind = len(self[key])
                for x in [m.start()+startind for m in re.finditer('[a-z.]', seq)]: 
                    insertions.add(x)
                
                self[key] += seq.replace(".","-")
        self._reset_col_names() 
        return insertions

    def build_extended_alignment(self, base_alignment, path_to_sto_extension):
        '''
        Given a base alignment, and the path to and .sto file, this 
        populates self with an extended alignment. This is equivalent of
        reading the base_alignment firs, and then merging in the extended 
        alignment. 
        '''
        if isinstance(base_alignment, ReadOnlyAlignment):
            self.set_alignment(copy.deepcopy(base_alignment))            
        elif isinstance(base_alignment, str):
            self.read_filepath(base_alignment, "FASTA")
        
        ext = ExtendedAlignment(self.fragments)
        ext.read_extended_alignment(path_to_sto_extension)
        
        self.merge_in(ext)
    
    def read_extended_alignment(self, path, aformat = "stockholm"):
        ''' Reads alignment from given path and figures out "insertion"
        columns. Labels insertion columns with special labels and labels the 
        rest of columns (i.e. original columns) sequentially. 
        '''         
        handle = open(path,'r')
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
            k=""
            assert all([self[k][c] != "-" for k in self.get_base_seq_names()]), (
                            "Insertion column has sequence among original"
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
        for i in xrange(0,self.get_length()):
            if not self._is_insertion_label(self.col_labels[i]):
                self.col_labels[i] = original_labels[j]
                j += 1
        assert j == len(original_labels), ("Some of original labels are unused."
                           " Some columns from original alignment went missing? %d %d" %(j,len(original_labels)))    
        
    def merge_in(self, other):
        '''
        Merges another alignment in with the current alignment. The other
        alignment needs to be an ExtendedAlignment as well. Since both 
        alignments are extended alignments, we know column labels, and we also
        know which columns are insertions. Columns with the same labels are
        merged together. In other cases gaps are introduced as necessary to
        merge the two alignments.
        '''
        assert isinstance(other, ExtendedAlignment)
        me = 0
        she = 0 # important alignments are assumed to be female!
        me_len = self.get_length() if not self.is_empty() else 0 
        she_len = other.get_length()   
        insertion = -1
        
        for f in other.fragments:
            self.fragments.add(f)
        for k in other.keys():
            assert k not in self.keys(), "Merging overlapping alignments not implemented"
            self[k] = ""
            
        while me < me_len or she < she_len:
            #print me, she, me_len, she_len
            if me == me_len and she == she_len:
                break
            if she != she_len and other.is_insertion_column(she):
                self.add_column(me, char = _AlignmentLookupHelper(she, other))                
                self.col_labels[me] = insertion
                insertion -= 1
                she += 1
                me += 1
                me_len += 1
            elif me != me_len and self.is_insertion_column(me):
                for k in other.keys():
                    self[k] += "-"
                self.col_labels[me] = insertion
                insertion -= 1                
                me += 1
            elif she == she_len or (me != me_len and self.col_labels[me] < other.col_labels[she]):
                ''' My column is not present (i.e. was allgap) in the "other"'''
                for k in other.keys():
                    self[k] += "-"
                me += 1                
            elif me == me_len or (she != she_len and self.col_labels[me] > other.col_labels[she]):
                ''' Her column is not present (i.e. was allgap) in "me"'''
                self.add_column(me, char= _AlignmentLookupHelper(she, other))
                self.col_labels[me] = other.col_labels[she]
                she += 1
                me += 1
                me_len += 1
            elif self.col_labels[me] == other.col_labels[she]:
                ''' A shared column'''                
                for k in other.keys():                    
                    self[k] += other[k][she]
                she += 1
                me += 1
            else:
                raise "hmmm, we thought this should be impossible? %d %d" %(me, she)
