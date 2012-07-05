#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple classes for reading and manipulating sequence data matrices
"""
# This file is part of SATe

# SATe is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Jiaye Yu and Mark Holder, University of Kansas

import re, sys
#from sepp import #_LOG
from sepp.filemgr import open_with_intermediates

_INDEL = re.compile(r"[-]")
_DANGEROUS_NAME_CHARS = re.compile(r"[^a-zA-Z0-9]")

DATASET_TAXA_ATTR = "taxon_sets"
DATASET_CHAR_ATTR = "char_matrices"

def is_sequence_legal(seq):
    """Check for illegal characters -- TODO, currently returns True"""
    return True

def read_fasta(src):
    """generator that returns (name, sequence) tuples from either a FASTA
    formatted file or file object.
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

def read_nexus(src):
    "TODO: use dendropy to do this."
    raise NotImplementedError('Input of NEXUS file format is not supported yet.')

def read_phylip(src):
    "TODO: use dendropy to do this."
    raise NotImplementedError('Input of PHYLIP file format is not supported yet.')

def write_fasta(alignment, dest):
    """Writes the `alignment` in FASTA format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    for name, seq in alignment.items():
        file_obj.write('>%s\n%s\n' % (name, seq) )
    if isinstance(dest, str):
        file_obj.close()

def write_phylip(alignment, dest):
    """Writes the `alignment` in relaxed PHYLIP format to either a file object or file"""
    file_obj = None
    if isinstance(dest, str):
        file_obj = open(dest, "w")
    else:
        file_obj = dest
    names = alignment.get_sequence_names()
    assert(names)
    ntax = len(names)
    seq = alignment[names[0]]
    nchar = len(seq)
    file_obj.write('%s\t%s\n' % (ntax, nchar) )
    for k in names:
        assert len(k.split()) == 1
        seq = alignment[k]
        assert len(seq) == nchar
        file_obj.write('%s\n%s\n' % (k, seq))
    if isinstance(dest, str):
        file_obj.close()

def write_nexus(alignment, file_obj):
    "TODO use dendropy"
    raise NotImplementedError('Output of NEXUS file format is not supported yet.')

def is_taxon_a_full_length_seq(x):
    if isinstance(x,str):
        return False if "_" in x else True
    return False if " " in x.label else True

def get_taxon_label(taxon):
    return taxon.label.replace(" ","_")

class Alignment(dict, object):
    """A simple class that maps taxa names to sequences.
    TODO: switch to dendropy character_matrix
    """
    def __init__(self):
        "creates an empty matrix"
        dict.__init__(self)
        self.datatype = None

    def get_datatype(self):
        return self._datatype

    def set_datatype(self, d):
        if d is None:
            self._datatype = None
        else:
            self._datatype = d.upper()

    datatype = property(get_datatype, set_datatype)
    
    def set_alignment(self, alignment):
	for name in alignment.keys():
	    self[name] = alignment[name].upper()

    def get_sequence_names(self):
        "returns a list of sequence names"
        return self.keys()

    def get_num_taxa(self):
        "returns the number sequences"
        return len(self.get_sequence_names())

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
        if ( file_format.upper() == 'FASTA' ):
            read_func = read_fasta
        elif ( file_format.upper() == 'NEXUS' ):
            read_func = read_nexus
        elif ( file_format.upper() == 'PHYLIP' ):
            read_func = read_phylip
        else:
            raise NotImplementedError("Unknown file format (%s) is not supported" % file_format)
        for name, seq in read_func(file_obj):
            self[name] = seq.upper()

    def write_to_path(self, filename, schema='FASTA'):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        file_obj = open_with_intermediates(filename,'w')
        self.write(file_obj, file_format=schema)
        file_obj.close()

    def write(self, file_obj, file_format):
        """Writes the sequence data in the specified `file_format` to `file_obj`"""
        if ( file_format.upper() == 'FASTA' ):
            write_func = write_fasta
        elif ( file_format.upper() == 'NEXUS' ):
            write_func = write_nexus
        elif ( file_format.upper() == 'PHYLIP' ):
            write_func = write_phylip
        else:
            write_func = write_fasta
        write_func(self, file_obj)

    def write_unaligned_fasta(self, filename):
        """Writes the sequence data without gaps as FASTA, but note that the
        lines may bet "ragged".
        """
        file_obj = open_with_intermediates(filename, 'w')
        for name, seq in self.items():
            file_obj.write('>%s\n%s\n' % (name, re.sub(_INDEL, '', seq)))
        file_obj.close()

    def sub_alignment(self, sub_keys):
        "Creates an new alignment with a subset of the taxa."
        new_alignment = Alignment()
        new_alignment.datatype = self.datatype
        for key in sub_keys:
            new_alignment[key] = self[key]
        return new_alignment

    def is_empty(self):
        return self.__len__() < 1

    def is_aligned(self):
        if self.is_empty():
            raise ValueError("The alignment is empty.\n")
        else:
            return all( [len(i)==len(self.values()[0]) for i in self.values()] )

    def partition_info(self, base=0):
        return (self.datatype, 1+base, self.sequence_length()+base)

    def sequence_length(self):
        if self.is_aligned():
            return len(self.values()[0])
        
    def add_column(self, pos, char='-'):
        #_LOG.debug("Added a column to reference alignment at position %d" %pos)
        for name in self.keys():
            seq = self.get(name)
            seq = seq[:pos] + char + seq[pos:]
            self[name] = seq


    def add_gap_at(self, short_reads, read_i):
        for seq in short_reads.keys():
            read = short_reads.get(seq)
            read = read[:read_i] + '-' + read[read_i:]
            short_reads[seq] = read            

    def merge_alignment_in(self, to_merge_alignment):
        short_reads = {}
        piv_seq = None
        failed = True
        for row in to_merge_alignment.items():
            if (not row[0] in self):
                short_reads[row[0]] = row[1] 
            elif piv_seq is None:
		if (row[1].replace('-','') == self[row[0]].replace('-','')):	      
		    piv_name, piv_seq = row
		    failed = False
		else:
		    sys.stderr.write("Sequence %s does not match up, look for another pivot\n" % row[0])		    
		    
            #else:
		#if (row[1].replace('-','') !=self[row[0]].replace('-','')):
		      #sys.stderr.write("ERROR " + str(row[0]) + "\n")
		      #failed = True
        if (failed):
	      sys.stderr.write("Failed to to find pivot, quitting\n")
	      exit()
        #_LOG.debug("Selected pivot is:\n%s" % piv_seq)
        ## The pivot is added to the short reads to enable
        ## easily testing the correctness of results
        short_reads[piv_name] = piv_seq
        new_seq = self[piv_name]                         
        
        pivot_i = 0 #Idx of character from original alignment
        read_i = 0  #current position in to merge in piv sequence
        for new_i, ch in enumerate(new_seq):
            # This case happens if reference alignment has - at the end
            # but pivot lacks them 
            if pivot_i >= len(piv_seq):
                if ch == '-':		    
		    #print "Changed 1 %s %s" % (ch, str(len(short_reads[piv_name])))
                    self.add_gap_at(short_reads, read_i)
                else:
                    raise RuntimeError("unalignable:\n%s\n%s" %(new_seq,piv_seq))
            # Pivot and reference are aligned. everything is fine. move on!
            elif (ch == piv_seq[pivot_i]):
		#print "Aligned %s %s" % (piv_seq[pivot_i], str(len(short_reads[piv_name])))
                pivot_i += 1
            # Reference has a gap that pivot lacks. Add a gap to short reads!
            elif ch == '-':			
		#print "Changed 2 %s %s" % (piv_seq[pivot_i], str(len(short_reads[piv_name])))
                self.add_gap_at(short_reads, read_i)
            # Pivot has a gap that reference lacks it. 
            # Add all such gaps to the reference alignment! 
            elif piv_seq[pivot_i] == '-':
                pc  = piv_seq[pivot_i]
                while (ch != piv_seq[pivot_i]):
		    #print "Changed me"
                    self.add_column(read_i)
                    pivot_i += 1
                    read_i += 1
                    pc = piv_seq[pivot_i]
                pivot_i += 1
            else:
                    raise "error"
            read_i += 1
            #print "(%s %s %s %s)" % (str(new_i), str(ch), str(pivot_i), str(read_i))

	#Take care of last case where adding gaps made short reads longer than original alignment	
	while (len(short_reads[piv_name]) > len(self[piv_name])):
	      self.add_column(len(short_reads[piv_name]))

        # Test that pivot is now aligned to reference
        #_LOG.debug("Ref: %s" %self[piv_name])
        #_LOG.debug("Piv: %s" %short_reads[piv_name])                    
        assert self[piv_name] == short_reads[piv_name]
        
        # Add short reads to the reference
        for seq, read in short_reads.items():
            self[seq]=read
