#!/usr/bin/env python
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

import re
import sys
from copy import deepcopy
from satelib import get_logger, log_exception
from satelib.filemgr import open_with_intermediates

_LOG = get_logger(__name__)
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
            self[name] = seq

    def write_filepath(self, filename, file_format='FASTA'):
        """Writes the sequence data in the specified `file_format` to `filename`"""
        file_obj = open_with_intermediates(filename,'w')
        self.write(file_obj, file_format=file_format)
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

def concatenate_alignments(alignments):
    combined_alignment = deepcopy(alignments[0])
    base = 0
    partitions = [ alignments[0].partition_info(base) ]
    base += alignments[0].sequence_length()
    if len(alignments) > 1:
        for a in alignments[1:]:
            partitions.append( a.partition_info(base) )
            for k in a.keys():
                if combined_alignment.has_key(k):
                    combined_alignment[k] += a[k]
                else:
                    for i in combined_alignment.keys():
                        if not a.has_key(i):
                            combined_alignment[i] += 'N'*len(a[k])
                    combined_alignment[k] = 'N'*base + a[k]
    if len(set([a.datatype for a in alignments])) == 1:
        combined_alignment.datatype = alignments[0].datatype
    else:
        combined_alignment.datatype = "MIXED"
    return (combined_alignment, partitions)

def disassemble_alignment(alignment, partition):
    alignments = []
    for p in partition:
        start = p[1] - 1
        end = p[2]
        segment = Alignment()
        for k in alignment.keys():
            if not all([i == 'N' for i in alignment[k][start:end]]):
                segment[k] = alignment[k][start:end]
            if alignment.datatype != 'MIXED' or 'UNKNOWN':
                segment.datatype = alignment.datatype
        alignments.append(segment)
    return alignments

class SequenceDataset(object):
    """Class for creating a dendropy reader, validating the input, and
    keeping mapping of real taxa names to "safe" versions that will not
    cause problems for our alignment and tree tools.

    The general order of calls should be:

    ############################################################################
    # Initialization
    ############################################################################
    sd = SequenceDataset()
    sd.read(file_obj, file_format='FASTA', datatype=datatype)

    ############################################################################
    # Check matrix
    ############################################################################
    assert sd.sequences_are_valid(remap_missing, map_missing_to)

    ############################################################################
    # read trees before changing taxa names
    ############################################################################
    tree_list = sd.dataset.read_trees(tree_f, 'NEWICK', encode_splits=True)

    ############################################################################
    # Go to safe labels
    ############################################################################
    alignment = sd.relabel_for_sate()

    ############################################################################
    # use the dataset object
    ############################################################################
    job = SateJob(alignment=alignment,
                    sate_team=sate_team,
                    name=options.jobname,
                    dataset=sd.dataset
                )
    job.tree = tree_list[0]

    job.run(tmp_dir_par=temporaries_dir)

    ############################################################################
    # restore the original names to change the dataset object held by the job
    ############################################################################
    sd.restore_taxon_names()

    ############################################################################
    # get the tree with the original labels
    ############################################################################
    tree_str = job.tree.as_newick_str()
    """

    def __init__(self):
        self.dataset = None
        self.alphabet = None
        self.real_to_safe_names = {}
        self.safe_to_real_names = {}
        self.datatype = None

    # def relabel_for_sate(self):

    def relabel_for_sate(self, make_names_safe=True):
        "Relabels the dataset taxa and returns an Alignment instance"
        try:
            taxa_block = self.taxa
            char_block = self.character_matrix
        except:
            log_exception(_LOG)
            return None

        self.safe_to_real_names = {}
        a = Alignment()
        a.datatype = self.datatype
        for taxon in taxa_block:
            char_vec = char_block[taxon]
            safe_name = self._register_safe_name(taxon.label, make_names_safe=make_names_safe)
            _LOG.debug("%s (%d) -> %s" % (taxon.label, id(taxon), safe_name))
            taxon.label = safe_name
            a[safe_name] = char_vec
        return a

    def get_character_matrix(self):
        """Returns the first character matrix or raises IndexError if no
        characters have been read."""
        return self.dataset.char_matrices[0]
    character_matrix = property(get_character_matrix)

    def get_taxa_block(self):
        """Returns the list of taxa."""
        return getattr(self.dataset, DATASET_TAXA_ATTR)[0]
    taxa = property(get_taxa_block)

    def set_alignment(self, aln):
        """replaces content of the current character matrix with the sequences
        from `aln`."""
        taxa = self.taxa
        char_block = self.character_matrix
        char_block.clear()
        for k, v in aln.iteritems():
            taxon = taxa.get_taxon(label=k)
            char_block[taxon] = v


    def _register_safe_name(self, name, make_names_safe=True):
        "Creates a unique entry in safe_to_real_names for name `n` (if needed)."
        ind = 0
        real_name = name
        if make_names_safe:
            safe_name_prefix = "".join(_DANGEROUS_NAME_CHARS.split(name))[:80].lower()
            safe_name = safe_name_prefix
        else:
            safe_name = name
        while True:
            if safe_name not in self.safe_to_real_names:
                self.safe_to_real_names[safe_name] = real_name
                return safe_name
            ind += 1
            safe_name = safe_name_prefix + str(ind)

    def restore_taxon_names(self):
        """Changes the labels in the contained dataset's first taxa_block."""
        try:
            taxa_block = getattr(self.dataset, DATASET_TAXA_ATTR)[0]
        except:
            return
        for taxon in taxa_block:
            safe_name = taxon.label
            real_name = self.safe_to_real_names[safe_name]
            taxon.label = real_name

        self.real_to_safe_names = {}
        self.safe_to_real_names = {}

    def read(self, file_obj, file_format='FASTA', datatype=None):
        """If the datatype is fasta (or some other type that does not
        specify the type of data, then the datatype arg should be DNA, RNA
        or 'PROTEIN'
        """
        fup = file_format.upper()
        amibig_formats = ['FASTA']
        if fup in amibig_formats:
            if not datatype:
                raise ValueError('datatype must be specified when the file_format is %s' % fup)
            dup = datatype.upper()
            datatype_list = ['DNA', 'RNA', 'PROTEIN']
            if not dup in datatype_list:
                raise ValueError('Expecting the datatype to be  DNA, RNA or PROTEIN')
            file_format = dup + file_format
        try:
            import dendropy
            self.dataset = dendropy.DataSet()
            self.dataset.read(file_obj, schema=file_format, row_type='str')
            n1 = len(self.dataset.taxon_sets[0].labels())
            n2 = len(set(self.dataset.taxon_sets[0].labels()))
            if n1 != n2:
                raise ValueError("There are redundant sequence names in your data set!")
                sys.exit(1)
        except:
            self.dataset = None
            raise
        try:
            tb = getattr(self.dataset, DATASET_TAXA_ATTR)[0]
            self.datatype = dup
        except:
            raise ValueError("No data was read from the file.")
        tb.lock()

    def sequences_are_valid(self, remap_missing=False, map_missing_to=None):
        """Check for ? in sequences"""
        try:
            taxa_block = getattr(self.dataset, DATASET_TAXA_ATTR)[0]
            char_block = getattr(self.dataset, DATASET_CHAR_ATTR)[0]
        except:
            raise ValueError("Data have not been read")
        try:
            sa_list = char_block.state_alphabets
            self.alphabet = sa_list[0]
            missing = self.alphabet.missing
        except:
            raise ValueError("Expecting a simple datatype with one state alphabet")
        if missing is None:
            raise ValueError("Expecting a DNA, RNA, or amino acid sequences")

        for taxon in taxa_block:
            char_vec = char_block[taxon]
            missing_inds = []
            for ind, s in enumerate(char_vec):
                if s == '?':
                    missing_inds.append(ind)
            if missing_inds:
                if remap_missing:
                    as_list = list(char_vec)
                    if map_missing_to:
                        for ind in missing_inds:
                            as_list[ind] = map_missing_to
                    else:
                        missing_inds.sort(reverse=True)
                        for ind in missing_inds:
                            as_list.pop(ind)
                    char_block[taxon] = ''.join(as_list)
                else:
                    return False
        return True

class MultiLocusDataset(list):
    pass
