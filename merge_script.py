'''
Created on Oct 10, 2012

@author: smirarab
'''
from sepp.alignment import MutableAlignment, ExtendedAlignment
from sepp.exhaustive import JoinAlignJobs
from sepp.problem import SeppProblem
import glob

job_joiner = JoinAlignJobs
original_backbone_file = ('/projects/sate8/namphuon/ultra_large/1000000/'
                          'sate.fasta')
original_backbone = MutableAlignment()
done = original_backbone.read_filepath(original_backbone_file)

original_frag_file = ('/projects/sate8/namphuon/ultra_large/1000000/'
                      'initial.fas.100')
original_frag = MutableAlignment()
done = original_frag.read_filepath(original_frag_file)

# First build extended alignment on entire fragment set
extendedAlignment = ExtendedAlignment(original_frag.get_sequence_names())

dirs = glob.glob('/projects/sate8/namphuon/ultra_large/1000000/upp_100_10_new/'
                 'temp/upp.1_HNlM/root/P_0/A_0_*/')

dirs.reverse()
for dir in dirs:
    print("Working on %s\n" % dir)
    aligned_files = glob.glob('%sFC_*/hmmalign.results.*' % dir)
    sequence_files = glob.glob('%sFC_*/hmmalign.frag.*' % dir)
    base_alignment_file = glob.glob('%s/*.fasta' % dir)
    base_alignment = MutableAlignment()
    done = base_alignment.read_filepath(base_alignment_file[0])
    subbackbone = original_backbone.get_soft_sub_alignment(
        base_alignment.get_sequence_names())
    frags = MutableAlignment()
    sequence_names = []
    for file in sequence_files:
        seq = MutableAlignment()
        done = seq.read_filepath(file)
        done = sequence_names.extend(seq.get_sequence_names())
        for name, seq in seq.items():
            frags[name] = seq.upper()
    problem = SeppProblem(sequence_names)
    problem.set_subalignment(subbackbone)

    mut_subalg = problem.subalignment.get_mutable_alignment()
    remaining_cols = mut_subalg.delete_all_gap()
    problem.annotations["ref.alignment.columns"] = remaining_cols
    problem.fragments = frags
    ap_alg = problem.read_extendend_alignment_and_relabel_columns(
        base_alignment_file, aligned_files)
    extendedAlignment.merge_in(ap_alg, convert_to_string=False)

extendedAlignment.write_to_path("/projects/sate8/namphuon/ultra_large/1000000/"
                                "upp_100_10_new/upp.unmasked.fasta")
extendedAlignment.remove_insertion_columns()
extendedAlignment.write_to_path("/projects/sate8/namphuon/ultra_large/1000000/"
                                "upp_100_10_new/upp.masked.fasta")
