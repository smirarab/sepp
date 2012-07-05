'''
Created on Jun 7, 2011

@author: smirarab
'''
from sepp.tools import ExternalTool
import os
from sepp.scheduler import DispatchableJob, jobq
from sepp.filemgr import temp_fs
from sepp.settings import settings
from sepp.utilities import convert_fasta_to_stockholm, read_sto_alignment
import copy
from sepp.alignment import is_taxon_a_full_length_seq

            
class ShortReadAlignment(object):
    '''
    Aligns many short reads to a subset of the reference alignment
    '''

    def __init__(self, fullReferenceAlignment):
        '''
        Constructor
        fullReferenceAlignment is the CharacterMatrix of all full length sequences        
        '''
        self.fullReferenceAlignment = fullReferenceAlignment;

    def alignShortReadToSubset(self, read_char_mat, subset_labels, tmp_dir_par, **kwargs):
        raise NotImplementedError("This method is not implement in abstract class ShortReadAlignment")

    def get_alignment_inducedby (self, subset_labels):
        return self.fullReferenceAlignment.sub_alignment(subset_labels)    
       
        
class HMMERShortReadAligner(ShortReadAlignment):

    def __init__(self, fullReferenceAlignment):
        ShortReadAlignment.__init__(self, fullReferenceAlignment)
        self.profilerBuilder = HMMERProfileBuilder() 
        self.aligner_to_profile = HMMERAlignerToProfile()
        self.name = 'alignshortreads'
        self.temp_dirs = []
        self.jobs = []
        self.indiviual_alignments = []
    
    def alignShortReadToSubset(self, read_char_mat, subset_labels, tmp_dir_par, **kwargs):        
                
        pref = 'temp'+ self.name
        scratch_dir = temp_fs.create_temp_subdir(parent=tmp_dir_par, prefix=pref)
        self.temp_dirs.append(scratch_dir)
        
        subset_matrix = self.get_alignment_inducedby(subset_labels);
         
        tbj = self.profilerBuilder.create_job(subset_matrix,
                                              scratch_dir,
                                              context_str=kwargs["context_str"] + self.name)
        jobq.put(tbj)
        profile_fn, seq_fn = tbj.get_results()
                        
        refstofn = os.path.join(scratch_dir, "ref.sto")
        convert_fasta_to_stockholm(seq_fn,refstofn)
        
        tbj = self.aligner_to_profile.create_job(read_char_mat,
                                              refstofn,
                                              profile_fn,
                                              scratch_dir,
                                              context_str=kwargs["context_str"] + self.name)
        jobq.put(tbj)
        self.jobs.append(tbj)

        #alignment = tbj.get_results()                

        return tbj
    
    def collect_results_from_jobs(self):
        for job in self.jobs:
            self.indiviual_alignments.append(job.get_results())
        if settings.get_setting("delete_temps"):
            for temp_dir in self.temp_dirs:
                temp_fs.remove_dir(temp_dir)    
        
    def merge_alignmets(self):
        new_alignment = copy.copy(self.fullReferenceAlignment)
        for alg in self.indiviual_alignments:
            new_alignment.merge_alignment_in(alg)
        return new_alignment
    
class HMMERProfileBuilder(ExternalTool):
    def __init__(self, **kwargs):
        self.name = 'hmmbuild'
        ExternalTool.__init__(self, self.name)        
        
    def _prepare_input(self, alignment_char_mat, scratch_dir):
        seqfn = os.path.join(scratch_dir, "ref.fasta")       
        alignment_char_mat.write_to_path(seqfn, schema='fasta')
        profilefn = os.path.join(scratch_dir, 'profile.hmm')
        return seqfn, profilefn

    def create_job(self, alignment_char_mat, scratch_dir, **kwargs):
        seqfn, outfn = self._prepare_input(alignment_char_mat, scratch_dir)
        #--informat afa ${d}.hmm ${d}.fasta
        invoc = [self.exe, '--informat', 'afa', outfn, seqfn]

        job_id = kwargs.get('context_str', '') + '_hmmbuild'


        def _finish_standard_job(outfn, seq_fn, invoc, scratch_dir, job_id):
            rpc = lambda: (outfn, seq_fn)
            job = DispatchableJob(invoc, result_processor=rpc,  cwd=scratch_dir, context_str=job_id)
            return job

        return _finish_standard_job(outfn=outfn, 
                                         seq_fn=seqfn,
                                         invoc=invoc,
                                         scratch_dir=scratch_dir,
                                         job_id=job_id)


class HMMERAlignerToProfile(ExternalTool):
    def __init__(self, **kwargs):
        self.name = 'hmmalign'
        ExternalTool.__init__(self, self.name)

    def _prepare_input(self, read_char_mat, scratch_dir):
        seqfn = os.path.join(scratch_dir, "reads.fasta")
        read_char_mat.write_to_path(seqfn, "fasta")
        alignedfn = os.path.join(scratch_dir, 'combined.sto')
        return seqfn, alignedfn

    def create_job(self, read_char_mat, ref_sto_fn, profile_fn, scratch_dir,**kwargs):
        seqfn, outfn = self._prepare_input(read_char_mat, scratch_dir)
        #hmmalign -o ${d}.combined.sto --mapali tmp.ref.sto ${d}.hmm tmp.fasta
        invoc = [self.exe, '--trim', '--allcol', '-o', outfn, '--mapali', ref_sto_fn, profile_fn, seqfn]

        job_id = kwargs.get('context_str', '') + '_hmmalign'
            
        def _finish_standard_job(outfn, invoc, scratch_dir, job_id):
            rpc = lambda: read_sto_alignment(outfn)
            job = DispatchableJob(invoc, result_processor=rpc,  cwd=scratch_dir, context_str=job_id)
            return job

        return _finish_standard_job(outfn=outfn, 
                                         invoc=invoc,
                                         scratch_dir=scratch_dir,
                                         job_id=job_id)

