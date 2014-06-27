from sepp.exhaustive_tipp import *
import os,tempfile,re
from sepp.alignment import MutableAlignment
from sepp.alignment import _write_fasta
from sepp.config import options
'''
Collection of functions for metagenomic pipeline for taxonomic classification
Created on June 3, 2014

@author: namphuon
'''
global character_map,taxon_map,level_map,key_map,marker_genes
character_map = {'A':'T', 'a':'t', 'C':'G', 'c':'g', 'T':'A', 't':'a', 'G':'C', 'g':'c', '-':'-'}
global levels
levels = ["species", "genus", "family", "order", "class", "phylum"]
marker_genes = ["nusA","rplB","rplK","rplS","rpsE","rpsS", "pgk","rplC","rplL","rplT","rpsI","smpB","dnaG","pyrg", "rplD","rplM","rpmA","rpsJ","frr","pyrG","rplE","rplN", "rpsB","rpsK","infC","rplA","rplF","rplP","rpsC","rpsM"]

#TODO Fix parameter passing
#TODO Make taxonomy loading a class
def load_taxonomy(taxonomy_file, lower=True):
  f = open(taxonomy_file, 'r')
  
  #First line is the keywords for the taxonomy, need to map the keyword to the positional 
  #index of each keyword
  results = f.readline().lower().replace('"','').strip().split(',')
  key_map = dict([(results[i],i) for i in xrange(0,len(results))])
    
  #Now fill up taxonomy, level maps keep track of what taxa exist at each level, taxon_map  
  #keep track of entire taxonomy
  taxon_map = {}
  level_map = {"species":{}, "genus":{}, "family":{}, "order":{}, "class":{}, "phylum":{}}

  for line in f:
    results = line.replace('"','').strip()
    if (lower):
      results.lower()
    results = results.split(',')
    #insert into taxon map
    taxon_map[results[0]] = results
    
    #insert into level map
    for level in levels:      
      if (results[key_map[level]] == ''):
        continue
      else:
        if (results[key_map[level]] not in level_map[level]):
          level_map[level][results[key_map[level]]] = {}
        level_map[level][results[key_map[level]]][results[0]]=results[0]
  return (taxon_map, level_map, key_map)
  

def build_profile(input,output_directory):  
  global taxon_map,level_map,key_map,levels
  temp_dir=tempfile.mkdtemp(dir=options().__getattribute__('tempdir'))
  if (options().bin == 'blast'):
    binned_fragments=blast_to_markers(input,temp_dir)
  else:
    binned_fragments=hmmer_to_markers(input,temp_dir)
  
  if binned_fragments:
    print "Finished binning"
  else:
    print "Unable to bin any fragments!\n"
    return
  
  #load up taxonomy for 30 marker genes
  (taxon_map, level_map, key_map) = load_taxonomy(os.path.join(options().__getattribute__('reference').path, 'refpkg/rpsB.refpkg/all_taxon.taxonomy'))
    
  #all classifications stored here  
  classifications = {}

  #Now run TIPP on each fragment    
  for (gene,frags) in binned_fragments.items():    
    #Get size of each marker
    total_taxa = 0
    with open(os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/sate.size'%gene), 'r') as f:
      total_taxa = int(f.readline().strip())
    decomp_size = options().alignment_size
    if (decomp_size > total_taxa):
      decomp_size = int(total_taxa/2)
    cpus = options().cpu
    if (len(frags.keys()) < cpus):
      cpus = len(frags.keys())
    os.system('run_tipp.py -c %s --cpu %s -m %s -f %s -t %s -adt %s -a %s -r %s -tx %s -txm %s -at %0.2f -pt %0.2f -A %d -P %d -p %s -o %s -d %s' % (options().config_file.name, cpus, options().molecule, temp_dir+"/%s.frags.fas.fixed" % gene,os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/sate.taxonomy'%gene),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/sate.tree'%gene),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/sate.fasta'%gene),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/sate.taxonomy.RAxML_info'%gene),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/all_taxon.taxonomy'%gene),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/species.mapping'%gene),options().alignment_threshold,options().placement_threshold,decomp_size,total_taxa,temp_dir+"/temp_file","tipp_%s" % gene,output_directory+"/markers/"))
    if (not os.path.exists(output_directory+"/markers/tipp_%s_classification.txt" % gene)):
      continue

    gene_classification = generate_classification(output_directory+"/markers/tipp_%s_classification.txt" % gene,options().placement_threshold)

    #Now write individual classification and also pool classifications    
    write_classification(gene_classification, output_directory+"/markers/tipp_%s.classification" % gene)    
    classifications.update(gene_classification)    
  remove_unclassified_level(classifications)
  write_classification(classifications, output_directory+"/markers/all.classification")
  write_abundance(classifications,output_directory)

def remove_unclassified_level(classifications,level=6):
  global taxon_map,level_map,key_map,levels
  frags = classifications.keys()
  for frag in frags:
    if classifications[frag][level] == 'NA':
      del classifications[frag]
  
def write_classification(class_input, output):
  '''Writes a classification file
  '''
  class_out = open(output, 'w')    
  class_out.write("fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\n");
  keys = class_input.keys()
  keys.sort()
  for frag in keys:
    class_out.write("%s\n" % "\t".join(class_input[frag]));
  class_out.close()    
  
def write_abundance(classifications,output_dir,labels=True,remove_unclassified=True):
  global taxon_map,level_map,key_map,levels  
  
  level_abundance = {1:{'total':0}, 2:{'total':0}, 3:{'total':0}, 4:{'total':0}, 5:{'total':0}, 6:{'total':0}}

  level_names = {1:'species', 2:'genus', 3:'family', 4:'order', 5:'class', 6:'phylum'}
  for lineage in classifications.values():
    #insert into level map
    for level in xrange(1,7):
      if (lineage[level] == 'NA'):
        continue
      else:
        if (lineage[level] not in level_abundance[level]):
          level_abundance[level][lineage[level]] = 0
        level_abundance[level][lineage[level]]+=1        
        level_abundance[level]['total']+=1    
  for level in level_names:
    f = open(output_dir + "/abundance.%s.csv" % level_names[level],'w');
    f.write('taxa\tabundance\n')
    lines = []
    for clade in level_abundance[level]:
      if clade == 'total':
        continue
      name = clade
      if labels:
        name = taxon_map[clade][key_map['tax_name']]
      lines.append('%s\t%0.4f\n' % (name,float(level_abundance[level][clade])/level_abundance[level]['total']))
    lines.sort()
    f.write(''.join(lines))  
    f.close()
    
def generate_classification(class_input,threshold):    
  global taxon_map,level_map,key_map,levels
  class_in = open(class_input, 'r')
  level_map_hierarchy = {"species":0, "genus":1, "family":2, "order":3, "class":4, "phylum":5,"root":6}
  #Need to keep track of last line so we can determine when we switch to new classification
  old_name = "";
  old_probability = 1;
  old_id = ""; 
  old_rank = ""; 
  
  #keep track of all fragment names
  names = {}
  classification = {}
  for line in class_in:    
    results = line.strip().split(',')
    if (len(results) > 5):
      results = [results[0], results[1], results[2], results[-2], results[-1]]
    (name, id, rank, probability) = (results[0], results[1], results[3], float(results[4]));
    names[name] = name;
    if (name != old_name):
      #when we switch to new fragment, output last classification for old fragment
      if (old_name != ""):
        lineage = taxon_map[old_id];
        output_line = [old_name]        
        for level in levels:
          clade = lineage[key_map[level]];
          if (clade == ""):
            clade = "NA"        
          output_line.append(clade)        
        classification[old_name]=output_line      
      old_name = name;
      old_rank = "root";
      old_probability = 1;
      old_id = '1';
    
    #Switch to new rank if the new probability is higher than threshold 
    #and our rank is more specific than our original rank
    if (rank in level_map_hierarchy and (level_map_hierarchy[old_rank] > level_map_hierarchy[rank]) and (probability > threshold)):
      old_rank = rank
      old_probability = probability
      old_id = id
    #Switch to new rank if the new rank matches old rank but has higher probability
    elif (rank in level_map_hierarchy and (level_map_hierarchy[old_rank] == level_map_hierarchy[rank]) and (probability > old_probability)):
      old_rank = rank
      old_probability = probability
      old_id = id               
      
  lineage = taxon_map[old_id];
  output_line = [old_name]        
  for level in levels:
    clade = lineage[key_map[level]];
    if (clade == ""):
      clade = "NA"        
    output_line.append(clade)        
  classification[name]=output_line
  return classification

def hmmer_to_markers(input,temp_dir):
  global marker_genes
  fragments = MutableAlignment()
  fragments.read_filepath(input)
  
  reverse = dict([(name+'_rev',reverse_sequence(seq)) for (name,seq) in fragments.items()])
  all_frags = MutableAlignment()
  all_frags.set_alignment(fragments)
  all_frags.set_alignment(reverse)
  frag_file=temp_dir+"/frags.fas"
  _write_fasta(all_frags,frag_file)
  
  #Now bin the fragments
  frag_scores = dict([(name,[-10000,'NA','NA']) for name in fragments.keys()])
  for gene in marker_genes:    
    #Now run HMMER search
    hmmer_search(frag_file,os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/sate.profile'%gene),temp_dir+"/%s.out" % gene)
    results=read_hmmsearch_results(temp_dir+"/%s.out" % gene)
    
    #Now select best direction for each frag
    for name in results.keys():
      bitscore = results[name][1]
      direction = 'forward'
      true_name = name
      if (name.find('_rev') != -1):
        true_name=true_name.replace('_rev','')
        direction = 'reverse'
      if frag_scores[true_name][0] < bitscore:
        frag_scores[true_name] = [bitscore,gene,direction]
    
  #Now bin the fragments
  genes = dict([])
  for name in frag_scores.keys():
    if (frag_scores[name][1] not in genes):
      genes[frag_scores[name][1]] = {}
    if (frag_scores[name][2] == 'forward'):
      genes[frag_scores[name][1]][name] = fragments[name]
    else:
      genes[frag_scores[name][1]][name] = reverse_sequence(fragments[name])    
  genes.pop("NA", None)      
  for gene in genes.keys():
    gene_file=temp_dir+"/%s.frags.fas" % gene
    _write_fasta(genes[gene],gene_file+".fixed")
  return genes  
  
  
def blast_to_markers(input,temp_dir):
  fragments = MutableAlignment()
  fragments.read_filepath(input)

  if (options().gene == None):    
    #First blast sequences against all markers    
    blast_results=temp_dir+"/blast.out"
    if (options().blast_file == None):    
      print "Blasting fragments against marker dataset\n"
      blast_fragments(input,blast_results)
    else:
      blast_results=options().blast_file
    #Next bin the blast hits to the best gene    
    gene_binning = bin_blast_results(blast_results)
  else:
    gene_binning = {options().gene:fragments.keys()}
  #Now figure out direction of fragments
  binned_fragments = dict([(gene,dict([(seq_name,fragments[seq_name]) for seq_name in gene_binning[gene]])) for gene in gene_binning])
  print "Finding best orientation of reads\n"
  for (gene,frags) in binned_fragments.items():
    #Add reverse complement sequence
    frags_rev = dict([(name+'_rev',reverse_sequence(seq)) for (name,seq) in frags.items()])
    gene_frags = MutableAlignment()
    gene_frags.set_alignment(frags)
    gene_frags.set_alignment(frags_rev)
    gene_file=temp_dir+"/%s.frags.fas" % gene
    _write_fasta(gene_frags,gene_file)
    
    #Now run HMMER search
    hmmer_search(gene_file,os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/sate.profile'%gene),temp_dir+"/%s.out" % gene)
    results=read_hmmsearch_results(temp_dir+"/%s.out" % gene)
    
    #Now select best direction for each frag
    for key in frags:
      forward_score = -10000
      backward_score = -10000
      if (key in results):
        forward_score = results[key][1]
      if (key+"_rev" in results):
        backward_score = results[key+"_rev"][1]
      if (backward_score > forward_score):
        frags[key]=gene_frags[key+"_rev"]
    
    #Now write to file
    _write_fasta(frags,gene_file+".fixed")
    binned_fragments[gene]=frags
  return binned_fragments  
    
def read_hmmsearch_results(input):
  #Group 1 (e-value) 2 (bitscore) and 9 (taxon name) contain the relevant information, other ones can be ignored unless we plan to do something later
  pattern = re.compile(r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
  start_reading = False
  infile=open(input)
  results = {}
  for line in infile:            
      line = line.strip()            
      if (not start_reading and line.startswith("E-value") == True):
          start_reading = True
      elif (start_reading and line == ""):
          start_reading = False
          break
      elif (start_reading):             
          matches = pattern.search(line)
          if (matches is not None and matches.group(0).find("--") == -1):
              results[matches.group(9).strip()] = (float(matches.group(1).strip()), float(matches.group(2).strip()))
  return results
    
def read_mapping(input, header=False, delimiter='\t'):
  '''Read a mapping file
  '''
  d = {}
  with open(input) as f:
      for line in f:
        if (header == True):
          next
        results = line.strip().split(delimiter)        
        d[results[0]] = results
  return d

def bin_blast_results(input):
  #Map the blast results to the markers
  gene_mapping = read_mapping(os.path.join(options().__getattribute__('reference').path, 'blast/markers/seq2marker.tab'))
  
  genes = {}
  with open(input) as f:
    for line in f:
      results = line.split('\t')        
      gene = gene_mapping[results[1]][1];  
      if gene in genes:
        genes[gene].append(results[0])
      else:
        genes[gene] = [results[0]]
  return genes

def hmmer_search(input, hmmer,output):
  '''Blast the fragments against all marker genes+16S sequences, return output
  '''    
  os.system('%s --noali -E 10000 --cpu %d -o %s %s %s' % (options().__getattribute__('hmmsearch').path,options().cpu, output, hmmer, input))  
  
  
def blast_fragments(input, output):
  '''Blast the fragments against all marker genes+16S sequences, return output
  '''
  os.system('%s -db %s -outfmt 6 -query %s -out %s -num_threads %d -max_target_seqs 1 ' % (options().__getattribute__('blast').path, os.path.join(options().__getattribute__('reference').path, "blast/markers/alignment.fasta.db"), input, output,options().cpu))
    
def reverse_sequence(sequence):
  global character_map 
  #Reverse complement the sequence
  return "".join([character_map.get(a,a) for a in sequence[::-1]])
  
def augment_parser():
    #default_settings['DEF_P'] = (100 , "Number of taxa (i.e. no decomposition)")
    parser = sepp.config.get_parser()
    
    tippGroup = parser.add_argument_group("TIPP Options".upper(), 
                         "These arguments set settings specific to TIPP")                                 
    
    tippGroup.add_argument("-at", "--alignmentThreshold", type = float, 
                      dest = "alignment_threshold", metavar = "N", 
                      default = 0.0,
                      help = "Enough alignment subsets are selected to reach a commulative probability of N. "
                             "This should be a number between 0 and 1 [default: 0.95]")                            

    tippGroup.add_argument("-pt", "--placementThreshold", type = float, 
                      dest = "placement_threshold", metavar = "N", 
                      default = 0.0,
                      help = "Enough placements are selected to reach a commulative probability of N. "
                             "This should be a number between 0 and 1 [default: 0.95]")    
    tippGroup.add_argument("-g", "--gene", type = str, 
                      dest = "gene", metavar = "N", 
                      default = None,
                      help = "Classify on only the specified gene. ")    
                      
    tippGroup.add_argument("-b", "--blast_file", type = str, 
                      dest = "blast_file", metavar = "N", 
                      default = None,
                      help = "Blast file with fragments already binned. ")    

    tippGroup.add_argument("-bin", "--bin_using", type = str, 
                      dest = "bin", metavar = "N", 
                      default = "blast",
                      help = "Tool for binning")    
                      
                      
def main():
    augment_parser() 
    sepp.config._options_singelton = sepp.config._parse_options()            
    if (options().alignment_size is None):
      options().alignment_size = 100
    input = options().fragment_file.name
    output_directory=options().outdir
    build_profile(input,output_directory)

if __name__ == '__main__':   
    main()
