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
global character_map,taxon_map,level_map,key_map,marker_genes,cog_genes
character_map = {'A':'T', 'a':'t', 'C':'G', 'c':'g', 'T':'A', 't':'a', 'G':'C', 'g':'c', '-':'-'}
global levels
levels = ["species", "genus", "family", "order", "class", "phylum"]
marker_genes = ["nusA","rplB","rplK","rplS","rpsE","rpsS", "pgk","rplC","rplL","rplT","rpsI","smpB","dnaG","pyrg", "rplD","rplM","rpmA","rpsJ","frr","pyrG1","rplE","rplN", "rpsB","rpsK","infC","rplA","rplF","rplP","rpsC","rpsM"]
cog_genes = ["COG0049", "COG0088", "COG0094", "COG0100", "COG0184", "COG0201", "COG0522", "COG0012", "COG0052", "COG0090", "COG0096", "COG0102", "COG0185", "COG0202", "COG0525", "COG0016", "COG0080", "COG0091", "COG0097", "COG0103", "COG0186", "COG0215", "COG0533", "COG0018", "COG0081", "COG0092", "COG0098", "COG0124", "COG0197", "COG0256", "COG0541", "COG0048", "COG0087", "COG0093", "COG0099", "COG0172", "COG0200", "COG0495", "COG0552"]

#TODO Fix parameter passing
#TODO Make taxonomy loading a class
def load_taxonomy(taxonomy_file, lower=True):
  global taxon_map, level_map, key_map
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
  if (options().genes == 'markers'):
    (taxon_map, level_map, key_map) = load_taxonomy(os.path.join(options().reference.path, 'refpkg/rpsB.refpkg/all_taxon.taxonomy'))
  else:
    (taxon_map, level_map, key_map) = load_taxonomy(os.path.join(options().reference.path, 'refpkg/COG0012.refpkg/all_taxon.taxonomy'))
    
  #all classifications stored here  
  classifications = {}
  classification_files = []
  #Now run TIPP on each fragment    
  gene_name = 'sate'
  if (options().genes == 'cogs'):
    gene_name = 'pasta'
  for (gene,frags) in binned_fragments.items():    
    #Get size of each marker
    total_taxa = 0
    with open(os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.size'%(gene,gene_name)), 'r') as f:
      total_taxa = int(f.readline().strip())
    decomp_size = options().alignment_size
    if (decomp_size > total_taxa):
      decomp_size = int(total_taxa/2)
    cpus = options().cpu
    if (len(frags.keys()) < cpus):
      cpus = len(frags.keys())
    extra = ''
    if options().dist == True:
      extra = '-D'
      
    if options().cutoff != 0:
      extra = extra+" -C %f" % options().cutoff
    print 'Cmd:\nrun_tipp.py -c %s --cpu %s -m %s -f %s -t %s -adt %s -a %s -r %s -tx %s -txm %s -at %0.2f -pt %0.2f -A %d -P %d -p %s -o %s -d %s %s' % (options().config_file.name, cpus, options().molecule, temp_dir+"/%s.frags.fas.fixed" % gene,os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.taxonomy'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.tree'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.fasta'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.taxonomy.RAxML_info'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/all_taxon.taxonomy'%gene),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/species.mapping'%gene),options().alignment_threshold,0,decomp_size,total_taxa,temp_dir+"/temp_file","tipp_%s" % gene,output_directory+"/markers/", extra)
    
    os.system('run_tipp.py -c %s --cpu %s -m %s -f %s -t %s -adt %s -a %s -r %s -tx %s -txm %s -at %0.2f -pt %0.2f -A %d -P %d -p %s -o %s -d %s %s' % (options().config_file.name, cpus, options().molecule, temp_dir+"/%s.frags.fas.fixed" % gene,os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.taxonomy'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.tree'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.fasta'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.taxonomy.RAxML_info'%(gene,gene_name)),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/all_taxon.taxonomy'%gene),os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/species.mapping'%gene),options().alignment_threshold,0,decomp_size,total_taxa,temp_dir+"/temp_file","tipp_%s" % gene,output_directory+"/markers/", extra))
    if (not os.path.exists(output_directory+"/markers/tipp_%s_classification.txt" % gene)):
      continue

    gene_classification = generate_classification(output_directory+"/markers/tipp_%s_classification.txt" % gene,options().placement_threshold)
    classification_files.append(output_directory+"/markers/tipp_%s_classification.txt" % gene)
    #Now write individual classification and also pool classifications    
    write_classification(gene_classification, output_directory+"/markers/tipp_%s.classification" % gene)    
    classifications.update(gene_classification)    
  remove_unclassified_level(classifications)
  write_classification(classifications, output_directory+"/markers/all.classification")
  write_abundance(classifications,output_directory)
  
  if (options().dist == True):
    distribution(classification_files,output_directory)

def distribution(classification_files,output_dir):  
  global taxon_map,level_map,key_map,levels,level_names
  distribution = {"species":{}, "genus":{}, "family":{}, "order":{}, "class":{}, "phylum":{}}
  total_frags = 0  
  for class_input in classification_files:
    class_in = open(class_input, 'r')
    frag_info = {"species":{'unclassified':1}, "genus":{'unclassified':1}, "family":{'unclassified':1}, "order":{'unclassified':1}, "class":{'unclassified':1}, "phylum":{'unclassified':1}}
    (old_name,old_rank,old_probability,old_line) = ("",None,-1,"")
    for line in class_in:          
      results = line.strip().split(',')
      if (len(results) > 5):
        results = [results[0], results[1], results[2], results[-2], results[-1]]
      (name, id, rank, probability) = (results[0], results[1], results[3], float(results[4]));
      if (rank not in distribution):
        continue      
      if (old_name == ""):
        (old_name,old_rank,old_probability,old_line) = (name,rank,probability,line)
      if (name != old_name):
        total_frags+=1
        assert frag_info['phylum']['unclassified'] != 1
        for clade in frag_info.keys():
          for clade_name in frag_info[clade].keys():
            if (clade_name not in distribution[clade]):
              distribution[clade][clade_name] = 0
            distribution[clade][clade_name]+=frag_info[clade][clade_name]
        frag_info = {"species":{'unclassified':1}, "genus":{'unclassified':1}, "family":{'unclassified':1}, "order":{'unclassified':1}, "class":{'unclassified':1}, "phylum":{'unclassified':1}} 
        (old_name,old_rank,old_probability,old_line) = (name,rank,probability,line)
      if (id not in frag_info[rank]):
        frag_info[rank][id] = 0
      frag_info[rank][id]+=probability
      frag_info[rank]['unclassified']-=probability        
    total_frags+=1
    assert frag_info['phylum']['unclassified'] != 1
    for clade in frag_info.keys():
      for clade_name in frag_info[clade].keys():
        if (clade_name not in distribution[clade]):
          distribution[clade][clade_name] = 0
        distribution[clade][clade_name]+=frag_info[clade][clade_name]
  
  level_names = {1:'species', 2:'genus', 3:'family', 4:'order', 5:'class', 6:'phylum'}
  for level in level_names:
    f = open(output_dir + "/abundance.distribution.%s.csv" % level_names[level],'w');
    f.write('taxa\tabundance\n')
    lines = []
    for clade in distribution[level_names[level]].keys():
      name = clade
      if (name != 'unclassified'):
        name = taxon_map[clade][key_map['tax_name']]      
      lines.append('%s\t%0.4f\n' % (name,float(distribution[level_names[level]][clade])/total_frags))
    lines.sort()
    f.write(''.join(lines))  
    f.close()  
  return distribution
  
    
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
  
#Fix problem with NA being unclassified  
def write_abundance(classifications,output_dir,labels=True,remove_unclassified=True):
  global taxon_map,level_map,key_map,levels  
  
  level_abundance = {1:{'total':0}, 2:{'total':0}, 3:{'total':0}, 4:{'total':0}, 5:{'total':0}, 6:{'total':0}}

  level_names = {1:'species', 2:'genus', 3:'family', 4:'order', 5:'class', 6:'phylum'}
  for lineage in classifications.values():
    #insert into level map
    for level in xrange(1,7):
      if (lineage[level] == 'NA'):        
        if ('unclassified' not in level_abundance[level]):
          level_abundance[level]['unclassified'] = 0      
        level_abundance[level]['unclassified']+=1        
        level_abundance[level]['total']+=1
        #continue
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
      if labels and name != 'unclassified':        
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
  
  if old_id in taxon_map:   
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
  gene_set = marker_genes
  align_name = 'sate'
  if (options().genes == 'cogs'):
    gene_set = cog_genes
    align_name = 'pasta'
  for gene in gene_set:    
    #Now run HMMER search
    hmmer_search(frag_file,os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%.profile'%(gene,align_name)),temp_dir+"/%s.out" % gene)
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
  align_name = 'sate'
  if (options().genes == 'cogs'):
    align_name = 'pasta'
  for (gene,frags) in binned_fragments.items():
    #Add reverse complement sequence
    frags_rev = dict([(name+'_rev',reverse_sequence(seq)) for (name,seq) in frags.items()])
    gene_frags = MutableAlignment()
    gene_frags.set_alignment(frags)
    gene_frags.set_alignment(frags_rev)
    gene_file=temp_dir+"/%s.frags.fas" % gene
    _write_fasta(gene_frags,gene_file)
    
    #Now run HMMER search
    hmmer_search(gene_file,os.path.join(options().__getattribute__('reference').path, 'refpkg/%s.refpkg/%s.hmm'%(gene,align_name)),temp_dir+"/%s.out" % gene)
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
  gene_mapping = read_mapping(os.path.join(options().__getattribute__('reference').path, 'blast/%s/seq2marker.tab' % options().genes))
  
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
  os.system('%s -db %s -outfmt 6 -query %s -out %s -num_threads %d -max_target_seqs 1 ' % (options().__getattribute__('blast').path, os.path.join(options().__getattribute__('reference').path, "blast/%s/alignment.fasta.db" % options().genes), input, output,options().cpu))
    
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
                      
    tippGroup.add_argument("-D", "--dist", 
                      dest = "dist", action='store_true', 
                      default = False,
                      help = "Treat fragments as distribution")    
                      
    tippGroup.add_argument("-C", "--cutoff", type = float, 
                      dest = "cutoff", metavar = "N", 
                      default = 0.0,
                      help = "Placement probability requirement to count toward the distribution. "
                             "This should be a number between 0 and 1 [default: 0.0]")    
                      
    tippGroup.add_argument("-G", "--genes", type = str, 
                      dest = "genes", metavar = "GENES", 
                      default = 'markers',
                      help = "Use markers or cogs genes [default: markers]")                   
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
