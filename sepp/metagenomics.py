import os
import tempfile
import re
import sepp
from sepp.alignment import MutableAlignment
from sepp.alignment import _write_fasta
from sepp.config import options
'''
Collection of functions for metagenomic pipeline for taxonomic classification
Created on June 3, 2014

@author: namphuon
'''
global character_map, taxon_map, level_map, key_map, refpkg
character_map = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'T': 'A',
                 't': 'a', 'G': 'C', 'g': 'c', '-': '-'}
global levels
levels = ["species", "genus", "family", "order", "class", "phylum"]


def load_reference_package():
    global refpkg

    refpkg = {}

    path = os.path.join(options().__getattribute__('reference').path,
                        options().genes)
    input = os.path.join(path, "file-map-for-tipp.txt")

    refpkg["genes"] = []
    with open(input) as f:
        for line in f.readlines():
            [key, val] = line.split('=')

            [key1, key2] = key.strip().split(':')
            val = os.path.join(path, val.strip())

            try:
                refpkg[key1][key2] = val
            except KeyError:
                refpkg[key1] = {}
                refpkg[key1][key2] = val

            if (key1 != "blast") and (key1 != "taxonomy"):
                refpkg["genes"].append(key1)


# TODO Fix parameter passing
# TODO Make taxonomy loading a class
def load_taxonomy(taxonomy_file, lower=True):
    global taxon_map, level_map, key_map
    f = open(taxonomy_file, 'r')

    # First line is the keywords for the taxonomy, need to map the keyword to
    # the positional index of each keyword
    results = f.readline().lower().replace('"', '').strip().split(',')
    key_map = dict([(results[i], i) for i in range(0, len(results))])

    # Now fill up taxonomy, level maps keep track of what taxa exist at each
    # level, taxon_map keep track of entire taxonomy
    taxon_map = {}
    level_map = {"species": {}, "genus": {}, "family": {}, "order": {},
                 "class": {}, "phylum": {}}

    for line in f:
        results = line.replace('"', '').strip()
        if (lower):
            results.lower()
        results = results.split(',')
        # insert into taxon map
        taxon_map[results[0]] = results

        # insert into level map
        for level in levels:
            if (results[key_map[level]] == ''):
                continue
            else:
                if (results[key_map[level]] not in level_map[level]):
                    level_map[level][results[key_map[level]]] = {}
                level_map[level][results[key_map[level]]][results[0]] = \
                    results[0]
    return (taxon_map, level_map, key_map)


def build_profile(input, output_directory):
    global taxon_map, level_map, key_map, levels, refpkg
    temp_dir = tempfile.mkdtemp(dir=options().__getattribute__('tempdir'))
    if (options().bin == 'blast'):
        binned_fragments = blast_to_markers(input, temp_dir)
    else:
        binned_fragments = hmmer_to_markers(input, temp_dir)

    if binned_fragments:
        print("Finished binning")
    else:
        print("Unable to bin any fragments!\n")
        return

    # Load up taxonomy for 30 marker genes
    (taxon_map, level_map, key_map) = load_taxonomy(refpkg["taxonomy"]["taxonomy"])

    # Store all classifications here
    classifications = {}
    classification_files = []

    # Run TIPP on each fragment
    for (gene, frags) in binned_fragments.items():
        # Set placement subset size to equal the size of each marker
        total_taxa = 0
        with open(refpkg[gene]["size"], 'r') as f:
            total_taxa = int(f.readline().strip())

        # Set alignment subset decomposition size
        decomp_size = options().alignment_size
        if (decomp_size > total_taxa):
            decomp_size = int(total_taxa / 2)

        # Set number of CPUS
        cpus = options().cpu
        if (len(frags) < cpus):
            cpus = len(frags)

        # Set extra arguments
        extra = ''
        if options().dist is True:
            extra = "-D"
        if options().max_chunk_size is not None:
            extra = extra + "-F %d" % options().max_chunk_size
        if options().cutoff != 0:
            extra = extra + " -C %f" % options().cutoff
 
        cmd = "run_tipp.py " \
                  + " -c " + options().config_file.name \
                  + " --cpu " + str("%d" % cpus) \
                  + " -m " + options().molecule \
                  + " -f " + temp_dir + '/' + gene + ".frags.fas.fixed" \
                  + " -t " + refpkg[gene]["placement-tree"] \
                  + " -adt " + refpkg[gene]["alignment-decomposition-tree"] \
                  + " -a " + refpkg[gene]["alignment"] \
                  + " -r " + refpkg[gene]["raxml-info-for-placement-tree"] \
                  + " -tx " + refpkg["taxonomy"]["taxonomy"] \
                  + " -txm " + refpkg[gene]["seq-to-taxid-map"] \
                  + " -at " + str("%0.2f" % options().alignment_threshold) \
                  + " -pt 0.0" \
                  + " -A " + str("%d" % decomp_size) \
                  + " -P " + str("%d" % total_taxa) \
                  + " -p " + temp_dir + "/temp_file" \
                  + " -o tipp_" + gene \
                  + " -d " + output_directory + "/markers/ " \
                  + extra

        print(cmd)
        os.system(cmd)

        tipp_output = output_directory \
                          + "/markers/tipp_" \
                          + gene \
                          + "_classification.txt"

        if (not os.path.exists(tipp_output)):
            continue

        classification_files.append(tipp_output)

        gene_classification = generate_classification(
            tipp_output,
            options().placement_threshold)

        # Apply placement threshold to classification data
        gene_classification_output =  output_directory \
                                          + "/markers/tipp_" \
                                          + gene \
                                          + "_classification_" \
                                          + str("%0.2f" % options().placement_threshold) \
                                          + ".txt"
        gene_classification = generate_classification(
            tipp_output,
            options().placement_threshold)
        write_classification(
            gene_classification,
            gene_classification_output)

        # Pool classification
        classifications.update(gene_classification)

    remove_unclassified_level(classifications)
    write_classification(
        classifications,
        output_directory + "/markers/all.classification")
    write_abundance(classifications, output_directory)

    if (options().dist is True):
        distribution(classification_files, output_directory)


def distribution(classification_files, output_dir):
    global taxon_map, level_map, key_map, levels, level_names
    distribution = {"species": {}, "genus": {}, "family": {}, "order": {},
                    "class": {}, "phylum": {}}
    total_frags = 0
    for class_input in classification_files:
        class_in = open(class_input, 'r')
        frag_info = {"species": {'unclassified': 1},
                     "genus": {'unclassified': 1},
                     "family": {'unclassified': 1},
                     "order": {'unclassified': 1},
                     "class": {'unclassified': 1},
                     "phylum": {'unclassified': 1}}
        old_name = ""
        for line in class_in:
            results = line.strip().split(',')
            if (len(results) > 5):
                results = [results[0], results[1], results[2],
                           results[-2], results[-1]]
            (name, id, rank, probability) = (
                results[0], results[1], results[3], float(results[4]))
            if (rank not in distribution):
                continue
            if (old_name == ""):
                old_name = name
            if (name != old_name):
                total_frags += 1
                assert frag_info['phylum']['unclassified'] != 1
                for clade, cladeval in frag_info.items():
                    for clade_name, cnc in cladeval.items():
                        if (clade_name, cnc not in distribution[clade]):
                            distribution[clade][clade_name] = 0
                        distribution[clade][clade_name] += cnc
                frag_info = {"species": {'unclassified': 1},
                             "genus": {'unclassified': 1},
                             "family": {'unclassified': 1},
                             "order": {'unclassified': 1},
                             "class": {'unclassified': 1},
                             "phylum": {'unclassified': 1}}
                old_name = name
            if (id not in frag_info[rank]):
                frag_info[rank][id] = 0
            frag_info[rank][id] += probability
            frag_info[rank]['unclassified'] -= probability
        total_frags += 1
        assert frag_info['phylum']['unclassified'] != 1
        for clade, cladeval in frag_info.items():
            for clade_name, cnc in cladeval.items():
                if (clade_name not in distribution[clade]):
                    distribution[clade][clade_name] = 0
                distribution[clade][clade_name] += cnc

    level_names = {1: 'species', 2: 'genus', 3: 'family', 4: 'order',
                   5: 'class', 6: 'phylum'}
    for level in level_names:
        f = open(output_dir + "/abundance.distribution.%s.csv" %
                 level_names[level], 'w')
        f.write('taxa\tabundance\n')
        lines = []
        for clade, value in distribution[level_names[level]].items():
            name = clade
            if (name != 'unclassified'):
                name = taxon_map[clade][key_map['tax_name']]
            lines.append('%s\t%0.4f\n' % (name, float(value) / total_frags))
        lines.sort()
        f.write(''.join(lines))
        f.close()
    return distribution


def remove_unclassified_level(classifications, level=6):
    global taxon_map, level_map, key_map, levels
    frags = list(classifications.keys())
    for frag in frags:
        if classifications[frag][level] == 'NA':
            del classifications[frag]


def write_classification(class_input, output):
    '''Writes a classification file
    '''
    class_out = open(output, 'w')
    class_out.write("fragment\tspecies\tgenus\tfamily\torder\tclass\tphylum\n")
    keys = list(class_input.keys())
    keys.sort()
    for frag in keys:
        class_out.write("%s\n" % "\t".join(class_input[frag]))
    class_out.close()


# Fix problem with NA being unclassified
def write_abundance(classifications, output_dir, labels=True,
                    remove_unclassified=True):
    global taxon_map, level_map, key_map, levels

    level_abundance = {
        1: {'total': 0}, 2: {'total': 0}, 3: {'total': 0}, 4: {'total': 0},
        5: {'total': 0}, 6: {'total': 0}}

    level_names = {1: 'species', 2: 'genus', 3: 'family', 4: 'order',
                   5: 'class', 6: 'phylum'}
    for lineage in classifications.values():
        # insert into level map
        for level in range(1, 7):
            if (lineage[level] == 'NA'):
                if ('unclassified' not in level_abundance[level]):
                    level_abundance[level]['unclassified'] = 0
                level_abundance[level]['unclassified'] += 1
                level_abundance[level]['total'] += 1
                # continue
            else:
                if (lineage[level] not in level_abundance[level]):
                    level_abundance[level][lineage[level]] = 0
                level_abundance[level][lineage[level]] += 1
                level_abundance[level]['total'] += 1
    for level in level_names:
        f = open(output_dir + "/abundance.%s.csv" % level_names[level], 'w')
        f.write('taxa\tabundance\n')
        lines = []
        for clade in level_abundance[level]:
            if clade == 'total':
                continue
            name = clade
            if labels and name != 'unclassified':
                name = taxon_map[clade][key_map['tax_name']]
            lines.append('%s\t%0.4f\n' % (
                name, float(level_abundance[level][clade]) / level_abundance[
                    level]['total']))
        lines.sort()
        f.write(''.join(lines))
        f.close()


def generate_classification(class_input, threshold):
    global taxon_map, level_map, key_map, levels
    class_in = open(class_input, 'r')
    level_map_hierarchy = {"species": 0, "genus": 1, "family": 2, "order": 3,
                           "class": 4, "phylum": 5, "root": 6}
    # Need to keep track of last line so we can determine when we switch to
    # new classification
    old_name = ""
    old_probability = 1
    old_id = ""
    old_rank = ""

    # keep track of all fragment names
    names = {}
    classification = {}
    for line in class_in:
        results = line.strip().split(',')
        if (len(results) > 5):
            results = [results[0], results[1], results[2],
                       results[-2], results[-1]]
        (name, id, rank, probability) = (
            results[0], results[1], results[3], float(results[4]))
        names[name] = name
        if (name != old_name):
            # when we switch to new fragment, output last classification for
            # old fragment
            if (old_name != ""):
                lineage = taxon_map[old_id]
                output_line = [old_name]
                for level in levels:
                    clade = lineage[key_map[level]]
                    if (clade == ""):
                        clade = "NA"
                    output_line.append(clade)
                classification[old_name] = output_line
            old_name = name
            old_rank = "root"
            old_probability = 1
            old_id = '1'

        # Switch to new rank if the new probability is higher than threshold
        # and our rank is more specific than our original rank
        if (
            rank in level_map_hierarchy and
            (level_map_hierarchy[old_rank] > level_map_hierarchy[rank]) and
            (probability > threshold)
           ):
            old_rank = rank
            old_probability = probability
            old_id = id
        # Switch to new rank if the new rank matches old rank but has higher
        # probability
        elif (
              rank in level_map_hierarchy and
              (level_map_hierarchy[old_rank] == level_map_hierarchy[rank]) and
              (probability > old_probability)
             ):
            old_rank = rank
            old_probability = probability
            old_id = id

    if old_id in taxon_map:
        lineage = taxon_map[old_id]
        output_line = [old_name]
        for level in levels:
            clade = lineage[key_map[level]]
            if (clade == ""):
                clade = "NA"
            output_line.append(clade)
        classification[name] = output_line
    return classification


def hmmer_to_markers(input, temp_dir):
    global refpkg
 
    fragments = MutableAlignment()
    fragments.read_filepath(input)

    reverse = dict([(name+'_rev', reverse_sequence(seq))
                    for (name, seq) in fragments.items()])
    all_frags = MutableAlignment()
    all_frags.set_alignment(fragments)
    all_frags.set_alignment(reverse)
    frag_file = temp_dir+"/frags.fas"
    _write_fasta(all_frags, frag_file)

    # Now bin the fragments
    frag_scores = dict([(name, [-10000, 'NA', 'NA'])
                        for name in fragments.keys()])

    for gene in refpkg["genes"]:
        # Now run HMMER search
        hmmer_output = temp_dir + '/' + gene + ".out"
        hmmer_search(
            frag_file,
            refpkg[gene]["hmm"],
            hmmer_output)
        results = read_hmmsearch_results(hmmer_output)

        # Now select best direction for each frag
        for name, value in results.items():
            bitscore = value[1]
            direction = 'forward'
            true_name = name
            if (name.find('_rev') != -1):
                true_name = true_name.replace('_rev', '')
                direction = 'reverse'
            if frag_scores[true_name][0] < bitscore:
                frag_scores[true_name] = [bitscore, gene, direction]

    # Now bin the fragments
    genes = dict([])
    for name, val in frag_scores.items():
        if (val[1] not in genes):
            genes[val[1]] = {}
        if (val[2] == 'forward'):
            genes[val[1]][name] = fragments[name]
        else:
            genes[val[1]][name] = reverse_sequence(fragments[name])
    genes.pop("NA", None)
    for gene, seq in genes.items():
        gene_file = temp_dir + "/%s.frags.fas" % gene
        _write_fasta(seq, gene_file + ".fixed")

    return genes


def blast_to_markers(input, temp_dir):
    global refpkg

    fragments = MutableAlignment()
    fragments.read_filepath(input)

    if (options().gene is None):
        # First blast sequences against all markers
        blast_results = temp_dir + "/blast.out"
        if (options().blast_file is None):
            print("Blasting fragments against marker dataset\n")
            blast_fragments(input, blast_results)
        else:
            blast_results = options().blast_file
        # Next bin the blast hits to the best gene
        gene_binning = bin_blast_results(blast_results)
    else:
        gene_binning = {options().gene: list(fragments.keys())}
    
    # Now figure out direction of fragments
    binned_fragments = dict([
        (gene, dict([(seq_name, fragments[seq_name])
                     for seq_name in gene_binning[gene]]))
        for gene in gene_binning])
    
    print("Finding best orientation of reads\n")
    for (gene, frags) in binned_fragments.items():
        # Add reverse complement sequence
        frags_rev = dict([(name + '_rev', reverse_sequence(seq))
                          for (name, seq) in frags.items()])

        gene_frags = MutableAlignment()
        gene_frags.set_alignment(frags)
        gene_frags.set_alignment(frags_rev)
        gene_file = temp_dir + "/%s.frags.fas" % gene
        _write_fasta(gene_frags, gene_file)

        # Now run HMMER search
        hmmer_output = temp_dir + '/' + gene + ".out"
        hmmer_search(
            gene_file,
            refpkg[gene]["hmm"],
            hmmer_output)
        results = read_hmmsearch_results(hmmer_output)

        # Now select best direction for each frag
        for key in frags:
            forward_score = -10000
            backward_score = -10000
            if (key in results):
                forward_score = results[key][1]
            if (key+"_rev" in results):
                backward_score = results[key + "_rev"][1]
            if (backward_score > forward_score):
                frags[key] = gene_frags[key + "_rev"]

        # Now write to file
        _write_fasta(frags, gene_file + ".fixed")
        binned_fragments[gene] = frags
    return binned_fragments


def read_hmmsearch_results(input):
    # Group 1 (e-value) 2 (bitscore) and 9 (taxon name) contain the
    # relevant information, other ones can be ignored unless we plan to do
    # something later
    pattern = re.compile(
        r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)"
        r"\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)")
    start_reading = False
    infile = open(input)
    results = {}
    for line in infile:
        line = line.strip()
        if (not start_reading and line.startswith("E-value") is True):
            start_reading = True
        elif (start_reading and line == ""):
            start_reading = False
            break
        elif (start_reading):
            matches = pattern.search(line)
            if (matches is not None and matches.group(0).find("--") == -1):
                results[matches.group(9).strip()] = (
                    float(matches.group(1).strip()),
                    float(matches.group(2).strip()))
    return results


def read_mapping(input, header=False, delimiter='\t'):
    '''Read a mapping file
    '''
    d = {}
    with open(input) as f:
        for line in f:
            if (header is True):
                next
            results = line.strip().split(delimiter)
            d[results[0]] = results
    return d


def bin_blast_results(input):
    global refpkg

    # Map the blast results to the markers
    gene_mapping = read_mapping(refpkg["blast"]["seq-to-marker-map"])

    genes = {}
    with open(input) as f:
        for line in f:
            results = line.split('\t')
            gene = gene_mapping[results[1]][1]
            if gene in genes:
                genes[gene].append(results[0])
            else:
                genes[gene] = [results[0]]
    return genes


def hmmer_search(input, hmmer, output):
    cmd = options().__getattribute__('hmmsearch').path \
              + " --noali -E 10000 " \
              + " --cpu " + str("%d" % options().cpu) \
              + " -o " + output \
              + " " + hmmer \
              + " " + input
    
    print(cmd)
    os.system(cmd)


def blast_fragments(input, output):
    '''Blast the fragments against all marker genes+16S sequences, return
    output'''
    global refpkg

    cmd = options().__getattribute__('blast').path \
              + " -db " + refpkg["blast"]["database"] \
              + " -outfmt 6 " \
              + " -query " + input \
              + " -out " + output \
              + " -num_threads " + str("%d" % options().cpu) \
              + " -max_target_seqs 1"

    print(cmd)
    os.system(cmd)


def reverse_sequence(sequence):
    global character_map
    #  Reverse complement the sequence
    return "".join([character_map.get(a, a) for a in sequence[::-1]])


def augment_parser():
    parser = sepp.config.get_parser()

    tippGroup = parser.add_argument_group(
        "TIPP Options".upper(),
        "These arguments set settings specific to TIPP")

    tippGroup.add_argument(
        "-at", "--alignmentThreshold", type=float,
        dest="alignment_threshold", metavar="N",
        default=0.0,
        help="Enough alignment subsets are selected to reach a commulative "
             "probability of N."
             "This should be a number between 0 and 1 [default: 0.0]")

    tippGroup.add_argument(
        "-pt", "--placementThreshold", type=float,
        dest="placement_threshold", metavar="N",
        default=0.95,
        help="Enough placements are selected to reach a commulative "
             "probability of N."
             "This should be a number between 0 and 1 [default: 0.95]")
    
    tippGroup.add_argument(
        "-g", "--gene", type=str,
        dest="gene", metavar="N",
        default=None,
        help="Classify on only the specified gene. ")

    tippGroup.add_argument(
        "-b", "--blast_file", type=str,
        dest="blast_file", metavar="N",
        default=None,
        help="Blast file with fragments already binned. ")

    tippGroup.add_argument(
        "-bin", "--bin_using", type=str,
        dest="bin", metavar="N",
        default="blast",
        help="Tool for binning")

    tippGroup.add_argument(
        "-D", "--dist",
        dest="dist", action='store_true',
        default=False,
        help="Treat fragments as distribution")

    tippGroup.add_argument(
        "-C", "--cutoff", type=float,
        dest="cutoff", metavar="N",
        default=0.0,
        help="Placement probability requirement to count toward the "
             "distribution. "
             "This should be a number between 0 and 1 [default: 0.0]")

    tippGroup.add_argument(
        "-G", "--genes", type=str,
        dest="genes", metavar="GENES",
        default="markers-v3",
        help="Set of markers to use [default: markers-v3]")


def main():
    augment_parser()

    sepp.config._options_singelton = sepp.config._parse_options()

    if (options().alignment_size is None):
        print("WARNING: Alignment subset size set to NONE, re-setting to 100.")
        options().alignment_size = 100

    input = options().fragment_file.name

    output_directory = options().outdir

    load_reference_package()

    build_profile(input, output_directory)


if __name__ == '__main__':
    main()
