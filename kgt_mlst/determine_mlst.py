import os
import sys
import subprocess
import gzip

from kgt_mlst import identify_species
from kgt_mlst import kmergenetyper

def determine_mlst(arguments):
    os.system("mkdir -p {0}".format(arguments.output))
    if arguments.species == None:
        header, genome_size = identify_species.auto_identifiy_species(arguments)
        arguments.species = header.split(" ")[1][0].lower() + header.split(" ")[2].lower()
    input_string = " ".join(arguments.input)
    name = arguments.input[0].split('/')[-1].split('.')[0]
    if not os.path.exists(arguments.db_dir + '/mlst_db/{0}/{0}.fsa'.format(arguments.species)):
        print (name, "Species not found in mlst database", header, sep='\t')
        sys.exit()
    kmergenetyper.kmergenetyperRunner(input_string,
                        arguments.db_dir + '/mlst_db/{0}/{0}'.format(arguments.species),
                        arguments.min_depth, #Insert relative min depth
                        arguments.output + '/mlst').run()
    mlst_type = get_mlst_type(arguments, arguments.output + '/mlst/{}.res'.format(name))
    print (name, mlst_type, header, sep='\t')

def get_mlst_type(arguments, res_file):
    """Returns the mlst results"""
    mlst_genes, mlst_genes_depths, mlst_genes_coverage = parse_kma_res_and_depth(res_file)
    mlst_type, expected_genes, st_included_mlst_genes = derive_mlst(
        arguments.species, mlst_genes, mlst_genes_depths, arguments.db_dir, mlst_genes_coverage)

    return mlst_type

def check_allele_template_coverage(mlst_genes, template_coverage, found_genes):
    """
    Checks if the allele depth is above 100 else returns false
    :param mlst_genes:
    :param template_coverage:
    :param found_genes:
    :return:
    """
    flag = True
    for i in range(len(found_genes)):
        if found_genes[i] in mlst_genes:
            if float(template_coverage[i]) < 100:
                flag = False
    return flag

def derive_mlst(species, found_genes, template_depth, db_dir, mlst_genes_coverage):
    """Returns the mlst results"""
    with open(db_dir + '/mlst_db/config', 'r') as fd:
        for line in fd:
            if line[0] != "#":
                line = line.rstrip().split("\t")
                if line[0] == species:
                    expected_genes = line[2].split(",")

    multiple_alelles, multiple_allele_list, mlst_bool, mlst_genes = check_muliple_alelles(found_genes, expected_genes)

    if mlst_bool:
        if multiple_alelles:
            #Simply selecting highest depth hits but gives warning!
            if template_depth != 'skip':
                mlst_genes = select_highest_depth_alleles(found_genes, template_depth, expected_genes, multiple_allele_list)
            mlst_type = look_up_mlst("{0}/mlst_db/{1}/{1}.tsv".format(db_dir, species), mlst_genes, expected_genes)
            mlst_type += '+'
        else:
            mlst_type = look_up_mlst("{0}/mlst_db/{1}/{1}.tsv".format(db_dir, species), mlst_genes, expected_genes)
        if not check_allele_template_coverage(mlst_genes, mlst_genes_coverage, found_genes):
            mlst_type += '*'
        return mlst_type, expected_genes, list(mlst_genes)
    else:
        return 'Unknown ST', expected_genes, []

def check_muliple_alelles(found_genes, expected_genes):
    flag = False
    hits = set()
    multiples = set()
    mlst_genes = set()
    for item in found_genes:
        gene = item.split("_")[0]
        if gene in expected_genes:
            if gene not in hits:
                hits.add(gene)
                mlst_genes.add(item)
            else:
                flag = True
                multiples.add(gene)
    if len(hits) == len(expected_genes):
        return flag, multiples, True, mlst_genes
    return flag, multiples, False, mlst_genes

def look_up_mlst(file, mlst_genes, expected_genes):
    if len(mlst_genes) != len(expected_genes):
        return "Not all alleles were found"
    with open(file, 'r') as infile:
        for line in infile:
            if not line.startswith('ST'):
                line = line.rstrip().split("\t")
                potential_mlst_set = set()
                for i in range(len(expected_genes)):
                    potential_mlst_set.add(expected_genes[i] + "_" + line[i + 1])
                if potential_mlst_set == mlst_genes:
                    return line[0]
    return "Unknown ST"

def select_highest_depth_alleles(found_genes, template_depth, expected_genes, multiple_allele_list):
    final_genes = set()
    for i in range(len(template_depth)):
        gene = found_genes[i].split("_")[0]
        if gene not in multiple_allele_list and gene in expected_genes:
            final_genes.add(found_genes[i])

    for item in multiple_allele_list:
        high_score = 0
        name = ''
        for i in range(len(found_genes)):
            if found_genes[i].startswith(item):
                if float(template_depth[i]) > high_score:
                    high_score = float(template_depth[i])
                    name = found_genes[i]
        final_genes.add(name)
    return final_genes


def derive_prefix(file):
    return os.path.basename(file).split('.')[0]

def check_for_kma():
    """Checks if kma is installed"""
    try:
        subprocess.call(["kma"], stdout=open(os.devnull, 'wb'))
    except Exception:
        sys.exit("kma is not installed correctly directly in the PATH.")

def parse_kma_res_and_depth(file):
    genes = []
    template_depth = []
    coverage = []
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                genes.append(line.strip().split('\t')[0])
                template_depth.append(line.strip().split('\t')[-3])
                coverage.append(line.strip().split('\t')[5])
    return genes, template_depth, coverage

def number_of_bases_in_file(filename):
    """Returns the number of bases in a file"""
    if filename.endswith('.gz'):
        type = 'fastq.gz'
    elif filename.endswith('.fastq') or filename.endswith('.fq'):
        type = 'fastq'
    elif filename.endswith('.fasta') or filename.endswith('.fa') or filename.endswith('.fna') or filename.endswith('.fsa'):
        type = 'fasta'
    #determine type#
    #TBD - add support for gzipped files
    sum = 0
    if type == 'fasta':
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sum += len(line.strip())

    elif type == 'fastq':
        line_count = 1
        with open(filename, 'r') as f:
            for line in f:
                if line_count == 2:
                    sum += len(line.strip())
                line_count += 1
                if line_count == 5:
                    line_count = 1
    elif type == 'fastq.gz':
        line_count = 1
        with gzip.open(filename, 'r') as f:
            for line in f:
                if line_count == 2:
                    sum += len(line.strip())
                line_count += 1
                if line_count == 5:
                    line_count = 1

    return sum
