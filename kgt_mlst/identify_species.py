import os
import sys

from kgt_mlst import kma

def auto_identifiy_species(arguments):
    input_string = " ".join(arguments.input)
    kma.KMARunner(input_string,
                  arguments.output + "/bac_species",
                  arguments.db_dir + "/bac_species_db/bac_species_db",
                  "-mem_mode -1t1 -t 8 -Sparse").run()
    with open(arguments.output + '/bac_species.spa', 'r') as f:
        best_score = 0
        genome_size = 0
        for line in f:
            if not line.startswith('#'):
                score = float(line.split('\t')[2])
                if score > best_score:
                    best_score = score
                    best_line = line.split('\t')[0]
                    genome_size = int(line.split('\t')[4])

    return best_line, genome_size





