import os
import sys

from kgt_mlst import kma

def auto_identifiy_species(arguments):
    input_string = " ".join(arguments.input)
    print (input_string)
    kma.KMARunner(input_string,
                  arguments.output + "/bac_species",
                  arguments.db_dir + "/bac_species_db/bac_species_db",
                  "-mem_mode -1t1 -t 8 -Sparse").run()
    print (arguments.input)
    name = arguments.input[0].split('/')[-1].split('.')[0]
    print (arguments.output + '/' + name + '.res')
    sys.exit()
    with open(arguments.output + '/' + name + '.res', 'r') as f:
        best_score = 0
        for line in f:
            if not line.startswith('#'):
                score = float(line.split('\t')[1])
                if score > best_score:
                    best_score = score
                    best_line = line.split('\t')[0]

    return best_line




