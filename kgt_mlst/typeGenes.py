import os
import sys
import subprocess

from kmergenetyper import realignConsensus

def type_genes(args):
    check_for_kma()
    os.system('mkdir -p {}'.format(args.output))
    collections = check_prefix_collisions(args)
    if args.illumina != []:
        args.illumina.sort()
        for i in range(0, len(args.illumina), 2):
            prefix = derive_prefix(args.illumina[i])
            if prefix in collections:
                prefix = prefix + '_illumina'
            input_seqs = ' '.join(args.illumina[i:i+2])
            os.system('kma -ipe {} -o {}/{} -t_db {} -ill -md {} -ID 95 -eq {} -t {}'.format(input_seqs, args.output, prefix, args.t_db, args.md, args.q_score, args.threads))
            realignConsensus.realign_consensus(args.output, prefix, args.t_db, args.keep, args.threads)
    if args.nanopore != []:
        for item in args.nanopore:
            prefix = derive_prefix(item)
            if prefix in collections:
                prefix = prefix + '_nanopore'
            os.system('kma -i {} -o {}/{} -t_db {} -ont -md {} -ID 95 -eq {} -t {}'.format(item, args.output, prefix, args.t_db, args.md, args.q_score, args.threads))
            realignConsensus.realign_consensus(args.output, prefix, args.t_db, args.keep, args.threads)
    if args.fasta != []:
        for item in args.fasta:
            prefix = derive_prefix(item)
            if prefix in collections:
                prefix = prefix + '_fasta'
            os.system('kma -i {} -o {}/{} -t_db {} -md 0.1 -ID 95 -t {}'.format(item, args.output, prefix, args.t_db, args.threads))
            realignConsensus.realign_consensus(args.output, prefix, args.t_db, args.keep, args.threads)

def check_prefix_collisions(args):
    prefixes = []
    collisions = []
    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            prefix = derive_prefix(args.illumina[i])
            if prefix in prefixes:
                print ('Prefix {} is not unique. Output name has been modified.'.format(prefix))
                collisions.append(prefix)
            else:
                prefixes.append(prefix)
    if args.nanopore != []:
        for item in args.nanopore:
            prefix = derive_prefix(item)
            if prefix in prefixes:
                print ('Prefix {} is not unique. Output name has been modified.'.format(prefix))
                collisions.append(prefix)
            else:
                prefixes.append(prefix)
    return collisions
def derive_prefix(file):
    return os.path.basename(file).split('.')[0]

def check_for_kma():
    """Checks if kma is installed"""
    try:
        subprocess.call(["kma"], stdout=open(os.devnull, 'wb'))
    except Exception:
        sys.exit("kma is not installed correctly directly in the PATH.")
