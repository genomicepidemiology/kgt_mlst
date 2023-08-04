import os
import sys
import subprocess

from kgt_mlst import identify_species

def determine_mlst(arguments):
    specie = identify_species.auto_identifiy_species(arguments)
    print (specie)


def derive_prefix(file):
    return os.path.basename(file).split('.')[0]

def check_for_kma():
    """Checks if kma is installed"""
    try:
        subprocess.call(["kma"], stdout=open(os.devnull, 'wb'))
    except Exception:
        sys.exit("kma is not installed correctly directly in the PATH.")
