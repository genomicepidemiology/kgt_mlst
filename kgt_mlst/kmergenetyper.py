import logging
import os
import sys
import subprocess
from pathlib import Path

class kmergenetyperRunner():
    def __init__(self, input, database, md, output):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.input = input
        self.database = database
        self.md = md
        self.output = output


    def run(self):
        """runs kmergenetyper"""
        if len(self.input.split()) > 1:
            self.nanopore = None
            self.illumina = self.input
        else:
            self.nanopore = self.input
            self.illumina = None
        kmergenetyper_cmd = "kmergenetyper -t_db {} -o {} -md {}".format(self.database, self.output, self.md)
        if self.nanopore:
            kmergenetyper_cmd += " -nanopore {}".format(self.nanopore)
        if self.illumina:
            kmergenetyper_cmd += " -illumina {}".format(self.illumina)
        self.logger.info("Running kmergenetyper with the following command: {}".format(kmergenetyper_cmd))
        os.system(kmergenetyper_cmd)