#!/usr/local/Anaconda/envs/py3.4.3/bin/python

"""
Fills out read group information for NISC-provided bam files
"""
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("file", help="Input bam file to generate read group information")
args = parser.parse_args()
file = args.file

samtools_input = "samtools view -h " + file + "| head -n 100 | grep ^@RG"
samtools_view = subprocess.check_output(samtools_input, shell=True)
