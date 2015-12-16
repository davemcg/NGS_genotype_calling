#!/usr/local/Anaconda/envs/py3.4.3/bin/python

"""
Fills out read group information for NISC-provided bam files
"""
import argparse
import subprocess

parser.add_argument("file", help="Input bam file to generate read group information")
args = parser.parse_args()
file = args.file

samtools_view = subprocess.check_output("samtools view -h CCGO_800067.bam  | head -n 100 | grep ^@RG", shell=True)
print(samtools_view)
print(samtools_view)
print(file)
