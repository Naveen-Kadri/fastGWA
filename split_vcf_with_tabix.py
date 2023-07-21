#!/usr/local/bin/python3.10

""" Use tabix to split a vcf into chunks (of size = chunk_size). Stop  when the new chunk has no variant
The last chunk with zero variants is deleted

Might not work for smaller chunks.. chunking stops when no variants are found. Should be ok for chunks >1mb

This is much faster than the split in the snakefile : should be replaced

"""

import gzip
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(prog="split vcf")
# args without the  -- are not optional
parser.add_argument(
    "infile", help="Input vcf [chromosome-wise]", type=str
)  # is there a file type ? what is argparse.FileType('w') ?
parser.add_argument(
    "outdir",
    help="Output directory for the  output files, suffixes 1 - n (number of chunks) will be added",
    type=str,
)


parser.add_argument(
    "-chr", "--chromosome", help="Chromosome number", type=int, required=True
)

parser.add_argument(
    "-c",
    "--chunk_size",
    help="Chunk size in bp, default = 5_000_000",
    default=5000000,
    type=int,
    required=False,
)

parser.add_argument(
    "-o",
    "--overlap",
    help="overlap size in bp, default = 1_000_000",
    default=1000000,
    type=int,
    required=False,
)


args = parser.parse_args()
infile = args.infile
mychr = args.chromosome
chunk_size = args.chunk_size
overlap = args.overlap
outdir = args.outdir

if not os.path.exists(outdir):
    os.makedirs(outdir)


def check_if_var(infile):
    with gzip.open(infile, "rt") as inf:
        for line in inf:
            if line[0] != "#":
                return 1
        return 0


var_present = 1
chunk_number = 0
start = 0
while var_present == 1:
    chunk_number += 1
    end = start + chunk_size
    out_file = f'{outdir}/chunk{chunk_number}.vcf.gz'
    mycommand = f"tabix -h {infile}  {mychr}:{start}-{end} | bgzip > {out_file}"
    print(mycommand)
    myprocess = subprocess.run(mycommand, shell=True)
    if myprocess.returncode != 0:
        print("something went wrong")
        exit()
    var_present = check_if_var(out_file)
    print(f"chunk {chunk_number} {var_present}")
    #start = end + 1
    start = end - overlap


print(
    f"vcf was split into {chunk_number} chunks the last chunk is with zero variants and will be deleted !"
)
os.remove(out_file)
