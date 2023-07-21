infiles = snakemake.input.infiles
outfile = snakemake.output.mbfile


with open(outfile, "w") as out:
    for myfile in infiles:
        out.write(f'{myfile[:-4]}\n')
