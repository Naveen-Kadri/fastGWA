infiles = snakemake.input.infiles
outfile = snakemake.output.outfile
out = open (outfile, "w")

print (infiles)

for infile in infiles:
    if "id" in infile:
        out.write (f'{infile [:-7]}\n')

out.close ()
