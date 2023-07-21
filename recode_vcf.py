'''
need to recode only recessive and dominant!
in recessive : heterozygotes should be recoded to 00 : het behaves as ref/ref
in dominant  : heterozygotes should be recoded to 11 : het behaves as hom/hom
'''
import gzip

infile = snakemake.input.infile
outfile = snakemake.output.outfile
inheri = snakemake.wildcards.inheri

# infile = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/GWAS/GENOTYPES/fv/female/CHR1/SPLIT/chunk1.vcf.gz'
# outfile = '/cluster/work/pausch/naveen/RECOMBINATION/todel.vcf'
out = open(outfile, "w")



#either the ref allele is recessive or the alt allele is recessive
#in the first case refref =0|1 and refalt,altalt =0|0
#in the second case altalt=0|1 and refalt, refref = 0|0
#so the modes are better named as ref and alt

if inheri == "ref":
    recode={
        "0|0" : "0|1",
        "0|1" : "0|0",
        "1|0" : "0|0",
        "1|1" : "0|0"
    }
elif inheri =="alt":
        recode={
        "0|0" : "0|0",
        "0|1" : "0|0",
        "1|0" : "0|0",
        "1|1" : "0|1"
    }
else:
    exit('inheritance not recognized')


with gzip.open(infile, 'rt') as inf:
    for line in inf:
        if line[0] != '#':
            spl = line.rstrip().split()
            if int(spl[1]) % 1000_000 == 0:
                print(spl[1])
            if ',' not in spl[3] and ',' not in spl[4]:
                info = spl[:9]
                gts = spl[9:]
                gts = [recode.get(gt.split(':')[0], ".|.") for gt in gts]
                tw = '\t'.join(info + gts)
                out.write(f'{tw}\n')

        else:
            out.write(line)

out.close()
