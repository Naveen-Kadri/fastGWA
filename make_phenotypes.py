from collections import defaultdict
# key_file = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/HIDDENPHASE/key.txt'
# pheno_file = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/LINKPHASE_RUN2/fv/cleaned_GRR.txt'
# out_file = 'test.txt'
# breed = 'fv'
# noff_thresh = 6
# zscore_thresh = 5
# sex = 'female'
# min_gam = 1


key_file = snakemake.input[0]
pheno_file = snakemake.input[1]

out_file = snakemake.output.pheno
breed = snakemake.wildcards.breed
sex = snakemake.wildcards.sex
min_gam = int(snakemake.wildcards.ngam)
zscore_thresh = float(snakemake.wildcards.std)
noff_thresh = snakemake.params.noff

out = open(out_file, 'w')
sexes = {
    'male': '1',
    'female': '2'
}

sex = sexes[sex]

print('reading the key file')
key = {}
with open(key_file) as inf:
    for line in inf:
        oid, nid = line.rstrip().split()
        key[oid] = nid


grr = defaultdict(list)
print('reading the phenotype file')
with open(pheno_file) as inf:
    for lnum, line in enumerate(inf):
        spl = line.rstrip().split()
        if lnum == 0:
            header = spl
        else:
            info = dict(zip(header, spl))
            if info['sex'] == sex and float(info['zscore']) <= zscore_thresh:
                if info['gp'] == "1" or int(info['noff']) >= noff_thresh:
                    grr[breed+info['par']].append(int(info['GRR']))


for par in grr:
    if par not in key:
        continue
    mygrr = grr[par]
    ngam = len(mygrr)
    if ngam >= min_gam:
        mgrr = sum(mygrr) / ngam
        out.write(f'{key[par]}\t{key[par]}\t{mgrr:.4f}\t{ngam}\n')
out.close()
