# from os.path import dirname
import os
localrules: make_mbfile, cat_mbfiles
#ruleorder : makebed2_additive, makebed2
# tools
GCTA = '/cluster/home/nkadri/PROGRAMS/GCTA/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 '

split_vcf = '/cluster/home/nkadri/RECOMBINATION/REFALT/IMPUTE/GWAS/split_vcf_with_tabix.py '

# files
OUT_DIR = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/GWAS'
vcf_hd = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/TOHD/{breed}/CHR{chr}/imputed.vcf.gz'
vcf_seq = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/IMPUTE/TOSEQ/{breed}/CHR{chr}/imputed_{sex}.vcf.gz'
# vep = '/cluster/work/pausch/vcf_UCD/2022_03/VEP_filtered_dv/{chr}_vep.vcf.gz'
vep = '/cluster/work/pausch/vcf_UCD/2022_03/vep_beagle/{chr}_vep.vcf.gz'
r2_file = '/cluster/work/pausch/naveen/SNPDATA/TOSEQ/CHR{chr}/r2.txt'


# files required to make phenotype files
# key mapping states file and the ori ids [ids with breed prefix]
states_key = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/HIDDENPHASE/key.txt'
phenotypes = {
    'GRR': '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/LINKPHASE_RUN2/{breed}/cleaned_GRR.txt'
}


# gametes with gp or noff >=6
# good_gams = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/LINKPHASE_RUN2/{breed}/good_gams.txt'


# wildcards
traits = ['GRR']
breeds = ['bv', 'fv']
inheritances = ['additive','ref', 'alt'] 
#inheritances = ['ref', 'alt']
chromosomes = range(1, 30)
sexes = ['female', 'male']
# ngams = [1, 3]  # min number of gametes to get the mean GRR per parent
ngams = [1]
# number of std deviations to filter out GRR outliers (before taking the mean per parent)
stds = [5]


#TEST
breed = ['bv']
sexes = ['male']
inheritances = ['additive']
chromosomes = [25]



rule all:
    input:
        expand(OUT_DIR +'/fastGWA/{trait}/std{std}/gam{ngam}/{inheri}/{breed}/{sex}/CHR{chr}/result.fastGWA.gz', trait=traits, inheri= inheritances, breed=breeds, sex=sexes, chr=chromosomes, std=stds, ngam=ngams)


rule make_phenotypes:
    input:
        states_key,
        lambda wc: phenotypes[wc.trait]
    output:
        pheno = OUT_DIR + \
            '/{trait}/{breed}/{sex}/std{std}/gam{ngam}/phenotypes.txt'
    resources:
        walltime = '00:30',
        mem_mb = 8000
    params:
        noff = 6,  # number of offspring when the grand parent is not genotyped
    script:
        'make_phenotypes.py'


rule makebed:
    input:
        vcf = vcf_hd,
        samples = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/HIDDENPHASE/GWAS/{trait}/{breed}/{sex}/phenotypes.txt'
    output:
        outfiles = temp(expand(
            OUT_DIR + '/GRM/{{breed}}/{{sex}}/chr{{chr}}.{ext}', ext=['bim', 'fam', 'bed']))
    params:
        lambda wc, output: output.outfiles[0][:-4]
    resources:
        mem_mb = 16000,
        walltime = "01:00"
    shell:
        '''
        module load gcc/8.2.0 plink/1.9-beta6.18
        plink --cow \
        --vcf {input.vcf} \
        --make-bed \
        --out {params[0]} \
        --double-id \
        --maf 0.05 \
        --keep {input.samples}
        '''

rule make_grm:
    input:
        infiles = rules.makebed.output.outfiles
    output:
        outfiles = temp(expand(
            OUT_DIR + '/GRM/{{breed}}/{{sex}}/chr{{chr}}.{ext}', ext=['grm.id', 'grm.bin', 'grm.N.bin']))
    params:
        lambda wc, input: input.infiles[0][:-4],
        lambda wc, output: output.outfiles[0][:-7]
    resources:
        mem_mb = 12000,
        walltime = "04:00"
    threads:
        10
    shell:
        '''
        {GCTA} --autosome-num 29 \
        --make-grm-bin \
        --bfile {params[0]} \
         --out {params[1]} \
        --threads {threads}
        '''

rule make_grm_list:
    input:
        infiles = expand(OUT_DIR + '/GRM/{{breed}}/{{sex}}/chr{chr}.{ext}', ext=[
            'grm.id', 'grm.bin', 'grm.N.bin'], chr=chromosomes)
    output:
        outfile = OUT_DIR + \
            "/GRM/{breed}/{sex}/grm_list.txt"
    script:
        "make_grm_list.py"

rule merge_grm:
    input:
        infiles = expand(OUT_DIR + '/GRM/{{breed}}/{{sex}}/chr{chr}.{ext}', ext=[
            'grm.id', 'grm.bin', 'grm.N.bin'], chr=chromosomes),
        grm_list = OUT_DIR + \
            "/GRM/{breed}/{sex}/grm_list.txt"
    output:
        outfiles = expand(
            OUT_DIR + '/GRM/{{breed}}/{{sex}}/genome_cleaned.{ext}', ext=['grm.id', 'grm.bin', 'grm.N.bin'])
    resources:
        mem_mb = 16000,
        walltime = "04:00"
    threads:
        10
    params:
        lambda wc, output: output.outfiles[0][:-7]
    shell:
        '''
        {GCTA} --autosome-num 29 \
         --mgrm {input.grm_list} \
         --make-grm-bin \
         --out {params[0]} \
         --threads {threads}
        '''

rule make_sparse_grm:
    '''
    writes grm.id and grm.sp
    '''
    input:
        grm = rules.merge_grm.output.outfiles
    output:
        outfiles = expand(
            OUT_DIR + '/GRM/{{breed}}/{{sex}}/sparse.{ext}', ext=['grm.id', 'grm.sp'])
    params:
        lambda wc, input: input[0][:-7],
        lambda wc, output: output[0][:-7]
    shell:
        '''
        {GCTA} \
        --grm {params[0]} \
        --make-bK-sparse 0.05 \
        --out {params[1]}
        '''

rule pca:
    input:
        infiles = rules.merge_grm.output.outfiles
    output:
        outfiles = expand(
            OUT_DIR + '/GRM/{{breed}}/{{sex}}/genome.{ext}', ext=['eigenval', 'eigenvec'])
    params:
        lambda wc, input: input.infiles[0][:-7],
        lambda wc, output: output.outfiles[0][:-9],
        num_pc = 4  # number of pcs to include in the association model
    resources:
        mem_mb = 8000,
        walltime = "24:00"
    threads:
        10
    shell:
        '''
        {GCTA} --autosome-num 29 \
        --pca {params.num_pc} \
        --grm {params[0]} \
        --out {params[1]} \
        --threads {threads} \
        '''

rule plot_pca:
    input:
        infile = OUT_DIR + \
            '/GRM/{breed}/{sex}/genome.eigenvec'
    output:
        plotfile = OUT_DIR + \
            '/GRM/{breed}/{sex}/pca_{breed}_{sex}.pdf'
    script:
        "plot_pca.R"


rule index_merged:
    input:
        vcf = vcf_seq
    output:
        index = vcf_seq + '.tbi'
    shell:
        '''
        module load gcc/8.2.0 htslib/1.10.2 \n
        tabix -p vcf  {input.vcf}
        '''

checkpoint split_vcf:
    '''split files have the format chunk_{number}.vcf.gz'''
    input:
        vcf = vcf_seq,
        index = rules.index_merged.output.index
    output:
        dir = directory(OUT_DIR +
                        '/GENOTYPES/{breed}/{sex}/CHR{chr}/SPLIT/')
    params:
        chunk_size = 3000000,
        overlap = 0
    shell:
        '''
        module load gcc/8.2.0 bcftools/1.6 htslib/1.10.2 \n
        python {split_vcf} {input.vcf} {output.dir} --chromosome {wildcards.chr} --chunk_size {params.chunk_size} --overlap {params.overlap}
        '''


fastGWA_mem = {

    'scan': {

        'bv': {
            'male': [2000, 1],
            'female': [1000, 6]
        },

        'fv': {
            'male': [2000, 6],
            'female': [2000, 10]
        }
    },

    'model': {

        'bv': {
            'male': [1000, 5],
            'female': [3000, 5]
        },
        'fv': {
            'male': [1000, 6],
            'female': [32000, 20]
        }
    }
}

rule recode_vcf:
    input:
        infile = OUT_DIR +
        '/GENOTYPES/{breed}/{sex}/CHR{chr}/SPLIT/chunk{chunk}.vcf.gz'
    output:
        outfile = OUT_DIR +
        '/GENOTYPES/{breed}/{sex}/CHR{chr}/{inheri}/chunk{chunk}.vcf'
    wildcard_constraints:
        inheri = "ref|alt" 
    resources:
        mem_mb = 2000,
        walltime = '04:00'
    script:
        'recode_vcf.py'

rule zip_recoded:
    input:
        infile = rules.recode_vcf.output.outfile
    output:
        outfile = rules.recode_vcf.output.outfile + ".gz"
    resources:
        mem_mb = 4000,
        walltime = '01:00'
    shell:
        '''
        module load  gcc/8.2.0 htslib/1.15.1
        bgzip {input.infile}
        '''

rule makebed2_additive:
    input:
        vcf = OUT_DIR +
        '/GENOTYPES/{breed}/{sex}/CHR{chr}/SPLIT/chunk{chunk}.vcf.gz'
    output:
        outfile = OUT_DIR +
        '/GENOTYPES/{breed}/{sex}/CHR{chr}/additive/chunk{chunk}.bed'
    params:
        lambda wc, output: output.outfile[:-4]
    resources:
        mem_mb = 2000,
        walltime = '01:00'
    shell:
        '''
        module load gcc/8.2.0 plink/1.9-beta6.18 \n
        plink --cow \
        --vcf {input.vcf} \
        --make-bed \
        --out {params[0]} \
        --double-id \
        --set-missing-var-ids @_#
        '''


rule makebed2:
    '''
    dirname(str(rules.zip_recoded.output)) + '/chunk{chunk}.bed'
    '''
    input:
        vcf = rules.zip_recoded.output
    output:
        outfile = OUT_DIR +
        '/GENOTYPES/{breed}/{sex}/CHR{chr}/{inheri}/chunk{chunk}.bed'
    wildcard_constraints:
        inheri = "ref|alt" 
    params:
        outfile=lambda wc, output: output.outfile[:-4],
    resources:
        mem_mb = 2000,
        walltime = '01:00'
    shell:
        '''
        module load gcc/8.2.0 plink/1.9-beta6.18 \n
        plink --cow \
        --vcf {input.vcf} \
        --make-bed \
        --out {params.outfile} \
        --double-id \
        --set-missing-var-ids @_#
        '''


def list_binaries(wc):
    out_dir = checkpoints.split_vcf.get(
        breed=wc.breed, sex=wc.sex, chr=wc.chr).output[0]
    chunks = glob_wildcards(out_dir + "/chunk{chunks}.vcf.gz").chunks
    chunks = sorted([int(el) for el in chunks])
    print(f'chunks are  {chunks}')
    outfiles = [
        OUT_DIR + f'/GENOTYPES/{wc.breed}/{wc.sex}/CHR{wc.chr}/{wc.inheri}/chunk{chunk}.bed' for chunk in chunks]
    return outfiles


rule make_mbfile:
    input:
        infiles = list_binaries
    output:
        mbfile = OUT_DIR +
        '/fastGWA/{trait}/std{std}/gam{ngam}/{inheri}/{breed}/{sex}/mbfiles_chr{chr}.txt'
    script:
        'make_mbfile.py'

rule cat_mbfiles:
    input:
        infiles = expand(OUT_DIR +
                         '/fastGWA/{{trait}}/std{{std}}/gam{{ngam}}/{{inheri}}/{{breed}}/{{sex}}/mbfiles_chr{chr}.txt', chr=chromosomes)
    output:
        outfile = OUT_DIR +
        '/fastGWA/{trait}/std{std}/gam{ngam}/{inheri}/{breed}/{sex}/mbfiles.txt'
    shell:
        '''
        cat {input.infiles} > {output.outfile}
        '''

rule mean_center:
    input:
        pheno = rules.make_phenotypes.output.pheno,
    output:
        pheno = OUT_DIR + '/{trait}/{breed}/{sex}/std{std}/gam{ngam}/phenotypes_centered.txt'
    script:
        'mean_center.R'
        
rule fastGWAS_model_only:
    '''
    --fastGWA-mlm-exact or --fastGWA-mlm
    both work only with sparse!
    geno : either whole genome or only one chromosomes .. markers are randomly sampled
    '''
    input:
        pheno=rules.mean_center.output.pheno,
        pc = rules.pca.output.outfiles[1],
        grm = rules.make_sparse_grm.output.outfiles,
        geno = OUT_DIR +'/fastGWA/{trait}/std{std}/gam{ngam}/{inheri}/{breed}/{sex}/mbfiles_chr25.txt'
        #geno = rules.cat_mbfiles.output.outfile
    output:
        id = OUT_DIR + '/fastGWA/{trait}/std{std}/gam{ngam}/{inheri}/{breed}/{sex}/model.fastGWA.mdl.id',
        bin = OUT_DIR + '/fastGWA/{trait}/std{std}/gam{ngam}/{inheri}/{breed}/{sex}/model.fastGWA.mdl.bin',
    params:
        grm=lambda wc, input: input.grm[0][:-7],
        out=lambda wc, output: output[0][:-15]
    resources:
        mem_mb = lambda wc :fastGWA_mem['model'][wc.breed][wc.sex][0],
        walltime = '120:00:00'
    threads:
        lambda wc :fastGWA_mem['model'][wc.breed][wc.sex][1],
    shell:
        '''
        {GCTA} --autosome-num 29 \
        --mbfile {input.geno} \
        --grm-sparse {params.grm} \
        --fastGWA-mlm-exact \
        --model-only \
        --pheno {input.pheno} \
        --qcovar {input.pc} \
        --thread-num {threads} \
        --out {params.out}
        '''
# --bfile {params[1]} \
rule fastGWA:
    input:
        model = rules.fastGWAS_model_only.output,
        mbfile = rules.make_mbfile.output.mbfile,
        pheno = rules.make_phenotypes.output.pheno
    output:
        outfile = OUT_DIR +
        '/fastGWA/{trait}/std{std}/gam{ngam}/{inheri}/{breed}/{sex}/CHR{chr}/result.fastGWA'
    params:
        lambda wc, input: input.model[0][:-7],
        lambda wc, output: output.outfile[:-8]
    resources:
        mem_mb = lambda wc :fastGWA_mem['scan'][wc.breed][wc.sex][0],
        walltime = '24:00'
    threads:
        lambda wc :fastGWA_mem['scan'][wc.breed][wc.sex][1]
    shell:
        '''
        {GCTA} --autosome-num 29 \
        --load-model {params[0]} \
        --mbfile {input.mbfile} \
        --thread-num {threads} \
        --out {params[1]}
        '''
rule zip:
    input:
        infile=rules.fastGWA.output.outfile
    output:
        outfile=rules.fastGWA.output.outfile + ".gz"
    shell:
        'gzip {input.infile}'
