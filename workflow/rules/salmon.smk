from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()


rule genome:
    input:
        FTP.remote(config['salmon']['genome'])
    output:
        'resources/genome.fa.gz'
    run:
        shell('mv {input} {output}')


rule transcriptome:
    input:
        FTP.remote(config['salmon']['transcriptome'])
    output:
        'resources/transcriptome.fa.gz'
    run:
        shell('mv {input} {output}')


rule genesets:
    input:
        FTP.remote(config['salmon']['genesets'])
    output:
        'resources/genesets.gtf.gz'
    run:
        shell('mv {input} {output}')


rule decoys:
    input:
        'resources/genome.fa.gz'
    output:
        'resources/decoys.txt'
    shell:
        '''
        grep "^>" <(gunzip -c {input}) |
        cut -d " " -f 1 |
        sed -e 's/>//g' > {output}
        '''


rule concatenate_genome_transcriptome:
    input:
        genome='resources/genome.fa.gz',
        transcriptome='resources/transcriptome.fa.gz'
    output:
        'resources/gentrome.fa.gz'
    shell:
        '''
        cat {input.transcriptome} {input.genome} > {output}
        '''


rule salmon_index:
    input:
        gentrome='resources/gentrome.fa.gz',
        decoys='resources/decoys.txt'
    output:
        directory('resources/salmon_index')
    params:
        threads=config['salmon']['threads']
    conda:
        "../envs/salmon.yaml"
    threads: config['salmon']['threads']
    resources:
        mem_free_gb=f"{config['salmon']['memory_per_cpu']}"
    log: 'resources/salmon_index.log'
    shell:
        '''
        salmon index -t {input.gentrome} -d {input.decoys} -p {params.threads} -i {output} 2> {log}
        '''


rule transcript_gene_map:
    input:
        gtf='resources/genesets.gtf.gz'
    output:
        tgmap='resources/txp2gene.tsv'
    params:
        renv=config['renv']
    log: script="results/logs/transcript_gene_map.log"
    script:
        "../scripts/txp2gene.R"
