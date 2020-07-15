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
        rules.genome.output
    output:
        'resources/decoys.txt'
    shell:
        '''
        grep "^>" <(gunzip -c {input}) |
        cut -d " " -f 1 |
        sed -e 's/>//g' > {output}
        '''


rule gentrome:
    input:
        genome=rules.genome.output,
        transcriptome=rules.transcriptome.output
    output:
        'resources/gentrome.fa.gz'
    shell:
        '''
        cat {input.transcriptome} {input.genome} > {output}
        '''


rule salmon_index:
    input:
        gentrome=rules.gentrome.output,
        decoys=rules.decoys.output
    output:
        directory('resources/salmon_index')
    params:
        threads=config['salmon']['threads']
    conda:
        "../envs/salmon.yaml"
    threads: config['salmon']['threads']
    resources:
        mem_free_gb=f"{config['salmon']['memory_per_cpu']}"
    log:
        out='results/logs/salmon_index/out',
        err='results/logs/salmon_index/err',
        time='results/logs/time/salmon_index'
    shell:
        '''
        {DATETIME} > {log.time} &&
        salmon index -t {input.gentrome} -d {input.decoys} -p {params.threads} -i {output} 2> {log} &&
        {DATETIME} >> {log.time}
        '''


rule transcript_gene_map:
    input:
        gtf=rules.genesets.output
    output:
        tgmap='resources/txp2gene.tsv'
    conda:
        "../envs/r.yaml"
    log: script="results/logs/transcript_gene_map.log"
    script:
        "../scripts/txp2gene.R"
