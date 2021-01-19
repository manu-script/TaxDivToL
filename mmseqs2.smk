import pandas as pd
import numpy as np
import re
import os

meta = pd.read_csv("metadata.tsv", sep='\t', header=0, index_col="sample")
split = 100000000

rule all:
    input:
        "db/nt",
        "db/nt_raw.stats",
        "db/nt_rep.stats",
        "db/nt_mapping",
        "db/nt.idx",
        expand("fastp/{sample}", sample=meta.index.unique()),
        expand("mmseqs/{sample}", sample=meta.index.unique()),


rule dwnld_nt:
    output: temp("db/nt.gz")
    threads: 1
    log:
        stderr = "logs/dwnld_nt.stderr",
        stdout = "logs/dwnld_nt.stdout"
    shell: "aria2c -x 16 -j 16 -s 16 -o {output} https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz 1> {log.stdout} 2> {log.stderr}"


rule nt_stats:
    input: rules.dwnld_nt.output
    output: "db/nt_raw.stats"
    threads: 2
    log: "logs/nt_stats.stderr"
    shell: "seqkit stats -j {threads} -a -o {output} {input} 2> {log}"


rule nt_seqs:
    input: rules.dwnld_nt.output
    output:
        fa = temp("db/nt_seqs.fasta"),
        stats = "db/nt_seqs.stats"
    threads: 2
    log: "logs/nt_seqs.stderr"
    shell: "seqkit seq -m 100 -M 1000000 -j {threads} {input} 2> {log} | seqkit grep -n -p 'satellite,repeat,repetitive,transpos,intersperse,tandem' -r -v -j {threads} | dustmasker -outfmt fasta | tee {output.fa} | seqkit stats -a -j {threads} -o {output.stats}"


rule easy_linclust:
    input: rules.nt_seqs.output.fa
    output:
        tmp = temp(directory("db/tmp")),
        all_seqs = temp("db/nt_all_seqs.fasta"),
        cluster = temp("db/nt_cluster.tsv"),
        rep_seq = temp("db/nt_rep_seq.fasta")
    threads: 40
    log:
       stderr = "logs/easy_linclust.stderr",
       stdout = "logs/easy_linclust.stdout"
    shell: "mmseqs easy-linclust {input} db/nt {output.tmp} --mask 1 --mask-lower-case 1 --min-seq-id 0.99 --cluster-mode 2 --cov-mode 1 -c 0.8 --threads {threads} 1> {log.stdout} 2> {log.stderr}"


rule rep_stats:
    input: rules.easy_linclust.output.rep_seq
    output: "db/nt_rep.stats"
    threads: 2
    log: "logs/rep_stats.stderr"
    shell: "seqkit stats -j {threads} -a -o {output} {input} 2> {log}"


rule createdb:
    input: rules.easy_linclust.output.rep_seq
    output: "db/nt"
    threads: 1
    log:
        stderr = "logs/createdb.stderr",
        stdout = "logs/createdb.stdout"
    shell: "mmseqs createdb {input} {output} 1> {log.stdout} 2> {log.stderr}"


rule createindex:
    input: rules.createdb.output
    output:
        dbindex = "db/nt.idx",
        tmpdir = temp(directory("db/tmpdir"))
    threads: 40
    log:
        stderr = "logs/createindex.stderr",
        stdout = "logs/createindex.stdout"
    shell: "mmseqs createindex {input} {output.tmpdir} --mask 1 --mask-lower-case 1 --search-type 3 --threads {threads} 1> {log.stdout} 2> {log.stderr}"


rule dwnld_tax:
    output:
        targz = temp("db/taxonomy/new_taxdump.tar.gz"),
        taxdump = temp(directory("db/taxonomy"))
    threads: 1
    log:
        stderr = "logs/dwnld_tax.stderr",
        stdout = "logs/dwnld_tax.stdout"
    shell: "aria2c -x 16 -j 16 -s 16 -o {output.targz} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz 1> {log.stdout} 2> {log.stderr} && tar -xvzf {output.targz} -C {output.taxdump}"


rule dwnld_taxmap:
    output:
          acc2tax = temp("db/nucl_gb.accession2taxid.gz"),
          taxmap = temp("db/nt.taxmap")
    threads: 1
    log:
        stderr = "logs/dwnld_taxmap.stderr",
        stdout = "logs/dwnld_taxmap.stdout"
    shell: "aria2c -x 16 -j 16 -s 16 -o {output.acc2tax} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz 1> {log.stdout} 2> {log.stderr} && gzip -cd {output.acc2tax} | cut -f 2,3 --output-delimiter ' ' 1> {output.taxmap}"


rule createtaxdb:
    input:
        db = rules.createdb.output,
        taxdump = rules.dwnld_tax.output.taxdump,
        taxmap = rules.dwnld_taxmap.output.taxmap
    output:
        tmp = temp(directory("db/tmp")),
        db = "db/nt_mapping"
    threads: 40
    log:
        stderr = "logs/createtaxdb.stderr",
        stdout = "logs/createtaxdb.stdout"
    shell: "mmseqs createtaxdb {input.db} {output.tmp} --ncbi-tax-dump {input.taxdump} --tax-mapping-file {input.taxmap} --threads {threads} 1> {log.stdout} 2> {log.stderr}"


rule fastq_qc:
    input:
        read1 = lambda wildcards: meta.loc[wildcards.sample, "read1"],
        read2 = lambda wildcards: meta.loc[wildcards.sample, "read2"]
    output:
        reads = temp(directory("fastp/{sample}")),
        html = "fastp/{sample}.html",
        json = "fastp/{sample}.json"
    threads: 8
    params:
        size = split
    log:
        fastp = "logs/fastq_qc/{sample}.log",
        seqkit = "logs/fastq_split/{sample}.log"
    shell: "fastp -i {input.read1} -I {input.read2} -w {threads} -j {output.json} -h {output.html} -l 35 -g -y -5 -3 --detect_adapter_for_pe --stdout 2> {log.fastp} | seqkit split2 -s {params.size} -O {output.reads} -j {threads} 2> {log.seqkit}"


rule easy_taxonomy:
    input:
        db = rules.createdb.output,
        dbindex = rules.createindex.output.dbindex,
        taxdb = rules.createtaxdb.output.db,
        reads = rules.fastq_qc.output.reads
    output:
        out = directory("mmseqs/{sample}"),
        tmp = temp(directory("mmseqs/{sample}/tmp")),
    threads: 40
    log:
        stderr = "logs/easy_taxonomy/{sample}/stderr",
        stdout = "logs/easy_taxonomy/{sample}/stdout"
    run:
        parts={}
        for f in os.listdir(input.reads):
            if f.endswith('.fastq'):
                parts[int(re.findall(r'part_(\d+).fastq', f)[0])] = os.path.join(input.reads, f)
        for p in sorted(parts)[:-1]:
            p = str(p)
            shell("mmseqs easy-taxonomy "+parts[int(p)]+" {input.db} {output.out}/part"+p+" {output.tmp} --exact-kmer-matching 1 -s 1.0 --max-seqs 10 --mask 1  --mask-lower-case 1 -e 0.00001 --min-seq-id 0.6 -c 0.8 --cov-mode 2 --tax-lineage 1 --search-type 3 --lca-mode 2 --threads {threads} 1> {log.stdout}.part"+p+" 2> {log.stderr}.part"+p)
            shell("awk '$2 != unclassified' {output.out}/part"+p+"_lca.tsv | gzip -c 1> {output.out}/part"+p+"_lca.tsv.gz")
            shell("rm {output.out}/part"+p+"_lca.tsv {output.out}/part"+p+"_report {output.out}/part"+p+"_tophit_report {output.out}/part"+p+"_tophit_aln "+parts[int(p)])

