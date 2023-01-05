from Bio import SeqIO
import gzip
from ete3 import NCBITaxa
from snakemake.utils import min_version
min_version("7.6.2")

rule all:
    input:
         "db/taxonomy",
         "db/ete3/taxa.sqlite",
         "db/uniref/uniref100.fasta.gz",
         "db/kaiju/uniref100.fmi",
         "db/kaiju/uniref100.ktaxonomy",


rule taxonomy:
    output:
          db = directory("db/taxonomy"),
          taxdump = "db/taxonomy/new_taxdump.tar.gz",
          nodes = "db/taxonomy/nodes.dmp",
          names = "db/taxonomy/names.dmp",
          taxidlineage = "db/taxonomy/taxidlineage.dmp",
          taxids = "db/taxonomy/ncbi.valid.taxids",
    threads: 1
    params:
        taxdump = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz",
    log: "logs/makedb/taxonomy.log"
    shell: """
            aria2c -x 16 -s 16 -j 16 -o {output.taxdump} {params.taxdump} &> {log}
            tar -C {output.db} -xvzf {output.taxdump} &>> {log}
            grep -w '2157\|2\|2759\|10239' {output.taxidlineage} | cut -f 1 > {output.taxids}
    """

rule ete3:
    input:
        taxdump=rules.taxonomy.output.taxdump
    output:
        sqlite = "db/ete3/taxa.sqlite",
        pkl = "db/ete3/taxa.sqlite.traverse.pkl"
    threads: 1
    run: NCBITaxa(dbfile=output.sqlite, taxdump_file=input.taxdump)

rule uniref:
    output:
        fa = "db/uniref/uniref100.fasta.gz",
        idmapping = "db/uniref/idmapping.dat.gz",
        version = "db/uniref/uniref100.version",
    threads: 1
    params:
        version = "https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.release_note",
        fa = "https://ftp.expasy.org/databases/uniprot/current_release/uniref/uniref100/uniref100.fasta.gz",
        idmapping = "https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz",
    log: "logs/makedb/uniref100.log"
    shell: """
            aria2c -o {output.version} {params.version} &> {log}
            aria2c -x 16 -j 16 -s 16 -o {output.idmapping} {params.idmapping} &>> {log}
            aria2c -x 16 -j 16 -s 16 -o {output.fa} {params.fa} &>> {log}
           """

rule uniref_acc2taxid:
    input:
        taxids = rules.taxonomy.output.taxids,
        idmapping = rules.uniref.output.idmapping,
    output:
        acc2tax = "db/kaiju/uniref100.acc2tax"
    threads: 8
    shell: """
            zgrep -Fwf <(zgrep -w "UniRef100" {input.idmapping} | cut -f 3 | sort -u -S 32G --parallel={threads} | cut -f 2 -d "_") {input.idmapping} | \
            grep -w 'NCBI_TaxID' | grep -Fwf {input.taxids} | awk -v FS='\t' -v OFS='\t' '{{print "UniRef100_"$1,$3}}' > {output.acc2tax}
    """

rule uniref_reformat:
    input:
        fa = rules.uniref.output.fa,
        acc2tax = rules.uniref_acc2taxid.output.acc2tax
    output:
        fa = temp("db/kaiju/uniref100.fa"),
    threads: 1
    run:
        with open(input.acc2tax, "r") as a2t:
            acc2tax = dict(line.split("_")[1].split() for line in a2t)
        with gzip.open(input.fa, "rt") as infa, open(output.fa, "w") as outfa:
            for rec in SeqIO.parse(infa, "fasta"):
                rec.id = rec.id.split("_")[1]
                taxid = acc2tax.get(rec.id)
                if taxid is not None:
                    rec.id = rec.id + "_" + taxid
                    rec.description = ""
                    SeqIO.write(rec, outfa, "fasta")

rule kaiju:
    input:
        fa = rules.uniref_reformat.output.fa,
    output:
        bwt = temp("db/kaiju/uniref100.bwt"),
        sa = temp("db/kaiju/uniref100.sa"),
        fmi = "db/kaiju/uniref100.fmi"
    threads: 40
    conda: "envs/kaiju.yaml"
    params:
        db = "db/kaiju/uniref100"
    log: "logs/makedb/kaiju_uniref.log"
    shell: """
            kaiju-mkbwt -o {params.db} -n {threads} {input.fa} &> {log}
            kaiju-mkfmi {params.db} &>> {log}
           """

rule kraken_taxonomy:
    input:
        names = rules.taxonomy.output.names,
        nodes = rules.taxonomy.output.nodes,
        acc2tax = rules.uniref_acc2taxid.output.acc2tax
    output:
        ktaxonomy = "db/kaiju/uniref100.ktaxonomy"
    threads: 1
    conda: "envs/krakentools.yaml"
    log: "logs/makedb/kraken_taxonomy.out"
    shell: """
            make_ktaxonomy.py --nodes {input.nodes} --names {input.names} --seqid2taxid {input.acc2tax} --output {output.ktaxonomy} &> {log}
           """