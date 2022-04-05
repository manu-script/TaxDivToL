import pandas as pd
import os
import itertools
import random
from snakemake.utils import min_version
min_version("6.10")

meta = pd.read_csv("metadata.tsv", sep='\t', header=0, index_col="sample", comment='#')
taxa = ["Tree of Life", "Archaea", "Bacteria", "Eukaryota", "Viruses", "Protists", "Metazoa", "Fungi", "Viridiplantae", "Streptophyta", "Arthropoda", "Actinopteri", "Mammalia"]

rule all:
    input:
        expand("data/qc/{sample}_R1.fq.gz", sample=meta.index.unique()),
        expand("diversity/centrifuge/{sample}.classification.out.gz", sample=meta.index.unique()),
        expand("diversity/centrifuge/{sample}.classification.kreport", sample=meta.index.unique()),
        "diversity/counts/combined.classification.mpa",
        "diversity/counts/combined.filtered.mpa",
        "diversity/counts/combined.filtered.counts.tsv",
        expand("diversity/counts/combined.filtered.{taxon}.tsv", taxon=["Actinopteri"]),
        "diversity/counts/combined.filtered.richness.krona",
        "diversity/counts/combined.filtered.richness.html",
        expand("diversity/counts/{sample}.filtered.mpa", sample=meta.index.unique()),
        expand("diversity/richness/{sample}.richness.tsv", sample=meta.index.unique()),
        expand("diversity/richness/{sample}.richness.png", sample=["S29DEC19"]),
        "diversity/richness/combined.richness.cv.tsv"


rule clumpify:
    input:
         r1 = lambda wildcards: meta.loc[wildcards.sample, "read1"],
         r2 = lambda wildcards: meta.loc[wildcards.sample, "read2"]
    output:
         r1 = temp("data/qc/{sample}_R1_clumped.fq.gz"),
         r2 = temp("data/qc/{sample}_R2_clumped.fq.gz")
    threads: 40
    log: "logs/diversity/clumpify/{sample}.log"
    conda: "envs/bbtools.yaml"
    shell: "clumpify.sh in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} reorder=p overwrite=t shortname=t dedupe=t ziplevel=1 -Xmx800g -eoom &"
    "> {log}"

rule duk:
    input:
        r1 = rules.clumpify.output.r1,
        r2 = rules.clumpify.output.r2
    output:
        r1 = "data/qc/{sample}_R1.fq.gz",
        r2 = "data/qc/{sample}_R2.fq.gz",
        unpaired = "data/qc/{sample}_unpaired.fq.gz"
    threads: 6
    log: "logs/diversity/duk/{sample}.log"
    conda: "envs/bbtools.yaml"
    shell: "bbduk.sh in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outs={output.unpaired} ref=adapters,artifacts showspeed=f ziplevel=2 k=31 ktrim=r qtrim=rl trimq=10 mink=11 minlength=51 maxns=3 minavgquality=10 trimpolygright=5 forcetrimmod=5 -Xmx128g -eoom &> {log}"

rule centrifuge:
    input:
         db = "db/centrifuge",
         r1 = rules.duk.output.r1,
         r2 = rules.duk.output.r2,
         unp = rules.duk.output.unpaired
    output:
         out = "diversity/centrifuge/{sample}.classification.out.gz",
         report = "diversity/centrifuge/{sample}.report.tsv",
         r1 = "diversity/centrifuge/{sample}.classified.1.fq.gz",
         r2 = "diversity/centrifuge/{sample}.classified.2.fq.gz"
    threads: 40
    params:
        db = "db/centrifuge/nt",
        reads = "diversity/centrifuge/{sample}.classified.fq",
        r1 = "diversity/centrifuge/{sample}.classified.1.fq",
        r2 = "diversity/centrifuge/{sample}.classified.2.fq"
    log: "logs/diversity/centrifuge/{sample}.out"
    shell: """
            set -eo pipefail
            tools/centrifuge-1.0.4/centrifuge -k 1 -p {threads} --al-conc-gz {params.reads} --min-hitlen 35 -x {params.db} -1 {input.r1} -2 {input.r2} --report-file {output.report} 2> {log} | bgzip -c -@ 6 > {output.out}
            mv {params.r1} {output.r1}
            mv {params.r2} {output.r2}
    """

rule report:
    input:
        db = "db/centrifuge",
        classification = rules.centrifuge.output.out
    output:
        kreport = "diversity/centrifuge/{sample}.classification.kreport",
        mpa = "diversity/centrifuge/{sample}.classification.mpa",
        splitdir = directory("diversity/centrifuge/{sample}.classification.split"),
    threads: 1
    params:
        split = "diversity/centrifuge/{sample}.classification.split/part.",
        db = "db/centrifuge/nt"
    conda: "envs/krakentools.yaml"
    log: "logs/diversity/report/{sample}.out"
    shell: """
            set +e
            header=$(zcat {input.classification} | head -1)
            mkdir -p {output.splitdir}
            bgzip -cd -@ 6 {input.classification} | tail -n +2 | \
            split -d -l 100000000 --additional-suffix=.out - {params.split} &> {log}

            for file in {output.splitdir}/part.*.out
            do
                tools/centrifuge-1.0.4/centrifuge-kreport -x {params.db} <(echo -e "$header"; cat $file) 1> {output.splitdir}/$(basename ${{file/.out/.kreport}}) 2>> {log}
                kreport2mpa.py -r {output.splitdir}/$(basename ${{file/.out/.kreport}}) -o {output.splitdir}/$(basename ${{file/.out/.mpa}}) &>> {log}
                sed -i 's/^k__/d__/' {output.splitdir}/$(basename ${{file/.out/.mpa}})
                rm $file
            done
            
            combine_kreports.py -r {output.splitdir}/part.*.kreport --no-headers --only-combined -o {output.kreport} &>> {log}
            kreport2mpa.py -r {output.kreport} -o {output.mpa} &>> {log}
            sed -i 's/^k__/d__/' {output.mpa}
           """

rule combine_reports:
    input:
        kreport = expand("diversity/centrifuge/{sample}.classification.kreport", sample=meta.index.unique())
    output:
        kreport = "diversity/counts/combined.classification.kreport",
        mpa = "diversity/counts/combined.classification.mpa"
    threads: 1
    log: "logs/diversity/combine_reports.log"
    conda: "envs/krakentools.yaml"
    shell: """
            combine_kreports.py -r {input.kreport} --no-headers --only-combined -o {output.kreport} &> {log}
            kreport2mpa.py -r {output.kreport} -o {output.mpa} &>> {log}
            sed -i 's/^k__/d__/' {output.mpa}
           """

rule filter_combined:
    input:
        mpa = rules.combine_reports.output.mpa,
    output:
        mpa = "diversity/counts/combined.filtered.mpa"
    threads: 1
    run:
        df = pd.read_csv(input.mpa, sep='\t', header=None, index_col=0, names=["lineage","count"])
        df = df[df.index.str.startswith("d__")]
        df["parent1_count"] = df.index.map(lambda x: df.loc[x,"count"] if "|".join(x.split("|")[:-1]) == '' else df.loc["|".join(x.split("|")[:-1]),"count"])
        df["parent2_count"] = df.index.map(lambda x: df.loc[x,"count"] if "|".join(x.split("|")[:-2]) == '' else df.loc["|".join(x.split("|")[:-2]),"count"])
        df["parent1_relabund"] = (df["count"]/df["parent1_count"])*100
        df["parent2_relabund"] = (df["count"]/df["parent2_count"])*100
        df = df[df.index.map(lambda x: x.split("|")[-1].startswith("f__"))]
        df = df[(df["parent2_relabund"]>=0.1) & (df["parent1_relabund"]>=1) & (df["count"]>=10)]
        df = df.drop(columns=["parent1_count", "parent1_relabund", "parent2_count", "parent2_relabund"])
        df = df.reset_index()
        df.to_csv(output.mpa, columns=["lineage", "count"], sep="\t", header=False, index=False)

rule prep_krona:
    input:
        mpa = rules.filter_combined.output.mpa
    output:
        richness = "diversity/counts/combined.filtered.richness.krona",
        abundance = "diversity/counts/combined.filtered.abundance.krona"
    threads: 1
    run:
        df = pd.read_csv(input.mpa, sep='\t', header=None, index_col=1, names=["lineage","total"])
        df = df["lineage"].str.split("|", expand=True).fillna("")
        df.to_csv(output.richness, header=False, sep= "\t", index=False)
        df.to_csv(output.abundance, header=False, sep= "\t", index=True)

rule krona:
    input:
        richness = rules.prep_krona.output.richness,
        abundance = rules.prep_krona.output.abundance
    output:
        richness = "diversity/counts/combined.filtered.richness.html",
        abundance = "diversity/counts/combined.filtered.abundance.html"
    threads: 1
    conda: "envs/krakentools.yaml"
    log: "logs/diversity/taxonomy/krona.log"
    shell: """
            ktImportText -q -n 'Tree of Life' -o {output.richness} {input.richness} &> {log}
            ktImportText -n 'Tree of Life' -o {output.abundance} {input.abundance} &>> {log}
           """

rule combine_filtered_counts:
    input:
        S27 = "diversity/centrifuge/S27DEC19.classification.mpa",
        S28 = "diversity/centrifuge/S28DEC19.classification.mpa",
        S29 = "diversity/centrifuge/S29DEC19.classification.mpa",
        combined = rules.filter_combined.output.mpa
    output:
        counts = "diversity/counts/combined.filtered.counts.tsv"
    threads: 1
    run:
        S27 = pd.read_csv(input.S27, sep='\t', header=None, names=["lineage","S27DEC19"])
        S28 = pd.read_csv(input.S28, sep='\t', header=None, names=["lineage","S28DEC19"])
        S29 = pd.read_csv(input.S29, sep='\t', header=None, names=["lineage","S29DEC19"])
        combined = pd.read_csv(input.combined, sep='\t', header=None, names=["lineage","Total count"])
        combined = combined.merge(S27, on="lineage", how="left").merge(S28, on="lineage", how="left").merge(S29, on="lineage", how="left").fillna(0)
        combined.to_csv(output.counts, header=True, index=False, sep="\t")

rule compare_checklist:
    input:
        counts = rules.combine_filtered_counts.output.counts,
        checklist = "data/{taxon}_families.tsv",
        genomelist = "data/{taxon}_annotated_genomes.tsv"
    output:
        tsv = "diversity/counts/combined.filtered.{taxon}.tsv"
    threads: 1
    run:
        checklist = "|".join(pd.read_csv(input.checklist, header=None, squeeze=True).tolist())
        genomelist = "|".join(pd.read_csv(input.genomelist, sep="\t", header=0, usecols=["family"], squeeze=True).tolist())
        df = pd.read_csv(input.counts, sep='\t', header=0)
        df = df[df["lineage"].str.contains(wildcards.taxon)]
        df["known_taxon"] = df["lineage"].str.contains(checklist)
        df["annotated_genome"] = df["lineage"].str.contains(genomelist)
        df.to_csv(output.tsv, sep="\t", header=True, index=False)

rule filter_sample:
    input:
        combined = rules.filter_combined.output.mpa,
        sample = rules.report.output.mpa
    output:
        mpa = "diversity/counts/{sample}.filtered.mpa"
    threads: 1
    run:
        combined = pd.read_csv(input.combined, sep='\t', header=None, usecols=[0], names=["lineage"])
        sample = pd.read_csv(input.sample, sep='\t', header=None,names=["lineage", "count"])
        sample = sample.merge(combined, on="lineage", how="inner")
        sample.to_csv(output.mpa, header=False, sep="\t", index=False)

def filter_taxon(dat, taxon):
    if taxon == "Protists":
        dat = dat[(dat["lineage"].str.contains("Eukaryota")) & (~dat["lineage"].str.contains("Fungi|Viridiplantae|Metazoa"))]
    elif taxon == "Tree of Life":
        dat = dat[(dat["lineage"].str.contains("Viruses|Archaea|Bacteria|Eukaryota"))]
    else:
        dat = dat[dat["lineage"].str.contains(taxon)]
    return dat

rule richness:
    input:
        dir = rules.report.output.splitdir,
        filtered = rules.filter_sample.output.mpa
    output:
        stats = "diversity/richness/{sample}.richness.tsv"
    threads: 1
    run:
        filtered = pd.read_csv(input.filtered, sep='\t', header=None, names=["lineage", "total_count"])
        parts = sorted([file for file in os.listdir(input.dir) if file.endswith(".mpa")])[:-1]
        stats = []
        for taxon in taxa:
            filt = filter_taxon(filtered, taxon)
            for depth in range(len(parts)):
                r = []
                cnt = 0
                if depth!=0:
                    for replicate in random.sample(list(itertools.combinations(parts, depth)), k=3):
                        cnt += 1
                        print(taxon, depth, cnt)
                        dfs = []
                        for part in replicate:
                            dfs.append(pd.read_csv(os.path.join(input.dir, part), sep='\t', header=None, names=["lineage", "count"]))
                        df = pd.concat(dfs).groupby("lineage").sum()
                        df = df.reset_index()
                        df = df[df["count"] >= 10]
                        df = df.merge(filt, on="lineage", how="inner")
                        r.append(len(df))
                    r = pd.Series(r)
                    depth *= 100
                    richness = round(r.mean(),2)
                    richness_std = round(r.std(),2)
                    total = len(filt)
                    coverage = round((richness/total)*100, 2)
                    coverage_std = round((richness_std/total)*100, 2)
                    stats.append((wildcards.sample, taxon, depth, richness, richness_std, coverage, coverage_std))
                else:
                    stats.append((wildcards.sample, taxon, depth, 0, 0, 0, 0))
        stats = pd.DataFrame(stats, columns=["sample", "taxon", "depth", "richness", "richness_std", "coverage", "coverage_std"])
        stats.to_csv(output.stats, sep="\t", header=True, index=False)

rule plot_richness_coverage:
    input:
        stats = rules.richness.output.stats
    output:
        rich = "diversity/richness/{sample}.richness.png",
        cov = "diversity/richness/{sample}.coverage.png"
    threads: 1
    conda: "envs/ggplot.yaml"
    log: "logs/diversity/plot_richness/{sample}.log"
    shell: "Rscript scripts/plot.richness.coverage.R --stats {input.stats} --rich {output.rich} --cov {output.cov} &> {log}"

rule richness_cv:
    input:
        S27DEC19 = "diversity/richness/S27DEC19.richness.tsv",
        S28DEC19 = "diversity/richness/S28DEC19.richness.tsv",
        S29DEC19 = "diversity/richness/S29DEC19.richness.tsv",
    output:
        cv =  "diversity/richness/combined.richness.cv.tsv"
    threads: 1
    conda: "envs/ggplot.yaml"
    log: "logs/diversity/richness_cv.log"
    shell: "Rscript scripts/calc.richness.cv.R --S27DEC19 {input.S27DEC19} --S28DEC19 {input.S28DEC19} --S29DEC19 {input.S29DEC19} --cv {output.cv} &> {log}"
