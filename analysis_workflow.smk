import pandas as pd
from ete3 import NCBITaxa
from functools import reduce
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
sns.set_theme(style="whitegrid")
import itertools
from snakemake.utils import min_version
min_version("7.6.2")

meta = pd.read_csv("metadata.tsv", sep='\t', header=0, index_col="sample", comment='#')
pilot = ["S27DEC19", "S28DEC19", "S29DEC19"]
seasonal = ["S27DEC19", "S28DEC19", "S29DEC19", "S1MAR20", "S6MAR20", "S26MAR20", "S1JUL20", "S6JUL20", "S14JUL20", "S17JUL20", "S26JUL20", "S29JUL20", "S30JUL20", "S6NOV20", "S14NOV20", "S17NOV20"]
taxa = ["ToL","Viruses","Archaea","Bacteria","Eukaryota","Metazoa","Fungi","Viridiplantae","Arthropoda","Chordata","Actinopteri","Mammalia","Aves","Amphibia"]


rule all:
    input:
        "diversity/taxonomy/supplementary/uniref100.families.png",
        "diversity/taxonomy/supplementary/samples.metaplot.png",
        expand("diversity/taxonomy/khist/{sample}.hist.tsv", sample=pilot),
        expand("diversity/taxonomy/preseq/{sample}.extrap.tsv", sample=pilot),
        "diversity/taxonomy/preseq/kmer.saturation.curve.png",
        expand("data/qc/{sample}_R1.fq.gz", sample=meta.index.unique()),
        expand("diversity/taxonomy/kaiju/{sample}.classification.R1.tsv.gz", sample=meta.index.unique()),
        expand("diversity/taxonomy/kaiju/{sample}.lca.tsv.gz", sample=meta.index.unique()),
        expand("diversity/taxonomy/counts/{sample}.mpa", sample=meta.index.unique()),
        "diversity/taxonomy/supplementary/classified.read.proportions.png",
        "diversity/taxonomy/counts/combined.mpa",
        "diversity/taxonomy/supplementary/classified.rank.proportions.png",
        "diversity/taxonomy/counts/combined.filtered.mpa",
        "diversity/taxonomy/counts/combined.filtered.richness.html",
        expand("diversity/taxonomy/counts/{sample}.filtered.mpa", sample=meta.index.unique()),
        "diversity/taxonomy/counts/combined.filtered.counts.matrix.tsv",
        expand("diversity/taxonomy/checklist/{taxon}.checklist.tsv", taxon="Actinopteri"),
        expand("diversity/taxonomy/checklist/{taxon}.relabund.png", taxon="Actinopteri"),
        expand("diversity/taxonomy/incidences/{taxon}.freq.mpa", taxon=taxa),
        expand("diversity/taxonomy/inext/{taxon}.inext.tsv", taxon=taxa),
        "diversity/taxonomy/inext/plots/ToL.accum.png",
        "diversity/taxonomy/inext/plots/asyest.png",
        expand("diversity/taxonomy/incidences/{taxon}.freq.counts.tsv", taxon=["ToL"]),
        expand("diversity/taxonomy/spader/jaccard.dist.matrix.{taxon}.tsv", taxon=["ToL"]),
        expand("diversity/taxonomy/spader/jaccard.spatiotemporal.dist.{taxon}.tsv", taxon=["ToL"]),
        "diversity/taxonomy/counts/combined.filtered.taxonomy.table.tsv",
        expand("diversity/taxonomy/phyloseq/NMDS.samples.{taxon}.png", taxon=["ToL","Viruses","Archaea","Bacteria","Eukaryota"]),
        "diversity/taxonomy/phyloseq/betadiv.png",

rule db_composition:
    input:
        acc2taxid = "db/kaiju/uniref100.acc2tax",
        taxdb = "db/ete3/taxa.sqlite"
    output:
        png = "diversity/taxonomy/supplementary/uniref100.families.png",
    threads: 1
    run:
        df = pd.read_csv(input.acc2taxid, sep='\t', header=None, names=["acc", "taxid"], usecols=["taxid"])
        df = df.taxid.value_counts().to_frame().rename(columns={"taxid": "count"})
        ncbi = NCBITaxa(dbfile=input.taxdb)
        df["lineage"] = df.index.map(lambda x: ncbi.get_rank(ncbi.get_lineage(x)))
        df["family"] = df.lineage.map(lambda x: [k for k,v in x.items() if v=="family"])
        df["domain"] = df.lineage.map(lambda x: [k for k,v in x.items() if v=="superkingdom"])
        df["family"] = df.family.map(lambda x: ncbi.get_taxid_translator(x).get(x[0]) if x else "NA")
        df["domain"] = df.domain.map(lambda x: ncbi.get_taxid_translator(x).get(x[0]) if x else "NA")
        df = df[(df["family"] != "NA") & (df["domain"] != "NA")]
        families = df.groupby("domain").apply(lambda x: len(x.family.unique())).reset_index().rename(columns={0: "families"})
        sequences = df.groupby("domain").agg({"count": "sum"}).reset_index().rename(columns={"count": "sequences"})
        dat = families.merge(sequences, on="domain")
        dat["families_labels"] = dat.apply(lambda x: x["domain"] + " (" + str(round((x["families"] / dat["families"].sum()) * 100, 2)) + "%)", axis=1)
        dat["sequences_labels"] = dat.apply(lambda x: x["domain"] + " (" + str(round((x["sequences"] / dat["sequences"].sum()) * 100, 2)) + "%)", axis=1)
        fig, axs = plt.subplots(ncols=2, figsize=(7.08, 4.72))
        ax1, ax2 = axs
        wedges1, texts1 = ax1.pie(dat.families.tolist(), wedgeprops=dict(width=0.3))
        wedges2, texts2 = ax2.pie(dat.sequences.tolist(), wedgeprops=dict(width=0.3))
        ax1.legend(wedges1, dat["families_labels"].tolist(), loc="upper center", ncol=1, bbox_to_anchor=(0.5, 0))
        ax2.legend(wedges2, dat["sequences_labels"].tolist(), loc="upper center", ncol=1, bbox_to_anchor=(0.5, 0))
        ax1.axis('equal')
        ax2.axis('equal')
        ax1.text(0., 0., "Families\nn={n:,}".format(n=dat.families.sum()), horizontalalignment='center', verticalalignment='center')
        ax2.text(0., 0., "Sequences\nn={n}M".format(n=round(dat.sequences.sum()/1000000, 2)), horizontalalignment='center', verticalalignment='center')
        plt.tight_layout()
        plt.rcParams.update({'font.size': 12})
        plt.savefig("diversity/taxonomy/supplementary/uniref100.families.png", dpi=300)

rule samples_metaplot:
    input:
         metadata = "metadata.tsv"
    output:
          plot = "diversity/taxonomy/supplementary/samples.metaplot.png"
    threads: 1
    run:
        df = pd.read_csv(input.metadata, sep='\t', header=0)
        df = df[df["year"]==2020]
        plt.figure(figsize=(7.08, 4.72))
        sectors = {"Northern": "s", "Central": "X", "Southern": "^"}
        seasons = ["Summer", "Monsoon", "Winter"]
        ax = sns.scatterplot(data=df, x="salinity", y="temp", hue="season", style="sector", palette="deep", markers=sectors, hue_order=seasons, s=50)
        plt.xlabel("Salinity (ppt)")
        plt.ylabel(u"Temperature (\u00B0C)")
        plt.xlim(0,16)
        plt.ylim(22,34)
        plt.tight_layout()
        plt.rcParams.update({'font.size': 12})
        plt.savefig(output.plot, dpi=300)

rule khist:
    input:
         r1 = "data/fastq/{sample}_R1.fastq.gz",
    output:
        hist = "diversity/taxonomy/khist/{sample}.hist.tsv",
    threads: 12
    conda: "envs/bbtools.yaml"
    log: "logs/diversity/taxonomy/khist/{sample}.log"
    shell: """
            kmercountexact.sh in={input.r1} khist={output.hist} threads={threads} k=31 minprob=0.95 reads=10000000 &> {log}
    """

rule preseq:
    input:
        hist = rules.khist.output.hist
    output:
        extrap = "diversity/taxonomy/preseq/{sample}.extrap.tsv"
    threads: 1
    conda: "envs/preseq.yaml"
    log: "logs/diversity/taxonomy/preseq/{sample}.log"
    shell: """
            preseq lc_extrap -e 100000000000 -s 100000000 -o {output.extrap} -H -v <(grep -v '#' {input.hist} | cut -f 1,2) &> {log}
    """

rule kmer_saturation:
    input:
        extraps = expand("diversity/taxonomy/preseq/{sample}.extrap.tsv", sample=pilot)
    output:
        plot = "diversity/taxonomy/preseq/kmer.saturation.curve.png"
    threads: 1
    run:
        plt.figure(figsize=(7.02, 4.33))
        dfs = []
        for extrap in input.extraps:
            df = pd.read_csv(extrap, sep="\t", header=0, names=["kmers", "uniq", "LCI", "UCI"])
            df["uniq"] = df["uniq"] / df["kmers"]
            df["LCI"] = df["LCI"] / df["kmers"]
            df["UCI"] = df["UCI"] / df["kmers"]
            df["kmers"] /= 1000000000
            dfs.append(df)
        df = pd.concat(dfs, ignore_index=True).fillna(1.0)
        ax = sns.lineplot(data=df, x="kmers", y="uniq", estimator="mean", errorbar=("sd", 1))
        ax.axhline(y=0.05, color="gray", linestyle='--')
        ax.text(0, 0.06, "0.05", color="gray")
        ax.axvline(x=50, color="gray", linestyle='--')
        ax.text(51, 0.06, "50", color="gray")
        plt.xlim(-1,100)
        plt.ylim(0,1)
        plt.xlabel("Observed kmers (billions)")
        plt.ylabel(u"Fraction of unique kmers")
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tight_layout()
        plt.savefig(output.plot, dpi=300)

rule clumpify:
    input:
         r1 = "data/fastq/{sample}_R1.fastq.gz",
         r2 = "data/fastq/{sample}_R1.fastq.gz"
    output:
         r1 = temp("data/qc/{sample}_R1_clumped.fq.gz"),
         r2 = temp("data/qc/{sample}_R2_clumped.fq.gz")
    threads: 40
    log: "logs/diversity/taxonomy/clumpify/{sample}.log"
    conda: "envs/bbtools.yaml"
    shell: "clumpify.sh in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} reorder=p overwrite=t shortname=t dedupe=t ziplevel=1 -Xmx800g -eoom &> {log}"

rule duk:
    input:
        r1 = rules.clumpify.output.r1,
        r2 = rules.clumpify.output.r2
    output:
        r1 = "data/qc/{sample}_R1.fq.gz",
        r2 = "data/qc/{sample}_R2.fq.gz",
        unpaired = "assembly/qc/{sample}_unpaired.fq.gz"
    threads: 6
    log: "logs/diversity/taxonomy/duk/{sample}.log"
    conda: "envs/bbtools.yaml"
    shell: "bbduk.sh in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outs={output.unpaired} ref=adapters,artifacts showspeed=f ziplevel=2 k=31 ktrim=r qtrim=rl trimq=10 mink=11 minlength=51 maxns=3 minavgquality=10 trimpolygright=5 forcetrimmod=5 -Xmx128g -eoom &> {log}"

rule kaiju:
    input:
         db = "db/kaiju/uniref100.fmi",
         nodes = "db/taxonomy/nodes.dmp",
         r1 = rules.duk.output.r1,
         r2 = rules.duk.output.r2,
    output:
         r1 = "diversity/taxonomy/kaiju/{sample}.classification.R1.tsv.gz",
         r2 = "diversity/taxonomy/kaiju/{sample}.classification.R2.tsv.gz"
    threads: 40
    conda: "envs/kaiju.yaml"
    log: "logs/diversity/taxonomy/kaiju/{sample}.log"
    shell: """
            kaiju -t {input.nodes} -f {input.db} -i {input.r1} -a mem -z {threads} -v 2> {log} | bgzip -@ 8 -c > {output.r1}
            kaiju -t {input.nodes} -f {input.db} -i {input.r2} -a mem -z {threads} -v 2>> {log} | bgzip -@ 8 -c > {output.r2}
    """

rule lca:
    input:
        r1 = rules.kaiju.output.r1,
        r2 = rules.kaiju.output.r2,
        nodes = "db/taxonomy/nodes.dmp"
    output:
        unique = temp("diversity/taxonomy/kaiju/{sample}.classification.unique.tsv"),
        duplicated = temp("diversity/taxonomy/kaiju/{sample}.classification.duplicated.tsv"),
        lca = "diversity/taxonomy/kaiju/{sample}.lca.tsv.gz"
    threads: 8
    conda: "envs/kaiju.yaml"
    params:
        tmp = "diversity/taxonomy/kaiju",
    log: "logs/diversity/taxonomy/lca/{sample}.log"
    shell: """
            export LC_ALL=C
            zcat {input.r1} {input.r2} | cut -f 1,2,3,7 | sort -S 128G -k 2,2 -T {params.tmp} --parallel={threads} | \
            tee >(uniq -u | cut -f 1,2,3 > {output.unique}) | uniq -D | awk -v FS='\t' -v OFS='\t' '{{print "U",$2,0}}' > {output.duplicated}

            kaiju-mergeOutputs -t {input.nodes} \
            -i <(cat {output.unique} {output.duplicated} | sed -n 'p;n' | sort -S 64G -k2,2 -T {params.tmp} --parallel={threads}) \
            -j <(cat {output.unique} {output.duplicated} | sed -n 'n;p' | sort -S 64G -k2,2 -T {params.tmp} --parallel={threads}) \
            -c lca 2> {log} | awk -v FS='\t' -v OFS='\t' '{{if ($3 == 1 || $3 == 131567) {{print "U",$2,0}} else {{print $1,$2,$3}}}}' | bgzip -@ {threads} -c > {output.lca}
    """

rule make_report:
    input:
        lca = rules.lca.output.lca,
        tax = "db/kaiju/uniref100.ktaxonomy"
    output:
        kreport = "diversity/taxonomy/counts/{sample}.kreport",
        mpa = "diversity/taxonomy/counts/{sample}.mpa"
    threads: 2
    log: "logs/diversity/taxonomy/make_kreport/{sample}.log"
    shell: """
            python tools/KrakenTools/make_kreport.py -t {input.tax} -o {output.kreport} -i <(bgzip -@ {threads} -cd {input.lca}) &> {log}
            python tools/KrakenTools/kreport2mpa.py -r {output.kreport} -o {output.mpa} &>> {log}
    """

rule classified_reads:
    input:
        kreports = expand("diversity/taxonomy/counts/{sample}.kreport", sample=meta.index.unique())
    output:
        plot = "diversity/taxonomy/supplementary/classified.read.proportions.png",
        tsv = "diversity/taxonomy/supplementary/classified.read.proportions.tsv"
    threads: 1
    run:
        plt.figure(figsize=(4.33, 4.33))
        dat = []
        for sample in meta.index.unique():
            kreport = f"diversity/taxonomy/counts/{sample}.kreport"
            df = pd.read_csv(kreport, sep='\t', header=None, names=["perc", "clade", "node", "rank", "taxid", "name"])
            total = df["node"].sum()
            classified = (total - df.loc[df["name"] == "unclassified", "node"].values[0]) / total
            dat.append((sample, total, classified))
        dat = pd.DataFrame(dat, columns=["sample", "total", "classified"])
        dat.to_csv(output.tsv, sep="\t", header=True, index=False)
        dat["total"] /= 1000000
        g = sns.JointGrid(data=dat, x="total", y="classified", xlim=(0,1300), ylim=(0,1), marginal_ticks=False)
        g.plot(sns.scatterplot, sns.boxplot)
        g.refline(x=416.6)
        g.ax_joint.text(425, 0.01, "416.6", color="gray")
        g.set_axis_labels(xlabel="PE150 reads passing QC (millions)", ylabel="Fraction of taxonomically classified reads")
        plt.tight_layout()
        plt.rcParams.update({'font.size': 12})
        plt.savefig(output.plot, dpi=300)

rule combine_reports:
    input:
        kreport = expand("diversity/taxonomy/counts/{sample}.kreport", sample=meta.index.unique())
    output:
        kreport = "diversity/taxonomy/counts/combined.kreport",
        mpa = "diversity/taxonomy/counts/combined.mpa"
    threads: 1
    log: "logs/diversity/taxonomy/combine_reports.log"
    shell: """
            python tools/KrakenTools/combine_kreports.py -r {input.kreport} --no-headers --only-combined -o {output.kreport} &> {log}
            python tools/KrakenTools/kreport2mpa.py -r {output.kreport} -o {output.mpa} &>> {log}
            sed -i 's/^k__/d__/' {output.mpa}
           """

rule classified_ranks:
    input:
        mpa = rules.combine_reports.output.mpa
    output:
        png = "diversity/taxonomy/supplementary/classified.rank.proportions.png"
    threads: 1
    run:
        df = pd.read_csv(input.mpa, sep='\t', header=None, names=["lineage", "count"])
        df["rank"] = df.lineage.map(lambda x: x.split("|")[-1].split("_")[0])
        ranks = {"Phylum": "p", "Class": "c", "Order": "o", "Family": "f", "Genus": "g", "Species": "s"}
        classified = df.loc[df["rank"]=="d", "count"].sum()
        dat = []
        for rank, r in ranks.items():
            if rank != "Species":
                reads = df.loc[df["rank"]==r, "count"].sum()
            else:
                spec = df[df["rank"]=="s"]
                spec = spec[spec.lineage.map(lambda x: x.split("|")[-2].split("_")[0] == "g")]
                reads = spec["count"].sum()
            dat.append((rank, round(reads/1e9, 2), str(round((reads/classified)*100, 2)) + "%"))
        dat = pd.DataFrame(dat, columns=["rank", "reads", "proportion"])
        plt.figure(figsize=(7.02, 4.33))
        pal = sns.color_palette("rocket", len(dat))
        ax = sns.barplot(x="rank", y="reads", data=dat, palette=pal)
        ax.bar_label(ax.containers[0], labels=dat["proportion"].tolist())
        plt.xlabel("Taxonomic rank")
        plt.ylabel("Classified reads (billions)")
        plt.tight_layout()
        plt.rcParams.update({'font.size': 12})
        plt.savefig(output.png, dpi=300)

rule filter_combined:
    input:
        mpa = rules.combine_reports.output.mpa
    output:
        mpa = "diversity/taxonomy/counts/combined.filtered.mpa",
        richness = "diversity/taxonomy/counts/combined.filtered.richness.krona",
        abundance = "diversity/taxonomy/counts/combined.filtered.abundance.krona"
    threads: 1
    run:
        df = pd.read_csv(input.mpa, sep='\t', header=None, names=["lineage","count"])
        df["rank"] = df.lineage.map(lambda x: x.split("|")[-1].split("_")[0])
        df["parent"] = df.lineage.map(lambda x: "|".join(x.split("|")[:-1]))
        df.index = df.lineage
        df["parent_count"] = df.parent.map(lambda x: 0 if x == '' else df.loc[x, "count"])
        df = df.merge(df.parent.value_counts().to_frame().reset_index().rename(columns={"index":"parent", "parent":"clade_taxa"}), on="parent")
        df["background_count"] = df["parent_count"]/df["clade_taxa"]
        df = df[(df["count"] > df["background_count"]) & (df["count"] >= 1000) & (df["rank"]=="f") &
                (~df["lineage"].str.contains("k__Orthornavirae|k__Pararnavirae|k__Shotokuvirae|k__Loebvirae|f__Aspergillaceae"))]
        df.to_csv(output.mpa, columns=["lineage", "count"], sep="\t", header=False, index=False)
        df.index = df["count"]
        df = df["lineage"].str.split("|", expand=True).fillna("")
        df.to_csv(output.richness, header=False, sep= "\t", index=False)
        df.to_csv(output.abundance, header=False, sep= "\t", index=True)

rule krona:
    input:
        richness = rules.filter_combined.output.richness,
        abundance = rules.filter_combined.output.abundance
    output:
        richness = "diversity/taxonomy/counts/combined.filtered.richness.html",
        abundance = "diversity/taxonomy/counts/combined.filtered.abundance.html"
    threads: 1
    conda: "envs/krakentools.yaml"
    log: "logs/diversity/taxonomy/krona.log"
    shell: """
            ktImportText -q -n 'Tree of Life' -o {output.richness} {input.richness} &> {log}
            ktImportText -n 'Tree of Life' -o {output.abundance} {input.abundance} &>> {log}
           """

rule filter_reports:
    input:
        raw = rules.make_report.output.mpa,
        filtered = rules.filter_combined.output.mpa
    output:
        mpa = "diversity/taxonomy/counts/{sample}.filtered.mpa"
    threads: 1
    run:
        raw = pd.read_csv(input.raw, sep='\t', header=None, names=["lineage","count"])
        filtered = pd.read_csv(input.filtered, sep='\t', header=None, names=["lineage","count"], usecols=["lineage"])
        filtered = filtered.merge(raw, on="lineage", how="left").fillna(0).astype({"count":int})
        filtered.to_csv(output.mpa, columns=["lineage", "count"], sep="\t", header=False, index=False)

rule count_matrix:
    input:
        mpas = expand("diversity/taxonomy/counts/{sample}.filtered.mpa", sample=meta.index.unique())
    output:
        counts = "diversity/taxonomy/counts/combined.filtered.counts.matrix.tsv",
    threads: 1
    run:
        mpas = []
        for mpa in input.mpas:
            sample = mpa.split("/")[-1].split(".")[0]
            mpas.append(pd.read_csv(mpa, sep="\t", header=None, names=["Lineage", sample]))
        df = reduce(lambda i, j: pd.merge(i, j, on='Lineage'), mpas)
        df.to_csv(output.counts, sep='\t', header=True, index=False)

rule checklist:
    input:
        mpa = rules.filter_combined.output.mpa,
        checklist = "data/{taxon}_families_chilika.tsv",
        proteomes = "data/{taxon}_uniprot-proteomes.tsv"
    output:
        tsv = "diversity/taxonomy/checklist/{taxon}.checklist.tsv"
    threads: 1
    run:
        checklist = pd.read_csv(input.checklist, header=0, sep='\t')
        checklist["checklist"] = True
        proteomes = pd.read_csv(input.proteomes, header=0, sep='\t')
        proteomes_families = []
        fams = []
        for idx, row in proteomes.iterrows():
            for taxon in row["Taxonomic lineage"].split(","):
                taxon = taxon.strip()
                if not taxon in fams and taxon.endswith("dae"):
                    fams.append(taxon)
                    proteomes_families.append((taxon, row["Proteome Id"]))
        proteomes_families = pd.DataFrame(proteomes_families, columns=["family", "Representative Proteome Id"])
        filtered_lineages = pd.read_csv(input.mpa, sep='\t', header=None, names=["lineage", "count"])
        filtered_lineages["family"] = filtered_lineages.lineage.map(lambda x: x.split("|")[-1].replace("f__",""))
        filtered_lineages = filtered_lineages[["family", "count"]][filtered_lineages.lineage.str.contains("c__Actinopteri")]
        df = checklist.merge(filtered_lineages, on="family", how="outer").fillna(False)
        df.loc[df["count"] == False, "count"] = 0
        df["uniprot_proteome"] = df["family"].map(lambda x: x in proteomes_families["family"].tolist())
        df = df.merge(proteomes_families, on="family", how="left")
        df.to_csv(output.tsv, sep="\t", header=True, index=False)

rule checklist_relabund:
    input:
        checklist = rules.checklist.output.tsv,
        counts = rules.count_matrix.output.counts
    output:
        plot = "diversity/taxonomy/checklist/{taxon}.relabund.png",
        tsv = "diversity/taxonomy/checklist/{taxon}.relabund.tsv"
    threads: 1
    run:
        checklist = pd.read_csv(input.checklist, sep='\t', header=0)
        checklist = checklist[["family"]][(checklist["checklist"]) & (checklist["uniprot_proteome"])]
        counts = pd.read_csv(input.counts, sep='\t', header=0)
        counts["family"] = counts["Lineage"].map(lambda x: x.split("|")[-1].replace("f__", ""))
        counts = checklist.merge(counts, on="family", how="left").fillna(0)
        counts.index = counts["family"].tolist()
        counts = counts.drop(columns=["Lineage", "family"])
        counts["total"] = counts.apply(lambda x: x.sum(), axis=1)
        counts = counts.sort_values(by="total", ascending=False)
        counts = counts.drop(columns=["total"])
        counts = counts.apply(lambda x: x/x.sum())
        counts = counts[seasonal]
        counts.to_csv(output.tsv, sep="\t", header=True, index=True)
        plt.figure(figsize=(4.33, 4.33))
        sns.set(font_scale=0.5)
        sns.heatmap(counts, xticklabels=1, yticklabels=1, cmap="rainbow", norm=LogNorm())
        plt.tick_params(axis='both', which='major', labelsize=8)
        plt.tight_layout()
        plt.savefig(output.plot, dpi=300)

rule incidence_freq:
    input:
        kreport = rules.combine_reports.output.kreport,
        mpa = rules.filter_combined.output.mpa
    output:
        inc = "diversity/taxonomy/incidences/{taxon}.freq.mpa"
    threads: 1
    run:
        kreport = pd.read_csv(input.kreport, sep='\t', header=None, names=["perc", "clade", "node", "rank", "taxid", "name"], usecols=["node"])
        df = pd.read_csv(input.mpa, sep='\t', header=None, names=["lineage","count"], index_col=["lineage"])
        reads = kreport["node"].sum()
        size = 100000000
        mincount = 1000
        units = reads//size
        if wildcards.taxon != "ToL":
            df = df[df.index.str.contains(wildcards.taxon)]
        df["count"] = df["count"].map(lambda x: units if (x//mincount) > units else x//mincount)
        df = pd.concat([pd.DataFrame([units],columns=["count"],index=["sampling_units"]), df])
        df.to_csv(output.inc, sep='\t', header=False, index=True)

rule inext:
    input:
        inc = rules.incidence_freq.output.inc,
        script = "scripts/inext.R"
    output:
        inext = "diversity/taxonomy/inext/{taxon}.inext.tsv",
        asyest = "diversity/taxonomy/inext/{taxon}.asyest.tsv",
    conda: "envs/inext.yaml"
    threads: 1
    log: "logs/diversity/taxonomy/inext/{taxon}.log"
    shell: """
            Rscript {input.script} --taxon {wildcards.taxon} --inc {input.inc} --inext {output.inext} --asyest {output.asyest} &> {log}
    """

rule accum_curve:
    input:
        tol_inext = expand("diversity/taxonomy/inext/{taxon}.inext.tsv", taxon="ToL"),
        domains_inext = expand("diversity/taxonomy/inext/{taxon}.inext.tsv", taxon=["Viruses","Archaea","Bacteria","Eukaryota"]),
        tol_script = "scripts/plot.tol.accum.R",
        domains_script = "scripts/plot.domains.accum.R"
    output:
        tol_plot = "diversity/taxonomy/inext/plots/ToL.accum.png",
        domains_plot = "diversity/taxonomy/inext/plots/domains.accum.png",
        domains_inext = temp("diversity/taxonomy/inext/plots/domains.inext.tsv")
    threads: 1
    log: "logs/diversity/taxonomy/plot_accum.log"
    run:
        shell("Rscript {input.tol_script} --inext {input.tol_inext} --plot {output.tol_plot} &> {log}")
        df = pd.concat([pd.read_csv(i, sep='\t', header=0) for i in input.domains_inext], ignore_index=True)
        df.to_csv(output.domains_inext, sep='\t', header=True, index=False)
        shell("Rscript {input.domains_script} --inext {output.domains_inext} --plot {output.domains_plot} &> {log}")

rule asyest:
    input:
        asyest = expand("diversity/taxonomy/inext/{taxon}.asyest.tsv", taxon=["Viruses","Archaea","Bacteria","Eukaryota","Metazoa","Fungi","Viridiplantae","Arthropoda","Chordata","Actinopteri","Mammalia","Aves","Amphibia"]),
        script = "scripts/plot.asyest.R"
    output:
        plot = "diversity/taxonomy/inext/plots/asyest.png",
        asyst = "diversity/taxonomy/inext/plots/asyest.tsv"
    threads: 1
    log: "logs/diversity/taxonomy/plot_asyest.log"
    run:
        df = pd.concat([pd.read_csv(i, sep='\t', header=0) for i in input.asyest], ignore_index=True)
        df = df[df["Diversity"]=="Species richness"]
        df.to_csv(output.asyst, sep='\t', header=True, index=False)
        shell("Rscript {input.script} --asyest {output.asyst} --plot {output.plot} &> {log}")

rule samples_inc_freq:
    input:
        kreport = expand("diversity/taxonomy/counts/{sample}.kreport", sample=seasonal),
        counts = rules.count_matrix.output.counts
    output:
        inc = "diversity/taxonomy/incidences/{taxon}.freq.counts.tsv"
    run:
        df = pd.read_csv(input.counts, sep='\t', header=0, index_col=["Lineage"])
        if wildcards.taxon != "ToL":
            df = df[df.index.str.contains(wildcards.taxon)]
        df = pd.concat([pd.DataFrame([], columns=df.columns, index=["sampling_units"]).fillna(0), df])
        for sample in seasonal:
            kreport = pd.read_csv(f"diversity/taxonomy/counts/{sample}.kreport", sep='\t', header=None, names=["perc", "clade", "node", "rank", "taxid", "name"], usecols=["node"])
            reads = kreport["node"].sum()
            size = 10000000
            mincount = 1000
            units = reads//size
            df[sample] = df[sample].map(lambda x: units if (x//mincount) > units else x//mincount)
            df.loc["sampling_units", sample] = units
        df.to_csv(output.inc, sep='\t', header=True, index=True)

rule spader:
    input:
        inc = rules.samples_inc_freq.output.inc,
        script = "scripts/spader.R"
    output:
        jaccard = "diversity/taxonomy/spader/jaccard.dist.matrix.{taxon}.tsv",
    conda: "envs/spader.yaml"
    threads: 1
    log: "logs/diversity/taxonomy/spader/{taxon}.log"
    shell: """
            Rscript {input.script} --inc {input.inc} --jaccard {output.jaccard} &> {log}
    """

rule richness_betadiv:
    input:
        jaccard = expand("diversity/taxonomy/spader/jaccard.dist.matrix.{taxon}.tsv", taxon=["ToL"]),
        meta = "metadata.tsv"
    output:
        dist = "diversity/taxonomy/spader/jaccard.spatiotemporal.dist.{taxon}.tsv"
    run:
        dist = []
        for taxon in ["ToL"]:
            distmat = pd.read_csv(f"diversity/taxonomy/spader/jaccard.dist.matrix.{taxon}.tsv", sep="\t", header=0, index_col=0)
            meta = pd.read_csv(input.meta, sep='\t', header=0, index_col="sample")
            meta = meta[meta["year"]==2020]
            if taxon == "ToL":
                taxon = "Tree of Life"
            for season,rows in meta.groupby("season"):
                for i,j in itertools.combinations(rows.index, 2):
                    jaccard = round(distmat.loc[i,j], 2)
                    dist.append((i,j,jaccard, "Spatial beta diversity", taxon))
            for station,rows in meta.groupby("station"):
                for i,j in itertools.combinations(rows.index, 2):
                    jaccard = round(distmat.loc[i,j], 2)
                    dist.append((i,j,jaccard, "Temporal beta diversity", taxon))
        dist = pd.DataFrame(dist, columns=["sample1", "sample2", "jaccard", "type", "domain"])
        dist.to_csv(output.dist, sep='\t', header=True, index=False)

rule taxonomy_table:
    input:
        counts = rules.count_matrix.output.counts
    output:
        tsv = "diversity/taxonomy/counts/combined.filtered.taxonomy.table.tsv"
    threads: 1
    run:
        counts = pd.read_csv(input.counts, sep='\t', header=0, usecols=["Lineage"])
        counts[["Domain", "Phylum", "Class", "Order", "Family"]] = "NA"
        for idx, row in counts.iterrows():
            levels = {i.split("_")[0]:i for i in row["Lineage"].split("|")}
            counts.loc[idx, "Domain"] = levels.get("d")
            counts.loc[idx, "Phylum"] = levels.get("p")
            counts.loc[idx, "Class"] = levels.get("c")
            counts.loc[idx, "Order"] = levels.get("o")
            counts.loc[idx, "Family"] = levels.get("f")
        counts.to_csv(output.tsv, sep='\t', header=True, index=False)

rule phyloseq:
    input:
        counts = rules.count_matrix.output.counts,
        tax = rules.taxonomy_table.output.tsv,
        meta = "metadata.tsv",
        script = "scripts/phyloseq.R"
    output:
        nmds = "diversity/taxonomy/phyloseq/NMDS.samples.{taxon}.png",
        distmat = "diversity/taxonomy/phyloseq/dist.matrix.{taxon}.tsv"
    threads: 1
    conda: "envs/phyloseq.yaml"
    log: "logs/diversity/taxonomy/phyloseq/{taxon}.log"
    shell: """
            Rscript {input.script} --otu {input.counts} --tax {input.tax} --meta {input.meta} --domain {wildcards.taxon} --nmds {output.nmds} --dist {output.distmat} &> {log}
    """

rule betadiv:
    input:
        distmat = expand("diversity/taxonomy/phyloseq/dist.matrix.{taxon}.tsv", taxon=["ToL","Viruses","Archaea","Bacteria","Eukaryota"]),
        meta = "metadata.tsv",
        script = "scripts/plot.betadiv.R"
    output:
        plot = "diversity/taxonomy/phyloseq/betadiv.png",
        dist = "diversity/taxonomy/phyloseq/betadiv.tsv",
    threads: 1
    log: "logs/diversity/taxonomy/betadiv.log"
    run:
        dist = []
        for domain in ["ToL","Viruses","Archaea","Bacteria","Eukaryota"]:
            distmat = pd.read_csv(f"diversity/taxonomy/phyloseq/dist.matrix.{domain}.tsv", sep="\t", header=0, index_col=0)
            meta = pd.read_csv(input.meta, sep='\t', header=0, index_col="sample")
            meta = meta[meta["year"]==2020]
            if domain == "ToL":
                domain = "Tree of Life"
            for season,rows in meta.groupby("season"):
                for i,j in itertools.combinations(rows.index, 2):
                    bray = round(distmat.loc[i,j], 2)
                    dist.append((i,j,bray, "Spatial beta diversity", domain))
            for station,rows in meta.groupby("station"):
                for i,j in itertools.combinations(rows.index, 2):
                    bray = round(distmat.loc[i,j], 2)
                    dist.append((i,j,bray, "Temporal beta diversity", domain))
        dist = pd.DataFrame(dist, columns=["sample1", "sample2", "bray", "type", "domain"])
        dist.to_csv(output.dist, sep='\t', header=True, index=False)
        shell("Rscript {input.script} --dist {output.dist} --plot {output.plot} &> {log}")
