import json
import pprint


configfile: "snek_config.yaml"


pangolin = [x + "_EPI_ISL_" + y for x, y in zip(config["vocs"], config["lineages"])]


wildcard_constraints:
    dataset="[^/]+",


rule create_benchmark:
    input:
        fasta="genome_data/sequences.fasta",
        metadata="genome_data/metadata.tsv",
        voc=expand("genome_data/{voc}.fasta", voc=pangolin),
    output:
        expand(
            "benchmarks/{{dataset}}/wwsim_{voc}_ab{ab}_1.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
        expand(
            "benchmarks/{{dataset}}/wwsim_{voc}_ab{ab}_2.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
        touch("benchmarks/{dataset}/snek"),
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        percs=lambda wildcards, input: ",".join(
            [str(ab) for ab in config["abundances"]]
        ),
        spike=lambda wildcards, input: "--spike_only" if config["spike_only"] else "",
    shell:
        "python benchmarking/create_benchmarks.py "
        "--voc_perc {params.percs} "
        "-m {input.metadata} -fr {input.fasta} "
        "-fv {params.vocs} "
        "-o benchmarks/{wildcards.dataset} "
        "--total_cov {config[tot_cov]} {params.spike} "


rule run_kallisto:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
        wwsim1=expand(
            "benchmarks/{{dataset}}/wwsim_{voc}_ab{ab}_1.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
        wwsim2=expand(
            "benchmarks/{{dataset}}/wwsim_{voc}_ab{ab}_2.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
    output:
        expand(
            "benchmarks/{{dataset}}/out/{voc}_ab{ab}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=config["abundances"],
            min_ab=config["min_ab"],
        ),
        dir=directory("benchmarks/{dataset}/out"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
    run:
        outdir = output.dir
        shell("mkdir -p {outdir}")
        for voc in pangolin:
            for ab in config["abundances"]:
                shell(
                    "kallisto quant -t {threads} -b {config[bootstraps]} "
                    "-i {input.idx} -o {outdir}/{voc}_ab{ab} "
                    "benchmarks/{wildcards.dataset}/wwsim_{voc}_ab{ab}_1.fastq "
                    "benchmarks/{wildcards.dataset}/wwsim_{voc}_ab{ab}_2.fastq "
                    "| tee {outdir}/{voc}_ab{ab}.log"
                )
                shell(
                    "python pipeline/output_abundances.py "
                    "-m {config[min_ab]} "
                    "-o {outdir}/{voc}_ab{ab}/predictions_m{config[min_ab]}.tsv "
                    "--metadata {config[ref]}/metadata.tsv "
                    "--voc {params.vocs} "
                    "{outdir}/{voc}_ab{ab}/abundance.tsv "
                )


rule create_figs:
    input:
        expand(
            "benchmarks/{{dataset}}/out/{voc}_ab{ab}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=config["abundances"],
            min_ab=config["min_ab"],
        ),
    output:
        "benchmarks/figs/{dataset}/freq_error_plot_logscale.png",  #add sub
        "benchmarks/figs/{dataset}/freq_error_plot.png",
        "benchmarks/figs/{dataset}/freq_scatter_loglog.png",
        snek="benchmarks/figs/{dataset}/info.snek",
        dir=directory("benchmarks/figs/{dataset}"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
        json=lambda wildcards, input: json.dumps(config),
    shell:
        "python benchmarking/evaluate_abundances.py "
        "--voc {params.vocs} "
        "-m {config[min_ab]} "
        "-o {output.dir} "
        "benchmarks/{wildcards.dataset}/out/*/predictions_m{config[min_ab]}.tsv && "
        "echo {params.json} > {output.snek} "
