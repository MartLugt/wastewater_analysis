import json
import pprint


configfile: "snek_config.yaml"


pangolin = [x + "_EPI_ISL_" + y for x, y in zip(config["vocs"], config["lineages"])]


# dataset should be in format:
# [datasetname]_se_ie_de
wildcard_constraints:
    dataset="[^/]+",


rule create_benchmark:
    input:
        fasta="genome_data/sequences.fasta",
        metadata="genome_data/metadata.tsv",
        voc=expand("genome_data/{voc}.fasta", voc=pangolin),
    output:
        expand(
            "benchmarks/{{dataset}}_s{{sub_err}}_i{{ins_err}}_d{{del_err}}/wwsim_{voc}_ab{ab}_1.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
        expand(
            "benchmarks/{{dataset}}_s{{sub_err}}_i{{ins_err}}_d{{del_err}}/wwsim_{voc}_ab{ab}_2.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
        touch("benchmarks/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/snek"),
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        percs=lambda wildcards: ",".join([str(ab) for ab in config["abundances"]]),
        spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
    run:
        shell(
            "python benchmarking/create_benchmarks.py "
            "--voc_perc {params.percs} "
            "-m {input.metadata} -fr {input.fasta} "
            "-fv {params.vocs} "
            "-o benchmarks/{wildcards.dataset}_s{wildcards.sub_err}_i{wildcards.ins_err}_d{wildcards.del_err} "
            "--total_cov {config[tot_cov]} "
            "{params.spike} "
            "--sub_error_rate {wildcards.sub_err} "
            "--ins_error_rate {wildcards.ins_err} "
            "--del_error_rate {wildcards.del_err} "
        )


rule run_kallisto:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
        wwsim1=expand(
            "benchmarks/{{dataset}}_s{{sub_err}}_i{{ins_err}}_d{{del_err}}/wwsim_{voc}_ab{ab}_1.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
        wwsim2=expand(
            "benchmarks/{{dataset}}_s{{sub_err}}_i{{ins_err}}_d{{del_err}}/wwsim_{voc}_ab{ab}_2.fastq",
            voc=pangolin,
            ab=config["abundances"],
        ),
    output:
        expand(
            "benchmarks/{{dataset}}_s{{sub_err}}_i{{ins_err}}_d{{del_err}}/out/{voc}_ab{ab}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=config["abundances"],
            min_ab=config["min_ab"],
        ),
        dir=directory("benchmarks/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/out"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
    threads: 12
    run:
        outdir = output.dir
        shell("mkdir -p {outdir}")
        for voc in pangolin:
            for ab in config["abundances"]:
                shell(
                    "kallisto quant -t {threads} -b {config[bootstraps]} "
                    "-i {input.idx} -o {outdir}/{voc}_ab{ab} "
                    "benchmarks/{wildcards.dataset}_s{wildcards.sub_err}_i{wildcards.ins_err}_d{wildcards.del_err}/wwsim_{voc}_ab{ab}_1.fastq "
                    "benchmarks/{wildcards.dataset}_s{wildcards.sub_err}_i{wildcards.ins_err}_d{wildcards.del_err}/wwsim_{voc}_ab{ab}_2.fastq "
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
            "benchmarks/{{dataset}}_s{{sub_err}}_i{{ins_err}}_d{{del_err}}/out/{voc}_ab{ab}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=config["abundances"],
            min_ab=config["min_ab"],
        ),
    output:
        "benchmarks/figs/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/freq_error_plot_logscale.png",
        "benchmarks/figs/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/freq_error_plot.png",
        "benchmarks/figs/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/freq_scatter_loglog.png",
        snek="benchmarks/figs/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/info.snek",
        dir=directory("benchmarks/figs/{dataset}_s{sub_err}_i{ins_err}_d{del_err}"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
        json=lambda wildcards, input: json.dumps(config),
    shell:
        "python benchmarking/evaluate_abundances.py "
        "--voc {params.vocs} "
        "-m {config[min_ab]} "
        "-o {output.dir} "
        "benchmarks/{wildcards.dataset}_s{wildcards.sub_err}_i{wildcards.ins_err}_d{wildcards.del_err}/out/*/predictions_m{config[min_ab]}.tsv && "
        "echo {params.json} > {output.snek} "


# rule create_figs_compare_error:
#     input:
#         lambda wildcards: expand(
#             "benchmarks/{ds}/out/{voc}_ab{ab}/predictions_m{min_ab}.tsv",
#             voc=pangolin,
#             ab=config["abundances"],
#             min_ab=config["min_ab"],
#             ds=wildcards.list.split(","),
#         ),
#     output:
#         touch("benchmarks/figs/error_compare/{dataset}/{list}"),
