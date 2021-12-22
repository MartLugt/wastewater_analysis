import json
import pprint
import numpy as np
from math import log10, floor


configfile: "snek_config.yaml"


pangolin = [x + "_EPI_ISL_" + y for x, y in zip(config["vocs"], config["lineages"])]


def round_sig(x, sig=2):
    if x == 0:
        return x
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


abus = [
    round(float(x), 2)
    for x in (
        config["abundances"]
        if not config["abundances_scale_log"]
        else np.geomspace(
            config["abundances_log"][0],
            config["abundances_log"][1],
            config["abundances_log"][2],
        )
    )
]
errors = [
    round_sig(float(x), 3)
    for x in (
        config["errors"]
        if not config["errors_scale_log"]
        else np.geomspace(
            config["errors_log"][1],
            config["errors_log"][0],
            config["errors_log"][2],
        )
    )
]


wildcard_constraints:
    dataset="[^/]+",
    # format="(?!chimeric)[^_]+",
    format="[^_]+",


rule create_benchmark_error_compare:
    input:
        fasta="genome_data/sequences.fasta",
        metadata="genome_data/metadata.tsv",
        voc=expand("genome_data/{voc}.fasta", voc=pangolin),
    output:
        expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_1.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
        expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_2.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
        snek=touch("benchmarks/{dataset}_{format}/snek"),
    threads: 2
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        errs=lambda wildcards: ",".join([str(er) for er in errors]),
        percs=lambda wildcards: ",".join([str(ab) for ab in abus]),
        spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
        json=lambda wildcards: json.dumps(config),
    run:
        if "chimeric" in wildcards.format:
            shell(
                "python benchmarking/create_chimeric_benchmarks.py "
                "--voc_perc {params.percs} "
                "--chim_perc {params.errs} "
                "-m {input.metadata} -fr {input.fasta} "
                "-fv {params.vocs} "
                "-o benchmarks/{wildcards.dataset}_{wildcards.format} "
                "--total_cov {config[tot_cov]} "
                "{params.spike} "
            )
        else:
            sid = "--sub_error " if "s" in wildcards.format else ""
            sid += "--ins_error " if "i" in wildcards.format else ""
            sid += "--del_error " if "d" in wildcards.format else ""
            shell(
                "python benchmarking/create_error_benchmarks.py "
                "--voc_perc {params.percs} "
                "--err_perc {params.errs} "
                "-m {input.metadata} -fr {input.fasta} "
                "-fv {params.vocs} "
                "-o benchmarks/{wildcards.dataset}_{wildcards.format} "
                "--total_cov {config[tot_cov]} "
                "{params.spike} "
                "{sid} "
            )
        shell("echo {params.json} > {output.snek}")


rule run_kallisto_batch_jobs:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
        wwsim1=expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{{er}}_1.fastq",
            voc=pangolin,
            ab=abus,
        ),
        wwsim2=expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{{er}}_2.fastq",
            voc=pangolin,
            ab=abus,
        ),
    output:
        preds=expand(
            "benchmarks/{{dataset}}_{{format}}/out/{voc}_ab{ab}_er{{er}}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            min_ab=config["min_ab"],
        ),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
    threads: 2
    run:
        outdir = str(output.preds).split("/")[0:-2]
        shell("mkdir -p {outdir}")
        for voc in pangolin:
            for ab in abus:
                shell(
                    "srun --ntasks=1 --cpus-per-task={threads} "
                    "kallisto quant -t {threads} -b {config[bootstraps]} "
                    "-i {input.idx} -o {outdir}/{voc}_ab{ab}_er{wildcards.er} "
                    "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{wildcards.er}_1.fastq "
                    "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{wildcards.er}_2.fastq "
                )
                shell(
                    "srun --ntasks=1 --cpus-per-task=1 "
                    "python pipeline/output_abundances.py "
                    "-m {config[min_ab]} "
                    "-o {outdir}/{voc}_ab{ab}_er{wildcards.er}/predictions_m{config[min_ab]}.tsv "
                    "--metadata {config[ref]}/metadata.tsv "
                    "--voc {params.vocs} "
                    "{outdir}/{voc}_ab{ab}_er{wildcards.er}/abundance.tsv "
                )


# rule run_kallisto_error_compare:
#     input:
#         idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
#         wwsim1=expand(
#             "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_1.fastq",
#             voc=pangolin,
#             ab=abus,
#             er=errors,
#         ),
#         wwsim2=expand(
#             "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_2.fastq",
#             voc=pangolin,
#             ab=abus,
#             er=errors,
#         ),
#     output:
#         expand(
#             "benchmarks/{{dataset}}_{{format}}/out/{voc}_ab{ab}_er{er}/predictions_m{min_ab}.tsv",
#             voc=pangolin,
#             ab=abus,
#             er=errors,
#             min_ab=config["min_ab"],
#         ),
#         dir=directory("benchmarks/{dataset}_{format}/out"),
#     params:
#         vocs=lambda wildcards, input: ",".join(config["vocs"]),
#     threads: 12
#     run:
#         outdir = output.dir
#         shell("mkdir -p {outdir}")
#         for voc in pangolin:
#             for ab in abus:
#                 for er in errors:
#                     shell(
#                         "kallisto quant -t {threads} -b {config[bootstraps]} "
#                         "-i {input.idx} -o {outdir}/{voc}_ab{ab}_er{er} "
#                         "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{er}_1.fastq "
#                         "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{er}_2.fastq "
#                         "| tee {outdir}/{voc}_ab{ab}_er{er}.log"
#                     )
#                     shell(
#                         "python pipeline/output_abundances.py "
#                         "-m {config[min_ab]} "
#                         "-o {outdir}/{voc}_ab{ab}_er{er}/predictions_m{config[min_ab]}.tsv "
#                         "--metadata {config[ref]}/metadata.tsv "
#                         "--voc {params.vocs} "
#                         "{outdir}/{voc}_ab{ab}_er{er}/abundance.tsv "
#                         "| tee -a {outdir}/{voc}_ab{ab}_er{er}.log"
#                     )


rule create_figs_compare_error:
    input:
        expand(
            "benchmarks/{{dataset}}_{{format}}/out/{voc}_ab{ab}_er{er}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            er=errors,
            min_ab=config["min_ab"],
        ),
    output:
        expand(
            "benchmarks/figs/{{dataset}}_{{format}}/{file}.{ext}",
            ext=config["plot_exts"],
            file=[
                "freq_error_plot",
                "freq_error_plot_logscale",
                "freq_scatter_loglog",
                "error_error_plot",
                "error_error_plot_logscale",
            ],
        ),
        snek="benchmarks/figs/{dataset}_{format}/snek",
        dir=directory("benchmarks/figs/{dataset}_{format}"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
        exts=lambda wildcards, input: ",".join(config["plot_exts"]),
        json=lambda wildcards, input: json.dumps(config),
    shell:
        "python benchmarking/evaluate_error.py "
        "--voc {params.vocs} "
        "--plot_abundance_value {config[plot_abundance_value]} "
        "--plot_error_value {config[plot_error_value]} "
        "-o {output.dir} "
        "--output_format {params.exts} "
        "benchmarks/{wildcards.dataset}_{wildcards.format}/out/*/predictions_m{config[min_ab]}.tsv "
        "&& echo {params.json} > {output.snek}"
        # "-m {config[min_ab]} "
