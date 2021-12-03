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
    format="(?!chimeric)[^_]+",


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
        snek=touch("benchmarks/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/snek"),
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        percs=lambda wildcards: ",".join([str(ab) for ab in config["abundances"]]),
        spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
        json=lambda wildcards: json.dumps(config),
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
        shell("cat {params.json} > {output.snek}")


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
    threads: 12
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        errs=lambda wildcards: ",".join([str(er) for er in errors]),
        percs=lambda wildcards: ",".join([str(ab) for ab in abus]),
        spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
        json=lambda wildcards: json.dumps(config),
    run:
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


rule create_benchmark_chimeric_compare:
    input:
        fasta="genome_data/sequences.fasta",
        metadata="genome_data/metadata.tsv",
        voc=expand("genome_data/{voc}.fasta", voc=pangolin),
    output:
        expand(
            "benchmarks/{{dataset}}_chimeric/wwsim_{voc}_ab{ab}_ch{er}.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
        # expand(
        #     "benchmarks/{{dataset}}_chimeric/wwsim_{voc}_ab{ab}_ch{er}_2.fastq",
        #     voc=pangolin,
        #     ab=abus,
        #     er=errors,
        # ),
        snek=touch("benchmarks/{dataset}_chimeric/snek"),
    threads: 12
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        errs=lambda wildcards: ",".join([str(er) for er in errors]),
        percs=lambda wildcards: ",".join([str(ab) for ab in abus]),
        spike=lambda wildcards: "--spike_only" if config["spike_only"] else "",
        json=lambda wildcards: json.dumps(config),
    run:
        shell(
            "python benchmarking/create_chimeric_benchmarks.py "
            "--voc_perc {params.percs} "
            "--chim_perc {params.errs} "
            "-m {input.metadata} -fr {input.fasta} "
            "-fv {params.vocs} "
            "-o benchmarks/{wildcards.dataset}_chimeric "
            "--total_cov {config[tot_cov]} "
            "{params.spike} "
        )
        shell("echo {params.json} > {output.snek}")


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


rule run_kallisto_error_compare:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
        wwsim1=expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_1.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
        wwsim2=expand(
            "benchmarks/{{dataset}}_{{format}}/wwsim_{voc}_ab{ab}_er{er}_2.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
    output:
        expand(
            "benchmarks/{{dataset}}_{{format}}/out/{voc}_ab{ab}_er{er}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            er=errors,
            min_ab=config["min_ab"],
        ),
        dir=directory("benchmarks/{dataset}_{format}/out"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
    threads: 12
    run:
        outdir = output.dir
        shell("mkdir -p {outdir}")
        for voc in pangolin:
            for ab in abus:
                for er in errors:
                    shell(
                        "kallisto quant -t {threads} -b {config[bootstraps]} "
                        "-i {input.idx} -o {outdir}/{voc}_ab{ab}_er{er} "
                        "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{er}_1.fastq "
                        "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_er{er}_2.fastq "
                        "| tee {outdir}/{voc}_ab{ab}_er{er}.log"
                    )
                    shell(
                        "python pipeline/output_abundances.py "
                        "-m {config[min_ab]} "
                        "-o {outdir}/{voc}_ab{ab}_er{er}/predictions_m{config[min_ab]}.tsv "
                        "--metadata {config[ref]}/metadata.tsv "
                        "--voc {params.vocs} "
                        "{outdir}/{voc}_ab{ab}_er{er}/abundance.tsv "
                        "| tee -a {outdir}/{voc}_ab{ab}_er{er}.log"
                    )


rule run_kallisto_chimeric_compare:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=config["ref"]),
        wwsim=expand(
            "benchmarks/{{dataset}}_chimeric/wwsim_{voc}_ab{ab}_ch{er}.fastq",
            voc=pangolin,
            ab=abus,
            er=errors,
        ),
        # wwsim2=expand(
        #     "benchmarks/{{dataset}}_chimeric/wwsim_{voc}_ab{ab}_ch{er}_2.fastq",
        #     voc=pangolin,
        #     ab=abus,
        #     er=errors,
        # ),
    output:
        expand(
            "benchmarks/{{dataset}}_chimeric/out/{voc}_ab{ab}_ch{er}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            er=errors,
            min_ab=config["min_ab"],
        ),
        dir=directory("benchmarks/{dataset}_chimeric/out"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
    threads: 12
    run:
        outdir = output.dir
        shell("mkdir -p {outdir}")
        for voc in pangolin:
            for ab in abus:
                for er in errors:
                    shell(
                        "kallisto quant -t {threads} -b {config[bootstraps]} "
                        "-i {input.idx} -o {outdir}/{voc}_ab{ab}_ch{er} "
                        # Sadly, badread doesnt seem to support paired end reads
                        "--single -l 150 -s 1 "
                        "benchmarks/{wildcards.dataset}_chimeric/wwsim_{voc}_ab{ab}_ch{er}.fastq "
                        # "benchmarks/{wildcards.dataset}_{wildcards.format}/wwsim_{voc}_ab{ab}_ch{er}_2.fastq "
                        "| tee {outdir}/{voc}_ab{ab}_ch{er}.log"
                    )
                    shell(
                        "python pipeline/output_abundances.py "
                        "-m {config[min_ab]} "
                        "-o {outdir}/{voc}_ab{ab}_ch{er}/predictions_m{config[min_ab]}.tsv "
                        "--metadata {config[ref]}/metadata.tsv "
                        "--voc {params.vocs} "
                        "{outdir}/{voc}_ab{ab}_ch{er}/abundance.tsv "
                        "| tee -a {outdir}/{voc}_ab{ab}_ch{er}.log"
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
        snek="benchmarks/figs/{dataset}_s{sub_err}_i{ins_err}_d{del_err}/snek",
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
        "--plot_abundance_value 10 "
        "-o {output.dir} "
        "--output_format {params.exts} "
        "benchmarks/{wildcards.dataset}_{wildcards.format}/out/*/predictions_m{config[min_ab]}.tsv "
        "&& echo {params.json} > {output.snek}"
        # "-m {config[min_ab]} "


rule create_figs_compare_chimeric:
    input:
        expand(
            "benchmarks/{{dataset}}_chimeric/out/{voc}_ab{ab}_ch{er}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=abus,
            er=errors,
            min_ab=config["min_ab"],
        ),
    output:
        expand(
            "benchmarks/figs/{{dataset}}_chimeric/{file}.{ext}",
            ext=config["plot_exts"],
            file=[
                "freq_error_plot",
                "freq_error_plot_logscale",
                "freq_scatter_loglog",
                "error_error_plot",
                "error_error_plot_logscale",
            ],
        ),
        snek="benchmarks/figs/{dataset}_chimeric/snek",
        dir=directory("benchmarks/figs/{dataset}_chimeric"),
    params:
        vocs=lambda wildcards, input: ",".join(config["vocs"]),
        exts=lambda wildcards, input: ",".join(config["plot_exts"]),
        json=lambda wildcards, input: json.dumps(config),
    shell:
        "python benchmarking/evaluate_error.py "
        "--voc {params.vocs} "
        "--plot_abundance_value 10 "
        "-o {output.dir} "
        "--output_format {params.exts} "
        "--chimeric "
        "benchmarks/{wildcards.dataset}_chimeric/out/*/predictions_m{config[min_ab]}.tsv "
        "&& echo {params.json} > {output.snek}"
        # "-m {config[min_ab]} "
