from glob import glob

VOC = ["B.1.1.7", "B.1.351", "B.1.617.2", "P.1"]
LIN = ["889440", "1001460", "1924762", "1194849"]
ABUNDANCES = [
    0.05,
    0.06,
    0.07,
    0.08,
    0.09,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    20,
    30,
    40,
    50,
    60,
    70,
    80,
    90,
    100,
]
REF = "reference_sets/USA"
MIN_AB = 1

pangolin = [x + "_EPI_ISL_" + y for x, y in zip(VOC, LIN)]

bootstraps = 0


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
            ab=ABUNDANCES,
        ),
        expand(
            "benchmarks/{{dataset}}/wwsim_{voc}_ab{ab}_2.fastq",
            voc=pangolin,
            ab=ABUNDANCES,
        ),
        touch("benchmarks/{dataset}/snek"),
    params:
        vocs=lambda wildcards, input: ",".join(input.voc),
        percs=lambda wildcards, input: ",".join([str(ab) for ab in ABUNDANCES]),
    shell:
        "python benchmarking/create_benchmarks.py "
        "--voc_perc {params.percs} "
        "-m {input.metadata} -fr {input.fasta} "
        "-fv {params.vocs} "
        "-o benchmarks/{wildcards.dataset} "
        "--total_cov 100 --spike_only "


rule run_kallisto:
    input:
        idx=expand("{ref}/sequences.kallisto_idx", ref=REF),
        wwsim1=expand(
            "benchmarks/{{dataset}}/wwsim_{voc}_ab{ab}_1.fastq",
            voc=pangolin,
            ab=ABUNDANCES,
        ),
        wwsim2=expand(
            "benchmarks/{{dataset}}/wwsim_{voc}_ab{ab}_2.fastq",
            voc=pangolin,
            ab=ABUNDANCES,
        ),
    output:
        expand(
            "benchmarks/{{dataset}}/out/{voc}_ab{ab}/predictions_m{min_ab}.tsv",
            voc=pangolin,
            ab=ABUNDANCES,
            min_ab=MIN_AB,
        ),
        dir=directory("benchmarks/{dataset}/out"),
    # script:
    #     "benchmarking/run_kallisto_snek.py"
    params:
        vocs=lambda wildcards, input: ",".join(VOC),
    run:
        outdir = output.dir
        shell("mkdir -p {outdir}")
        for voc in pangolin:
            for ab in ABUNDANCES:
                shell(
                    "kallisto quant -t {threads} -b {bootstraps} "
                    "-i {input.idx} -o {outdir}/{voc}_ab{ab} "
                    "benchmarks/{wildcards.dataset}/wwsim_{voc}_ab{ab}_1.fastq "
                    "benchmarks/{wildcards.dataset}/wwsim_{voc}_ab{ab}_2.fastq "
                    "| tee {outdir}/{voc}_ab{ab}.log"
                )
                shell(
                    "python pipeline/output_abuncances.py "
                    "-m {MIN_AB} "
                    "-o {outdir}/{voc}_ab{ab}/predictions_m{MIN_AB}.tsv "
                    "--metadata {REF}/metadata.tsv "
                    "--voc {params.vocs} "
                    "{outdir}/{voc}_ab{ab}/abundance.tsv "
                )
