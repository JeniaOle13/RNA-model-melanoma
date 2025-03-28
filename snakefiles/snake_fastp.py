import glob
import numpy as np

work_dir='/path/to/work_dir/'

samples = [s.split('/')[-1].split('_R1.fastq.gz')[0] for s in glob.glob(work_dir + '/reads/raw/*_R1.fastq.gz')]
samples.sort()
print(samples)

rule run_fastp:
    input: expand("{work_dir}/reads/raw_filtered/{sample}_R1.fastq.gz", sample = samples, work_dir = work_dir)

rule fastp:
    input:
        R1 = "{work_dir}/reads/raw/{sample}_R1.fastq.gz",
        R2 = "{work_dir}/reads/raw/{sample}_R2.fastq.gz"
    threads: 10
    conda: "/data12/bio/runs-jeniaole/snake-files/envs/fastp.yaml"
    output:
        R1 = "{work_dir}/reads/raw_filtered/{sample}_R1.fastq.gz",
        R2 = "{work_dir}/reads/raw_filtered/{sample}_R2.fastq.gz",
        html = "{work_dir}/reports/fastp/{sample}.html"
    shell: 'fastp -w {threads} --detect_adapter_for_pe --dedup -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --html {output.html}'
