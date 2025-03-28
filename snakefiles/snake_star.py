import glob
import numpy as np

work_dir = "work_dir='/path/to/work_dir/"
index_dir = "/path/to/STAR/GRCm39"

samples = [s.split('/')[-1].split('_R')[0] for s in glob.glob(work_dir + '/reads/raw_filtered/*_R1.fastq.gz')]
samples.sort()

print(samples)

rule run_star:
    input: expand("{work_dir}/star/{sample}_Aligned.sortedByCoord.out.bam", sample = samples, work_dir = work_dir, inedx_dir = index_dir)

rule star:
    input:
        index = index_dir,
        R1 = '{work_dir}/reads/raw_filtered/{sample}_R1.fastq.gz',
        R2 = '{work_dir}/reads/raw_filtered/{sample}_R2.fastq.gz'
    threads: 10
    conda: "/data12/bio/runs-jeniaole/snake-files/envs/star.yaml"
    params: prefix = '{work_dir}/star/{sample}_'
    output: '{work_dir}/star/{sample}_Aligned.sortedByCoord.out.bam'
    shell: 'STAR --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --genomeDir {input.index} --runThreadN {threads} --readFilesCommand zcat --readFilesIn {input.R1} {input.R2} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard'
