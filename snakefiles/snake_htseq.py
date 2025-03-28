import glob
import numpy as np

work_dir = "work_dir='/path/to/work_dir/"
annot_dir = "/path/to/STAR/GRCm39"

samples = [s.split('/')[-1].split('_R1.fastq.gz')[0] for s in glob.glob(work_dir + '/reads/raw_filtered/*R1.fastq.gz')]
samples.sort()

print(samples)

rule run_htseq:
    input: expand('{work_dir}/htseq/{sample}.counts', sample = samples, work_dir = work_dir, annot_dir = annot_dir)

rule htseq:
    input: 
        bam = "{work_dir}/star/{sample}_Aligned.sortedByCoord.out.bam", 
        annotation = annot_dir + "/gencode.vM36.chr_patch_hapl_scaff.annotation.gtf"
    output: "{work_dir}/htseq/{sample}.counts"
    conda: "/data12/bio/runs-jeniaole/snake-files/envs/htseq.yaml"
    threads: 1
    shell: "htseq-count -n {threads} --idattr gene_id -s no --add-chromosome-info --additional-attr=ID {input.bam} {input.annotation} > {output}"

