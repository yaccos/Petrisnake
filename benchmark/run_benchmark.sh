script_dir=$(pwd)/../scripts
sample=EC_3h_deep
sample_full=EC_3h_deep_S5
bc_cutoff=60000
annotation_dir=GCF_000005845.2

# Makes sure the cutadapt jobs can run in parallel
cutadapt -g NNNNNNNNAGAATACACGACGCTCTTCCGATCT --cores 3 -o ${sample}/${sample_full}_L001_R1_001.fastq.gz ${sample}/SRR28148450_1.fastq &
cutadapt -g NNNNNNNNAGAATACACGACGCTCTTCCGATCT --cores 3 -o ${sample}/${sample_full}_L001_R2_001.fastq.gz ${sample}/SRR28148450_2.fastq &
wait

python $script_dir/sc_pipeline_15_generic_v2.py ${sample_full} 1
chmod +x $script_dir/pipeline_v2_generic.sh
$script_dir/pipeline_v2_generic.sh ${sample} $bc_cutoff GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna  GCF_000005845.2/genomic.gtf ${sample}
