script_dir=$(pwd)/../scripts
sample=SA_9h_deep
sample_full=SA_9h_deep_S5
bc_cutoff=8000
annotation_dir=GCF_000756205.1

# Makes sure the cutadapt jobs can run in parallel
cutadapt -g NNNNNNNNAGAATACACGACGCTCTTCCGATCT --cores 3 -o ${sample}/${sample_full}_L001_R1_001.fastq.gz ${sample}/SRR28148452_1.fastq &
cutadapt -g NNNNNNNNAGAATACACGACGCTCTTCCGATCT --cores 3 -o ${sample}/${sample_full}_L001_R2_001.fastq.gz ${sample}/SRR28148452_2.fastq &
wait

python $script_dir/sc_pipeline_15_generic_v2.py ${sample_full} 1
chmod +x $script_dir/pipeline_v2_generic.sh
$script_dir/pipeline_v2_generic.sh ${sample} $bc_cutoff GCF_000756205.1/GCF_000756205.1_ASM75620v1_genomic.fna  GCF_000756205.1/genomic.gff ${sample}
