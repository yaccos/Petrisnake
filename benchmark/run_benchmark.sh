script_dir=$(pwd)/../scripts
sample=EC_3h_deep
sample_full=EC_3h_deep_S5
bc_cutoff=60000
annotation_dir=1_06SaEc_ref

python $script_dir/sc_pipeline_15_generic_v2.py ${sample_full} 1
chmod +x $script_dir/pipeline_v2_generic.sh
$script_dir/pipeline_v2_generic.sh ${sample} $bc_cutoff GCF_000005845.2/1_06SaEc_sequence.fa  GCF_000005845.2/1_06SaEc_annotation.fa ${sample}
