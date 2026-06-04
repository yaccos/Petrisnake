## edited May 26 2020 so logs folder is made if not already there

dir=$(echo $0 | sed "s/pipeline_v2_generic.sh//")
sample=$1
n_BCs=$2
fasta=$3
gff=$4
custom_name=$5

if [ ! -d "${sample}_logs" ]
then
	mkdir ${sample}_logs
fi

if [ -d "${sample}_bc1" ]
then
	python $dir/remove_cells_v1.py ${sample} ${n_BCs} # Select cells based on threshold (knee plot); n_BCs is the number of cells to use for additional analysis (typically 10000-40000)
	echo "selecting ${n_BCs} cells"
	num=1
	msg="selected ${n_BCs} cells"
elif [ -d "${sample}_selected_cells" ]
then
	echo "Directory ${sample}_bc1 does not exist but ${sample}_selected_cells already exists; using previously selected cells"
	num=1	
	msg="PLEASE NOTE: used previously selected cells, n_BCs of ${n_BCs} ignored!"
else
	echo "Directory ${sample}_bc1 does not exist. Maybe need to run sc_pipeline_11.py?"
	num=0
fi
if [ $num -gt 0 ]
then
	python $dir/trim_R2_v4.py ${sample} # trim bc1 and linker sequences from read 2; discard reads less than 17 bp
	python $dir/hairpin_2nd_trim.py ${sample} # trim revcomp bc1 and linker sequences (from hairpins)
	if [ ! -f "${fasta}.bwt" ]
	then
		bwa index ${fasta}
	fi
	python $dir/align_HiSeq_after_hairpin_trim_v2.py ${sample} ${fasta} ${sample} # align read 2 to genomes
	python $dir/featureCounts_directional_5.py ${sample} ${custom_name} ${sample} ${gff} # annotate features from gff and identify UMI groups
	python $dir/sc_sam_processor_11_generic.py 0 ${custom_name} ${sample} # generates a single file of collapsed UMIs (output suffix is _filtered_mapped_UMIs.txt})
	python $dir/make_matrix_mixed_species.py ${custom_name}_v11_threshold_0 # make matrix from filtered UMIs
	## Below are cleanup steps - comment out to see intermediate files
	# rm -r ${sample}/${sample}*QF*.fastq
	# rm -r ${sample}_R2_trimmed #fastq files after trim_R2_v4
	# rm -r ${sample}_2nd_trim
	# rm -r ${sample}_bwa_sam # sam files after align_v4
	# rm -r ${custom_name}_no_XT # sam and bam files without XT tag (generated during featureCounts_directional_5)
	# rm -r ${custom_name}_FC # bam files after feature counts (featureCounts_directional_5)
	# rm -r ${custom_name}_FC_directional_grouped_2 # bam files after umi tools group (featureCounts_directional_5)
	echo $msg
fi
