## edited May 26 2020 (after upload) to handle different feature label in gff

import sam_edit_tools
import sys
import os
import multiprocessing

def pipeline(sample,ID,name,gff):
    R2_file_name = ID + '_R2' + name[name.find('_bc1_'):len(name)]
    R2_file_name = R2_file_name.replace(' ' , '')
    R2_short_name = R2_file_name.replace(ID+'_','')
    sam_edit_tools.no_xt_new(old_sample + '_bwa_sam/1/' + R2_short_name + '/stdout', sample + '_no_XT/' + R2_file_name + '_no_XT.sam')
    samtools_view_command = (
        f'samtools view -bS {sample}_no_XT/{R2_file_name}_no_XT.sam | '
        f'samtools sort - > {sample}_no_XT/{R2_file_name}_no_XT_sorted.bam'
    )
    os.system(samtools_view_command)
    
    samtools_index_command_1 = (
        f'samtools index {sample}_no_XT/{R2_file_name}_no_XT_sorted.bam'
    )
    os.system(samtools_index_command_1)
    
    featureCounts_command = (
        f"featureCounts -t 'gene' -g 'Name' -s 1 -a {gff} -o {sample}_FC/{R2_file_name}_FC "
        f"-R BAM {sample}_no_XT/{R2_file_name}_no_XT_sorted.bam"
    )
    os.system(featureCounts_command)
    
    samtools_index_command_2 = (
        f'samtools index {sample}_FC/{R2_file_name}_no_XT_sorted.bam.featureCounts.bam'
    )
    os.system(samtools_index_command_2)
    
    umi_tools_command = (
        f'umi_tools group --per-gene --gene-tag=XT -I {sample}_FC/{R2_file_name}_no_XT_sorted.bam.featureCounts.bam '
        f'--group-out={sample}_FC_directional_grouped_2/{R2_file_name}_UMI_counts.tsv '
        f'--method=directional --output-bam -S {sample}_FC_directional_grouped_2/{R2_file_name}_group_FC.bam '
        f'>> {old_sample}_logs/featureCounts_directional_5/{sample}_umi_group.log'
    )
    os.system(umi_tools_command)

ID = sys.argv[1] 
if len(sys.argv) > 2:
    sample = sys.argv[2]
else:
    sample = ID
if len(sys.argv) > 3:
    old_sample = sys.argv[3]
os.system('mkdir ' + sample + '_no_XT/')
gff = sys.argv[4]
table = open(ID+'_selected_cumulative_frequency_table.txt')

os.system('mkdir ' + sample + '_FC')
os.system('mkdir ' + sample + '_FC_directional_grouped_2')
os.system('mkdir ' + old_sample + '_logs/featureCounts_directional_5')
jobs = []
pool =  multiprocessing.Pool(10)
for line in table:
    name = line.split('\t')[1]
    pool.apply_async(pipeline,args=(sample,ID,name,gff))
pool.close()
pool.join()
