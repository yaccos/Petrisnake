## Edited so maximum edit distance for bwa aln is 1

import sys
import os

ID = sys.argv[1]
if len(sys.argv) > 3:
    sample = sys.argv[3]
else:
    sample = ID

template = sys.argv[2]
table = open(ID+'_selected_cumulative_frequency_table.txt')

os.system(f'mkdir {sample}_bwa_sam')
R2_list = {}
R2_list[0] = ''
i = 0
group = 0
for line in table:
    i += 1
    name = line.split('\t')[1]
    R2_file_name = f'R2{name[name.find("_bc1_"):len(name)]}'
    R2_file_name = R2_file_name.replace(' ', '')
    if R2_list[group] == '':
        R2_list[group] = R2_file_name
    else:
        R2_list[group] = f'{R2_list[group]}\n{R2_file_name}'
    if  i%4000 == 0:
        group +=1
        R2_list[group] = ''
for group in R2_list:
    bwa_mem_command = (
        f'echo "{R2_list[group]}" | time parallel --bar --files --results '
        f'{sample}_bwa_sam -j10 bwa mem -a -k 19 -Y {template} '
        f'{sample}_2nd_trim/{ID}_{{}}_2trim.fastq'
    )
    os.system(bwa_mem_command)
#     bwa_samse_command = (
#         f'echo "{R2_list[group]}" | time parallel --bar --files --results '
#         f'{sample}_bwa_sam -j10 bwa samse -n 14 {template} '
#         f'{sample}_bwa_sai/1/{{}}/stdout {sample}_2nd_trim/{ID}_{{}}_2trim.fastq'
#     )
#     os.system(bwa_samse_command)

# os.system(f'rm -r {sample}_bwa_sai/')
