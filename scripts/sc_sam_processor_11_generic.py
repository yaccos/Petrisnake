#sam processor for directional UMI grouping
import sys
import os

sample = sys.argv[2]
if len(sys.argv) > 3:
    ID = sys.argv[3]
else:
    ID = sample
table = open(f'{ID}_selected_cumulative_frequency_table.txt')
version = 11
threshold = int(sys.argv[1])

output = open(f'{sample}_v{version}_threshold_{threshold}_filtered_mapped_UMIs.txt', 'w')
output.write('Cell Barcode\tUMI\tcontig:gene\ttotal_reads\n')

def is_ambiguous(line) -> bool:
    columns = line.split("\t")
    mapping_quality = int(columns[4])
    if mapping_quality == 0:
        return True
    elif "SA:Z" in line:
        return True
    else:
        return False

for cell in table:
    cell_barcode = cell.split('\t')[1]
    R2_file_name = f"{ID}_R2{cell_barcode[cell_barcode.find('_bc1_'):]}"
    R2_file_name = R2_file_name.replace(' ', '')
    
    samtools_command = f'samtools view {sample}_FC_directional_grouped_2/{R2_file_name}_group_FC.bam > {sample}_FC_directional_grouped_2/{R2_file_name}_group_FC.sam'
    os.system(samtools_command)
    
    sam = open(f'{sample}_FC_directional_grouped_2/{R2_file_name}_group_FC.sam')
    
    rm_command = f'rm {sample}_FC_directional_grouped_2/{R2_file_name}_group_FC.sam'
    os.system(rm_command)
    
    UMIgene_to_count = {}
    UMI_to_count = {}
    
    for line in sam:
        UMI = line[line.find('BX:Z:')+5:line.find('BX:Z:')+12]
        if line.split('\t')[1] == '4':
            contig = 'unaligned'
            gene = 'unaligned'
        else:
            edit_dist = line[line.find('NM:i:')+5:line.find('NM:i:')+6]
            contig = line.split('\t')[2]
            if 'XT:Z' in line:
                gene = line[line.find('XT:Z:')+5:len(line)+1]
                gene = gene.split('\t')[0].split('\n')[0]
                gene = gene.split('\t')[0]
                if ('rRNA' not in gene) and (is_ambiguous(line)):
                    gene = 'ambiguous'
            else:
                gene = 'no_feature'
        UMIgene = f'{UMI}:{contig}:{gene}'
        if UMIgene in UMIgene_to_count:
            UMIgene_to_count[UMIgene] += 1
        else:
            UMIgene_to_count[UMIgene] = 1

    for UMIgene in UMIgene_to_count:
        UMI = UMIgene.split(':')[0]
        contig_gene = f'{UMIgene.split(":")[1]}:{UMIgene.split(":")[2]}'
        if UMIgene_to_count[UMIgene] > threshold:
            output.write(f'{cell_barcode}\t{UMI}\t{contig_gene}\t{UMIgene_to_count[UMIgene]}\n')
