# long_read_path=/research/groups/ma1grp/home/zyu/work_2026/SV_2_Feb/sniffles2/fastq_fna/fastq_files/nanopore/sample_3_nanopore.fastq.gz
# fna_ref_path=/research/groups/ma1grp/home/zyu/work_2026/SV_2_Feb/sniffles2/fastq_fna/fna/sample_1.fna
# output_path=/research/groups/ma1grp/home/zyu/work_2026/SV_2_Feb/SV_pipeline/output
# sample_name="sample_3"
# long_quality=15
# stat_frequency_script=/research/groups/ma1grp/home/zyu/work_2026/SV_2_Feb/SV_pipeline/script/stat_breakpoints.py


# Parse command line arguments
while getopts "l:f:o:s:q:p:h" opt; do
    case $opt in
        l) long_read_path="$OPTARG" ;;
        f) fna_ref_path="$OPTARG" ;;
        o) output_path="$OPTARG" ;;
        s) sample_name="$OPTARG" ;;
        q) long_quality="$OPTARG" ;;
        p) stat_frequency_script="$OPTARG" ;;
        h) echo "Usage: $0 -l long_read_path -f fna_ref_path -o output_path -s sample_name [-q quality] [-p script_path]"; exit 0 ;;
        *) echo "Invalid option: -$OPTARG"; exit 1 ;;
    esac
done

# Validate required arguments
if [ -z "$long_read_path" ] || [ -z "$fna_ref_path" ] || [ -z "$output_path" ] || [ -z "$sample_name" ]; then
    echo "Error: Missing required arguments. Use -h for help."
    exit 1
fi


output_path_dir=${output_path}/${sample_name}
if [ ! -d ${output_path_dir} ]; then
  mkdir -p ${output_path_dir}
fi

#=================================================================
#+++++++++++++++++++++++Step 1 mapping Long ++++++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/1_mapping_Long
cd ${output_path}/1_mapping_Long

ln -s ${long_read_path} .

zcat ${long_read_path} | NanoFilt -q ${long_quality} > ${sample_name}_qc.fastq

seqkit stat *fastq *.gz > Long_read_stats.txt

minimap2 -ax map-ont ${fna_ref_path} ${sample_name}_qc.fastq | samtools sort -o ${sample_name}.bam

samtools index ${sample_name}.bam


#=================================================================
#+++++++++++++++++++++++Step 2 SV calling ++++++++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/2_SV_calling
cd ${output_path}/2_SV_calling

sniffles --input ${output_path}/1_mapping_Long/${sample_name}.bam --vcf ${sample_name}_SV_all.vcf --reference ${fna_ref_path} --all-contigs
sniffles --input ${output_path}/1_mapping_Long/${sample_name}.bam --vcf ${sample_name}_SV_chromosome.vcf --reference ${fna_ref_path} 

#=================================================================
#+++++++++++++++++++++++Step 3 Chimeric reads ++++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/3_Chimeric_reads
cd ${output_path}/3_Chimeric_reads


samtools view -h ${output_path}/1_mapping_Long/${sample_name}.bam | awk 'BEGIN {OFS="\t"} /^@/ || $0 ~ /SA:Z:/ {print $0}' > ${sample_name}_chimeric.sam
samtools view -bS ${sample_name}_chimeric.sam > ${sample_name}_chimeric.bam
samtools index ${sample_name}_chimeric.bam

#=================================================================
#+++++++++++++++++++++++Step 4 SV calling by chimeric reads ++++++
#=================================================================
mkdir -p ${output_path}/4_SV_calling_chimeric
cd ${output_path}/4_SV_calling_chimeric
sniffles --input ${output_path}/3_Chimeric_reads/${sample_name}_chimeric.bam --vcf ${sample_name}_SV_chimeric_all.vcf --reference ${fna_ref_path} --all-contigs
sniffles --input ${output_path}/3_Chimeric_reads/${sample_name}_chimeric.bam --vcf ${sample_name}_SV_chimeric_chromosome.vcf --reference ${fna_ref_path}

#=================================================================
#+++++++++++++++++++++++Step 5 stat frequency ++++++++++++++++++++
#=================================================================
mkdir -p ${output_path}/5_SV_frequency
cd ${output_path}/5_SV_frequency


cat ${output_path}/4_SV_calling_chimeric/${sample_name}_SV_chimeric_chromosome.vcf |grep -v "^#" |grep "PRECISE" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$8}' | sed 's/PRECISE.*END=//g' | sed 's/;.*//g' > ${sample_name}_SV_chromosome.bed

count=1
cat ${sample_name}_SV_chromosome.bed | while read line; do
    contig_a=$(echo $line | cut -f1 -d' ')
    breakpoint_a=$(echo $line | cut -f2 -d' ')
    sv_type=$(echo $line | cut -f3 -d' '| cut -f2 -d'.')
    forth_field=$(echo $line | cut -f4 -d' ')
    fifth_field=$(echo $line | cut -f5 -d' ')
    echo "sv_type: ${sv_type}"
    if [[ ${sv_type} == "BND" ]]; then
        contig_b=$(echo ${forth_field} | cut -f1 -d':'|sed 's/^.*\]//g')
        breakpoint_b=$(echo ${forth_field} | cut -f2 -d':' | sed 's/\].*//g')
    else
        contig_b=${contig_a}
        breakpoint_b=${fifth_field}
    fi
    python stat_breakpoints.py --bam ${output_path}/1_mapping_Long/${sample_name}.bam --contig-a ${contig_a} --breakpoint-a ${breakpoint_a} --contig-b ${contig_b} --breakpoint-b ${breakpoint_b} --slop 5 --span-cutoff 15 --sample-name ${sample_name}-${count}-${sv_type}
    count=$((count + 1))
done 





python stat_breakpoints.py \
  --bam ${sample_name}.bam --contig-a contig_1 --breakpoint-a 573076 \
  --contig-b contig_3 --breakpoint-b 103809 \
  --slop 5 --span-cutoff 15 \
  --sample-name test