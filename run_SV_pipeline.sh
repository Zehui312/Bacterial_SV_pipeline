meta_data=./meta_data.csv

cat $meta_data |tail -n +2 | while read line; do
  sample_name=$(echo $line | cut -d ',' -f 1)
  fastq_path=$(echo $line | cut -d ',' -f 2)
  ref_path=$(echo $line | cut -d ',' -f 3)
  output_path=$(echo $line | cut -d ',' -f 4)
  pipline_path=$(echo $line | cut -d ',' -f 5)
  QC_quality=$(echo $line | cut -d ',' -f 6)
  span_cutoff=$(echo $line | cut -d ',' -f 7)
  breakpoint_distance=$(echo $line | cut -d ',' -f 8)
  echo "sh $pipline_path/script/sniffles2_analysis.sh -l $fastq_path -f $ref_path -o $output_path -s $sample_name -q $QC_quality -p $pipline_path/script/stat_breakpoints.py -a $span_cutoff -b $breakpoint_distance"
done > run_pipeline.sh

mkdir -p ./logs/
count=1
while read runcode; do
    bsub -P SV_${count} -J SV_${count} -n 2 -R "rusage[mem=8GB]" -eo ./logs/SV_${count}.err -oo ./logs/SV_${count}.out $runcode
    count=$((count + 1))
done < run_pipeline.sh