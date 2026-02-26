# Bacterial_SV_pipeline
Bacterial_SV_pipeline is a workflow for analyzing bacterial structure variation based on reads-level.

# 1 Installing Conda 
We recommend installing miniconda to creat a enviroment for this pipeline.

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
conda install mamba==2.0.5
```

# 2 Creat enviroment for SV_pipeline

```
mamba env create -f sv_pipeline_env.yml
```

# 3 Running test 

```
sh sniffles2_analysis.sh \
-l /research/groups/ma1grp/home/zyu/work_2026/Pipeline_data/SV_pipeline_data/sample_3_nanopore.fastq.gz \
-f /research/groups/ma1grp/home/zyu/work_2026/Pipeline_data/SV_pipeline_data/sample_1.fna \
-o /research/groups/ma1grp/home/zyu/work_2026/SV_2_Feb/SV_pipeline/output \
-s F2 \
-q 15 \
-p /research/groups/ma1grp/home/zyu/work_2026/SV_2_Feb/SV_pipeline/script/stat_breakpoints.py
```
