# Bacterial_SV_pipeline
Bacterial_SV_pipeline is a workflow for analyzing bacterial structure variation based on reads-level.

# 1 Installing Conda 
We recommend installing miniconda to creat a enviroment for this pipeline.
## 1-1 install conda
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_25.3.1-1-Linux-x86_64.sh /miniconda.sh
bash miniconda.sh -p ~/miniconda3

```

## 1-2 install mamba
```
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

**Command-line Arguments**

* `-l` : Path to the Nanopore reads FASTQ file
* `-f` : Reference genome file
* `-o` : Output directory path
* `-s` : Sample name
* `-q` : Quality score cutoff for QC filtering
* `-p` : Path to `stat_breakpoints.py` script
