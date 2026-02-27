# Bacterial_SV_pipeline

Bacterial_SV_pipeline is a workflow for analyzing bacterial structural variations at the read level.

The workflow diagram is shown below:

![workflow](/image/workflow.png)

---

## 1. Installing Conda

We recommend installing Miniconda to create an isolated environment for this pipeline.

### 1.1 Install Conda

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_25.3.1-1-Linux-x86_64.sh -O miniconda.sh  
bash miniconda.sh -p ~/miniconda3
```


### 1.2 Install Mamba

Mamba is a fast reimplementation of Conda that significantly improves package installation speed.
```bash
conda install mamba==2.0.5 -c conda-forge
```

---

## 2. Create Environment for SV_pipeline

All required software dependencies are listed in the YAML file. You can create the environment and install all tools with:
```bash
mamba env create -f sv_pipeline_env.yml
```
After creating the environment, activate it:
```bash
conda activate sv_pipeline
```
---

## 3. Fill meta_data.csv
You just need to enter your specific parameters into **meta_data.csv**, and then run the pipeline.
- Sample_name : Sample name  
- Fastq_path : Path to Nanopore reads FASTQ file  
- Ref_path: Path to Reference genome file  
- Output_path: Path to Output directory  
- QC_quality: Quality score cutoff for QC filtering  
- span_cutoff: Minimum bp on EACH side of the breakpoint to count as spanning
- breakpoint_distance: Distance for breakpoint detection

---

## 4. Running SV_pipeline Analysis

```bash
sh run_pipeline.sh
```


