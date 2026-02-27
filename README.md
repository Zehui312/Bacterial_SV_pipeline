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

## 3. Running a Test Analysis

You just need to enter your specific parameters into **meta_data.csv**, and then run the pipeline.

```bash
sh run_pipeline.sh
```
---

## Command-line Arguments

- -l : Path to Nanopore reads FASTQ file  
- -f : Reference genome file  
- -o : Output directory  
- -s : Sample name  
- -q : Quality score cutoff for QC filtering  
- -p : Path to stat_breakpoints.py script