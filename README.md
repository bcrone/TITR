# TITR
TITR (Trans-ancestral Iterative Tissue Refinement)
This repository encapsulates two PRS models prioritizing functional genomic data aiming to improve trans-ancestral PRS portability.
## Installation
git clone https://github.com/bcrone/TITR
## Tool Requirements
LDSC: https://github.com/bulik/ldsc  
PLINK v1.9: https://www.cog-genomics.org/plink/1.9/  
PLINK v2.0: https://www.cog-genomics.org/plink/2.0/  
Python 3.9: https://www.python.org/downloads/release/python-390/  
Snakemake 8.2.1: https://snakemake.readthedocs.io/  
## Data Requirements
1. LDSC-formatted GWAS summary statistics: columns rsID, A1, A2, N, Z
2. PLINK-formatted GWAS summary statistics: columns SNP, A1, BETA, P
3. LDscore estimates per annotation (see LDSC wiki for tutorial)
4. 1000 Genomes Phase 3 European Reference Panel (Availble from 1000 Genomes FTP site or from LDSC tutorial)
5. UKBiobank training genotypes/phenotypes (White British European) + testing genotypes/phenotypes (African) (UKBiobank access required - https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access)
6. RegulomeDB-derived functional probablistic models (SURF/TURF/TLand) - resource TBD
## Computational Requirements
Both single tissue and TITR models were developed in a SLURM-managed HPC environment   
HPC nodes were installed with Redhat Enterprise Linux 8 operating system  
Processes were run on nodes with following hardware:  
- Processor: 2x 3.0 GHz Intel Xeon Gold 6154 (36 cores/node)
- Memory: 1.5TB RAM (65GB max requested)
Time requested varies from 5-10hrs/job - trait dependent
## Directory Structure
```
ROOT_PATH
├── data
│   ├── 1KG
│   │   ├── 1000G_EUR_Phase3_plink
│   │   ├── 1000G_Phase3_frq
│   │   ├── 1000G_Phase3_weights_hm3_no_MHC
│   │   ├── hapmap3_snps
│   │   ├── RegDB_custom_baselineLD_v1.2
│   │   │   └── serial
│   │   │       └── quantiles
│   │   └── TLand_custom_baselineLD_v1.2
│   │       └── serial
│   │           └── quantiles
│   ├── GWAS
│   │   ├── <traits>
│   ├── IMPACT
│   │   ├── <traits>
│   ├── LDSC
│   ├── RegDB
│   │   ├── <tissues>
│   ├── TLand
│   │   ├── <tissues>
│   └── UKB
│       ├── AFR
│       ├── EUR
│       └── phenos
│           ├── AFR
│           └── EUR
├── repo
│   └── TITR
│       ├── profile
│       │   └── slurm
│       ├── single-tissue-model
│       │   ├── logs
│       │   └── slurm
│       ├── TITR
│       └── snakemake
│           ├── config
│           └── workflow
│               └── logs
└── results
    ├── <tissues>
    ├── coding
    │   ├── <traits>
    ├── h2
    │   ├── <traits>
    ├── single_tissue
    │   └── <traits>
    │       ├── IMPACT
    │       ├── standard
    │       ├── SURF
    │       ├── TURF
    │       └── TLand
    └── traits
        └── <traits>
```
## Single Tissue Model
The single tissue functional prioritization model is implemented in a Snakemake workflow, developed for use on SLURM-managed HPC environment. 
### Workflow overview
The Snakemake workflow is as follows:  
- ```calculate_LDscores```: Calculate LDscores for each tissue model per chromsome (this is only done once)  
- ```calculate_h2```: Calculate partitioned heritability by TURF/TLand tissue models for each trait  
- ```calculate_tau_star```: Calculate tau star estimates by TURF/TLand tissue models for each trait  
- ```select_lead_tissue```: Identify TURF/TLand tissue model with greatest significant non-zero tau star estimate for each trait  
- ```run_PLINK_scores```: Calculate PRS for standard model and top 1/5/10/20/50% of functional models (IMPACT/SURF/TURF/TLand)  
- ```bootstrap_PLINK_scores```: Generate bootstraps (n=1000) for best p-value partition for all PRS models  
### Configuration
- ```snakemake/config/cluster_config.yml```: update account and partitions to system  
- ```snakemake/workflow/config.yml```: update paths to local installs/data paths  
- ```snakemake/workflow/Snakefile```: update traits lists for target traits  
- ```profile/slurm```: copy directory to ~/.config/snakemake  
### Usage 
For a dry-run of Snakemake:  
```
snakemake -np --profile slurm
```
To run pipeline:  
```
snakemake --profile slurm --latency-wait 300
```
## Multiple Tissue Model (TITR)
The TITR algorthim is implemented in a custom-developed workflow to allow for iterative modelling and analysis in a SLURM-managed HPC environment.  
### Workflow overview
The TITR algorithm is as follows:  
- ```build_run_ldsc.sh/batch_run_ldsc.sh```: calculate partitioned trait heritability by tissue model
- ```build_run_summarize_enrichment_tau_star.sh/run_summarize_enrichment_tau_star.sh/summarize_enrichment_tau_star.py```: summarize h2 enrichments for all tissues, calculate tau star estimates, identify significant non-zero tissue estimates as candidate tissues for PRS calculation
- ```build_run_score.sh/batch_run_score.sh```: calculate PRS models for candidate tissues
- ```build_run_optimize_TURF_PRS.sh/run_optimize_TURF_PRS.sh/optimize_TURF_PRS.py```: select lead PRS model on each iteration, test for terminating condition, integrate model SNPs into master PRS model
- ```clean_iteration.sh```: remove intermediate files from iteration, determine if terminating condition met, begin next iteration
- ```build_validate_score.sh/batch_validate_score.sh```: calculate PRS model for validation samples
- ```build_run_validate_TURF_PRS.sh/run_validate_TURF_PRS.sh/validate_TURF_PRS.py```: construct TITR-optimized PRS models for validation samples
### Configuration
- SLURM account information: update parameters in .sh scripts
- GWAS binary vs. quantitative traits: update ISQUANT params in ```build_run_optimize_TURF_PRS.sh``` and ```build_run_validate_TURF_PRS.sh```
- Max iterations: update ```clean_iteration.sh``` to terminate TITR after specified max iteration count
### Usage
To run TITR model on training data:  
```bash build_run_ldsc.sh <trait> 1 <root_path>``` (1 indicates 1st iteration of TITR)  
To construct scores of TITR-optimized PRS model for validation cohort:  
```bash build_validate_score.sh <ancestry> <trait> <iteration> <root_path>``` (ancestry: ancestry of validation cohort; iteration: max iteration of TITR)  
To evaluate TITR-optimized PRS model for validation scores:  
```bash build_run_validate_TURF_score.sh <ancestry> <trait> <iteration> <root_path>```
