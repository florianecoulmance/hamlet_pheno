#!/bin/bash
#SBATCH --job-name=smk_{rule}
#SBATCH --output=/fs/dss/work/doau0129/hamlet_pheno/logs/{rule}.%j.out
#SBATCH --error=/fs/dss/work/doau0129/hamlet_pheno/logs/{rule}.%j.err
#SBATCH --runtime={resources.runtime}
#SBATCH --mem={resources.mem_mb}
#SBATCH --cpus-per-task={threads}
#SBATCH --slurm_partition=rosa.p
#SBATCH --slurm_account=agfisheco

# Execute job
{exec_job}