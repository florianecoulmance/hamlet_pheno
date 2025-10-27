#!/bin/bash
#SBATCH --job-name=smk_{rule}
#SBATCH --output=logs/{rule}.{wildcards}.out
#SBATCH --error=logs/{rule}.{wildcards}.err
#SBATCH --time={resources.time}
#SBATCH --mem={resources.mem_mb}
#SBATCH --cpus-per-task={threads}
#SBATCH --partition=rosa.p


# Execute job
{exec_job}
