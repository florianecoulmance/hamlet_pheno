#!/bin/bash
#SBATCH --job-name=smk_{rule}
#SBATCH --output=/fs/dss/work/doau0129/hamlet_pheno/logs/{rule}.%j.out
#SBATCH --error=/fs/dss/work/doau0129/hamlet_pheno/logs/{rule}.%j.err
#SBATCH --time={resources.time}
#SBATCH --mem={resources.mem_mb}
#SBATCH --cpus-per-task={threads}
#SBATCH --partition=rosa.p
#SBATCH --account=agfisheco

sbatch \
  -A agfisheco \
  -p rosa.p \
  --cpus-per-task={threads} \
  --mem={resources.mem_mb} \
  --time={resources.time} \
  --wrap "{exec_job}"