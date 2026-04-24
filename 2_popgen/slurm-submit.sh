#!/bin/bash

sbatch \
  -A agfisheco \
  -p rosa.p \
  --cpus-per-task={threads} \
  --mem={resources.mem_mb} \
  --time={resources.time} \
  --wrap "{exec_job}"