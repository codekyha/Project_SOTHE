#!/bin/bash
# SOTHE-P2-RAD on a SLURM cluster (written for UHeM Altay; see slurm/README_SLURM.md).
# Replace ACCOUNT_NAME and PARTITION_NAME, then: sbatch slurm/rad_altay.sbatch  (from the folder SOTHE-P2-RAD)
#SBATCH -J sothe-rad
#SBATCH -A ACCOUNT_NAME
#SBATCH -p PARTITION_NAME
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 44
#SBATCH -t 0-01:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err
set -uo pipefail
if ! type module >/dev/null 2>&1; then
  for f in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /usr/share/lmod/lmod/init/bash /etc/profile.d/modules.sh; do
    [ -f "$f" ] && . "$f" && break
  done
fi
module load ANACONDA/Anaconda3-2024.06-1-python-3.12 || echo "WARNING: module load failed; using the python on PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
python rad.py verify || exit 1
python rad.py run --out results --set all --workers "${SLURM_CPUS_PER_TASK:-44}" || exit 1
R=$(ls -d results/*/ | sort | tail -1)
python rad.py analyse --results "$R" && python rad.py figures --results "$R" && python rad.py manifest --results "$R"
