# ---------------------------------------------------------------------------------------
# HPC settings. EDIT THIS FILE (once), everything in bin/hpc/ reads it.
# It's used on the Mac (to know where to upload) and on Ada (to know where things are).
# ---------------------------------------------------------------------------------------

# How you log in: ssh HPC_USER@HPC_HOST
HPC_USER="nehar"
HPC_HOST="ada.nri.bcm.edu"

# Work folder on Ada. Needs ~60 GB. No spaces in this path, please.
# Empty = "cln3_msi" in your home folder on Ada. If your home folder is small (check with
# `df -h ~` on Ada, see docs/ADA_RUN_GUIDE.md), put a scratch/project path here instead,
# e.g. HPC_WORK="/scratch/nehar/cln3_msi"
HPC_WORK=""

# Where the imzML data sit on Ada (the folder that CONTAINS "imzML files for biostats").
# Empty = HPC_WORK/data, which is where upload_to_hpc.sh puts it. A second user who reads
# the data from someone else's folder sets it here, e.g. HPC_DATA="/home/nehar/cln3_msi/data"
# (the data are only read, never written; see docs/ADA_RUN_GUIDE.md, "Running as another user").
HPC_DATA=""

# SLURM: leave empty to use the cluster defaults. Ask the HPC team / docs if unsure.
SLURM_PARTITION=""          # e.g. "defq", "compute"
SLURM_ACCOUNT=""            # only if Ada makes you charge jobs to an account / project
SLURM_QOS=""                # only if Ada uses QOS

# Only needed if Ada has conda as a module (otherwise setup installs Miniforge in HPC_WORK).
CONDA_MODULE=""             # e.g. "miniforge3" or "anaconda3"; check with: module avail conda

# On the Mac: the folder that contains "imzML files for biostats"
# (the Mac ignores upper/lower case, so "Brain - Pilar - Spatial" points to the same folder)
LOCAL_DATA_ROOT="/Users/neharao/Desktop/BRAIN - PILAR - SPATIAL"

# Which modes to run (one chain of jobs per mode, they run side by side)
MODES="lip_neg lip_pos met_neg"
