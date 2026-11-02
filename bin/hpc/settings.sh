# ---------------------------------------------------------------------------------------
# HPC settings. EDIT THIS FILE (once), everything in bin/hpc/ reads it.
# It's used on the Mac (to know where to upload) and on Ada (to know where things are).
# ---------------------------------------------------------------------------------------

# How you log in: ssh HPC_USER@HPC_HOST
HPC_USER="your_ada_username"
HPC_HOST="ada.your-institution.edu"

# Work folder on Ada: on the big storage (scratch / project space), NOT your home folder.
# Needs ~60 GB. No spaces in this path, please.
HPC_WORK="/scratch/your_ada_username/cln3_msi"

# SLURM: leave empty to use the cluster defaults. Ask the HPC team / docs if unsure.
SLURM_PARTITION=""          # e.g. "defq", "compute"
SLURM_ACCOUNT=""            # only if Ada makes you charge jobs to an account / project
SLURM_QOS=""                # only if Ada uses QOS

# Only needed if Ada has conda as a module (otherwise setup installs Miniforge in HPC_WORK).
CONDA_MODULE=""             # e.g. "miniforge3" or "anaconda3"; check with: module avail conda

# On the Mac: the folder that contains "imzML files for biostats"
LOCAL_DATA_ROOT="/Users/neharao/Desktop/BRAIN - PILAR - SPATIAL"

# Which modes to run (one chain of jobs per mode, they run side by side)
MODES="lip_neg lip_pos met_neg"
