# Running the CLN3 MSI pipeline on Ada: every command, in order

Commands are typed in the **Mac Terminal** unless a step says **on Ada**. You can tell where
you are from the prompt: on Ada it looks like `nehar@ada:~$`. Lines starting with `#` are
comments, you don't type them.

## Ada, as checked on your account

| | |
|---|---|
| Login | `ssh nehar@ada.nri.bcm.edu` (BCM VPN if off campus) |
| Work folder | `~/cln3_msi` = `/home/nehar/cln3_msi` (`/home` has 25 TB free, no personal quota shown) |
| SLURM | one partition `defq` (the default), one machine, no time limit, so nothing to set |
| Conda | shared Anaconda at `/home/anaconda/bin/conda`; our env lives in `~/cln3_msi/envs/cln3-msi` |
| Internet | yes (downloads work) |
| Code on the Mac | `/Users/neharao/Downloads/lipidomics-testbed-main` (unzipped GitHub download) |
| Data on the Mac | `/Users/neharao/Desktop/BRAIN - PILAR - SPATIAL` (same folder as "Brain - Pilar - Spatial") |

`bin/hpc/settings.sh` is already set up for all of this, so you don't need to edit it.

Layout on Ada:

```
~/cln3_msi/
  repo/                          the code (a mirror of the Mac folder; edit on the Mac, not here)
  data/imzML files for biostats/ the data (36 files, ~40 GB)
  envs/cln3-msi/                 the Python environment
  cln3_env.sh                    loads the environment (every job uses it)
  outputs/                       everything the pipeline writes (summaries/, figures/, report/, ...)
  slurm_logs/                    one log per job
```

**Only `imzML files for biostats/` gets uploaded.** These all stay on the Mac and aren't used:

- SCiLS projects (`.slx/.sbd`, `raw lip neg/pos`, `raw met neg/`, `Female Data - Pilar analysis/`)
- the `.srd` and `.mca` files
- the `export_All_*` files
- Pilar's results (`Results_CLN3_Analysis*/`, `PC 38.1 - depleted.png`, `.Rhistory`, `.RDataTmp`)
- the `09SEP26/` cerebellum run

---

## Where you are now

You have:

- uploaded the code: `~/cln3_msi/repo` exists on Ada
- run the setup: the prompt shows `(cln3-msi)`

Next come steps 1 to 7 below.

### Step 1. On Ada: confirm the setup finished

```bash
ls ~/cln3_msi/cln3_env.sh                           # must exist
ls ~/cln3_msi/outputs/cache                         # want: atlas  hmdb  lipidmaps
grep -A6 "Reference databases" ~/cln3_msi/outputs/summaries/00_setup.md
```

* If all three show up (the `grep` lists lipidmaps and hmdb as cached), go to step 2.
* If something is missing, run the setup (again). It's safe to rerun and quick when things
  are already there:

```bash
cd ~/cln3_msi/repo
bash bin/hpc/hpc_setup.sh          # want "20 passed" and "== setup OK" at the end
```

It warns "Data folder check failed" until the data is uploaded in step 4. That's expected.

### Step 1b. HMDB, the metabolite database (optional)

HMDB's website refuses automatic downloads from Ada ("HTTP Error 403"). Without HMDB the
metabolite mode uses the small curated list in the repo; the lipid modes use LIPID MAPS and
aren't affected. To add HMDB anyway:

1. In your Mac's browser go to https://hmdb.ca/downloads, find the metabolite structures file (SDF; labelled something like "Metabolite Structures")
   (SDF), and download `structures.zip` (it lands in Downloads).
2. Copy it to Ada and load it, from the Mac Terminal:

```bash
scp ~/Downloads/structures.zip nehar@ada.nri.bcm.edu:~/cln3_msi/
ssh nehar@ada.nri.bcm.edu "source ~/cln3_msi/cln3_env.sh && cd ~/cln3_msi/repo && python scripts/00_setup.py --config config.hpc.yaml --hmdb-file ~/cln3_msi/structures.zip | tail -n 5"
```

The summary should now list `hmdb: ... entries with formulas cached`.

### Step 2. On Ada: how big is the machine? (optional, 5 seconds)

```bash
nproc               # number of CPU cores
free -g | head -2   # memory in GB ("total" column)
```

Each mode's job asks for 4 cores and 32 GB:

* **At least 12 cores and ~100 GB:** all three modes run at the same time.
* **Less than that:** the jobs queue and run one after another. Slower, but nothing to change.

### Step 3. On Ada: end-to-end test on fake data (~10 min)

This runs the whole pipeline as a real SLURM job on small fake data, to prove the setup works
before the big upload.

```bash
cd ~/cln3_msi/repo
sbatch --output=$HOME/cln3_msi/slurm_logs/test_%j.out bin/hpc/job_test.sbatch
squeue -u nehar          # PD = waiting, R = running; when it's gone from the list, it finished
```

When `squeue -u nehar` no longer lists it:

```bash
tail -n 25 $(ls -t ~/cln3_msi/slurm_logs/test_*.out | head -1)      # want the last line: "19/19 checks passed"
exit                                            # leave Ada, back to the Mac
```

If it doesn't say 19/19, send me that log and stop here.

### Step 4. Mac: get the latest code and upload the data (the long one)

```bash
cd /Users/neharao/Downloads/lipidomics-testbed-main
bash bin/hpc/upload_to_hpc.sh
```

(If I've told you there's a newer version of the code, first replace the folder as described in
"Getting a newer version of the code" below.)

* **First** it prints `Local data: 36 imzML/ibd files (expected 36)` and the free space on Ada.
* **Then** it uploads the code (seconds), then the data with a progress line per file (~40 GB).
* **At the end** it compares Ada with the Mac and prints `OK: all files arrived`.
* **If it stops** (Ctrl+C, Wi-Fi drop, Mac asleep), run the same command again. It continues
  where it stopped.
* **Speed:** on the BCM network expect roughly 15 to 60 minutes; from home over VPN it can take
  hours.

### Step 5. Mac: start the analysis

```bash
bash bin/hpc/start_on_hpc.sh
```

This logs in, quickly re-checks the setup (the data check should now pass), and submits all
the jobs. It prints one job number per line:

```
01 inventory ─┬─ 02 QC slide E
              ├─ lip_neg: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┐
              ├─ lip_pos: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┼─ 10 report
              └─ met_neg: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┘
```

Each job waits for the one before it and only runs if that one succeeded. You can close the
Terminal: the jobs run on Ada on their own. Expect roughly 2 to 4 hours (longer if the three
modes have to queue, see step 2).

### Step 6. Mac: check progress (any time)

```bash
cd /Users/neharao/Downloads/lipidomics-testbed-main
bash bin/hpc/status.sh
```

It lists the jobs still queued or running, the phase summaries written so far, and the logs of
any job that hit an error.

More detail **on Ada** (`ssh nehar@ada.nri.bcm.edu` first):

```bash
squeue -u nehar                                                          # queued / running
sacct -u nehar -S today --format=JobID,JobName%22,State,Elapsed,MaxRSS   # finished + how it went
ls -t ~/cln3_msi/slurm_logs | head                                       # newest logs
tail -n 50 ~/cln3_msi/slurm_logs/<log name>.out                            # read one
scancel -u nehar                                                         # stop everything (only if needed)
```

* `COMPLETED` = fine.
* `FAILED` = look at its log.
* `CANCELLED` or `DependencyNeverSatisfied` = an earlier job in the same chain failed. Look for
  the first `FAILED` one.

### Step 7. Mac: get the results

```bash
cd /Users/neharao/Downloads/lipidomics-testbed-main
bash bin/hpc/fetch_results.sh
open ~/cln3_msi_outputs_from_ada/report/report.html
```

This copies the summaries, figures, tables, logs and the report to
`~/cln3_msi_outputs_from_ada`. The big data tables stay on Ada.

**Look at this first yourself:**
`~/cln3_msi_outputs_from_ada/figures/03_qc_all/lip_neg/sections_all.png`. On every slide the
**Control** animal should be on the **upper** section.

**Then send me:**

* `summaries/01_inventory_*.md`
* `summaries/02_qc_lip_neg.md`
* `summaries/03_qc_all_*.md`
* the `.out` log of any failed job (in `~/cln3_msi_outputs_from_ada/slurm_logs/`)

---

## Later: after I change the code, or to rerun part of it

```bash
# first get the newer code onto the Mac (see "Getting a newer version of the code" below), then:
cd /Users/neharao/Downloads/lipidomics-testbed-main
bash bin/hpc/upload_to_hpc.sh --code      # code only; the data is already on Ada
FROM=05 bash bin/hpc/start_on_hpc.sh      # rerun from phase 05, reusing 01-04
```

Use the phase you need: `FROM=03` reruns almost everything, `FROM=09` only the statistics and
the report.

Edit `config.yaml` **on the Mac**, then upload with `--code`. The upload makes Ada's copy
identical to the Mac's, so edits made directly on Ada get overwritten.

## Getting a newer version of the code (ZIP download, no git needed)

1. On GitHub, open https://github.com/N3haRao/lipidomics-testbed, make sure the branch menu
   says `main`, then Code > Download ZIP.
2. In Finder, delete (or rename) the old `Downloads/lipidomics-testbed-main` folder and unzip
   the new one in its place. It unzips to the same name.
3. Upload the code:

```bash
cd /Users/neharao/Downloads/lipidomics-testbed-main
bash bin/hpc/upload_to_hpc.sh --code
```

Your settings are part of the code (`bin/hpc/settings.sh`, `config.yaml`). If you changed
anything there, make the same change again in the new folder before uploading.

## From scratch (for a new computer or a fresh start)

```bash
# Mac: download the ZIP as above into Downloads and unzip it, then
cd /Users/neharao/Downloads/lipidomics-testbed-main
bash bin/hpc/upload_to_hpc.sh --code

# Ada
ssh nehar@ada.nri.bcm.edu       # first time: type "yes" to the fingerprint question, then your password
cd ~/cln3_msi/repo
bash bin/hpc/hpc_setup.sh
exit
```

Then steps 3 to 7 above.

Optional, so you don't type the Ada password every time (Mac):

```bash
ls ~/.ssh/id_ed25519.pub 2>/dev/null || ssh-keygen -t ed25519    # press Enter 3 times
ssh-copy-id nehar@ada.nri.bcm.edu                                 # Ada password, one last time
```

## Running as another user (for example, a colleague replicating the analysis)

A second person can run the whole pipeline under their own Ada account, reading the data
straight from `/home/nehar/cln3_msi/data` (no second 40 GB copy) and writing everything to
their own folder. The pipeline only ever reads the data folder.

Where things are on Ada:

| what | path |
|---|---|
| input data (the 36 imzML/ibd files) | `/home/nehar/cln3_msi/data/imzML files for biostats/` with subfolders `Lipid Negative/`, `Lipid Positive/`, `Metabolite Negative/` |
| Neha's copy of the code | `/home/nehar/cln3_msi/repo/` |
| Neha's outputs (report, tables, figures) | `/home/nehar/cln3_msi/outputs/` (report: `outputs/report/report.html`) |
| Neha's Python environment | `/home/nehar/cln3_msi/envs/cln3-msi/` |

**1. Neha gives read access (once, on Ada as nehar).** Replace `<user>` with the colleague's
Ada username. If Ada supports access-control lists, this gives access to that one person only:

```bash
setfacl -m u:<user>:x /home/nehar /home/nehar/cln3_msi
setfacl -R -m u:<user>:rX /home/nehar/cln3_msi/data /home/nehar/cln3_msi/outputs
```

If `setfacl` is not available, ask the Ada administrators to put both accounts in a shared
group, then `chgrp -R <group>` and `chmod -R g+rX` those folders (plus `chmod g+x` on
`/home/nehar` and `/home/nehar/cln3_msi`). Avoid `chmod o+r`: that would let every Ada user
read unpublished data. Check from the colleague's account with
`ls "/home/nehar/cln3_msi/data/imzML files for biostats/Lipid Negative" | head`.

**2. The colleague gets the code (on Ada, in their own account).**

```bash
ssh <user>@ada.nri.bcm.edu
mkdir -p ~/cln3_msi && cd ~/cln3_msi
git clone https://github.com/N3haRao/lipidomics-testbed.git repo
cd repo && git checkout claude/hopeful-ride-f4pcg6     # or the commit/branch to reproduce
```

Without git: `curl -L -o code.zip <ZIP link of the branch on GitHub>`, `unzip code.zip`, and
rename the folder to `repo`.

**3. Point the settings at the shared data.** In `bin/hpc/settings.sh` (open with `nano`):

```bash
HPC_USER="<user>"
HPC_WORK=""                                   # = ~/cln3_msi in the colleague's home
HPC_DATA="/home/nehar/cln3_msi/data"          # read-only use of Neha's data
```

**4. Set up and run (all on Ada).**

```bash
cd ~/cln3_msi/repo
bash bin/hpc/hpc_setup.sh       # own conda env + reference downloads, first time 10-20 min; ends "== setup OK"
bash bin/hpc/submit_all.sh      # submits every phase as SLURM jobs
squeue -u <user>                # watch; empty = finished
```

The setup stops with a clear message if the data folder isn't readable yet. Outputs go to
`~/cln3_msi/outputs/` of the colleague, so the two runs never overwrite each other. To check
that the runs match, compare the Reproducibility section of the two reports: the same code
fingerprint and config hash should give the same numbers.

**5. Look at the report** from a Mac or PC:

```bash
scp <user>@ada.nri.bcm.edu:~/cln3_msi/outputs/report/report.html .
```

To only look at Neha's existing results without rerunning anything, step 1 is enough, then
`scp <user>@ada.nri.bcm.edu:/home/nehar/cln3_msi/outputs/report/report.html .`

## If something goes wrong

| what you see | what to do |
|---|---|
| `ssh: ... timed out` | connect to the BCM VPN, try again |
| `Permission denied` on ssh | wrong password or account issue: BCM HPC support |
| `Not found: .../imzML files for biostats` | `LOCAL_DATA_ROOT` in `bin/hpc/settings.sh` is wrong |
| `hpc_setup.sh` stops at `conda create` or `pip install` | send me the last 30 lines |
| the fake-data test isn't 19/19 | send me `~/cln3_msi/slurm_logs/test_*.out` |
| a job ends in `TIMEOUT` | raise its value on the `TIME=` line of `bin/hpc/submit_all.sh`, upload `--code`, rerun with `FROM=` that phase |
| a job ends in `OUT_OF_MEMORY` | change `--mem=32G` to `64G` in `bin/hpc/job_phase.sbatch`, upload `--code`, rerun with `FROM=` |
| anything else | send me the job's `.out` log and the newest file in `~/cln3_msi_outputs_from_ada/logs/` |

## All remaining commands in one block

```bash
# ---- on Ada (you're here now) ----
ls ~/cln3_msi/cln3_env.sh ~/cln3_msi/outputs/cache
nproc; free -g | head -2
cd ~/cln3_msi/repo
sbatch --output=$HOME/cln3_msi/slurm_logs/test_%j.out bin/hpc/job_test.sbatch
squeue -u nehar
tail -n 25 $(ls -t ~/cln3_msi/slurm_logs/test_*.out | head -1)
exit

# ---- on the Mac ----
cd /Users/neharao/Downloads/lipidomics-testbed-main
bash bin/hpc/upload_to_hpc.sh
bash bin/hpc/start_on_hpc.sh
bash bin/hpc/status.sh
bash bin/hpc/fetch_results.sh
open ~/cln3_msi_outputs_from_ada/report/report.html
```
