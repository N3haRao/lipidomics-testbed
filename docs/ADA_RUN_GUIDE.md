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
| Code on the Mac | `/Users/neharao/Downloads/lipidomics-testbed` |
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
mkdir -p ~/cln3_msi/slurm_logs
sbatch --job-name=cln3_test --cpus-per-task=2 --mem=8G --time=00:45:00 \
  --output=$HOME/cln3_msi/slurm_logs/test_%j.out \
  --wrap "source $HOME/cln3_msi/cln3_env.sh && cd $HOME/cln3_msi/repo && bash bin/run_synthetic_test.sh"
squeue -u nehar          # PD = waiting, R = running; when it's gone from the list, it finished
```

When `squeue -u nehar` no longer lists it:

```bash
tail -n 25 ~/cln3_msi/slurm_logs/test_*.out      # want the last line: "19/19 checks passed"
exit                                            # leave Ada, back to the Mac
```

If it doesn't say 19/19, send me that log and stop here.

### Step 4. Mac: get the latest code and upload the data (the long one)

```bash
cd /Users/neharao/Downloads/lipidomics-testbed
git pull
bash bin/hpc/upload_to_hpc.sh
```

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
cd /Users/neharao/Downloads/lipidomics-testbed
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
cd /Users/neharao/Downloads/lipidomics-testbed
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
cd /Users/neharao/Downloads/lipidomics-testbed
git pull                                  # my changes
bash bin/hpc/upload_to_hpc.sh --code      # code only; the data is already on Ada
FROM=05 bash bin/hpc/start_on_hpc.sh      # rerun from phase 05, reusing 01-04
```

Use the phase you need: `FROM=03` reruns almost everything, `FROM=09` only the statistics and
the report.

Edit `config.yaml` **on the Mac**, then upload with `--code`. The upload makes Ada's copy
identical to the Mac's, so edits made directly on Ada get overwritten.

## From scratch (for a new computer or a fresh start)

```bash
# Mac
cd /Users/neharao/Downloads
git clone https://github.com/N3haRao/lipidomics-testbed.git     # login: N3haRao + a GitHub personal access token
cd lipidomics-testbed
git checkout claude/hopeful-ride-f4pcg6
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
mkdir -p ~/cln3_msi/slurm_logs
sbatch --job-name=cln3_test --cpus-per-task=2 --mem=8G --time=00:45:00 \
  --output=$HOME/cln3_msi/slurm_logs/test_%j.out \
  --wrap "source $HOME/cln3_msi/cln3_env.sh && cd $HOME/cln3_msi/repo && bash bin/run_synthetic_test.sh"
squeue -u nehar
tail -n 25 ~/cln3_msi/slurm_logs/test_*.out
exit

# ---- on the Mac ----
cd /Users/neharao/Downloads/lipidomics-testbed
git pull
bash bin/hpc/upload_to_hpc.sh
bash bin/hpc/start_on_hpc.sh
bash bin/hpc/status.sh
bash bin/hpc/fetch_results.sh
open ~/cln3_msi_outputs_from_ada/report/report.html
```
