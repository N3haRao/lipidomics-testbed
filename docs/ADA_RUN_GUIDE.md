# Running the CLN3 MSI pipeline on Ada: every command, in order

Everything below is typed in the **Mac Terminal** (Applications > Utilities > Terminal),
unless a step says "on Ada". Copy one block at a time. Lines starting with `#` are comments,
you don't need to type them.

| | |
|---|---|
| Ada login | `nehar@ada.nri.bcm.edu` |
| Code on the Mac | `/Users/neharao/Downloads/lipidomics-testbed` |
| Data on the Mac | `/Users/neharao/Desktop/BRAIN - PILAR - SPATIAL` |
| Work folder on Ada | `~/cln3_msi` by default (step 3 checks it's big enough) |

## What gets uploaded to Ada

Only one folder: **`imzML files for biostats/`** (36 files, ~40 GB). It holds 6 slides
(E, F, G, I, J, K) x 3 modes (Lipid Negative, Lipid Positive, Metabolite Negative), and each
slide comes as an `.imzML` + `.ibd` pair.

Everything else stays on the Mac. The analysis doesn't use it:

- the SCiLS projects (`.slx/.sbd`, `raw lip neg/`, `raw lip pos/`, `raw met neg/`, `Female Data - Pilar analysis/`)
- the region and feature files (`.srd`, `.mca`)
- the SCiLS exports (`export_All_*`)
- Pilar's results (`Results_CLN3_Analysis*/`, `PC 38.1 - depleted.png`, `.Rhistory`, `.RDataTmp`)
- the cerebellum run (`09SEP26/`)

The upload script picks the right folder by itself.

---

## Step 0. Before you start

* **Off campus?** Connect to the BCM VPN first. Ada usually only accepts connections from
  inside the BCM network.
* **Keep the Mac awake.** The data upload in step 8 takes a while. The script prevents sleep,
  but keep the laptop plugged in and don't close the lid.

## Step 1. Get the code into Downloads

```bash
cd /Users/neharao/Downloads
git clone https://github.com/N3haRao/lipidomics-testbed.git
cd lipidomics-testbed
git checkout claude/hopeful-ride-f4pcg6
ls
```

You should see `README.md`, `bin`, `config.yaml`, `scripts`, `src` and other folders.

* The repo is private, so git asks for a login. Username: `N3haRao`. Password: a GitHub
  **personal access token**, not your GitHub password. You can make one at GitHub > Settings >
  Developer settings > Personal access tokens. If you have the `gh` tool, run `gh auth login`
  first instead.
* No git at all? On GitHub, open the repo, switch to the branch
  `claude/hopeful-ride-f4pcg6`, click Code > Download ZIP, unzip it into Downloads, and rename
  the folder to `lipidomics-testbed`.

## Step 2. First login to Ada

```bash
ssh nehar@ada.nri.bcm.edu
```

* The first time, it asks "Are you sure you want to continue connecting (yes/no)?". Type `yes`.
* Then type your password. Nothing appears on screen while you type, that's normal. Press Enter.
* The prompt changes to something like `[nehar@ada ~]$`. You are now **on Ada**.

## Step 3. On Ada: look around (2 minutes)

Paste this whole block. It only reads, it changes nothing:

```bash
hostname
echo "home: $HOME"
df -h ~                                   # free space in your home folder
quota -s 2>/dev/null                      # your quota, if Ada uses quotas
ls -d /scratch/$USER /scratch /storage* /project* /data* 2>/dev/null   # other big storage
sinfo -s                                  # SLURM partitions (queues)
module avail 2>&1 | grep -i -E "conda|mamba|python" | head
which conda 2>/dev/null
curl -sI https://github.com | head -1     # does Ada have internet? (want "HTTP/2 200")
```

Note down:

1. **Space.** Is there ~60 GB free in your home folder (`Avail` column of `df -h ~`, or your
   quota)? If not, find a bigger folder in the `ls -d` line (e.g. `/scratch/nehar`).
2. **Partitions.** The names in the first column of `sinfo -s`. The one marked with `*` is
   the default.
3. **Conda module.** Any conda/miniforge/anaconda module listed?
4. **Internet.** Did `curl` print `HTTP/... 200`?

Then leave Ada:

```bash
exit
```

## Step 4. Back on the Mac: fill in the settings (only if step 3 needs it)

```bash
cd /Users/neharao/Downloads/lipidomics-testbed
open -e bin/hpc/settings.sh        # opens in TextEdit; or: nano bin/hpc/settings.sh
```

Already filled in for you: `HPC_USER="nehar"`, `HPC_HOST="ada.nri.bcm.edu"` and the Mac data
path. Change only what step 3 told you:

| if step 3 showed... | set |
|---|---|
| home folder has < 60 GB free | `HPC_WORK="/scratch/nehar/cln3_msi"` (or whichever big folder you found) |
| you must pick a partition (no `*` default) | `SLURM_PARTITION="name"` |
| a conda module exists | `CONDA_MODULE="miniforge3"` (the exact name from `module avail`) |
| Ada asks for an account/project when submitting jobs | `SLURM_ACCOUNT="..."` |

Save and close. In the rest of this guide, **`~/cln3_msi` means your `HPC_WORK`** if you
changed it.

## Step 5 (optional, recommended). Log in without typing the password every time

```bash
ls ~/.ssh/id_ed25519.pub 2>/dev/null || ssh-keygen -t ed25519      # press Enter 3 times
ssh-copy-id nehar@ada.nri.bcm.edu                                   # your Ada password, once
ssh nehar@ada.nri.bcm.edu "echo key login works"                    # should not ask anymore
```

If `ssh-copy-id` is refused, Ada doesn't allow keys. That's fine, skip this step. The
scripts reuse one connection, so you only type the password once per script run anyway.

## Step 6. Upload the code only (small, seconds)

```bash
cd /Users/neharao/Downloads/lipidomics-testbed
bash bin/hpc/upload_to_hpc.sh --code
```

## Step 7. On Ada: build the software and test it on fake data

```bash
ssh nehar@ada.nri.bcm.edu
```

Then, **on Ada**:

```bash
cd ~/cln3_msi/repo
bash bin/hpc/hpc_setup.sh
```

This takes 5 to 15 minutes the first time. It:

1. finds or installs conda
2. builds the `cln3-msi` environment in `~/cln3_msi/envs`
3. runs the unit tests (should say `20 passed`)
4. downloads LIPID MAPS, HMDB and the Allen atlas

It warns "Data folder check failed" because the data isn't there yet. That's expected at
this point. If it ends with `== setup OK`, continue.

Still on Ada: run the end-to-end test on fake data as a real SLURM job (about 10 minutes).
Add `-p <partition>` after `sbatch` if you set a partition in step 4.

```bash
mkdir -p ~/cln3_msi/slurm_logs
sbatch --job-name=cln3_test --cpus-per-task=2 --mem=8G --time=00:45:00 \
  --output=$HOME/cln3_msi/slurm_logs/test_%j.out \
  --wrap "source $HOME/cln3_msi/cln3_env.sh && cd $HOME/cln3_msi/repo && bash bin/run_synthetic_test.sh"
squeue -u nehar                       # PD = waiting, R = running; empty = finished
```

When it's no longer listed in `squeue`:

```bash
tail -25 ~/cln3_msi/slurm_logs/test_*.out      # want: "19/19 checks passed"
exit                                            # back to the Mac
```

## Step 8. Upload the data (the long one)

On the Mac:

```bash
cd /Users/neharao/Downloads/lipidomics-testbed
bash bin/hpc/upload_to_hpc.sh
```

* It first prints `Local data: 36 imzML/ibd files (expected 36)`.
* Then it shows a progress line per file (~40 GB in total).
* It ends by printing what arrived on Ada. The first number must be **36**.

If the connection drops or you stop it (Ctrl+C), run the same command again. It skips
everything that's already complete and continues where it stopped.

## Step 9. Start the analysis

On the Mac:

```bash
bash bin/hpc/start_on_hpc.sh
```

This logs in, re-checks the setup (fast now, and this time the data check should pass), and
submits the jobs:

```
01 inventory ─┬─ 02 QC slide E
              ├─ lip_neg: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┐
              ├─ lip_pos: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┼─ 10 report
              └─ met_neg: 03 → 04 → 05 → 06 → 07 → 08 → 09 ─┘
```

It prints a job number for each. You can close the Terminal now: the jobs run on Ada on
their own. Expect roughly 2 to 4 hours in total, longer if the queue is busy.

## Step 10. Check on it

From the Mac, any time:

```bash
cd /Users/neharao/Downloads/lipidomics-testbed
bash bin/hpc/status.sh
```

It shows:

* jobs still queued or running
* the phase summaries written so far
* the log files of any job that hit an error

Or, if you prefer to look on Ada (after `ssh nehar@ada.nri.bcm.edu`):

```bash
squeue -u nehar                                                     # what's queued/running
sacct -u nehar -S today --format=JobID,JobName%22,State,Elapsed,MaxRSS   # what finished, how it went
ls -t ~/cln3_msi/slurm_logs | head                                  # newest job logs
tail -50 ~/cln3_msi/slurm_logs/<name of a log>.out                  # read one
scancel -u nehar                                                    # STOP everything (only if needed)
```

A job shown as `CANCELLED` or `DependencyNeverSatisfied` means an earlier job in its chain
failed. Look at the log of the first `FAILED` one.

## Step 11. Get the results onto the Mac

```bash
cd /Users/neharao/Downloads/lipidomics-testbed
bash bin/hpc/fetch_results.sh
open ~/cln3_msi_outputs_from_ada/report/report.html
```

This copies the summaries, figures, tables, logs and the report into
`~/cln3_msi_outputs_from_ada`. It skips the big data tables, which stay on Ada.

**What to send me first:**

* `~/cln3_msi_outputs_from_ada/summaries/01_inventory_*.md`
* `~/cln3_msi_outputs_from_ada/summaries/02_qc_lip_neg.md`
* `~/cln3_msi_outputs_from_ada/summaries/03_qc_all_*.md`
* `~/cln3_msi_outputs_from_ada/figures/03_qc_all/lip_neg/sections_all.png`. Check this one
  yourself first: on every slide, the **Control** animal should be on the **upper** section.
* the `.out` log of any job that failed

## Step 12. Later: after I change the code, or to rerun part of it

```bash
cd /Users/neharao/Downloads/lipidomics-testbed
git pull                                  # get my changes
bash bin/hpc/upload_to_hpc.sh --code      # send them to Ada (the data is already there)
FROM=05 bash bin/hpc/start_on_hpc.sh      # rerun from phase 05 on, reusing 01-04
```

Use the phase number you need: `FROM=03` reruns almost everything, `FROM=09` only the
statistics and the report.

Edit `config.yaml` **on the Mac**, then upload with `--code`. The upload mirrors the Mac copy,
so anything edited directly on Ada gets overwritten.

---

## If something goes wrong

| what you see | what to do |
|---|---|
| `ssh: connect to host ada.nri.bcm.edu ... timed out` | connect to the BCM VPN, try again |
| `Permission denied (publickey,password)` | wrong password, or your Ada account isn't active yet: ask BCM HPC support |
| `Not found: .../imzML files for biostats` | the data path in `bin/hpc/settings.sh` (`LOCAL_DATA_ROOT`) is wrong |
| upload stops with `No space left on device` or `Disk quota exceeded` | home folder too small: set `HPC_WORK` to scratch (step 4), run step 6 onwards again |
| `hpc_setup.sh` fails at `curl` / `conda env create` / downloads | Ada has no internet: tell me, we'll build the env on the Mac side and copy it over |
| `sbatch: error: ... partition` or `... account` | set `SLURM_PARTITION` / `SLURM_ACCOUNT` in step 4, upload `--code`, rerun step 9 |
| a job `TIMEOUT` | raise its limit in `bin/hpc/submit_all.sh` (the `TIME=` line), upload `--code`, rerun with `FROM=` that phase |
| a job `OUT_OF_MEMORY` | change `--mem=32G` to `64G` in `bin/hpc/job_phase.sbatch`, upload `--code`, rerun with `FROM=` |
| anything else | send me the `.out` log and the newest file from `~/cln3_msi_outputs_from_ada/logs/` |

---

## All commands in one place (after steps 3 and 4 are done)

```bash
# --- Mac: code ---
cd /Users/neharao/Downloads
git clone https://github.com/N3haRao/lipidomics-testbed.git
cd lipidomics-testbed
git checkout claude/hopeful-ride-f4pcg6

# --- Mac: optional password-free login ---
ssh-keygen -t ed25519
ssh-copy-id nehar@ada.nri.bcm.edu

# --- Mac: code to Ada ---
bash bin/hpc/upload_to_hpc.sh --code

# --- Ada: setup + fake-data test ---
ssh nehar@ada.nri.bcm.edu
cd ~/cln3_msi/repo
bash bin/hpc/hpc_setup.sh
mkdir -p ~/cln3_msi/slurm_logs
sbatch --job-name=cln3_test --cpus-per-task=2 --mem=8G --time=00:45:00 \
  --output=$HOME/cln3_msi/slurm_logs/test_%j.out \
  --wrap "source $HOME/cln3_msi/cln3_env.sh && cd $HOME/cln3_msi/repo && bash bin/run_synthetic_test.sh"
squeue -u nehar
tail -25 ~/cln3_msi/slurm_logs/test_*.out
exit

# --- Mac: data, start, follow, results ---
cd /Users/neharao/Downloads/lipidomics-testbed
bash bin/hpc/upload_to_hpc.sh
bash bin/hpc/start_on_hpc.sh
bash bin/hpc/status.sh
bash bin/hpc/fetch_results.sh
open ~/cln3_msi_outputs_from_ada/report/report.html
```
