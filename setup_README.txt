# Setup Instructions

This document walks through everything needed to get my demultiplexing and genotyping
pipeline running from scratch: raw sequencing data on Rosalind, all the way through to
genotyped output on the HPC.

## Step 1: Software setup

### a) Create the conda environment

I've exported my environment to `ont-lmgseq.yml` using `conda env export --from-history`,
which only lists the packages I explicitly installed.

To recreate it run the following in your terminal:
```
conda env create -f ont-lmgseq.yml
conda activate ont-lmgseq
```
The environment has the following tools: cutadapt, seqkit, seqtk, blast, biopython, and porechop** (extra setup required for porechop - see below**).
usearch is separate - it's installed as an HPC module, not through conda. Load it with:
```
module load usearch
```

Both `demultiplexing.sh` and `genotyping.sh` check for their required tools on PATH before
doing anything else, so if you forget one of these steps you'll get an error message.

### b) Install porechop

Porechop needs to be installed from source, because it needs a custom adapter file (see
part c) swapped in before installing - a plain `pip install porechop` or `conda install
porechop` will give you the stock version looking for the wrong adapter sequences entirely.

The following installation instructions come from the porechop github: https://github.com/rrwick/porechop

Fetch the porechop tool and all files by cloning the porechop repo to your desired location:

First, make sure the `cutadapt` conda environment is active before you install it:
```
conda activate cutadapt
```
Porechop's executable gets installed directly into whichever environment is active at
install time. That's the reason `cutadapt` has to be active every time you run porechop
later, even though it isn't a conda package - it's a pip-style install that happened to
land inside that environment.

```
git clone https://github.com/rrwick/Porechop.git
cd Porechop
```

### c) Add the custom adapter sequences

The `adapters.py` file in the Porechop source defines what sequences porechop looks for
and trims. Our version only needs the PCR2 landing pad sequences (forward and
reverse-complement), not the default list of Nanopore kit adapters.

Copy the custom `adapters.py` from this repo over the one in the cloned source:
```
cp /path/to/repo/porechop_setup/adapters.py Porechop/porechop/adapters.py
```

Then install:
```
cd Porechop
python3 setup.py install
```

Confirm it worked:
```
porechop -h
python3 -c "from porechop.adapters import ADAPTERS; print([a.name for a in ADAPTERS])"
```
The second command should print `PCR2_landing_pads` and `PCR2_landing_pads_revcomp`, not a
long list of Nanopore kit names. If you ever need to change the adapter sequences again,
edit `adapters.py` and rerun `python3 setup.py install` from inside the Porechop folder with
`cutadapt` still active.

## Step 2: Basecalling raw data with Dorado (on Rosalind)

Raw sequencing data comes off the MinION as either pod5 or fast5 files. Most of the time
it'll already be pod5 - in that case you can skip straight to the basecalling command below.
If it's fast5, convert it to pod5 first.

### a) Convert fast5 to pod5 (skip this if you already have pod5s)

Data is usually split into `fast5_pass` and `fast5_fail` folders. Convert each separately:
```
pod5 convert fast5 --output converted_pod5s/ fast5_fail/*/*.fast5
pod5 convert fast5 --output converted_pass_pod5s.pod5 fast5_pass/*/*.fast5
```
Check the read count in each converted file matches what you'd expect with `pod5 inspect
reads <file>.pod5 | wc -l`.

### b) Merge into a single pod5 file

```
pod5 merge converted_pass_pod5s.pod5 output.pod5 --output <run_name>_all.pod5
```
Check the merged read count is the sum of the two inputs, then delete the intermediate pod5
files once you've confirmed that - they're large and no longer needed once merged.

### c) Download the basecalling model (only needed once per model)

```
~/dorado/bin/dorado download --model sup --data <run_name>_all.pod5
```
`sup` is the super-accuracy model (slowest, most accurate). `fast` and `hac` are faster,
lower-accuracy alternatives if you ever need them. This command will also download a few
related modification-calling models alongside the base model - that's normal, we don't use
those for this pipeline, they just come with it.

### d) Run basecalling

```
~/dorado/bin/dorado basecaller -v -r --emit-fastq --min-qscore 8 --disable-read-splitting \
  --no-trim --output-dir <output_dir> sup <pod5_folder>
```
What each flag does:
- `-v` - verbose logging
- `-r` - recursive, so it finds pod5 files in subfolders too
- `--emit-fastq` - output fastq instead of the default BAM. This is required - my pipeline
  needs a fastq to start, not a BAM.
- `--min-qscore 8` - only keep reads with an average quality score of 8 or higher (this is
  the standard ONT pass/fail quality cutoff)
- `--disable-read-splitting` - don't let dorado auto-split reads on suspected internal
  adapters. Read splitting is handled downstream by porechop instead.
- `--no-trim` - don't trim adapters during basecalling. Adapter/barcode trimming is handled
  downstream by porechop and cutadapt, so trimming here would interfere with that.
- `sup` - the model name (must match what you downloaded in step c)
- last argument - the folder containing your pod5 files

### e) Combine into one fastq

Dorado may write out more than one fastq file depending on version/settings. My pipeline
needs a single combined fastq of all quality-passed reads to start from:
```
cat <output_dir>/*.fastq > <run_name>_calls.fastq
```
That combined fastq is what you'll transfer to the HPC and use as the `-r` input to
`demultiplexing.sh`.

## Step 3: Understanding the scripts

- `demultiplexing.sh` - runs porechop then cutadapt to demultiplex reads by plate and well
- `genotyping.sh` - runs usearch (and optionally blastn) to assign a strain/segment identity
  to each demultiplexed well
- `sbatch_demultiplexing.sh` / `sbatch_genotyping.sh` - sbatch wrappers for running each step
  as a cluster job instead of interactively
- `run_pipeline_sbatch.sh` - runs both steps back to back in a single sbatch job

Each script prints its flags near the top if you're unsure what to pass in, and will prompt
you interactively for anything you leave out - except the sbatch scripts, which need every
flag passed explicitly up front, since there's no one there to answer a prompt once a job is
queued.

## Step 4: Get your data and supporting files

You'll need all of the following before running the pipeline:

### a) Basecalled reads
A single fastq file containing all quality-passed reads (see Step 2e).

### b) Barcode files
Two fasta files: one with plate barcode sequences, one with well barcode sequences. Each
fasta header is the barcode's name (e.g. `>plate01`, `>well01`) and the sequence below it is
the actual barcode sequence.

### c) Sample list CSV
A CSV listing every plate/well combination and what it is: sample type (positive, negative,
coinfection, or misassigned) and the parent strain(s) involved.

### d) Library prep reference database
A single fasta file containing every reference genome used anywhere in this library prep.
One reference per strain per genome segment, named `>Strain_Segment` - for example
`>PAN99_M`, `>CA09_NP`.

Once you have all four of these, you're ready to run `demultiplexing.sh` followed by
`genotyping.sh`.
