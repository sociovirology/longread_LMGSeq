# Setup Instructions

This document walks through everything needed to get my demultiplexing and genotyping
pipeline running from scratch: raw sequencing data on Rosalind, all the way through to
genotyped output on the HPC.


## Step 1: Basecalling raw data with Dorado (on Rosalind)

Raw sequencing data comes off the MinION as either pod5 or fast5 files. The best thing to do is set MinKNOW to generate .pod5 files, but it is possible to start with fast5s instead. If you're starting with pod5s, you can skip straight to the basecalling command below.

### a) Convert fast5 to pod5 (skip this if you already have pod5s)

fast5 data is usually split into `fast5_pass` and `fast5_fail` folders by fast basecalling algorithm run during sequencing. Convert each separately since we will want to re-run a more accurate basecalling algorithm on the complete data:
```
pod5 convert fast5 --output converted_fail_pod5s.pod5 fast5_fail/*/*.fast5
pod5 convert fast5 --output converted_pass_pod5s.pod5 fast5_pass/*/*.fast5
mkdir converted_pod5s
mv converted_*pod5s.pod5 converted_pod5s
cd converted_pod5s
```
If you want, you can check the read count in each converted file matches what you'd expect with `pod5 inspect
reads <file>.pod5 | wc -l`.

### b) Merge into a single pod5 file

```
pod5 merge converted_pass_pod5s.pod5 converted_fail_pod5s --output <run_name>_all.pod5
```
You can again check the merged read count is the sum of the two inputs, if you like. You can also delete the intermediate pod5 files since they're large and no longer needed after merging.

### c) Download the `sup` basecalling model (only needed once per model)

Run this command in the new directory with your pod5 data:
```
~/dorado/bin/dorado download --model sup --data <run_name>_all.pod5
```

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

## Step 2: Software setup

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
porechop` will give you the stock version with the default adapter sequences.

The following installation instructions are based on the porechop github: https://github.com/rrwick/porechop.

a) Important: First, make sure the `ont-lmgseq` conda environment is active before you install it:
```
conda activate ont-lmgseq
```
Porechop's executable gets installed directly into whichever environment is active at the time of installation. For this reason, the `ont-lmgseq` conda environment also has to be active every time you run porechop, even though it isn't a conda package. This installation method worked reliably for me but I will mention another option below that should work***, but I have not tested.

b) Fetch the porechop tool and all files by cloning the porechop repo to your desired location:
```
git clone https://github.com/rrwick/Porechop.git
cd Porechop
```

### c) Add the custom adapter sequences

The `adapters.py` file in the Porechop source defines what sequences porechop will look for and trim. Our version only needs the PCR2 landing pad sequences (forward and reverse-complement), not the default, which is a list of Nanopore kit adapters.

Copy the custom `adapters.py` from this repo over the one in the cloned source:
```
cp /path/to/repo/porechop_setup/adapters.py Porechop/porechop/adapters.py
```

Then continue with the porechop installation. Inside your cloned Porechop/ directory run:
```
python3 setup.py install
```

Confirm it worked:
```
porechop -h
python3 -c "from porechop.adapters import ADAPTERS; print([a.name for a in ADAPTERS])"
```
The second command should print `PCR2_landing_pads` and `PCR2_landing_pads_revcomp`, not a long list of Nanopore kit names. If you ever need to change the adapter sequences again, manually edit `adapters.py` and rerun `python3 setup.py install` from inside the Porechop folder with `ont-lmgseq` still active.

*** Alternative option: theoretically you should be able to use `pip3 install local/path/to/Porechop` AFTER you clone the repo and replace adapters.py, and would then avoid the need to use a conda environment, but I haven't tried it since the environment is needed for the other tools anyway.



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
