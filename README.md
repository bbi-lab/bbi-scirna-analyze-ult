# bbi-scirna-analyze-ult

## Intro

This *bbi-scirna-analyze-ult* pipeline reads the unaligned BAM files and hash read *.tsv* files from the *bbi-scirna-demux-ult* pipeline, and processes them
to make the required sample-specific output files.

## Installation

Install the following software

- Nextflow

You may also need to build and install

- process_hashes
- STAR aligner

The following sections have additional information.

### Install Nextflow

See the Nextflow installation instructions at

https://www.nextflow.io/docs/latest/install.html

### *process_hashes* program

We include a *process_hashes* executable in the *bbi-scirna-analyze-ult/bin* directory. It runs on the Shendure cluster nodes. If it does not run on your CPUs, you will need to build the executable from the source code. See the following sections on installing Rust and building and installing *process_hashes*.

### Install Rust

See the Rust installation instructions at

https://www.rust-lang.org/tools/install

### Build and install *process_hashes*

Run the following commands

```
cd bbi-scirna-analyze-ult/src/process_hashes
cargo build --release
cp target/release/process_hashes ../../bin
```

I recommend that you build *process_hashes* on a newer cluster node, for example, s020 on the Shendure cluster.

### *STAR* aligner program

We include a *STAR* aligner executable in the *bbi-scirna-analyze-ult/bin* directory. It runs on the Shendure cluster nodes. If it does not run on your CPUs, you will need to build the executable from the source code. See the following sections on building and installing the *STAR* aligner program.

### Load the *STAR* aligner on the Genome Sciences cluster

The *STAR* aligner can be loaded by Nextflow from a GS module. The module load command is in the file *bbi-scirna-analyze-ult/nextflow.config*.

I found that the most recent *STAR* aligner version, *STAR 2.7.11b*, can fail with a segmentation fault when run on input files with only a few reads. You can fix this problem by cloning the *STAR* git repository and editing the source-code file called *STAR-2.7.11b/source/serviceFuns.cpp* to add the lines

```
    if (N==0)
        return -1;
```

at line 300 so that the start of the edited function looks like

```
template <class argType>
inline int64 binarySearchExact(argType x, argType *X, uint64 N) {
    //binary search in the sorted list
    //check the boundaries first
    //returns -1 if no match found
    //if X are not all distinct, no guarantee which element is returned
    
    if (N==0)
        return -1;
    
    if (x>X[N-1] || x<X[0])
        return -1;
```

Then build the *STAR* executable using the command

```
make
```

in the *STAR-2.7.11b/source* directory. Copy the resulting *STAR* executable file to the *bbi-scirna-analyze-ult/bin* directory.

Notes:
- the STAR compilation failed when the samtools/1.19 module was loaded. I believe that the linker was using the samtools htslib rather than the STAR htslib.

### Edit the *experiment.config* file.

Use the *experiment.config* file that you prepared for the *bbi-scirna-demux-ult* pipeline. You can add Nextflow configuration parameters such as

- params.run_empty_drops
- params.run_scrublet
- params.cpuid_level

See *bbi-scirna-analyze-ult/main.nf* for additional parameters and their default values. You can include these parameters in *experiment.config* before you run the *bbi-scirna-demux-ult* pipeline.

## Run bbi-scirna-analyze-ult

Use the *run.analyze.sh* bash script to start the pipeline run.

The output files are in the directory *analyze_out*. They are organized by sample name and *process_group*.
