# Workflow management systems and tools





## Snamemake

Snakemake workflow management system is a tool to create reproducible and scalable data analyses. Workflows are described via a human readable, Python based language. 

They can be scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition.

### How to install Snamemake

The easiest way is to use the Mambaforge Python 3 distribution (Mambaforge is a Conda based distribution like Miniconda, which however uses Mamba a fast and more robust replacement for the Conda package manager). 

The tutorial assumes that you are using either Linux or MacOS X.

#### Step 1: Installing Mambaforge

First, please open a terminal or make sure you are logged into your Vagrant Linux VM. Assuming that you have a 64-bit system, on Linux, download and install Miniconda 3 with

```
$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
$ bash Mambaforge-Linux-x86_64.sh
```

On MacOS with x86_64 architecture, download and install with

```
$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
$ bash Mambaforge-MacOSX-x86_64.sh
```

On MacOS with ARM/M1 architecture, download and install with

```
$ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -o Mambaforge-MacOSX-arm64.sh
$ bash Mambaforge-MacOSX-arm64.sh
```

When you are asked the question

Do you wish the installer to prepend the install location to PATH ...? [yes|no]

answer with yes. 

#### Step 2: Preparing a working directory

First, create a new directory snakemake-first-steps at a place you can easily remember and change into that directory in your terminal:

```
$ mkdir snakemake-first-steps
$ cd snakemake-first-steps
```

#### Step 3: Creating an environment

First, make sure to activate the conda base environment with

```
$ conda activate base
```

#### Step 4: Activating the environment

To activate the snakemake-first-steps environment, execute

```
$ conda activate snakemake-first-steps
```

### Basics: An example workflow with snakemake

**Features:** 

- A Snakemake workflow is defined by specifying rules in a Snakefile. 
- Rules decompose the workflow into small steps (for example, the application of a single tool) by specifying how to create sets of output files from sets of input files. 
- Snakemake automatically determines the dependencies between the rules by matching file names.

The Snakemake language extends the Python language, adding syntactic structures for rule definition and additional controls. 

**Basic syntax:**

- All added syntactic structures begin with a `keyword` followed by a `code block` that is either in the same line or indented and consisting of multiple lines. 


#### Step 1 - Composing a simple workflow

Snakemake rule maps reads of a given data-sample or data-source. For this example, we will use the tool `bwa`, specifically the subcommand `bwa mem`. In the working directory, create a new file called `Snakefile` with an text editor.

In the Snakefile, define the following rule:

```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

A Snakemake rule has a name (here `bwa_map`) and a number of directives, here `input`, `output` and `shell`. 

The `input` and `output` directives are followed by lists of files that are expected to be used or created by the rule. 

The `shell` directive is followed by a Python string containing the `shell` command to execute. In the `shell` command string, we can refer to elements of the rule via **braces** notation (similar to the Python format function). Here, we refer to the `output` file by specifying `{output}` and to the `input` files by specifying `{input}`. 

Since the *rule* has multiple `input` files, Snakemake will concatenate them, separated by a whitespace. In other words, Snakemake will replace `{input}` with `data/genome.fa` `data/samples/A.fastq` before executing the command. 

The shell command invokes `bwa mem` with reference genome and reads, and pipes the output into samtools which creates a compressed BAM file containing the alignments. The output of samtools is redirected into the output file defined by the rule with `>`.






