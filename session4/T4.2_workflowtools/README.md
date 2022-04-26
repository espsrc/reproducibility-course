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

To execute this example, type the following:


```
snakemake -n
```

> NOTE: You are executing snakemake with the option dry-run mode:

It will produce the next:

```
Building DAG of jobs...
Job stats:
job        count    min threads    max threads
-------  -------  -------------  -------------
bwa_map        1              1              1
total          1              1              1


[Tue Apr 26 15:17:41 2022]
rule bwa_map:
    input: data/genome.fa, data/samples/A.fastq
    output: mapped_reads/A.bam
    jobid: 0
    resources: tmpdir=/tmp

Job stats:
job        count    min threads    max threads
-------  -------  -------------  -------------
bwa_map        1              1              1
total          1              1              1

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

Now we are going to execute the workflow without dry-run option:

```
snakemake
```

> NOTE: Snakemake needs to know what will be the amount of resources it will use for the workflow.

So you will need to run the command by typing:

```
snakemake -cores 1
```

```
snakemake -cores all
```

And this is the execution output:


```
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job stats:
job        count    min threads    max threads
-------  -------  -------------  -------------
bwa_map        1              1              1
total          1              1              1

Select jobs to execute...

[Tue Apr 26 15:24:06 2022]
rule bwa_map:
    input: data/genome.fa, data/samples/A.fastq
    output: mapped_reads/A.bam
    jobid: 0
    resources: tmpdir=/tmp

[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 0.919 CPU sec, 0.919 real sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem data/genome.fa data/samples/A.fastq
[main] Real time: 1.084 sec; CPU: 0.965 sec
[Tue Apr 26 15:24:07 2022]
Finished job 0.
1 of 1 steps (100%) done
Complete log: /home/ubuntu/snakemake/.snakemake/log/2022-04-26T152406.190684.snakemake.log
```

After that the output of the workflow will be in `./mapped_reads/A.bam`. To check the results you can use:

```
samtools view -h mapped_reads/A.bam | less
```

And you will see the next:

```
@SQ     SN:I    LN:230218
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem data/genome.fa data/samples/A.fastq
SRR800764.1     4       *       0       0       *       *       0       0  CATCTTTGGAGTAACTATTATTTCGCCCCTTTTGTTTGCTGCATATCGCCCCGCTCTCTGCATACACGATTGGATAATGACCAAAGCAAGGTTTAATACGC   ;:8:DDDDH>F+AEHHGIGECHJGIEHDGHGIHIGIIJJDHGJJGGGEBEFFHG>H?HEFBDBCCCCBDDDDDCD>CACC@CCD@??CC??((49CA:@<B   AS:i:0  XS:i:0
SRR800764.2     4       *       0       0       *       *       0       0       AGAATTTGACGTATATACTATCCAAAGTGGATAAATATTTGAGTTCCAATGAGGTAACTAAACTAGTTTTTTAAAACCCAGCAGAAAGAAAAAATACACTA   ?@@;DDF>FF?AABDEF9AEFHHB:CD<<CC<E<9E9EGIA<CDGDG>BD<?<FDBFHI4B>FFFG)=CHHECA;CE5;B;DE<@;@;;=;?B=9AC(:@>   AS:i:0  XS:i:0
SRR800764.3     4       *       0       0       *       *       0       0       ATTACTCATCAAATCATGTACGAACTTCGATGGCAGTGAACCAGTTGGAGGGCCATTGCTGGTTTCAGGGCCGGAAAAAGTTGAAGGGGGGCGTCCGGTTA   @@CFFFFDHHHHHIIIJJIAIHIJJJGJGFAHIIIJJJFEAH8BF@GGEHGGGGCHJJCHIIIBEE>?CHFFDD;=B?98:>C>CDC@05<&959@9>059   AS:i:0  XS:i:0
SRR800764.4     4       *       0       0       *       *       0       0       TATAATATTTTTTATTTTTAATTTATAATATATAATAATAAATTATAAATAAATTTTAATTAAAAGTAGTATTAACATATTATAAATAGACAAAAGAGTCT   <:+A:DADDAABH;CFEGH@BC,ABH>@HHCBBF9BFHI>CD@4CGCGGGH<4DBFG49BF??DGC3==FGGGA4CGCE4CG>DGCD>??EF>@BDB####   AS:i:0  XS:i:0

```

Now is time to show more deatils on this workflow with a report of this execution:

```
snakemake --report report.html
```

### Writing Workflows

In Snakemake, workflows are specified as Snakefiles. Snakefile contains rules that denote how to create output files from input files. Dependencies between rules are handled implicitly, by matching filenames of input files against output files. Thereby wildcards can be used to write general rules.


```
snakemake    = statement | rule | include | workdir | module | configfile | container
rule         = "rule" (identifier | "") ":" ruleparams
include      = "include:" stringliteral
workdir      = "workdir:" stringliteral
module       = "module" identifier ":" moduleparams
configfile   = "configfile" ":" stringliteral
userule      = "use" "rule" (identifier | "*") "from" identifier ["as" identifier] ["with" ":" norunparams]
ni           = NEWLINE INDENT
norunparams  = [ni input] [ni output] [ni params] [ni message] [ni threads] [ni resources] [ni log] [ni conda] [ni container] [ni benchmark] [ni cache]
ruleparams   = norunparams [ni (run | shell | script | notebook)] NEWLINE snakemake
input        = "input" ":" parameter_list
output       = "output" ":" parameter_list
params       = "params" ":" parameter_list
log          = "log" ":" parameter_list
benchmark    = "benchmark" ":" statement
cache        = "cache" ":" bool
message      = "message" ":" stringliteral
threads      = "threads" ":" integer
resources    = "resources" ":" parameter_list
version      = "version" ":" statement
conda        = "conda" ":" stringliteral
container    = "container" ":" stringliteral
run          = "run" ":" ni statement
shell        = "shell" ":" stringliteral
script       = "script" ":" stringliteral
notebook     = "notebook" ":" stringliteral
moduleparams = [ni snakefile] [ni metawrapper] [ni config] [ni skipval]
snakefile    = "snakefile" ":" stringliteral
metawrapper  = "meta_wrapper" ":" stringliteral
config       = "config" ":" stringliteral
skipval      = "skip_validation" ":" stringliteral
```




