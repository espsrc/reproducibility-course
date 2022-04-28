# Workflow management systems and tools


**Table of contents:**

<!-- vscode-markdown-toc -->
* 1. [Snakemake](#Snakemake)
	* 1.1. [How to install Snamemake](#HowtoinstallSnamemake)
		* 1.1.1. [Step 1: Installing Mambaforge](#Step1:InstallingMambaforge)
		* 1.1.2. [Step 2: Preparing a working directory](#Step2:Preparingaworkingdirectory)
		* 1.1.3. [Step 3: Creating an environment](#Step3:Creatinganenvironment)
		* 1.1.4. [Step 4: Activating the environment](#Step4:Activatingtheenvironment)
	* 1.2. [Basics: An example workflow with snakemake](#Basics:Anexampleworkflowwithsnakemake)
		* 1.2.1. [Step 1 - Composing a simple workflow](#Step1-Composingasimpleworkflow)
	* 1.3. [Writing Workflows](#WritingWorkflows)
	* 1.4. [A full example](#Afullexample)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

<HR>
	
**Aim of this session**
	
The aim of this session is to carry out a series of practical examples with which to practice the development of a workflow, through two examples: a) one related to genomics tools and b) another workflow related to an ETL process: downloading data, processing and obtaining results.


##  1. <a name='Snakemake'></a>Snakemake

Snakemake workflow management system is a tool to create reproducible and scalable data analyses. Workflows are described via a human readable, Python based language. 

They can be scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition.

###  1.1. <a name='HowtoinstallSnamemake'></a>How to install Snamemake

The easiest way is to use the Mambaforge Python 3 distribution (Mambaforge is a Conda based distribution like Miniconda, which however uses Mamba a fast and more robust replacement for the Conda package manager). 

The tutorial assumes that you are using either Linux or MacOS X.

####  1.1.1. <a name='Step1:InstallingMambaforge'></a>Step 1: Installing Mambaforge

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

####  1.1.2. <a name='Step2:Preparingaworkingdirectory'></a>Step 2: Preparing a working directory

First, create a new directory snakemake-first-steps at a place you can easily remember and change into that directory in your terminal:

```
$ mkdir snakemake-first-steps
$ cd snakemake-first-steps
```

 we download some example data on which the workflow shall be executed:

```
$ curl -L https://github.com/snakemake/snakemake-tutorial-data/archive/v5.24.1.tar.gz -o snakemake-tutorial-data.tar.gz
```

Next we extract the data. On Linux, run

```
$ tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"
```



####  1.1.3. <a name='Step3:Creatinganenvironment'></a>Step 3: Creating an environment

First, make sure to activate the conda base environment with

```
$ conda activate base
```

####  1.1.4. <a name='Step4:Activatingtheenvironment'></a>Step 4: Activating the environment


The `environment.yaml` file that you have obtained with the previous step can be used to install all required software into an isolated Conda environment with the name snakemake-tutorial via

```
$ mamba env create --name snakemake-tutorial --file environment.yaml
```

To activate the `snakemake-tutorial` environment, execute

```
$ conda activate snakemake-tutorial
````



###  1.2. <a name='Basics:Anexampleworkflowwithsnakemake'></a>Basics: An example workflow with snakemake

**Features:** 

- A Snakemake workflow is defined by specifying rules in a Snakefile. 
- Rules decompose the workflow into small steps (for example, the application of a single tool) by specifying how to create sets of output files from sets of input files. 
- Snakemake automatically determines the dependencies between the rules by matching file names.

The Snakemake language extends the Python language, adding syntactic structures for rule definition and additional controls. 

**Basic syntax:**

- All added syntactic structures begin with a `keyword` followed by a `code block` that is either in the same line or indented and consisting of multiple lines. 


####  1.2.1. <a name='Step1-Composingasimpleworkflow'></a>Step 1 - Composing a simple workflow

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

###  1.3. <a name='WritingWorkflows'></a>Writing Workflows

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

###  1.4. <a name='Afullexample'></a>A full example

Now we are going to create a complete example involving different operations, such as:

- Downloading various files from the web
- Performing operations on the files
- Apply a python function to plot the processed data on a graph.

The workflow we want to achieve will be the following:

- Retrieve the data from source A
- Retrieve the data from source B
- Merge the two data sources, adding to source A, source B as a new column.
- Sort the resulting file
- Draw a graph with the previous data.

You will need this `Snakefile`:


```
rule create_bio_plot:
        input:
                "output/sorted/full-bio-indicators.csv"
        output:
                "bio-habits.png"
        shell:
                "python create-bio-plot.py {input}"

rule sort_results:
        input:
                "output/merged/full-bio-indicators.csv"
        output:
                "output/sorted/full-bio-indicators.csv"
        shell:
                "cat {input} | sort > {output}"   
            
rule merge_bio_habits:
        input:
                a="output/output-biostats.csv",
                b="output/output-habits.csv"
        output:
                "output/merged/full-bio-indicators.csv"
        shell:
                "paste -d ' ' {input.a} <(awk '{{print $NF}}' {input.b}) > {output}"

rule get_biostats:
        output:
                "output/output-biostats.csv"
        params:
                biostats = "https://raw.githubusercontent.com/manuparra/reproducibility-course/main/session4/T4.2_workflowtools/data/output-biostats.csv"
        shell:
                "wget -O {output} {params.biostats}"

rule get_habits:
        output:
                "output/output-habits.csv"
        params:
                habits="https://raw.githubusercontent.com/manuparra/reproducibility-course/main/session4/T4.2_workflowtools/data/output-habits.csv"
        shell:
                "wget -O {output} {params.habits}"

```

And the external files that are called from the workflow:


```
import matplotlib.pyplot as plt
import csv
  
x = []
y = []
  
with open('output/sorted/full-bio-indicators.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ';')
      
    for row in plots:
        x.append(row[0])
        y.append(int(row[2]))
  
plt.bar(x, y, color = 'g', width = 0.72, label = "Age")
plt.xlabel('Names')
plt.ylabel('Ages')
plt.title('Ages of different persons')
plt.legend()
plt.savefig('bio-habits.png')
```

In order to execute this workflow type the following:

```
snakemake -c1
```



