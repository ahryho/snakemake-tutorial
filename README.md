# Snakemake 


- [About Snakemake](#About-Snakemake)
- [Setup](#setup)
- [Exercises](#Exercises)
   - [Getting started](#Exercise-ex00-Getting-started)
   - [On-the-fly coding](#Excercise-ex01-On-the-fly-coding)
   - [Parallel processing with wildcards](#Excercise-ex02-Parallel-processing-with-wildcards)
   - [Parametrizing rules](#Excercise-ex03-Parametrizing-rules)
   - [Parametrizing rules. Config files](#Excercise-ex04-Parametrizing-rules-Config-files)
   - [Reproducible environment](#Exercise-ex05-Reproducible-environment)
   - [Custom scripts](#Exercise-ex05-Custom-scripts)
   - [Parallelizing workflow. Locally](#Exercises-ex06-Parallelizing-workflow-Locally)
   - [Parallelizing workflow. On cluster](#Exercise-ex07-Parallelizing-workflow-On-cluster)
   - [Debugging the Snakefile](#Exercise-ex08-Debugging-the-Snakefile)
- [Useful command line arguments](#Useful-command-line-arguments)
- [Useful keywords and functions](#Useful-keywords-and-functions)
- [Useful resources](#Useful-resources)

# About Snakemake

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a Pythonic workflow engine developed by [Dr. Johannes Köster](https://johanneskoester.bitbucket.io/), also the founder of the Bioconda project for sustainably distributing bioinformatics software.

- Software to create reproducible and scalable workflows. The workflow is defined by specifying rules how to create sets of output files from sets of input files.

- Easy to prototype workflows due to embedded Python scripting

- Feature-rich and configurable: many different command-line options and multiple ways to configure rules

- __Snakemake__ doesn't (re-)run a job if a taget (output) file already exists and/or an output file is older than an input file. 

# Setup

Clone the repository: 

```sh
git clone \
    --depth 3 \
    --filter=blob:none \
    --sparse \ https://github.molgen.mpg.de/mpip/IT-Wiki/
;
cd IT-Wiki
git sparse-checkout init --cone
git sparse-checkout set 04\ MPIP\ cluster/snakemake-tutorial/
```

<!-- To setup a minimal snakemake environment, first install Miniconda with Bioconda channels by following the instructions [here](https://bioconda.github.io/user/install.html).  -->

Then, use either __Snakemake__ installed on the cluster, or, alternatively, install a __Snakemake__ environment. To install your own  environment, run:

```sh
module load anaconda/anaconda3
conda create --name py38 python=3.8 --channel conda-forge
conda install -n py38 -c conda-forge mamba
conda activate py38
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate /home/USERNAME/.conda/envs/py38/envs/snakemake
```

<!-- ## Installation via pip

Instead of conda, snakemake can be installed with pip. However, note that snakemake has non-python dependencies, such that the pip based installation has a limited functionality if those dependencies are not manually installed in addition. -->

# Exercises

Minimal examples corresponding to this README content can be found in the exercises folder. To run them, use: 

<!--- conda activate smktutorial --->
```console
user@slurmgate$ snakemake --core 1 --printshellcmds -- all
```

### __Exercise ex00:__ _Getting started_

```sh
rule hellomake:
    input: ["hellomake.c", "hellofunc.c"]
    output: "hellomake"
    shell:
        "gcc -o {output} {input} -I."
```

From __input__ files create some __output__ using the __shell__ command. Inside the shell command, all local and global variables, including input and output files can be accessed via their names in the [python format minilanguage](https://docs.python.org/py3k/library/string.html#formatspec). Instead of a __shell__ command, a rule can run some python code to generate the output. For this, use the ```run``` instead of ```shell```:

```sh
rule hellomake:
    input: ["hellomake.c", "hellofunc.c"]
    output: "hellomake"
    run:
        dummy python code to create output[0] from input[0] and input[1]
```

__{input}__ and __{output}__ are wildcards for the input and output files, and are parsed with Python string formatting rules

_Declarative language:_ Tell __Snakemake__ what to do, and it will make it for you:

```sh
rule all:
    input: "hellomake"    
```

__rule all__: the rule that contains all the target files that should be created by the current workflow. This rule doesn't contain any output file, any script or shell command. It allows to avoid mentioning the target in the command line while executing.

___
### __Excercise ex01:__ _On-the-fly coding_

```sh
import os

data = "data/{id}.txt"
ids, = glob_wildcards(data)

rule concatenate:
    input: expand(data, id = ids)
    output: "result/conc.txt"
    shell:
	"cat {input} > {output}"

rule all:
    input: rules.concatenate.output
```

It is possible to use arbitary Python outside the rules: ``` import os```

__glob_wildcards__: a special Snakemake function to look for all files matching a certain pattern. In this example, it returns __["f1", "f2"]__

__expand__: a special Snakemake function to create a list of strings from a pattern. In this example, it returns __["data/f1.txt", "data/f2.txt"]__

__rules.concatenate.output__: It is possible to access the output from other rules programmatically:
- __rules__: is a special variable that contains all the infromation that is present in the other rules of the workflow
- __concatenate__: is a name of the rule in the workflow
- __output__: is an output of the _concatenate_ rule

### __Example execution of a Snakefile__

To execute, run the following command:

```snakemake --core 1 --printshellcmds -- all```

__where__

__-- all__: the name of the rule to execute. It is a good practice to specify the target rule at the end, after two dashes to clarify that it is not an argument

__--printshellcmds__: see the executed commands

Automatic logging is placed in the .snakemake directory by default

This is what you see once run the command:

```js
[user@slurmgate 01]$ snakemake --core 1 --printshellcmds -- all
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	1	concatenate
	2

[Wed Apr 28 12:33:18 2021]
rule concatenate:
    input: data/f1.txt, data/f2.txt
    output: result/conc.txt
    jobid: 1

cat data/f1.txt data/f2.txt > result/conc.txt
[Wed Apr 28 12:33:18 2021]
Finished job 1.
1 of 2 steps (50%) done

[Wed Apr 28 12:33:18 2021]
localrule all:
    input: result/conc.txt
    jobid: 0

[Wed Apr 28 12:33:19 2021]
Finished job 0.
2 of 2 steps (100%) done
Complete log: /path_to/IT-Wiki/04 MPIP cluster/snakemake-tutorial/exercises/01/.snakemake/log/2021-04-28T123318.778561.snakemake.log
```
___

### __Excercise ex02:__ _Parallel processing with wildcards_

```py
data = "data/{id}.txt"
ids, = glob_wildcards(data)

rule sort:
    input: "data/{id}.txt"
    output: "result/{id}_sorted.txt"
    shell:
        "sort {input} > {output}"

rule all:
    input: expand(rules.sort.output, id=ids)
```


__{id}__: the wildcard that matches any partial filename that is a string.

The wildcards must be consistent between input and output files.

It is also possible to restrict the wildcards to a specific format by using regular expressions.
In this example, to constrain the __id__ to alphanumeric characters, use __{id,[A-Za-z0-9]+}.txt__:

```py
data = "data/{id,[A-Za-z0-9]+}.txt"
ids, = glob_wildcards(data)
```
___
### __Excercise ex03:__ _Parametrizing rules_
If you want to run a generic rule differently based on the type of sample being used, or, run the same rule several times, using different sets of parameters, it is possible :

1. Specifying input files using __input functions__
2. Specifying parameters using the __params__ directive
____
__Goal 1:__ Map fastq files to index files

__Goal 2:__ Subsample reads from both samples

```py
data = "data/{id}.fq"
ids, = glob_wildcards(data)
subsampling = [12, 20]
sample_dict = {"chr21.test1": "index/index1/",
               "chr21.test2": "index/index2/"}

rule all:
    input: 
        expand("result/{id}_subsample_{subsample}_Aligned.out.sam", id=ids, subsample=subsampling)

rule get_subsample:
    input:
        data
    output:
        temp("{id}_tmp_{subsample}.fq")
    params:
        "{subsample}"
    shell:
        "head -n {params} {input} > {output}"

rule mapping:
    input:
        fastq="{id}_tmp_{subsample}.fq",
        index=lambda wc: sample_dict[wc.id]
    output: 
        "result/{id}_subsample_{subsample}_Aligned.out.sam"
    params:
        subsample="{subsample}"
    shell:
        "STAR "
        "--genomeDir {input.index} "
        "--readFilesIn {input.fastq} "
        "--outFileNamePrefix result/{wildcards.id}_subsample_{params.subsample}_"
```

To execute the example, run:

```snakemake --core 1 -j 8 --printshellcmds -- all```

__Goal 1:__

Use __input functions__ to determine the input of a rule based on the matching __wildcard of the output__

__Input functions__ take a single variable: the __wildcards__ object

Create a dictionary telling which mapping index to use for which sample:
```py
sample_dict = {"chr21.test1": "index/index1/",
               "chr21.test2": "index/index2/"}
```
Look up the ID using the lambda function:
```py
index=lambda wc: sample_dict[wc.id]
```

__Goal 2:__

Use the __params__ keyword to add additional non-file parameters to the rule:
```py
params: subsample="{subsample}"
```
Add an additional subsample parameter and run each mapping twice using this parameter:
```py
output: 
    "result/{id}_subsample_{subsample}_Aligned.out.sam"
params: 
    subsample="{subsample}"
...

rule all:
    input: 
        expand("result/{id}_subsample_{subsample}_Aligned.out.sam", id=ids, subsample=subsampling)
```

The __params__ object gets passed to the shell directive and can be accessed in a way: __{params.subsample}__

___
__Lexical scoping and string formatting__

_Snakemake_ creates the following 'special' variables as objects:
__input, output, params, wildcards, log, threads, resources, config__.

These are then passed to the rule, shell directive in this example, and expanded using the __.format(**vars)__
string formatting function.

You can refer any variable within the scope of the rule, including external variables defined. 

___
__Workflow graph__

To create the _Snakemake_ workflow graph, execute:

```snakemake --dag -- all | dot -Tpdf > dag.pdf```

The results will be written out to the ```dag.pdf``` file.
___

### __Excercise ex04:__ _Parametrizing rules. Config files_

Alternatively, metadata can be stored in a __config__ file:

```yaml
sample_dict:
 chr21.test1: data/index
 chr21.test2: data/index
subsample:
 [12,20]
```

Then, the _Snakemake_ file looks as follow:

```py
data = "data/{id}.fq"
ids, = glob_wildcards(data)

configfile: "config.yaml"

rule all:
    input: 
        expand("result/{id}_subsample_{subsample}_Aligned.out.sam", id=ids, subsample=config["subsample"])

...

rule mapping:
    input:
        fastq="{id}_tmp_{subsample}.fq",
        index=lambda wc: config["sample_dict"][wc.id]
...      
```
___
### __Exercise ex05:__ _Reproducible environment_

It is possible to define isolated software environments per rule. 

Upon execution of a workflow, the Conda package manager is used to obtain and deploy the defined software packages in the specified versions. 

Packages will be installed into the working directory.

The environment specifications are stored in the __environmnet yaml file__:

```yaml
channels:
 - r
dependencies:
 - r=4.0.3
 - r-ggplot2=3.3.3
```

To install the required conda environments without running the full workflow, run:

```snakemake --use-conda --conda-create-envs-only```

To run the full workflow:

```snakemake --cores 1 --use-conda -- all```

The result of execution:

```js
[user@slurmgate 05]$ snakemake --core 1 -j 2 --use-conda --printshellcmds -- all
Building DAG of jobs...
Creating conda environment envs/plot.yaml...
Downloading and installing remote packages.
Environment for envs/plot.yaml created (location: .snakemake/conda/c00adc15fcb7ce78cd2de9774da7df74)
...
```
___
### __Exercise ex05:__ _Custom scripts_

A rule can also point to an external script instead of a shell command.

Inside the script, you have access to an object __snakemake__ that provides access to the same objects that are available in the run and shell directives: __input, output, params, wildcards, log, threads, resources, config__.

For example, in R, you can use __snakemake@input[[1]]__, or __snakemake@input[["filename"]]__, to get the first file. 

In Python, you would use __snakemake.input[0]__.

_Snakefile_:

```py
data = "data/{file}.csv"
ids, = glob_wildcards(data)

rule all:
    input: 
        expand("plots/plot_{file}.pdf", file = ids)

rule plot:
    input:
        data
    output:
        "plots/plot_{file}.pdf"
    conda:
        "envs/plot.yaml"
    script:
        "scripts/plot.R"
```

_plot.R_:

```r
library(ggplot2)

gene <- read.csv2(snakemake@input[[1]], dec = ",")

pdf(file = snakemake@output[[1]])
ggplot(gene, aes(x = id)) +
    geom_density(alpha=.2, fill="#FF6666") +
    ggtitle(paste0("Plot "), snakemake@wildcards[["file"]])     
dev.off()
```
where

__snakemake@input[[1]]__: the first input file of the above rule _plot_

__snakemake@output[[1]]__: the first output file of the above rule _plot_

__snakemake@wildcards[["file"]]__: the wildcard variale which itself contains the __file__ wildcard for the above rule _plot_

To execute the example, run:

```snakemake --core 1 -j 2 --use-conda --printshellcmds -- all```

The result of execution:

```js
[user@slurmgate 05]$ snakemake --core 1 -j 2 --use-conda --printshellcmds -- all
Building DAG of jobs...
Creating conda environment envs/plot.yaml...
Downloading and installing remote packages.
Environment for envs/plot.yaml created (location: .snakemake/conda/c00adc15fcb7ce78cd2de9774da7df74)
Using shell: /bin/bash
Provided cores: 2
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	2	plot
	3
Select jobs to execute...

[Thu May  6 18:09:20 2021]
rule plot:
    input: data/x.csv
    output: plots/plot_x.pdf
    jobid: 1
    wildcards: file=x

[Thu May  6 18:09:20 2021]
rule plot:
    input: data/y.csv
    output: plots/plot_y.pdf
    jobid: 2
    wildcards: file=y

Rscript --vanilla /05/.snakemake/scripts/tmp1imjv2r6.plot.R
Rscript --vanilla /05/.snakemake/scripts/tmpp_00n_lx.plot.R
Activating conda environment: /05/.snakemake/conda/c00adc15fcb7ce78cd2de9774da7df74
Activating conda environment: /05/.snakemake/conda/c00adc15fcb7ce78cd2de9774da7df74
Registered S3 methods overwritten by 'ggplot2':
  method         from
  [.quosures     rlang
  c.quosures     rlang
  print.quosures rlang

...
```
___

### __Exercises ex06:__ _Parallelizing workflow. Locally_

There are several options for parallelizing a workflow. 

On a single node with many cores, you can use:

```snakemake -j $NUM_THREADS```

By default, each rule is assumed to take 1 thread. You can use the __threads:__ keyword to specify how many threads each rule takes.

__For example:__

```py
rule sort_parallel:
    input:
        "data/big_data.txt"
    output:
        "results/big_data_sorted.txt"
    threads: 8
    shell:
        "sort --parallel={threads} {input} > {output}"
```
where

__threads:__ specifies how many cores a rule takes

The __{threads}__ keyword is accessible in the shell directive through lexical scoping.

To execute the example, run:

```snakemake --core 1 -j 8 --use-conda --printshellcmds```

__Note:__ If __-j__ is less than threads, threads will be reduced to __–j__.
___

### __Exercise ex07:__ _Parallelizing workflow. On cluster_

To run the rule on the cluster, use the __--cluster__ to provide specific keywords when submitting:

```sh
snakemake --cluster 'sbatch  --mem=2g -c 1' -j 4 -- all
```

__Snakemake__ creates a job script for each job to be run and uses this command to submit that job.

__-j 4__ argument limits to at most 4 jobs running at the same time.

__Note:__ Using __--cluster__ like shown above will cause all jobs to request same amount of resources from HPC. To make __Snakemake__ request different amount of resources based on each job, you can:

1. Define resources in rule __parameters__
2. Define resources in the __cluster configuration file__. e.g. __cluster.config__:

___
__Define resources in rule parameters__

The string following __--cluster__ can access all variables within a rule in the same way as in a shell command within a rule:

```py
rule mapping:
    input:
        fastq="{id}_tmp_{subsample}.fq",
        index="data/index"
    output: 
        "results/{id}_subsample_{subsample}_Aligned.out.sam"
    params:
        subsample = "{subsample}",
        mem = "2g"
    threads: 4
    shell:
        "STAR "
        "--genomeDir {input.index} "
        "--readFilesIn {input.fastq} "
        "--outFileNamePrefix result/{wildcards.id}_subsample_{params.subsample}_"
```

Then invoke workflow:

```snakemake --cluster 'sbatch -t --mem={params.mem} -c {threads}' -j 4```

__Note:__ These __params__ must be defined for all rules or you will get an error when Snakemake tries to submit a job from a rule that doesn’t have one or more of these params defined.
____
__Define resources in the cluster configuration file__

Instead of putting parameters for cluster submission into each individual rule, you can write a separate config file that contains these parameters.

```yaml
__default__:
  partition: "pe"
  cpus: "{threads}"
  mem: "2g"

get_index:
  mem: "8g"

mapping:
  mem: "4g"
```

By default, rules take resources corresponding to those annotated in __default__.

For specific rules that require more resources, you can annotate them specifically using __the rule name__ in the config file (e.g. get_index or mapping).

Now in the rule definition, you only need to specify threads (to allow __Snakemake__ utilize this paramter for multi-threading when running locally)

```py
rule mapping:
    input:
        fastq="{id}_tmp_{subsample}.fq",
        index="data/index"
    output: 
        "results/{id}_subsample_{subsample}_Aligned.out.sam"
    params:
        "{subsample}"
    threads: 4
    shell:
        "STAR "
        "--genomeDir {input.index} "
        "--readFilesIn {input.fastq} "
        "--outFileNamePrefix result/{wildcards.id}_subsample_{params}_"
```
___

__Redirect cluster output__

You can tell __Snakemake__ where to place __stdout__ and __stderr__ by specifiynf the output location in the cluster configuration file:

```yaml
__default__:
  partition: "hp"
  cpus: "{threads}"
  mem: "2g"
  output: "logs_slurm/{rule}_{wildcards.id}_{wildcards.subsample}.out"
  error: "logs_slurm/{rule}_{wildcards.id}_{wildcards.subsample}.err"
```

Then run the workflow:

```snakemake -d ~/path_to_ex/07 --cluster-config cluster.yaml --cluster 'sbatch --mem={cluster.mem} -c {cluster.cpus} -o {cluster.output} -e {cluster.error}' -j 4 --printshellcmds -- all```

where 

__d__: specifies the path to the working directory.

This will redirect __stdout__ and __stderr__ of every job to their corresponding __log__ files named as by rule and wildcards under __logs_slurm__ directory. 

### __Exercise ex08:__ _Debugging the Snakefile_

Use the __log:__ keyword to capture STDERR automatically, or you can print to it explicitly, e.g.:

```py
my_vars = ["A", "B", "C"]

rule incorrect:
    output: 
        "results/test.txt"
    log: "log.txt"
    run:
        with open(output[0], "w") as outh, open(log[0], "w") as logh:
            print("Running incorrect rule...", file=logh)
            print(my_vars[3], file=outh)
```

You can use __pdb__, a Python debugger, to step through rules and examine variables at each step:

```py
my_vars = ["A", "B", "C"]

rule incorrect:
    output: 
        "results/test.txt"
    run:
        import pdb; pdb.set_trace()
        with open(output[0], "w") as outh, open(log[0], "w") as logh:
            print("Running incorrect rule...", file=logh)
            print(my_vars[3], file=outh)
```
Use the __--debug__ flag to enable __pdb__ in ```run``` or ```script``` directives:

```snakemake -s Snakefile_incorrect --core 1 --debug```

The output:

```js
(py38) [ahryhorzhevska@slurmgate 08]$ snakemake -s Snakefile_incorrect --core 1 --debug
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	incorrect
	1

[Sat May  8 14:27:26 2021]
rule incorrect:
    output: results/test.txt
    log: log.txt
    jobid: 0

Job counts:
	count	jobs
	1	incorrect
	1
> /home/ahryhorzhevska/mpip/code/IT-Wiki/04 MPIP cluster/snakemake-tutorial/exercises/08/Snakefile_incorrect(16)__rule_incorrect()
(Pdb) my_vars
['A', 'B', 'C']
(Pdb) my_vars[1]
'B'
(Pdb) my_vars[3]
*** IndexError: list index out of range
```
The correct rule with the error fixed:

```py
my_vars = ["A", "B", "C"]

rule correct:
    output: "results/test.txt"
    log: "log.txt"
    run:
        with open(output[0], "w") as outh, open(log[0], "w") as logh:
            print("Running correct rule...", file=logh)
            print(my_vars[2], file=outh)
```
__Note:__ for debbuging execution of Snakefiles themselves, use the __--verbose__ option:

```snakemake -s Snakefile_incorrect --core 1 --verbose```

```json
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	incorrect
	1
Resources before job selection: {'_cores': 1, '_nodes': 9223372036854775807}
Ready jobs (1):
	incorrect
Selected jobs (1):
	incorrect
Resources after job selection: {'_cores': 0, '_nodes': 9223372036854775806}

[Sat May  8 14:33:36 2021]
rule incorrect:
    output: results/test.txt
    log: log.txt
    jobid: 0

Job counts:
	count	jobs
	1	incorrect
	1
[Sat May  8 14:33:36 2021]
Error in rule incorrect:
    jobid: 0
    output: results/test.txt
    log: log.txt (check log file(s) for error message)

RuleException:
IndexError in line 13 of /home/ahryhorzhevska/mpip/code/IT-Wiki/04 MPIP cluster/snakemake-tutorial/exercises/08/Snakefile_incorrect:
list index out of range
  File "/home/ahryhorzhevska/mpip/code/IT-Wiki/04 MPIP cluster/snakemake-tutorial/exercises/08/Snakefile_incorrect", line 13, in __rule_incorrect
  File "/usr/lib64/python3.6/concurrent/futures/thread.py", line 56, in run
Exiting because a job execution failed. Look above for error message
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: /home/ahryhorzhevska/mpip/code/IT-Wiki/04 MPIP cluster/snakemake-tutorial/exercises/08/.snakemake/log/2021-05-08T143336.287182.snakemake.log
unlocking
removing lock
removing lock
removed all locks
```

# Useful command line arguments

- Use ```--forceall``` to rerun an analysis from scratch
- Use ```--keep-going``` to prevent your workflow from stopping if a step fails
- Use ```--until``` or ```--omit-from``` to stop the pipeline at a certain step.
- Use ```--delete-all-output``` to remove all files generated by the workflow
- Use ```--rulegraph | dot –Tpdf > ruledag.pdf``` to print an overview of the workflow
- Use ```--archive``` to create a gzipped tarball of the entire workflow, including Conda packages

# Useful keywords and functions

- Use the __priority__ keyword or __--prioritize TARGET__ on the command line to increase the priority of a file, e.g. if your PI is bugging you about some particular figure)
- Use the __temp()__ function around files that you do not need, and they will be automatically deleted when no rules depend on them anymore
- Use __protected()__ around important files to protect them from accidental deletion
- Use __directory()__ to indicate the creation of a directory (by default, Snakemake only recognizes files)
- Use the __remote()__ function to download files via S3, Google, Dropbox, https, ftp, etc.
- Use the __include__ keyword to include sub-Snakefiles and better organize your code

# Useful resources

- [Snakemake on Stackoverflow](https://stackoverflow.com/questions/tagged/snakemake)
- [Snakemake Github](https://github.com/snakemake/snakemake)
- [Bioconda](https://bioconda.github.io/)
- [Cluster profiles](https://github.com/Snakemake-Profiles)
- [Snakemake wrappers](https://snakemake-wrappers.readthedocs.io/)

