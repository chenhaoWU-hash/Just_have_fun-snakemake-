# Just_have_fun-snakemake-

This is an exercise following the official Snakemake tutorial, along with some notes from the practice.(https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).
It is good to know that the exemple from the tutorial is actually in bioinformatics,which help me understanding faster this workflow.  

In setup parts,I chose pixi for the package management. I use vscode for daily coding, and need to connect to the selected terminal via the icon in the bottom-left corner. this is use to determine the working path.  As I'm working in my personal computer, and I already use linux before, setting up the enviroment are very smoothly by following the steps:  

```curl -fsSL https://pixi.sh/install.sh | bash```  
``` mkdir snakemake-tutorial ```  
```cd snakemake-tutorial```  

In this new directory,download and extract the data via:  
```curl -L https://api.github.com/repos/snakemake/snakemake-tutorial-data/tarball -o snakemake-tutorial-data.tar.gz```  

```tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"```  

This will create a folder data and a file environment.yaml in the working directory.  

Then create and activate the env:  
```pixi init --import environment.yaml```  
```pixi shell```  

Now it's time to start.

**Step 1 Mapping reads**  
Create a new file called Snakefile. The first Snakemake rule maps reads of a given sample to a given reference genome. For this we used tool bwa,specifically the subcommand bwa mem (already install from pixi when we prepare the enviroment). Then define the first rule bwa_map in Snakefile.  

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
excude the first workflow, Snakemake tries to generate given target files. Target files can be specified via the command line(A.bam files is provide by this tutorial in mapped_read directory).  
```snakemake -np mapped_reads/A.bam```  

**Step 2 Generalizing the read mapping rule**  
Snakemake allows generalizing rules by using named wildcards. Simply replace the A in the second input file and in the output file with the wildcard {sample}  

```rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```  
Snakemake infers the value of the wildcard based on the specific output filename you request (e.g. {sample}=A).This value will automatically substitute into all {sample} positions within the input/output, determining which inputs are required for this job. So to prevent multiple jobs from writing to the same file in parallel, all outputs of the same rule must contain an identical set of wildcards. For exemple: 
```snakemake -np mapped_reads/B.bam```  
Snakemake will determine that the rule bwa_map can be applied to generate the target file by replacing the wildcard {sample} with the value B.  

We can also specify multiple targets:  
```snakemake -np mapped_reads/A.bam mapped_reads/B.bam```  
```snakemake -np mapped_reads/{A,B}.bam```  

In both cases, Snakemake only proposes to create the output file mapped_reads/B.bam.  bcs we executed the workflow before and no input file is newer than the output file mapped_reads/A.bam. 

to modify update the file modification date of the input file, we can:  
```touch data/samples/A.fastq``` then rerun with ``` snakemake -np mapped_reads/A.bam mapped_reads/B.bam```  

**Step 3 Sorting read alignments**  
 Sorted the read alignments in the BAM files with the samtools sort command. We add the following rule beneath the bwa_map rule:  

 ```rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```
Snakemake automatically creates missing directories before jobs are executed. In addition to using the string format {sample} (input/output), Snakemake also provides an object named `wildcards` {wildcards.sample} for accessing wildcard values.

**Step 4 Indexing read alignments and visualizing the DAG of jobs**  
Next, we need to use samtools again to index the sorted read alignments so that we can quickly access reads by the genomic location they were mapped to. This can be done with the following rule:  

```rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
```  
Having three steps already, it is a good time to take a closer look at the resulting directed acyclic graph (DAG) of jobs. By executing:  

```snakemake sorted_reads/{A,B}.bam.bai --dag | dot -Tsvg > dag.svg```  

We create a visualization of the DAG using the dot command provided by Graphviz.The DAG contains a node for each job with the edges connecting them representing the dependencies. The frames of jobs that don’t need to be run (because their output is up-to-date) are dashed. For rules with wildcards, the value of the wildcard for the particular job is displayed in the job node. 


***exercise***  
Run parts of the workflow using different targets. Recreate the DAG and see how different rules’ frames become dashed because their output is present and up-to-date.  

For details, please refer to dag.svg and dag2.svg  by running the code with different status.  

**Step 5 Calling genomic variants**  
***exercise***  
```snakemake calls/all.vcf --dag | dot -Tsvg > dag_all.svg```  

**Step 6 Using custom scripts**  
**Step 7 Adding a target rule**  
***exercise***  


