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

**Step 4 Indexing read alignments and visualizing the DAG of jobs**  
***exercise***  
**Step 5 Calling genomic variants**  
***exercise***  
```snakemake calls/all.vcf --dag | dot -Tsvg > dag_all.svg```  

**Step 6 Using custom scripts**  
**Step 7 Adding a target rule**  
***exercise***  


