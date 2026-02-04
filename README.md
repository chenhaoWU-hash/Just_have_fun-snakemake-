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
The next step in our workflow will aggregate the mapped reads from all samples and jointly call genomic variants on them. For the variant calling, we will combine the two utilities samtools and bcftools.  
Samtools: Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format  
BCFtools: Reading/writing BCF2/VCF/gVCF files and calling/filtering/summarising SNP and short indel sequence variants  

Add the sample list at the very top of the Snakefile with:
```SAMPLES = ["A", "B"]```  
Then add the following rule:
```rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```
because for variant calling, we need to combined all sorted BAM files (A, B...) as input for the same rule, that's why we need expand() explicitly tell it which samples are available. Output: calls/all.vcf (only one file, without {sample}) It requires listing all input files for each sample, such as:

```sorted_reads/A.bam```, ```sorted_reads/B.bam```, along with their corresponding indexes: ```sorted_reads/A.bam.bai```, ```sorted_reads/B.bam.bai``` Generate a Python list: ```["sorted_reads/A.bam","sorted_reads/A.bam.bai","sorted_reads/B.bam","sorted_reads/B.bam.bai"]```  
Snakemake won't searching the require file by itself, so we need to tell it which input are we looking for, and the sample names contained with A & B in this tutorial.

***exercise***  
obtain the updated DAG of jobs for the target file calls/all.vcf, it's actually the process we have so far.  
```snakemake calls/all.vcf --dag | dot -Tsvg > dag_all.svg```  

**Step 6 Using custom scripts**  
Although Snakemake permits Python code to be written directly within rules, it is more advisable to place such logic in separate script files and invoke them within rules using the `script:` directive.  
```rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```
The actual Python code to generate the plot is hidden in the script scripts/plot-quals.py. Script paths are always relative to the referring Snakefile. Create the file scripts/plot-quals.py, with the following content:
```
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])
```
Here we want to generate a quals plot by using a python script, In Python scripts triggered by `script:`, `snakemake.input` and `snakemake.output`, only indicates the list of input and output files defined and parsed for the current rule (current job) within the Snakefile. It is also possible to use R scripts.  

**Step 7 Adding a target rule**  
Snakemake can use filenames or rule names without wildcards as run targets. If no target is specified on the command line, it defaults to executing the first rule in the Snakefile. Therefore, it is advisable to place a rule named `all` at the top, which aggregates all final outputs, as the default target.  
```rule all:
    input:
        "plots/quals.svg"
```  
Snakemake resolves dependencies by constructing a DAG based on input/output filename matching, rather than determining execution order through rule sequencing.

Now, we excuted with:  
``` snakemake -n```  
In case you have multiple reasonable sets of target files, you can add multiple target rules at the top of the Snakefile. While Snakemake will execute the first per default, you can target any of them via the command line (for example, snakemake -n mytarget).Therefore, when no target is specified, Snakemake treats the default target as "generating plots/quals.svg" and prints the execution plan detailing which rules/steps need to be executed for this purpose.  

***exercise***  
Create the DAG of jobs for the complete workflow.  
```snakemake --dag | dot -Tsvg > dag.svg```  

Execute the complete workflow and have a look at the resulting plots/quals.svg.  
```snakemake --cores 8```

Snakemake provides handy flags for forcing re-execution of parts of the workflow. Have a look at the command line help with snakemake --help and search for the flag --forcerun. Then, use this flag to re-execute the rule samtools_sort and see what happens.  
```snakemake -np --forcerun samtools_sort```  

Snakemake displays the reason for each job (under reason:). Perform a dry-run that forces some rules to be reexecuted (using the --forcerun flag in combination with some rulename) to understand the decisions of Snakemake.
```snakemake -np --forcerun samtools_sort,bcftools_call```  





