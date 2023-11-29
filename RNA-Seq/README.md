
RNA-Seq pipeline

This nextflow script does the following

1) Trim adapters and run FASTQC

2) Maps the reads to the genome using STAR

3) Removes PCR duplicates

4) Get read counts

5) Converts BAM to bigwig to visualise in IGV


Create a directory “Files” to store Fastq files. In the nextflow script pairedend Fastq files Paired end sequencing files have an extension *.s_1.r_1.fq.gz and *.s_1.r_2.fq.gz. If the fastq files have different naming conventions please make the changes in the Nextflow script accordingly

1) Keep nextflow script and nextflow.config in the same folder

2) If the nextflow is installed locally, please add the below text in .bashrc file. In the below example nextflow execultable is in /data/Softwares/. 

	export PATH=${PATH}:/data/Softwares/Nextflow

3) Output files are generated in the current working directory with sub folders named after each sample (For example: Output/Sample1, Output/Sample2……..)

4) Use the approproprite reference genome (hg38 or mm39). Prior to running this program reference genome should be generated keeping read length in mind.


5) How to run the pipeline:

	nextflow run RNASeq-pipeline.nf



Important notes: 

•	Please be mindful of the shared resources on the server. Nextflow is designed for parallel processing. So don’t run the pipeline for more than 10 samples at a time. Here we are using only one processor 

•	Nextflow also creates a directory called “work” for all the output files. As the output files are also stored in a directory called Output/, delete the content in directory “work” by typing below commands

	nextflow clean -f (OR) rm -rf work 

Notes for Differential gene expression analysis

https://gist.github.com/stephenturner/f60c1934405c127f09a6

