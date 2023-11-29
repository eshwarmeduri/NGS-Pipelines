eCLIP pipeline

Workflow is incorporated from ENCODE pipeline - https://www.encodeproject.org/pipelines/ENCPL357ADL/

This nextflow script does the following

1) Trim adapters

2) Dicards the reads mapped to Repetetive elements and maps to the genome

3) Removes PCR duplicates


Create a directory “Files” to store Fastq files. In the nextflow script pairedend Fastq files Paired end sequencing files have an extension *.s_1.r_1.fq.gz and *.s_1.r_2.fq.gz. If the fastq files have different naming conventions please make the changes in the Nextflow script accordingly

Keep nextflow script and nextflow.config in the same folder

If the nextflow is installed locally, please add the below text in .bashrc file. In the below example nextflow execultable is in /data/Softwares/.

export PATH=${PATH}:/data/Softwares/Nextflow

Output files are generated in the current working directory with sub folders named after each sample (For example: Output/Sample1, Output/Sample2……..)

Use the approproprite reference genome (hg38 or mm39). Prior to running this program reference genome should be generated keeping read length in mind.

How to run the pipeline:

nextflow run eCLIP-Pipeline.nf

Important notes:

• Please be mindful of the shared resources on the server. Nextflow is designed for parallel processing. So don’t run the pipeline for more than 10 samples at a time. Here we are using only one processor

• Nextflow also creates a directory called “work” for all the output files. As the output files are also stored in a directory called Output/, delete the content in directory “work” by typing below commands

nextflow clean -f (OR) rm -rf work 

Peak calling can be done using clipper - https://github.com/YeoLab/clipper
