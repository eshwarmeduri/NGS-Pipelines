
HiChIP-Seq pipeline


Create a directory “Files” to store Fastq/fq files. Paired end sequencing files from the CI are usually named as *.s_1.r_1.fq.gz and *.s_1.r_2.fq.gz. If the fastq files have different naming conventions make the necessary changes in the script.

>	Please run only one sample at a time

>	Keep nextflow script and nextflow.config in the same folder

>	Nextflow executable is in /data/Softwares/. Please add this text to .bashrc file. 

		export PATH=${PATH}:/data/Softwares/Nextflow

>	Output files are generated in the current working directory with sub folders named after each sample (For example: Output/Sample1, Output/Sample2……..)

>	Please note that the current reference genome is hg38 

How to run the pipeline:

	nextflow run HiChIP-Pipeline.nf



Important notes: 

•	Please be mindful of the shared resources on the server. Nextflow is conducive for parallel processiing of samples. Run one sample at a time

•	Nextflow also creates a directory called “work” for all the output files. As the output files are also stored in a directory called Output/, delete the content in directory “work” by typing below commands

o	nextflow clean -f
o	rm -rf work 


Notes

https://hichip.readthedocs.io/en/latest/index.html

