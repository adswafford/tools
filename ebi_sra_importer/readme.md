First start by creating a conda environment with the needed packages:  
conda create -name ebi-sra-importer -c bioconda xmltodict lxml pandas requests sra-tools entrez-direct

Then source activate ebi-sra-importer
cd to your directory of choice and then run the script. __Note, if you want to first check the metadata before downloading a large number of files, pass the -debug true parameter.__

__To run the python script:__  
	python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number]  
	python3 EBI_SRA_Downloader.py -sra [sraaccession_number]

__Available flags:__  
 -ebi [ebiaccession_number]  
 required flag for script to download EBI study  
 e.g. python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number]

 -sra [sraaccession_number]  
 required flag for script to download SRA study  
 e.g. python3 EBI_SRA_Downloader.py -sra [sraaccession_number]

 -sample [sample_file_name]  
 optional flag for user to customize sample file name  
 e.g. python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number] -sample [sample_file_name]

 -prep [prep_file_name]  
 optional flag for user to customize prep file name  
 e.g. python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number] -prep [prep_file_name]

 -study [study_info_file_name]  
 optional flag for user to customize study info file name  
 e.g. python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number] -study [study_info_file_name]

 -all-seqs true  
 optional flag for user to filter out non-metagenomic samples  
 without setting the flag, the non-metagenomic samples will be filtered out  
 e.g. python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number] -all-seqs true

 -all-platform true  
 optional flag for user to filter out non-illumina samples  
 without setting the flag, the non-illumina samples will be filtered out  
 e.g. python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number] -all-platform true


 -debug true  
 optional flag for developer to enter debug mode  
 debug mode is simply for developer to debug, which will not download massive size data files  
 e.g. python3 EBI_SRA_Downloader.py -ebi [ebiaccession_number] -debug true

__Generated files after running the script:__  
 [accession_number]_detail.txt  
 	The details for the corresponding study.  
 	The following information is listed in this file (in order):  
 		Library name, Sample accession, Run accession, Download path, Library source, Platform,  
 		Submitted format, Library strategy, Library layout

 [accession_number]_study_info.txt  
  The information about the corresponding study.  
  The following information is listed in this file (in order):  
 		Study title, Alias, Study abstract, Study description, Principal investigator, Environmental packages  
 	PI and EP are default value for now.

 [accession_number]_prep_info.txt  
 	The information about each run in each sample in the study.  
 	The following information is listed in this file (in order):  
 		sample_name, run_prefix, instrument_platform, library_strategy, library_source, library_layout  
 	run_prefix is made up by [sample accession].[run accession].

 [accession_number]_sample_info.txt  
 	The information about each sample in the study.  
 	The contained information deffers according to the information avalibility online.  
 	The generated sample_info.txt from SRA study contains the information about each run in the study.

 [accession_number]/[library_name]/[EBI_sample_accession].[EBI_run_accession].txt  
 [accession_number]/[library_name]/[SRA_run_accession].txt  
 	These files are in xml format and contain the information about the run accession.

 [accession_number]/[library_name]/[EBI_sample_accession].[EBI_run_accession].fastq  
 [accession_number]/[library_name]/[EBI_sample_accession].[EBI_run_accession].sff  
 	These files are the FASTQ or SFF files for each EBI run accession.

 FASTQ files from SRA runs are placed in the folder where the script is placed.  
