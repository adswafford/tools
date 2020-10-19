__First, start by creating a conda environment with the needed packages:__ 

```bash
conda create --name ebi-sra-importer -c bioconda xmltodict lxml pandas requests sra-tools entrez-direct pyyaml xlrd
```

Then activate the new conda environment (`source activate ebi-sra-importer`), cd to your directory of choice, and run the script. __Note, if you want to first check the metadata before downloading a large number of files, pass the --no_seqs flag.__

__To run the python script:__  
	python3 EBI_SRA_Downloader.py -project [ebi or sra accession_number]

__Available flags:__  
 -project [ebia or sra ccession_number]  
 required flag for script to download EBI/SRA study, can pass mulitple
 e.g. python3 EBI_SRA_Downloader.py -project [accession_number] [accession_number ...N]
 
 -output [directory where files will be saved]
 
 -mode [specifies which repository to use]                    
 
 -prefix [list of prefixes for sample and prep info files]
 
 --strategy [list of one or more library strategies to select]
 
 --sources [list of one or more library sources to select]
 
 --platforms [list of one or more sequencing platforms to select]
 
 --validators [list of one or more yaml files to use in validating]
 
 --no_seqs [skip downloading files]
 
 --verbose          
__Advanced use:__
--force
Note: requires a single file provided to the -yaml flag or a single file in the --yaml_dir directory
this will be used for validation regardless of the scientific_name found in the study


__Generated files after running the script:__  
 
 [accession_number]_prep_info.txt  
 	The information about each run in each sample in the study.  
 	The following information is listed in this file (in order):  
 		sample_name, run_prefix, platform, library_strategy, library_source, library_layout  
 	run_prefix is made up by [sample accession].[run accession].

 [accession_number]_sample_info.txt  
 	The information about each sample in the study.  
 	The contained information deffers according to the information avalibility online.  
 	The generated sample_info.txt from SRA study contains the information about each run in the study.
