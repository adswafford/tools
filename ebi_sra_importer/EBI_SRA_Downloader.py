#!/usr/bin/env python3
# conda command to install dependencies for amplicon and stool files (no host depletion needed):
#   conda create -n ebi_sra_importer pandas requests entrez-direct sra-tools
#   xmltodict lxml pyyaml xlrd -c bioconda -c conda-forge -y
# to enable host depletion use:
#   conda create -n ebi_sra_importer pandas requests entrez-direct sra-tools
#   bowtie2 minimap2 samtools bedtools xmltodict lxml pyyaml xlrd fastqc multiqc
#   -c bioconda-c conda-forge -y
# or to enable later after installation:
#   conda activate ebi_sra_importer
#   conda install bowtie2 minimap2 samtools bedtools fastqc multiqc
#
# pip commands to install all dependencies:
#   pip install csv glob requests subprocess xmltodict sys lxml os urllib pyyaml xlrd
#   pip install argparse pandas bioconda sra-tools entrez-direct
#   pip install bowtie2 minimap2 samtools bedtools
#   pip install fastqc multiqc
#
# Instruction:
#       python EBI_SRA_Downloader.py -project [accession] [accession ... N]
#                   Generate the study info, study detail, prep, and  sample
#                   files for the entered EBI accession, and download the
#                   FASTQ files.
#               Optional flags:
#                   -output [directory where files will be saved]
#                   -mode [specifies which repository to use]                    
#                   -prefix [list of prefixes for sample and prep info files]
#                   --strategy [list of one or more library strategies to select]
#                   --platforms [list of one or more sequencing platforms to select]
#                   --validators [list of one or more yaml files to use in validating]
#                   --no_seqs [skip downloading files]
#                   --verbose          
#

# libraries used
import csv
import sys
import glob
import logging
import requests
import subprocess
from xmltodict import parse
from lxml import etree
import os
from os import path, makedirs, remove
import shutil
from pathlib import Path
from urllib.request import urlretrieve
from argparse import ArgumentParser
import pandas as pd
from pandas import read_csv, DataFrame
import numpy as np
import re
import yaml

### START: GENERIC CLEANUP METHODS

def split_caps(cap_string):
    """Splits up a string based on capital letters

    Metadata from SRA tends to include multiple entries destinguished by their
    capitalization. This method splits up these strings and returns them
    with underscore delimiting.

    Parameters
    ----------
    cap_string : string
        Thestring to be parsed

    Returns
    -------
    string
        underscore delimited string
    """
    temp_list = []
    temp_list = re.findall('[A-Z][^A-Z]*', cap_string)
    if len(temp_list) > 0:
        return "_".join(temp_list)
    else:
        return cap_string
        
def scrub_special_chars(input_string,custom_dict={}):
    """Removes special characters from a string

    This method scrubs special characters from any input string and
    replaces them with underscores for punctuation and spelled-out 
    strings for other characters. Users may supply a custom dictionary
    which will be applied prior to the base rules to ensure that user
    preference for renaming is respected, e.g. "/" may be converted to
    '_or_' by the user to avoid conversion to '_per' by the default
    dictionary. Returned strings are also compatible with Qiita column
    requirements.   

    Parameters
    ----------
    input_string : string
        The string to be scrubbed

    Returns
    -------
    string
        returns the scrubbed string
    """
    
    replace_dict ={"__":"_",
                   " ":"_",
                   "-": "_",
                   "(": "_leftparen_",
                   ")": "_rightparen_",
                   "/": "_per_",
                   "-":"_",
                   "|":"_bar_",
                   "~":"_tilde_",
                   "`":"_",
                   "@":"_at_",
                   "#":"_number_",
                   "$":"dollar",
                   "%":"_perc_",
                   "^":"_carrot_",
                   "&":"_and_",
                   "*":"_star_",
                   "+":"_plus",
                   "=":"_equals",
                   "\\":"_per_",
                   "{":"_leftbracket_",
                   "}":"_rightbracket_",
                   "[":"_leftbracket_",
                   "]":"_rightbracket_",
                   "?":"_question",
                   "<":"_less_",
                   ">":"_more_",
                   ",":"_",
                   ".":"_",}
    for k in custom_dict.keys():
        input_string=input_string.replace(k,custom_dict[k])
    for k in replace_dict.keys():
        input_string=input_string.replace(k,replace_dict[k])
    
    return input_string

def sub_sample_study(sample_df,max,rand_sample):
    """Subsamples a dataframe for the specified number of samples

    Simple helper function to take a dataframe and return the top or
    a random subset based on the number requested

    Parameters
    ----------
    sample_df : pd.DataFrame
        The dataframe to be subsampled
    max: int
        The number of samples to take
    rand_sample: boolean
        Whether to perform a random subsample

    Returns
    -------
    sample_df
        subsampled dataframe
    """
    if rand_sample:
        return sample_df.sample(max)
    else:
        return sample_df.head(max)

### END: GENERIC CLEANUP METHODS

### START: NORMALIZATION STAGING METHODS

def qiimp_parser(filename):
    """Parses Qiimp-format yaml or xlsx files to obtain a dictionary of
    rules for metadata normalization

    This method will accept an input xlsx from Qiimp (qiita.ucsd.edu/qiimp)
    with a yaml formatted set of rules in cell A1 of the Sheet metadata_schema,
    or a yaml/yml file with the same format which may be created by extracting
    the rules from such a Qiimp template into a new yaml/yml file. Throws errors
    if unexpected .xlsx formatted file or corrupted yaml is found.
    
    N.B. this method can be used to parse any yaml/yml file into a Python dict.
    
    Parameters
    ----------
    filename : string
        The filename of the Qiimp .xlsx or yml to be parsed

    Returns
    -------
    parse_yml
        returns the scrubbed string
    """
    ext=Path(filename).suffix
    parsed_yml={}
    if ext == '.xlsx':
        try:
            temp_yml=pd.read_excel(filename,sheet_name='metadata_schema',header=None) #assume QIIMP-style excel
            parsed_yml = yaml.load(temp_yml.at[0,0])
        except:
            logger.warning("Invalid .xlsx file. Please ensure file is from QIIMP or contains a compliant yaml " + 
                          "in cell A1 of the sheet labelled 'metadata_schema'.")
    elif ext == '.yml' or ext == '.yaml':
        try:
            with open(filename) as file:
            # The FullLoader parameter handles the conversion from YAML
            # scalar values to Python dictionary format
                parsed_yml=yaml.load(file, Loader=yaml.FullLoader)
                if DEBUG: logger.info(parsed_yml)
        except:
            logger.warning("Could not load yaml from " + filename + "Contents:\n")
            logger.warning(subprocess.run(['cat',filename]))
    else:
        logger.warning("Invalid file extension for yaml parsing: " + str(ext))

    if len(parsed_yml) == 0:
        logger.warning("The file " + filename +" contains no yaml data. Please check contents and try again.")

    return parsed_yml

def set_yaml_validators(validators,fields=['scientific_name']):
    """Used to set up dictionary of validators for metadata normalization
    
    This method takes a set of supplied yaml, yml or xlsx files and looks
    for a scientific_name field to serve as the key for a dictionary of
    validators used to normalize the metadata. All files parsed are expected to
    be in Qiimp-format and validators that do not contain a default for 
    scientific_name will be omitted. Users may supply a custom field to look for
    instead of scientific_name for key creation.
    
    #TODO: this method currently cannot differentiate between sufficiently vague
    scientific names, e.g. human metagenome for multiple body sites
    skin, vaginal, etc.

    
    Parameters
    ----------
    validators : list
        The filenames of the Qiimp .xlsx or yml to be stored

    Returns
    -------
    yaml_dict
        returns the dictionary of validators
    """
    yaml_dict={}
    for v in validators:
        new_yaml={}
        #try:
        new_yaml = qiimp_parser(v)
        #except:
        #    logger.warning("Could not open yaml file " + v + " Please ensure the file exists and is a valid yaml file.")

        for field in fields:
            if field not in new_yaml.keys():
                logger.warning("Invalid validator yaml. Please ensure a default value is provided for '" + field +
                      "'. Available keys in file: " + str(new_yaml.keys()))
            elif 'default' not in new_yaml[field].keys():
                logger.warning("Invalid validator yaml. Please ensure a default value is provided for '" + field +
                      "'. Available keys in : " + field + str(new_yaml[field].keys()))
            else:
                yaml_dict[new_yaml[field]['default']] = v

    return yaml_dict

def set_empo_normalizers(sample_type_map):
    """Reads in an mapping file that for EMPO-related comparisons

    This method will accept an input tab or comma-delimited file that contains
    the column sample_type which will be used to normalize sample_type temrs
    and EMPO fields to support better metaanalysis across public studies. Returns
    a pandas DataFrame for use in normalization
    
    Parameters
    ----------
    sample_type_map : string
        The filename of the mapping file to be read

    Returns
    -------
    empo_type_df
        returns the DataFrame for mapping
    """
    empo_type_df =pd.DataFrame()
    empo_type_df=pd.read_csv(sample_type_map, sep='\t', dtype=str)
    if 'sample_type' not in empo_type_df: #try to load again as a csv
        empo_type_df = pd.read_csv(sample_type_map)

    if 'sample_type' not in empo_type_df:
        logger.warning("sample_type not found in sample type mapping file: " + sample_type_map + "\nAvailable columns: " + empo_type_df.columns)
    else:
        empo_type_df=empo_type_df.set_index('sample_type')
        logger.info("Sample type mapping file found. Adding normalization columns: " + empo_type_df.columns)

    return empo_type_df


### END: NORMALIZATION STAGING METHODS

### START: NORMALIZATION METHODS
def set_prep_type(expt_df,row):
    """Maps EBI library strategies to Qiita prep types

    Method to map  EBI library strategies to Qiita prep types. For amplicon
    and 'other' types, also looks for target_gene to set prep type. Any
    ambigious prep type is labeled as such and this dictionary can be updated
    as new prep types are supported. Also attempts to account for common typos
    and format changes for target_gene.

    Parameters
    ----------
    expt_df : pd.DataFrame
        The dataframe of experiment data from EBI to be normalized
    row: index value
        The specific line of the dataframe to look up

    Returns
    -------
    library_strat_to_qiita_dict[ebi_strat] OR
    amplicon_dict[tg] OR
    'AMBIGUOUS': string
        Qiita-normalized prep type term
    """
    amplicon_list = ['AMPLICON','OTHER']
    amplicon_dict= {'16S rRNA':'16S',
                    '16S':'16S',
                    '16S rDNA':'16S',
                    '16S RNA':'16S',
                    '16S DNA':'16S',
                    '16S ':'16S',
                    'ITS':'ITS',
                    'ITS1':'ITS',
                    'ITS2':'ITS',
                    'ITS ':'ITS',
                    '18S rRNA':'18S',
                    '18S rDNA':'18S',
                    '18S':'18S',
                    '18S RNA':'18S',
                    '18S DNA':'18S',
                    '18S ':'18S',
                    }
    library_strat_to_qiita_dict={'POOLCLONE':'Metagenomic',
        'CLONE':'Metagenomic',
        'CLONEEND':'AMBIGIOUS',
        'WGS':'Metagenomic',
        'WGA':'Metagenomic',
        'WCS':'AMBIGIOUS',
        'WXS':'Metatranscriptomic',
        'ChIP-Seq':'Metagenomic',
        'RNA-Seq':'Metatranscriptomic',
        'MRE-Seq':'AMBIGIOUS',
        'MeDIP-Seq':'AMBIGIOUS',
        'MBD-Seq':'AMBIGIOUS',
        'MNase-Seq':'AMBIGIOUS',
        'DNase-Hypersensitivity':'AMBIGIOUS',
        'Bisulfite-Seq':'Metgenomic',
        'EST':'AMBIGIOUS',
        'FL-cDNA':'AMBIGIOUS',
        'miRNA-Seq':'Metatranscriptomic',
        'ncRNA-Seq':'Metatranscriptomic',
        'FINISHING':'AMBIGIOUS',
        'CTS':'AMBIGIOUS',
        'Tn-Seq':'AMBIGIOUS',
        'VALIDATION':'AMBIGIOUS',
        'FAIRE-seq':'AMBIGIOUS',
        'SELEX':'AMBIGIOUS',
        'RIP-Seq':'AMBIGIOUS',
        'ChIA-PET':'AMBIGIOUS',
        'RAD-Seq':'AMBIGIOUS',}
    ebi_strat = expt_df.at[row,'library_strategy']

    if ebi_strat not in amplicon_list:
        try:
            return library_strat_to_qiita_dict[ebi_strat]
        except:
            logger.warning(ebi_strat + " not found in EBI-Qiita library strategy mapping. Setting to 'AMBIGUOUS'.")
            return 'AMBIGIOUS'
    else:
        #since this is amplicon data, there should be a target gene, if not return AMBIGUOUS
        try:
            tg=expt_df.at[row,'target_gene']
            return amplicon_dict[tg]
        except:
            logger.warning(ebi_strat + " not found in amplicon resolution mapping. Setting to 'AMBIGUOUS'.")
            return 'AMBIGIOUS'

def normalize_types(md,mapping):
    """Maps EBI library strategies to Qiita prep types

    Method to map  EBI library strategies to Qiita prep types. For amplicon
    and 'other' types, also looks for target_gene to set prep type. Any
    ambigious prep type is labeled as such and this dictionary can be updated
    as new prep types are supported. Also attempts to account for common typos
    and format changes for target_gene.

    Parameters
    ----------
    expt_df : pd.DataFrame
        The dataframe of experiment data from EBI to be normalized
    row: index value
        The specific line of the dataframe to look up

    Returns
    -------
    library_strat_to_qiita_dict[ebi_strat] OR
    amplicon_dict[tg] OR
    'AMBIGUOUS': string
        Qiita-normalized prep type term
    """
    #small modifications from Daniel McDonald normalization code
    simple_sample_type = mapping['simple_sample_type'].to_dict()
    empo_1 = mapping['empo_1'].to_dict()
    empo_2 = mapping['empo_2'].to_dict()
    empo_3 = mapping['empo_3'].to_dict()

    qsst = [simple_sample_type.get(v) for v in md['sample_type']]
    qemp1 = [empo_1.get(v) for v in md['sample_type']]
    qemp2 = [empo_2.get(v) for v in md['sample_type']]
    qemp3 = [empo_3.get(v) for v in md['sample_type']]

    md['qiita_empo_1'] = qemp1
    md['qiita_empo_2'] = qemp2
    md['qiita_empo_3'] = qemp3
    md['simple_sample_type'] = qsst
    
    #to check for empo values, need to ensure columns are there to avoid error
    empo_cols=['empo_1','empo_2','empo_3']
    for e in empo_cols:
        if e not in md.columns:
            md[e]= np.nan
    # if there is an existing empo value, then let's prefer it
    emp1 = md[~md['empo_1'].isnull()]
    emp2 = md[~md['empo_2'].isnull()]
    emp3 = md[~md['empo_3'].isnull()]
    md.loc[emp1.index, 'qiita_empo_1'] = emp1['empo_1']
    md.loc[emp2.index, 'qiita_empo_2'] = emp1['empo_2']
    md.loc[emp3.index, 'qiita_empo_3'] = emp1['empo_3']

    return md

def validate_samples(raw_df,sample_type_col,yaml_validator_dict,prefix):
    """Applies yaml validation rules to metadata

    #TODO fill in description

    Parameters
    ----------
    raw_df: pd.DataFrame
        un-normalized/validated dataframe of metadata    
    sample_type_col : string
        column for identifying sample type 
    yaml_validator_dict : dict
        dictionary of Qiimp-style validators to use for normalization
    prefix: string
        prefix to assign to validator error log, includes path
    
    Returns
    -------
    valid_df: pd.DataFrame
       Validated dataframe
    """
    
    #initialize variables
    st_list = []
    msg = ''
    validator_yaml={}

    for st in raw_df[sample_type_col].unique():
        df_to_validate = raw_df[raw_df[sample_type_col]==st]
        if force:
            validator_yaml = qiimp_parser(yaml_validator_dict[0])
            if 'scientific_name' in df_to_validate.keys():
                df_to_validate=df_to_validate.rename({'scientific_name':'orig_scientific_name'},axis=1)
        else:
            if st not in yaml_validator_dict.keys():
                logger.warning("No yaml file for validating " + st + " You may provide one or more custom files " +
                  " using the --validators flag.")
            else:
                validator_yaml = qiimp_parser(yaml_validator_dict[st])

        for k in validator_yaml.keys():
            if k not in df_to_validate.columns:
                msg= msg + k + ' not found in metadata.\n' 
                try:
                    df_to_validate[k]= validator_yaml[k]['default']
                    msg = msg + "Setting " + k + " to " + validator_yaml[k]['default'] + " for " + st + " samples\n"                     
                except:
                    df_to_validate[k]= 'not provided'
                    msg = msg + k + " has no default in yaml template. Encoding as 'not provided'\n"
            else:
                #construct rules
                uniq = df_to_validate[k].unique()
                allowed_list = []
                min_value= ''
                max_value = ''
                min_value_excl= ''
                max_value_excl = ''
                if 'anyof' in validator_yaml[k].keys():
                    anyof_list = validator_yaml[k]['anyof']            
                    for r in anyof_list:
                        if r['type'] == 'string':
                            for a in r['allowed']:
                                allowed_list.append(a)
                        elif r['type'] == 'number':
                            if 'min' in r.keys():
                                min_value = r['min']
                            if 'max' in r.keys():
                                max_value = r['max']
                            if 'min_exclusive' in r.keys():
                                min_value = r['min_exclusive']
                            if 'max_exclusive' in r.keys():
                                max_value = r['max_exclusive']
                elif validator_yaml[k]['type'] in validator_yaml[k].keys():
                    if validator_yaml[k]['type']== 'string':
                        allowed_list=validator_yaml[k]['allowed']
                    if validator_yaml[k]['type'] == 'number' or validator_yaml[k]['type'] =='integer':
                        if 'min' in validator_yaml[k].keys():
                            min_value = validator_yaml[k]['min']
                        if 'max' in validator_yaml[k].keys():
                            max_value = validator_yaml[k]['max']
                        if 'min_exclusive' in validator_yaml[k].keys():
                            min_value_excl = validator_yaml[k]['min']
                        if 'max_exclusive' in validator_yaml[k].keys():
                            max_value_excl = validator_yaml[k]['max']

                #alert user of issues
                for u in uniq.astype(str):
                    if not u.isnumeric():
                        if u not in allowed_list and len(allowed_list) > 0:
                            msg = msg + "Warning " + u + " found in column " + k + " but not allowed per Qiimp template." +\
                            "valid values: " + str(allowed_list) + "\n"          
                    else:
                        if u not in allowed_list: #assume it's actually a number
                            if min_value != '' and u < min_value:
                                msg = msg + "Warning " + u + " found in column " + k + " but less than min value per yaml: " +\
                                str(min_value) + "\n"
                            if max_value != '' and u > max_value:
                                msg = msg + "Warning " + u + " found in column " + k + " but more than max value per yaml: " +\
                                     str(max_value) + "\n"
                            if min_value_excl != '' and u <= min_value_excl:
                                msg = msg + "Warning " + u + " found in column " + k + " but less than min value per yaml: " +\
                                     str(min_value_excl) + "\n"
                            if max_value_excl != '' and u >= max_value_excl:    
                                lmsg = msg + "Warning " + u + " found in column " + k + " but not allowed per yaml: " +\
                                str(max_value_excl) + "\n"
        st_list.append(df_to_validate)
        if len(msg) > 0:
            logger.warning("Errors found during validation:")
            logger.warning(msg)
            if DEBUG:
                st_log_name=st.replace(' ','_')
                valid_log_filename = prefix + "_" + st_log_name + '_validation_errors.log'
                errors = open(valid_log_filename, "w")
                n = errors.write(msg)
                errors.close()
                logger.warning("Validation errors written to " + valid_log_filename)
    valid_df = pd.concat(st_list)
    return valid_df

### END: NORMALIZATION METHODS

### START: EBI/ENA AND NCBI/SRA DATA RETRIEVAL METHODS
def get_study_details(study_accession,mode='ebi',prefix=''):
    """Retrieve study information including list of samples

    #TODO fill in description

    Parameters
    ----------
    study_accession:string
        project or study name from EBI/ENA or NCBI/SRA
    mode : string
        'ebi' or 'sra'; repo to use for data retrieval
    prefix: string
        prefix to assign to validator error log, includes path
    
    Returns
    -------
    study_df: pd.DataFrame
       dataframe of study information as provided by the chosen repository
    """
    if mode == 'ebi':
        studyUrl = "http://www.ebi.ac.uk/ena/data/view/" + study_accession \
                    + "&display=xml"
        try:
            response = requests.get(studyUrl)
            xml_dict=parse(response.content)
        except:
            logger.error("Could not obtain study information for " + study_accession + " at URL: " + studyUrl
                             + " . Please check connection and study/project name and try again.")

        if 'STUDY_SET' in xml_dict.keys():
            logger.info(study_accession + " is study ID. Writing config file")
            id_tuple = write_config_file(xml_dict,prefix,mode)
        elif 'PROJECT_SET' in xml_dict.keys():
            try:
                secondary_accession = xml_dict['PROJECT_SET']['PROJECT']['IDENTIFIERS']['SECONDARY_ID']
                secondaryUrl = "http://www.ebi.ac.uk/ena/data/view/" + secondary_accession + "&display=xml"
                logger.warning(study_accession + " is project ID. Retrieved secondary ID: " + secondary_accession + " Writing config file.")
                logger.info("Alternate url" + secondaryUrl)
                response2 = requests.get(secondaryUrl)
                xml_study_dict=parse(response2.content)
                id_tuple = write_config_file(xml_study_dict,prefix,mode,xml_dict)

            except:
                logger.error("No matching study ID found for project " + study_accession)
        else:
            logger.warning("No study information found for " + study_accession)

        host = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
        read_type = "&result=read_run&"
        fields = "library_name,secondary_sample_accession,run_accession," + \
                 "experiment_accession,fastq_ftp,library_source," + \
                 "instrument_platform,submitted_format,library_strategy," +\
                 "library_layout,tax_id,scientific_name,instrument_model," + \
                "library_selection,center_name,experiment_title," +\
                "study_title,study_alias,experiment_alias,sample_alias,sample_title"
        url = ''.join([host, study_accession, read_type, "fields=", fields])
        if DEBUG: logger.info(url)
        study_df = pd.read_csv(url,sep='\t')
        study_df.dropna(axis=1,how='all',inplace=True)

        #add primary (project) and secondary (study) ID to dataframe
        study_df['primary_ebi_id']=id_tuple[0]
        study_df['secondary_ebi_id']=id_tuple[1]
        study_df['ebi_import']='TRUE'

    elif mode == 'sra':
        p1 = subprocess.Popen(['esearch', '-db', 'sra', '-query', study_accession],
                   stdout=subprocess.PIPE)
        for i in p1.stdout:
            if "<Count>0</Count>" in i.decode("utf-8"):
                p1.stdout.close()
                raise Exception(study_accession + " is not a valid SRA study ID")

        p1 = subprocess.Popen(['esearch', '-db', 'sra', '-query', study_accession],
                           stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['efetch', '-format', 'native'], stdin=p1.stdout,
                           stdout=subprocess.PIPE)
        p1.stdout.close()

        for i in p2.stdout:
            line = i.decode("utf-8")
            start = line.find("<STUDY ")
            if(start != -1):
                end = line.find("</STUDY>")
                info = line[start: end+8]
                doc = parse(info)
                break
        p2.stdout.close()
        write_config_file(doc,study_accession,mode)
        
        p1 = subprocess.Popen(['esearch', '-db', 'sra', '-query', study_accession],
                          stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['efetch', '-format', 'runinfo'], stdin=p1.stdout,
                          stdout=subprocess.PIPE)
        count=0
        headers=[]
        for i in p2.stdout:
            if count == 0:
                sra_headers = i.decode("utf-8").replace('\n','').split(',')
                for h in sra_headers:
                    headers.append(scrub_special_chars(split_caps(h).lower(),custom_dict={'i_d':'id',
                                                                                    's_r_a':'sra',
                                                                                    'e_b_i':'ebi',
                                                                                    's_r_r':'srr',
                                                                                    'm_b':'mb'}))
                    study_df=pd.DataFrame({},columns=headers)
            else:
                tmp=str(i.decode("utf-8"))

                if len(tmp) > 1:
                    a_series = pd.Series(list(csv.reader(tmp.splitlines(),delimiter=','))[0], index = study_df.columns)
                    study_df = study_df.append(a_series, ignore_index=True)
            count += 1
    else:
        raise Exception(mode + " is not a valid repository.")

    return study_df


#now feed in study_df above
def get_sample_info(input_df,mode='ebi',plat=[],strat=[],validator_files={},prefix='',names=[],
                    src=[],empo_mapping=pd.DataFrame(),sample_type_column='scientific_name'):
    """Retrieve information for each sample including experiment/prep information

    #TODO fill in description

    Parameters
    ----------
    input_df : pd.DataFrame
        dataframe of samples for retreival
    mode : string
        'ebi' or 'sra'; repo to use for data retrieval
    plat :  list
        list of valid plaforms to filter for
    strat: list
        list of valid library strategies to filter for
    validator_files: dict
        dictionary of available validation files for metadata normalization
    prefix: string
        prefix to assign to validator error log, includes path
    names : list
        list of scientific names to filter for
    src : list
        list of library sources to filter for
    empo_mapping : pd.DataFrame()
        datagframe of loaded EMPO and/or sample_type mappings
    sample_type_column : string
        column for identifying sample type 

  
    Returns
    -------
    study_df: pd.DataFrame
       dataframe of study information as provided by the chosen repository
    """
    
    #need to provide way to list additional prep_info columns, now function will return tuple
    add_prep_cols=[]

    #check to see if the .part file already exists:
    pre_validation_file = prefix + "_unvalidated_sample_info.part"
    if not path.isfile(pre_validation_file):
        if mode =='ebi':
            identifier = 'secondary_sample_accession'
            run_accession = 'run_accession'
            input_df['platform']=input_df['instrument_platform']
        elif mode == 'sra':
            identifier = 'sample'
            run_accession = 'run'
            input_df['instrument_model']=input_df['model']
        else:
            raise Exception(mode + " is not a valid repository.")

        #Note: for now this loop just uses the data in EBI since it is mirrored with NCBI
        sample_info_list=[]
        sample_count_dict = {}
        prep_df_dict={}

        #apply filters for platforms and strategies
        except_msg=''
        if len(plat) > 0:
            except_msg = except_msg + "Selected Platforms: " + str(plat) + "\n Available Platforms:" +\
                        str(input_df['platform'].str.lower().unique()) + "\n"
            input_df = input_df[input_df['platform'].str.lower().isin(plat)]

        if len(strat) > 0:
            except_msg = except_msg + "Selected Strategies: " + str(strat) + "\n Available Strategies:" +\
                        str(input_df['library_strategy'].str.lower().unique()) + "\n"
            input_df = input_df[input_df['library_strategy'].str.lower().isin(strat)]

        if len(names) > 0:
            except_msg = except_msg + "Selected Scientific Names: " + str(names) + "\n Available Scientific Names:" +\
                        str(input_df['scientific_name'].str.lower().unique()) + "\n"
            input_df = input_df[input_df['scientific_name'].str.lower().isin(names)]

        if len(src) > 0:
            except_msg = except_msg + "Selected Library Sources: " + str(src) + "\n Available Library Sources:" +\
                        str(input_df['library_source'].str.lower().unique()) + "\n"
            input_df = input_df[input_df['library_source'].str.lower().isin(src)]

        if len(input_df) == 0:
            raise Exception("No files after selection criteria:\n" + except_msg)

        #need to catch common issue where identifier is identical for unique samples
        #for now assume that in this case, library_name will be unique, and if it isn't combine sample and run names
        if len(input_df) > 1 and input_df[identifier].nunique() == 1:
            if input_df['library_name'].nunique() != 1:
                 #print(str(len(_input_df)) + " is length and input_df[identifier].nunique() = " + str(input_df[identifier].nunique()))
                 input_df['sample_name']=input_df['library_name'].apply(lambda x: scrub_special_chars(x).replace('_','.'))
            else:
                 input_df['sample_name']=input_df[identifier]+ '.' + input_df[run_accession]
        else:
            input_df['sample_name']=input_df[identifier]

        input_df['run_prefix']=input_df[run_accession]

        for index, row in input_df.iterrows():
            sample_accession = row[identifier]
#            ebi_lib_strat = row['library_strategy']
            sample_name = row['sample_name']
            expt_accession = row['experiment_accession']

            sampleUrl = "http://www.ebi.ac.uk/ena/data/view/" + sample_accession \
                      + "&display=xml"
            if DEBUG: logger.info(sampleUrl)

            sample_response = requests.get(sampleUrl)
            sample_xml_dict=parse(sample_response.content)

            exptUrl = "http://www.ebi.ac.uk/ena/data/view/" + expt_accession \
                      + "&display=xml"
            if DEBUG: logger.info(exptUrl) 
            expt_response = requests.get(exptUrl)
            expt_xml_dict=parse(expt_response.content)


            if 'EXPERIMENT_SET' not in expt_xml_dict.keys():
                logger.warning('No experimental metadata found for sample named: ' + sample_accession + ' with experiment_accession ' + expt_accession + ', omitting.')
            else:
                if 'SAMPLE_SET' not in sample_xml_dict.keys():
                    logger.warning('No metadata found for sample named: ' + sample_accession + ' omitting.')
                else:
                    input_df.at[index,'experiment_title_specific']=expt_xml_dict['EXPERIMENT_SET']['EXPERIMENT']['TITLE']
                    try:
                        ea=expt_xml_dict['EXPERIMENT_SET']['EXPERIMENT']['EXPERIMENT_ATTRIBUTES']['EXPERIMENT_ATTRIBUTE']
                        for e in ea:
                            col = scrub_special_chars(e['TAG']).lower()
                            #print(col)
                            if col not in add_prep_cols:
                                add_prep_cols.append(col)
                            try:
                                input_df.at[index,col]=e['VALUE']
                            except:
                                logger.warning('No value found for experiment attribute: ' + col + '. Setting to "not provided".') 
                    except:
                        logger.warning("No experiment attributes found for " + expt_accession + " corresponding to sample " + sample_accession)
                    #need to get prep_type and convert based on dictionary and additional information, e.g.target_gene for amplicon
                    prep_type = set_prep_type(input_df, index)

                    if prep_type not in sample_count_dict.keys() :
                        sample_count_dict[prep_type]= {sample_name:0}
                    elif sample_name not in sample_count_dict[prep_type].keys():
                        sample_count_dict[prep_type][sample_name]= 0
                    else:
                        sample_count_dict[prep_type][sample_name] = sample_count_dict[prep_type][sample_name] + 1

                    #now add sample specific information
                    input_df.at[index,'sample_title_specific']=sample_xml_dict['SAMPLE_SET']['SAMPLE']['TITLE']

                    sn=sample_xml_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_NAME']
                    for s in sn.keys():
                        col = scrub_special_chars(s).lower()
                        #print(col)
                        input_df.at[index,col]=sn[s]

                    sa=sample_xml_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']
                    for s in sa:
                        col = scrub_special_chars(s['TAG']).lower()
                        #print(col)
                        try:
                            input_df.a[index,col]=s['VALUE']
                        except:
                            logger.warning('No value found for sample attribute: ' + col + '. Setting to "not provided".')
                        input_df.at[index,col]='not provided'
                    #this sets the key used for splitting the files into prep_info templates
                    input_df.at[index,'prep_file']=prep_type + '_' + str(sample_count_dict[prep_type][sample_name])


        #the loops above takes the most time so write out a placeholder file for faster re-running if interrupted
        input_df.to_csv(pre_validation_file,sep='\t',index=False)
    else:
        input_df = pd.read_csv(pre_validation_file,sep='\t',dtype=str)

    #start by normalizing sample types and EMPO fields if sample_type is present
    if len(empo_mapping) > 0:
        if 'sample_type' in input_df.columns:
            input_df = normalize_types(input_df,empo_mapping)
        else:
            logger.warning("'sample_type' not found in metadata. Skipping EMPO normalization.")

    output_df=validate_samples(input_df,sample_type_column,validator_files,prefix)

    #tidy output before returning
    output_df.columns = [scrub_special_chars(col).lower() for col in output_df.columns]

    return output_df, add_prep_cols


def fetch_sequencing_data(download_df,output_dir="./",mode='ebi',host_deplete=False,host_db='/databases/minimap2/human-phix-db.mmi',cpus=4,hd_method='minimap2',qc=True,quality=False):
    """Fetch fastq files and run quality or host depletion per file (pair) as requested
    
    #TODO: complete description
    
    Parameters
    ----------
    download_df: pd.DataFrame
        dataframe with sample and path information for downloading
    output_dir : string
        the path to save the output files
    mode : string
        'ebi' or 'sra'; repo to use for data retrieval
    host_deplete: boolean
        whether to run host depletion on the files as they are downloaded
    host_db : string
        path to the databased to be used for filtering
    cpus : int
        number of threads to use for fastp and fastq
    hd_method : string
        method for host depletion
    qc : boolean
        whether to run fastqc on output
    quality: boolean
        whether to run fastp on downloaded files

    Returns
    -------
    valid_df : pd.DataFrame
        dataframe of samples that were successfully downloaded
    
    """

    logger.info("Downloading the fastqs")
    # Download the fastqs
    
    for index, row in download_df.iterrows():
        sequencer = row['instrument_model']
        failed_list=[]
        hd_file_list =[]
        if mode == 'ebi':
            try:
                files = row['fastq_ftp'].split(';')
            except:
                try:
                    files = row['download_path'].split(';')
                except:
                    logger.warning("Skipping sample:" + row['sample_name']
                                    + ", run: " + row['run_accession'] + "No fastq ftp found.")
            for f in files:
                fq_path = output_dir + "/" + f.split('/')[-1]
                #hd_file_list.append(fq_path)
                if DEBUG: logger.info(f)
                if type(f) != str:
                    logger.warning("Skipping sample:" + row['sample_name']
                                   + ", run: " + row['run_accession'] + "fastq ftp path is not string.")
                elif path.isfile(fq_path):
                    logger.warning("Skipping download of " + fq_path + " File exists.")

                elif path.isfile(fq_path.replace('.fastq.gz','.R1.fastp.fastq.gz')): #if the quality filtered file exists, also skip
                    logger.warning("Skipping download of " + fq_path.replace('.fastq.gz','.R1.fastp.fastq.gz') + " fastp file exists.")
                elif path.isfile(fq_path.replace('.fastq.gz','.R1.filtered.fastq.gz')) or path.isfile(fq_path.replace('.fastq.gz','.R2.filtered.fastq.gz')): #if the host depleted file exists, also skip
                    logger.warning("Skipping download of " + fq_path.replace('.fastq.gz','.R1.filtered.fastq.gz') + " Host depleted file exists.")
                elif row['sample_name'] in failed_list: #this should only happen if read 1 failes
                    logger.warning("Skipping download of read 2 for " + row['sample_name'] + " . Other read(s) failed to download.")
                else:
                #add catch in case there is an issue with the connection
                    try:
                        urlretrieve("ftp://" +f, fq_path)
                        hd_file_list.append(fq_path)
                    except:
                        logger.warning("Issue with urlretrieve for " + row['sample_name'] + "Skipping.")
                        failed_list.append(row['sample_name'])
        elif mode =='sra':
            subprocess.run(['fastq-dump', '-I', '--split-files', '--gzip', '--outdir', output_dir, row['run_accession']])
            sra_files = glob.glob(output_dir + row['run_accession'] + '*.gz')
            for f in sra_files:
                hd_file_list.append(f)

        else:
            raise Exception(mode + " is not a valid repository")

        if len(hd_file_list) > 0: #make sure one of the above worked
            #run fastqc
            if qc:
                logger.info('Running fastqc')
                run_qc(hd_file_list,output_dir,cpus)

            #run quality check if requested
            if quality and host_deplete:
                logger.warning("Quality check run as part of host depletion.")
            elif quality:
                run_quality_check(hd_file_list,cpus,output_dir,sequencer)

            #enter host depletion checks
            valid_hd_prep_types = ['WGS','WGA','WXS','ChIP-Seq','RNA-Seq']

            if host_deplete:
                if len(hd_file_list) > 2:
                    logger.warning("More than 2 files in the fastq download path. Skipping host depletion for "
                         + row['sample_name'])
                else:
                    if row['library_strategy'] in valid_hd_prep_types:
                        logger.info("Running host depletion on " + row['sample_name'])
                        run_host_depletion(hd_file_list,host_db,cpus,output_dir,hd_method,qc,sequencer)
                    elif force_hd:
                        logger.warning("Forcing host depletion for " + + " with library_strategy " + + ". This may break...")
                        run_host_depletion(hd_file_list,host_db,cpus,output_dir,hd_method,qc,sequencer)
                    else:
                        logger.warning("Skipping host depletion for " + row['sample_name'] + " Library Strategy is "
                                     + row['library_strategy'] + ". Valid formats are: "+ str(valid_hd_prep_types))

    #drop samples from dataframe if the download fails
    valid_df =download_df[~download_df['sample_name'].isin(failed_list)]

    return valid_df

### END: EBI/ENA AND NCBI/SRA DATA RETRIEVAL METHODS


### START: OUTPUT FILE WRITING METHODS

def write_config_file(xml_dict,prefix,mode,id_type='study',xml_proj_dict={}):
    """Writes the study_config and study_title files needed for Qiita loading

    #TODO write description
    
    Parameters
    ----------
    xml_dict: dict
        distionary of study information for parsing
    prefix: string
        prefix to assign to validator error log, includes path
    mode : string
        'ebi' or 'sra'; repo to use for data retrieval
    
    id_type : string
        type of id supplied, either study or project
    xml_proj_dict: dict
        dictionary of project information for parsing

    Returns
    -------
    empo_type_df : tuple
        returns project and study id
    """

    config_string='[required]\ntimeseries_type_id = 1\nmetadata_complete = True\nmixs_compliant = True'

    #for now at least, setting the PI to default to Qiita-Help
    config_string = config_string + '\nprincipal_investigator = Qiita-EBI Import, qiita.help@gmail.com, See study details'

    config_string = config_string + '\nreprocess = False'
    if mode == 'ebi':
        parse_dict = xml_dict['STUDY_SET']['STUDY']
    elif mode == 'sra':
        parse_dict = xml_dict['STUDY']
    else:
        logger.warning('Received invalid mode: ' + mode)

    title= 'XXEBIXX'
    alias= '\nstudy_alias = XXEBIXX'
    abstract ='\nstudy_abstract = XXEBIXX'
    description = '\nstudy_description = XXEBIXX'
    proj_id =parse_dict['IDENTIFIERS']['SECONDARY_ID']
    study_id =parse_dict['IDENTIFIERS']['PRIMARY_ID']


    if '@alias' in parse_dict.keys():
        alias=alias.replace('XXEBIXX',parse_dict['@alias'])
    elif '@alias' in xml_proj_dict.keys():
        alias=alias.replace('XXEBIXX',_dict['PROJECT_SET']['PROJECT']['@alias'])
    else:
        logger.warning("No alias found, using XXEBIXX for alias.")

    desc_dict={}
    if 'DESCRIPTOR' in parse_dict.keys():
        desc_dict = parse_dict['DESCRIPTOR']
    elif 'IDENTIFIERS' in xml_proj_dict.keys():
        desc_dict = xml_proj_dict['IDENTIFIERS']
    else:
        logger.warning("No DESCRIPTOR or IDENTIFIER values found. Using XXEBIXX for values.")

    if len(desc_dict) > 0:
        if 'STUDY_ABSTRACT' in desc_dict.keys():
            abstract=abstract.replace('XXEBIXX',desc_dict['STUDY_ABSTRACT'])
        elif 'ABSTRACT' in desc_dict.keys():
            abstract=abstract.replace('XXEBIXX',desc_dict['ABSTRACT'])
        else:
            logger.warning("No abstract found, using XXEBIXX for abstract")

        if 'STUDY_DESCRIPTION' in desc_dict.keys():
            abstract=abstract.replace('XXEBIXX',desc_dict['STUDY_DESCRIPTION'])
        elif 'DESCRIPTION' in desc_dict.keys():
            abstract=abstract.replace('XXEBIXX',desc_dict['DESCRIPTION'])
        else:
            logger.warning("No description found, using XXEBIXX for description")

        if 'STUDY_TITLE' in desc_dict.keys():
            title=title.replace('XXEBIXX',desc_dict['STUDY_TITLE'])
        elif 'TITLE' in desc_dict.keys():
            title=title.replace('XXEBIXX',desc_dict['TITLE'])
        else:
            logger.warning("No title found, using XXEBIXX for title")

    config_string= config_string + alias + description + abstract + '\nefo_ids = 1\n[optional]' #To add? + '\nproject_id = ' + proj_id + '\nstudy_id = ' + study_id

    study_config_file = prefix + '_study_config.txt'
    study_title_file = prefix + '_study_title.txt'

    # Write out files
    c_file = open(study_config_file, "w")
    c_file.write(config_string)
    c_file.close()

    t_file = open(study_title_file, "w")
    t_file.write(title)
    t_file.close()

    return proj_id, study_id

def create_details_file(study_details_df, study_accession,mode='ebi',prefix='',file_suffix="_detail"):
    """Returns the details of the EBI/SRA study

    If the accession ID is valid, generate a .details.txt, and return the
    detail file name of this EBI/SRA study. Else return None
    
    #TODO: consider dropping since no longer used?

    Parameters
    ----------
    study_details_df : pd.Dataframe
        dataframe of study information to write out
    study_accession:string
        project or study ID from EBI/ENA or NCBI/SRA
    mode : string
        'ebi' or 'sra'; repo to use for data retrieval
    prefix: string
        prefix to assign to validator error log, includes path
    file_suffix : string
        The suffix for the output study detail file

    Returns
    -------
    string
        study details file name
    """
    if len(prefix)==0:
        prefix = study_accession
    study_details = prefix + "_" + mode + file_suffix + ".txt"
    study_details_df.to_csv(study_details,sep='\t',header=True,index=False)
    return study_details
    
def write_info_files(final_tuple,max_prep,prefix=''):
    """Writes out the prep and sample information files

    This is the key method that creates the sample_information
    and preparation information files needed for Qiita.
    
    #TODO: complete description
    
    #TODO: consider dropping since no longer used?

    Parameters
    ----------
    final_tuple : tuple
        tuple of metadata (0) and prep_info columns (1)
    max_prep: int
        max number of samples to write into any prep info file
    prefix: string
        prefix to assign to validator error log, includes path

    Returns
    -------
    None
    
    """
    
    final_df = final_tuple[0]
    if max_prep > len(final_df):
        max_prep = len(final_df)

    prep_info_columns = ['run_prefix','experiment_accession','platform','instrument_model','library_strategy',
                         'library_source','library_layout','library_selection','fastq_ftp','ena_checklist',
                         'ena_spot_count','ena_base_count','ena_first_public','ena_last_update','instrument_platform',
                         'submitted_format','sequencing_method','target_gene','target_subfragment','primer']

    #add multiqc_columns
    multiqc_cols=['raw_reads','non_human_reads','post_qc_reads','frac_non_human_from_raw','frac_non_human_from_qc']
    prep_info_columns = prep_info_columns + multiqc_cols
    #add ebi/user supplied prep info columns
    ebi_prep_info_cols= final_tuple[1]
    for c in ebi_prep_info_cols:
        if c not in prep_info_columns:
            prep_info_columns.append(c)

    amplicon_min_prep_list=['target_gene','target_subfragment','primer']
    amplicon_type_preps = ['16S','ITS','18S']
    final_df.columns =[scrub_special_chars(col).lower() for col in final_df.columns]

    #write sample_info
    sample_df=final_df[final_df.columns[~final_df.columns.isin(prep_info_columns)]]

    #check for duplicates here. N.B. Need to retain previously to enable download of all runs without faffing around with prep files
    sample_df=sample_df.drop_duplicates('sample_name')
    sample_df=sample_df.set_index('sample_name').dropna(axis=1,how='all')
    sample_df.to_csv(prefix+'_sample_info.tsv',sep='\t',index=True,index_label='sample_name')

    #clean up pre-validation file assuming correct validation
    if not DEBUG:
        pre_validation_file = prefix + "_unvalidated_sample_info.part"
        if path.isfile(pre_validation_file):
            remove(pre_validation_file)

    prep_info_columns = ['sample_name'] + prep_info_columns #add to list for writing out prep files
    for prep_file in final_df['prep_file'].unique():

        #adding way to check for min essential target gene information where needed
        logger.info(prep_file.split('_')[0])
        prep_df = final_df[final_df['prep_file']==prep_file]
        if prep_file.split('_')[0] in amplicon_type_preps: #check to see if the prep is amplicon-style, specified by list above
            for min_prep in amplicon_min_prep_list: #if amplicon-style, enforce presence or null values for minimum prep info information
                if min_prep not in prep_df.columns:
                    prep_df[min_prep]='XXEBIXX' #will throw warning, but okay with current pandas

        #now write out the prep info files
        prep_df= prep_df[prep_df.columns[prep_df.columns.isin(prep_info_columns)]].set_index('sample_name')
        prep_df=prep_df.dropna(axis=1,how='all')
        prep_df_list = [prep_df[i:i+max_prep] for i in range(0,prep_df.shape[0],max_prep)]
        prep_count=0
        for prep in prep_df_list:
            prep.to_csv(prefix+'_prep_info_'+ prep_file + '_part' + str(prep_count) +'.tsv',sep='\t',index=True,index_label='sample_name')
            prep_count += 1
            
### END: OUTPUT FILE WRITING METHODS

### START: FASTQ FILE PROCESSING METHODS

def run_quality_check(raw_list,cpus=4,output_dir='./',model='',method='fastp',min_length=45,qc=True,keep=False):
    """Runs fastp quality filtering and optionally fastqc
    
    #TODO: complete description, consider dropping fastqc if fastp alternative can be used
    

    Parameters
    ----------
    raw_list: list
        list of files to quality filter
    cpus : int
        number of threads to use for fastp and fastq
    output_dir : string
        the path to save the output files
    model : string
       the model of instrument used for sequencing
    method : string
        method for qc, only fastp supported currently
    min_length : int
        minimum sequencing length for quality filtering; set to 45 to permit RNA-Seq data
    qc : boolean
        whether to run fastqc on output
    keep : boolean
        whether to keep the raw files after filtering

    Returns
    -------
    return_list : list
        list of quality-filtered files
    
    """
    valid_methods = ['fastp']
    polyG_model_list = ['Illumina MiniSeq', 'Illumina NovaSeq 6000', 'NextSeq 500', 'NextSeq 550']

    read1_fastq = raw_list[0]

    if not method in valid_methods:
        logger.warning("Method " + method + " not supported. Available quality check methods: "
                       + valid_methods + ". Returning raw files: " + raw_list)
        return raw_list

    elif method == 'fastp':
        qc_fastq1 = read1_fastq.replace('.fastq.gz','.R1.fastp.fastq.gz')
        qc_fastq2 = '' #dummy file to help with cleanup
        return_list=[qc_fastq1]

        #set fastp parameters
        fastp_args=['fastp','-l',str(min_length),'-i', read1_fastq,'-o',qc_fastq1]
        #to determine %microbial reads, need to write out and fastqc files. will assume this is default when host depleting

        if len(raw_list)==2:
            read2_fastq = raw_list[1]
            qc_fastq2 = read2_fastq.replace('.fastq.gz','.R2.fastp.fastq.gz')
            fastp_args = fastp_args + ['-I',read2_fastq,'-O',qc_fastq2]
            return_list.append(qc_fastq2)

        #add final flags for fastp after determining input file count
        #fastp_html =read1_fastq.replace('_1','')+'_fastp.html' #TODO evaluate if fastp can replace fastqc?
        fastp_args = fastp_args + ['-w',str(cpus)] #,'-h',fastp_html]
        if model in polyG_model_list:
            fastp_args + fastp_args + ['-g','--poly_g_min_len'] #polyG filtering, 10 is default

        #see if file already exists
        if path.isfile(qc_fastq1): #process should make both, so just check
            logger.warning("Output fastp file already exists for " + read1_fastq + " Skipping.")

        else:
            fastp_ps = subprocess.Popen(fastp_args)
            fastp_ps.wait()
            logger.info("past fastp")

        #clean up raw file if qc done
        if path.isfile(qc_fastq1):
            if qc: #very likely
                result=run_qc([qc_fastq1],output_dir,cpus)
            if not keep:
                subprocess.run(['rm', read1_fastq])
        if path.isfile(qc_fastq2):
            if qc: #very likely
                result=run_qc([qc_fastq2],output_dir,cpus)
            if not keep:
                subprocess.run(['rm', read2_fastq])

        if DEBUG: logger.info("Return list is: " + str(return_list))
        return return_list

def run_host_depletion(fastq_file_list,filter_db='',cpus=4,output_dir='./',method='minimap2',qc=True,model='',keep=False,min_length=45):
    """Quality filter and host deplete list of files
    
    #TODO: complete description
    
    Parameters
    ----------
    fastq_file_list: list
        list of raw files to quality filter and host deplete
    filter_db : string
        path to the databased to be used for filtering
    cpus : int
        number of threads to use for fastp and fastq
    output_dir : string
        the path to save the output files
    method : string
        method for host depletion
    model : string
       the model of instrument used for sequencing
    qc : boolean
        whether to run fastqc on output
    keep : boolean
        whether to keep the raw files after quality filtering, and fastp files after host filtering    
    min_length : int
        minimum sequencing length for quality filtering; set to 45 to permit RNA-Seq data

    Returns
    -------
    None
    
    """
    
    db_dict = {'bowtie2':'/databases/bowtie/Human_phiX174/Human_phix174',
                'minimap2':'/databases/minimap2/human-phix-db.mmi'
              }
    if filter_db == '':
        filter_db = db_dict[method]

    read1_fastq=str(fastq_file_list[0])
    filtered_fastq_1=read1_fastq.replace('.fastq.gz','.R1.filtered.fastq')
    output_fastq1=filtered_fastq_1 + '.gz'

    #new helper list for qc checks
    to_qc=[read1_fastq]

    if len(fastq_file_list) == 2:
        read2_fastq=str(fastq_file_list[1])
        filtered_fastq_2=read2_fastq.replace('.fastq.gz','.R2.filtered.fastq')
        output_fastq2 =filtered_fastq_2 + '.gz'
        to_qc.append(read2_fastq)
    else:
        read2_fastq = '' #to reduce code complexity for checking for paired vs unpaired, dummy file string here
        output_fastq2 = '' #dummy file to help with cleanup 


    #run some checks to help wtih resumed runs
    skip = False
    skip = path.isfile(output_fastq1)
    if output_fastq2 != '':
        skip = path.isfile(output_fastq2)

    if skip:
        logger.warning("All filtered files found for " + read1_fastq.replace('.fastq.gz','').split('/')[-1] + ". Skipping host depletion.")

    elif method == 'bowtie2':
        if DEBUG: logger.info("Starting bowtie2 depletion")
        bowtie2_args=['bowtie2', '-p', str(cpus), '-x', filter_db]
        stv_args_1 =['samtools', 'view','-F', '256']
        sts_args =['samtools', 'sort', '-@', str(cpus),'-n']
        stv_args_2 = ['samtools', 'view','-bS']
        btb_args = ['bedtools', 'bamtofastq', '-i', '-', '-fq', filtered_fastq_1]

        if len(fastq_file_list) == 2:
            bowtie2_args=bowtie2_args + ['-1',read1_fastq,'-2',read2_fastq]
            stv_args_1 = stv_args_1 + ['-f','12']
            btb_args = btb_args + ['-fq2',filtered_fastq_2]
        else:
            bowtie2_args=bowtie2_args + ['-U',read1_fastq]
            stv_args_1 = stv_args_1 + ['-f','4']

        bowtie2_args.append('--fast-local')

        #now run bowtie2 commands in chain
        bt2_ps = subprocess.Popen(bowtie2_args, stdout=subprocess.PIPE)
        stv_ps1 = subprocess.Popen(stv_args_1,stdin=bt2_ps.stdout, stdout=subprocess.PIPE)
        sts_ps = subprocess.Popen(sts_args,stdin=stv_ps1.stdout, stdout=subprocess.PIPE)
        stv_ps2 = subprocess.Popen(stv_args_2,stdin=sts_ps.stdout, stdout=subprocess.PIPE)
        btb_ps = subprocess.Popen(btb_args,stdin=stv_ps2.stdout, stdout=subprocess.PIPE)
        btb_ps.wait()

        if path.isfile(filtered_fastq_1):
            subprocess.run(['gzip', filtered_fastq_1])
            run_qc([output_fastq1],output_dir,cpus)
            if not keep:
                subprocess.run(['rm', read1_fastq])
        if path.isfile(filtered_fastq_2):
            subprocess.run(['gzip', filtered_fastq_2])
            run_qc([output_fastq2],output_dir,cpus)
            if not keep:
                subprocess.run(['rm', read2_fastq])

    elif method == 'minimap2':

        #separating out filtering for flexibility and % microbial read determination
        qced_fastq=run_quality_check(to_qc,cpus,output_dir,model)
        if DEBUG: logger.info("Starting minimap2 depletion.")
        minimap2_args = ['minimap2','-ax', 'sr', '-t', str(cpus), filter_db,'-a'] + qced_fastq
        stf_args=['samtools','fastq', '-@', str(cpus),'-F', '256', '-']

        if len(fastq_file_list)==2:
            stf_args= stf_args + ['-f', '12','-1',output_fastq1,'-2',output_fastq2]
        else: #different samtools parameter for unpaired
            stf_args= stf_args + ['-f', '4','-0',output_fastq1]

        minimap2_ps = subprocess.Popen(minimap2_args, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stf_ps = subprocess.Popen(stf_args,stdin=minimap2_ps.stdout, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stf_ps.wait()

        #remove source files if not needed
        if path.isfile(output_fastq1):
            run_qc([output_fastq1],output_dir,cpus)
            if not keep:
                #subprocess.run(['rm', read1_fastq]) #moved to post quality
                subprocess.run(['rm', qced_fastq[0]])

        if path.isfile(output_fastq2):
            run_qc([output_fastq2],output_dir,cpus)
            if not keep:
                #subprocess.run(['rm', read2_fastq]) #moved to post-quality
                subprocess.run(['rm', qced_fastq[1]])
    else:
        logger.warning("Selected depletion method '" + method + "' not currently supported. Please select either bowtie2 or minimap2.")

def run_qc(input_list,output_dir='./',cpus=4,multiqc=False):
    """Run fastqc and/or multiqc list of files
    
    #TODO: complete description
    
    Parameters
    ----------
    input_list: list
        list of files to qc
    output_dir : string
        the path to save the output files
    cpus : int
        number of threads to use for fastp and fastq

    multiqc : boolean
        whether to run multiqc on the output directory

    Returns
    -------
    'qc_done': string
        indicator that qc completed    
    
    """
    #check to see if fastqc has been run on the files in output, and if not, run it
    output_dir = output_dir + '/fastqc/'
    if not path.exists(output_dir):
        makedirs(output_dir)
    if not len(input_list) == 0: #bypass fastqc run by passing empty list
        fastqc_args = ['fastqc']

        for f in input_list:
            fastqc_name = output_dir+f.split('/')[-1].replace('.fastq.gz','_fastqc.html')

            if not path.isfile(fastqc_name):
                fastqc_args = fastqc_args + [f]

        #add a check to ensure there is something to run
        if len(fastqc_args) != 1:
            fastqc_args = fastqc_args + ['-t',str(cpus),'-o',output_dir]
            if not DEBUG:
                fastqc_args = fastqc_args + ['-q']

            fastqc_ps = subprocess.Popen(fastqc_args)
            fastqc_ps.wait()

    if multiqc:
        multiqc_args = ['multiqc',output_dir,'-o',output_dir,'-s','-f']
        subprocess.run(multiqc_args)

    return 'qc_done'

def add_qc_data(partial_df,output_dir='./'):
    """Automatically augment prep information with read data and % microbial reads
    
    #TODO: complete description
    
    Parameters
    ----------
    partial_df: pd.DataFrame
        dataframe of combined sample and prep information for augmenting
    
    output_dir : string
        the path to save the output files

    Returns
    -------
    partial_df: pd.DataFrame
        augmented dataframe with read information
    
    """
    
    #read in multiqc data
    multiqc_data_file = output_dir+'/multiqc_data/multiqc_fastqc.txt'
    try:
        mqc = pd.read_csv(multiqc_data_file,sep='\t')
    except:
        logger.warning("Multiqc data file " + multiqc_data_file + " missing. Please check path or re-run with -qc flag.")

    #add catch to skip with warning if the above fails
    if len(mqc) > 0:
        raw_files = mqc[~mqc['Sample'].str.contains('.fastp')] #drops .fastp
        raw_files = raw_files[~raw_files['Sample'].str.contains('.filtered')] #drops filtered
        fastp_files = mqc[mqc['Sample'].str.contains('.fastp')] #drops .fastp and .filtered
        filt_files = mqc[mqc['Sample'].str.contains('.filtered')] #drops .fastp and .filtered

        #to date read 1 length = read 2 length, so assume this for now. can update to sum with groupby and compare later
        raw_files['run_prefix'] = raw_files['Sample'].apply(lambda x: x.split('.')[0].split('_')[0]) #accounts for unpaired and paired files
        raw_files = raw_files.drop_duplicates(subset='run_prefix',keep = 'first')
        raw_files = raw_files.rename({'Total Sequences':'raw_reads'},axis=1)
        raw_files = raw_files[['run_prefix','raw_reads']]

        fastp_files['run_prefix'] = fastp_files['Sample'].apply(lambda x: x.split('.')[0].split('_')[0]) #accounts for unpaired and paired files
        fastp_files = fastp_files.rename({'Total Sequences':'post_qc_reads'},axis=1)
        fastp_files = fastp_files.drop_duplicates(subset='run_prefix',keep = 'first')
        fastp_files = fastp_files[['run_prefix','post_qc_reads']]

        filt_files['run_prefix'] = filt_files['Sample'].apply(lambda x: x.split('.')[0].split('_')[0]) #accounts for unpaired and paired files
        filt_files = filt_files.rename({'Total Sequences':'non_human_reads'},axis=1)
        filt_files = filt_files.drop_duplicates(subset='run_prefix',keep = 'first')
        filt_files = filt_files[['run_prefix','non_human_reads']]

        #now merge everything
        merge_raw=partial_df.merge(raw_files,how='left',left_on='run_prefix',right_on='run_prefix')
        merge_fastp=merge_raw.merge(fastp_files,how='left',left_on='run_prefix',right_on='run_prefix')
        merge_filt=merge_fastp.merge(filt_files,how='left',left_on='run_prefix',right_on='run_prefix')

        #last calculate fraction
        merge_filt['frac_non_human_from_raw']=merge_filt['non_human_reads']/merge_filt['raw_reads']
        merge_filt['frac_non_human_from_qc']=merge_filt['non_human_reads']/merge_filt['post_qc_reads']

        return merge_filt
    else:
        logger.warning("No multiqc file found. Check output directory and paths and try again. Returning unchanged dataframe.")
        return partial_df


### END: FASTQ FILE PROCESSING METHODS 

###START: MAIN FUNCTION

if __name__ == '__main__':

    #set up logging
    handler = logging.StreamHandler()
    fmt_str = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handler.setFormatter(logging.Formatter(fmt_str))
    logger = logging.getLogger(__name__)
    logger.addHandler(handler)
    
    #set up parameters options and defaults
    parser = ArgumentParser(description='Please note that the following ' +
                            'packages have to be installed for running this ' +
                            'script: 1)lxml 2)pandas 3)glob 4)csv 5)sys ' +
                            '6)urllib 7)argparse 8)requests 9)xmltodict ' +
                            '10)subprocess 11)bioconda 12)sra-tools 13)os ' +
                            '14)entrez-direct 15)pyyaml')
    parser.add_argument("-project","--project", nargs='*',
                        help="EBI/ENA project or study accession(s) " +
                        "to retrieve")
    parser.add_argument("-o","--output", default='./',
                        help='directory for output files. Default is working directory.')
    parser.add_argument("-mode", "--mode", default='ebi',
                        help="sra accession " +
                        "repository to be queried.", choices=['ebi','sra'])
    parser.add_argument("-prefix", "--prefix", nargs='*',
                        help="prefix(es) to prepend to output info files")
    parser.add_argument("-src","--sources",nargs='*',choices=['GENOMIC','GENOMIC SINGLE CELL','TRANSCRIPTOMIC',
                                                              'TRANSCRIPTOMIC SINGLE CELL','METAGENOMIC',
                                                              'METATRANSCRIPTOMIC','SYNTHETIC','VIRAL RNA','OTHER'],
                        help="list of one or more sources for restricting sample selection.")
    parser.add_argument("-strat","--strategies",nargs='*',choices=['POOLCLONE','CLONE','CLONEEND','WGS','WGA',
                                                       'WCS','WXS','AMPLICON','ChIP-Seq','RNA-Seq',
                                                       'MRE-Seq','MeDIP-Seq','MBD-Seq','MNase-Seq',
                                                       'DNase-Hypersensitivity','Bisulfite-Seq','EST',
                                                       'FL-cDNA','miRNA-Seq','ncRNA-Seq','FINISHING',
                                                       'TS','Tn-Seq','VALIDATION','FAIRE-seq','SELEX',
                                                       'RIP-Seq','ChIA-PET','RAD-Seq','Other'],
                        help="list of one or more libary strategies to restrict sample selection.")
    parser.add_argument("-plat", "--platforms", nargs='*', choices=['LS454','Illumina','Ion Torrent','PacBio_SMRT','OXFORD_NANOPORE'],
                        help="List of one or more platforms to restrict sample selection.")
    parser.add_argument("-name","--scientific_names",nargs='*',
                        help="List of scientific_names to restrict for selection.")
    parser.add_argument("-yaml","--validators",nargs='*',
                        help="one or more yaml files in QIIMP format for validation.")
    parser.add_argument("-yaml-dir","--yaml_dir", default ='./',
                        help="One or more yaml files in QIIMP format for validation. Loads yml files in ./ by default.")
    parser.add_argument("-no-seqs", "--no_seqs", default=False,action='store_true',
                        help="Omit download of fastq files.")
    parser.add_argument("-v", "--verbose", default=False, action='store_true',
                        help="Output additional messages.")
    parser.add_argument("-log", "--log",default='./output.log',
                        help="filename for logger. Defaults to [output_dir]/[ProjectID]_output.log")
    parser.add_argument("-prep-max","--prep_max",type=int, default=10000,
                        help="Max number of samples per prep info file.")
    parser.add_argument("-f","--force_yaml",default=False, action='store_true',
                        help="Advanced: force use of specified yaml for validation.")
    parser.add_argument("-hd","--host_deplete",default=False, action='store_true',
                        help="Advanced: host deplete using bowtie2. Uses human_PhiX db on barnacle by default.")
    parser.add_argument("-db","--host_db",
                        help="Advanced: specify the path to the host database for depletion.")
    parser.add_argument("-p","--cpus",type=int, default=4,
                        help="Number of processors to use during host depletion. Default is 4.")
    parser.add_argument("-fhd","--force_host_depletion",default=False, action='store_true',
                        help="Advanced: force host depletion for non-supported libary strategies. May break.")
    parser.add_argument("-map","--empo_mapping",
                        help="Advanced: .tsv with sample_type to EMPO mappings to use for normalization. " 
                                                    + "Should contain the following columns: " 
                                                    + " [sample_type, simple_sample_type, empo_1, empo_2,empo_3]")
    parser.add_argument("-method","--depletion_method",default='bowtie2',
                        help="Advanced: set host depletion method. bowtie2 by default. TODO implement minimap2 as second option.")
    parser.add_argument("-qc","--run_qc",default=False,action='store_true',
                        help="Advanced: run fastqc and generate multiqc report on downloaded and host-depleted files.")
    parser.add_argument("-sip","--max_samples",type=int,
                        help="Advanced: Max number of samples to grab from the study.")
    parser.add_argument("-rand","--random",default=False,action='store_true',
                        help="Advanced: when sampling, randomly select subset for processing. N.B. must supply a number"
                             + " with -sip (or --max_samples).")
    parser.add_argument("-qual","--quality_filter",default=False,action='store_true',
                        help="Advanced: run quality filtering of raw data. Automatically enabled when host depleting.")

    # parse the flags and initialize output file names
    args = parser.parse_args()

    if args.project is None:
        logger.warning("""
                python EBI_SRA_Downloader.py -project [accession] [accession ... N]
                    Generate the study info, study detail, prep, and  sample
                    files for the entered EBI accession, and download the
                    FASTQ files.
                Optional flags:
                    -output [directory where files will be saved]
                    -mode [specifies which repository to use]                    
                    -prefix [list of prefixes for sample and prep info files]
                    --strategy [list of one or more library strategies to select]
                    --sources [list of one or more library sources to select]
                    --platforms [list of one or more sequencing platforms to select]
                    --scientific_names [list of one or more scientific names to select]
                    --validators [list of one or more yaml files to use in validating]
                    --no_seqs [skip downloading files]
                    --prep_max [Max number of samples per prep info file: https://qiita.ucsd.edu/static/doc/html/faq.html? 
                      highlight=size#how-should-i-split-my-samples-within-preparations]
                    --verbose
               """)
        sys.exit(2)
    else:
        #settings
        mode=args.mode 
        DEBUG = args.verbose 
        omit_seqs = args.no_seqs 
        force = args.force_yaml
        max_prep =args.prep_max
        host_deplete = args.host_deplete
        proc= args.cpus
        force_hd = args.force_host_depletion
        hd_method=args.depletion_method
        qc = args.run_qc
        quality_filter=args.quality_filter
        random_sample = args.random

        if args.max_samples is not None:
            max_samples = args.max_samples
            subset = True
        else:
            subset = False

        if random_sample and subset:
            raise Exception("Must supply max_samples when requesting a random sample selection.")

        db_dict = {'bowtie2':'/databases/bowtie/Human_phiX174/Human_phix174',
                   'minimap2':'/databases/minimap2/human-phix-db.mmi'
                   }

        if args.host_db is not None:
            host_db = args.host_db
        else:
            try:
                host_db=db_dict[hd_method]
            except:
                raise Exception("No database found for host depletion method " + hd_method + " Please check spelling and try again.")

        if args.log is not None:
            fh=logging.FileHandler(args.log)
            logger.addHandler(fh)

        if DEBUG: logger.setLevel(logging.INFO)

        # Output directory
        output = args.output
        if not path.exists(output):
            makedirs(output)
        if list(output)[-1] != '/':
            output = output + '/'

        #set up validators
        yaml_validator_dict = {}
        yaml_list =  []

        if args.validators is not None:
            for y in args.validators:
                yaml_list.append(y)
        else:
            for file in os.listdir(args.yaml_dir):
                if file.endswith(".yml") or file.endswith(".yaml"):
                    yaml_list.append(os.path.join(args.yaml_dir, file))

        #since we're providing a way to override parsing, check assumption that a single validator is being passed
        if force:
            if len(yaml_list) > 1:
                raise Exception("Error: more than one yaml supplied with 'force' mode. Please supply only one yaml file." +\
                                "Note .yml and .yaml files in working directory (or specified by --yaml_dir) are loaded if " +\
                                " no -yaml flag is supplied. Loaded: " + str(yaml_list))
            elif len(yaml_list) == 0:
                raise Exception("Error: no yaml supplied with 'force' mode. Please supply only one yaml file.")
            else:
                yaml_validator_dict = yaml_list
        else:
            yaml_validator_dict = set_yaml_validators(yaml_list)

        if args.empo_mapping is not None:
            empo_mapping_dict = set_empo_normalizers(args.empo_mapping)
        else:
            empo_mapping_dict = {}

        platforms=[]
        if args.platforms is not None:
            for p in args.platforms:
                platforms.append(p.lower())

        strategies =[]
        if args.strategies is not None:
            for s in args.strategies:
                strategies.append(s.lower())

        names=[]
        if args.scientific_names is not None:
            for n in args.scientific_names:
                names.append(n.lower())

        sources=[]
        if args.sources is not None:
            for c in args.sources:
                sources.append(c.lower())

        # Retreive study information
        p_count=0
        for p in args.project:
            if args.prefix is not None:
                if len(args.prefix) == 1:
                    file_prefix = output + args.prefix[0] + '_' + p
                elif len(args.prefix) == len(args.project):
                    file_prefix = output + args.prefix[p_count] + '_' + p
                else:
                    raise Exception("Number of prefixes does not match number of projects. Set a single prefix or matched prefixes.")
            else:
                file_prefix = output + p
            p_count +=1    

            study = get_study_details(p,mode,file_prefix)

            if subset:
                study = sub_sample_study(study,max_samples,random_sample)

            #get and tidy sample and prep metadata
            md_tuple=get_sample_info(study,mode,platforms,strategies,yaml_validator_dict,file_prefix,names,sources,empo_mapping_dict)

            #write out sample and prep info files for reference while downloading
            write_info_files(md_tuple,max_prep,file_prefix)

            if not omit_seqs:
                valid_samples = fetch_sequencing_data(md_tuple[0],output,mode,host_deplete,host_db,proc,hd_method,qc,quality_filter)
                md_tuple = valid_samples, md_tuple[1]
                

            #run multiqc if requested
            if qc:
                multiqc_list = glob.glob(output + '/*_fastqc.html')
                run_qc(multiqc_list,output,proc,True)
                if host_deplete:
                    qc_df = add_qc_data(md_tuple[0],output+'/fastqc')
                    md_tuple = qc_df, md_tuple[1]

                    #re-write out sample and prep info files if qc+host_deplete performed since % non-human reads will be added
                    write_info_files(md_tuple,max_prep,file_prefix)

            #kludge to handle output log for now
            if path.isfile('./output.log'):
                shutil.move('./output.log',file_prefix + '_output.log')
