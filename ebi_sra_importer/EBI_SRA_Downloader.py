#!/usr/bin/env python3
# conda command to install all dependencies:
#   conda create -n ebi_sra_importer pandas requests entrez-direct sra-tools xmltodict lxml pyyaml xlrd -c bioconda -c conda-forge -y
# to enable host depletion use:
#   conda create -n ebi_sra_importer pandas requests entrez-direct sra-tools bowtie2 samtools bedtools xmltodict lxml pyyaml xlrd -c bioconda -c conda-forge -y
# or to enable later after installation:
#   conda activate ebi_sra_importer
#   conda install bowtie2 samtools bedtools
#
# pip command to install all dependencies:
#   pip install csv glob requests subprocess xmltodict sys lxml os urllib pyyaml xlrd
#   pip install argparse pandas bioconda sra-tools entrez-direct
#   pip install bowtie2 samtools bedtools
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
#                   --sep [provide a delimiter for parsing the description field]
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

def split_caps(cap_string):
    temp_list = []
    temp_list = re.findall('[A-Z][^A-Z]*', cap_string)
    if len(temp_list) > 0:
        return "_".join(temp_list)
    else:
        return cap_string
        
def scrub_special_chars(input_string,custom_dict={}):
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
    for k in replace_dict.keys():
        input_string=input_string.replace(k,replace_dict[k])
    for k in custom_dict.keys():
        input_string=input_string.replace(k,custom_dict[k])
    return input_string

def qiimp_parser(filename):
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
            logger.warning("Could not load yaml from " + filename + "Contents:")
            logger.warning(subprocess.run(['cat',file]))
    else:
        logger.warning("Invalid file extension for yaml parsing: " + str(ext))

    if len(parsed_yml) == 0:
        logger.warning("The file " + filename +" contains no yaml data. Please check contents and try again.")

    return parsed_yml

def set_yaml_validators(validators,fields=['scientific_name']):
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

#these functions create the study_details files for ebi and sra respectively
def write_config_file(xml_dict,prefix,mode,xml_proj_dict={}):
    config_string='[required\ntimeseries_type_id = 1\nmetadata_complete = True\nmixs_compliant = True\n principal_investigator = ebi-import\nreprocess = False'
    if mode == 'ebi':
        parse_dict = xml_dict['STUDY_SET']['STUDY']
    elif mode == 'sra':
        parse_dict = xml_dict['STUDY']
    else:
        logger.warning('Received invalid mode: ' + mode)

    title= '\n study_title = XXEBIXX'
    alias= '\nstudy_alias = XXEBIXX'
    abstract ='\nstudy_abstract = XXEBIXX'
    description = '\nstudy_description = XXEBIXX'

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

    config_string= config_string + alias + description + abstract + 'efo_ids = 1\[optional]' + title

    study_file = prefix + '_study_config.txt'

    logger.info("In write config")
    # Write out study_file
    file = open(study_file, "w")
    file.write(config_string)
    file.close()

def get_study_details(study_accession,mode='ebi',prefix=''):
    #need to add check for invalid url!
    if mode == 'ebi':
        studyUrl = "http://www.ebi.ac.uk/ena/data/view/" + study_accession \
                    + "&display=xml"
        response = requests.get(studyUrl)
        xml_dict=parse(response.content)

        if 'STUDY_SET' in xml_dict.keys():
            logger.info(study_accession + " is study ID. Writing config file")
            write_config_file(xml_dict,prefix,mode)
        elif 'PROJECT_SET' in xml_dict.keys():
            try:
                secondary_accession = xml_dict['PROJECT_SET']['PROJECT']['IDENTIFIERS']['SECONDARY_ID']
                secondaryUrl = "http://www.ebi.ac.uk/ena/data/view/" + secondary_accession + "&display=xml"
                logger.warning(study_accession + " is project ID. Retrieved secondary ID: " + secondary_accession + " Writing config file.")
                logger.info("Alternate url" + secondaryUrl)
                response2 = requests.get(secondaryUrl)
                xml_study_dict=parse(response2.content)
                write_config_file(xml_study_dict,prefix,mode,xml_dict)

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

    #create_details_file(study_df,study_accession,mode,prefix)
    return study_df


#now feed in study_df above
def get_sample_info(input_df,mode='ebi',plat=[],strat=[],validator_files={},prefix='',names=[],src=[],empo_mapping=pd.DataFrame(),sample_type_column='scientific_name'):

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
        
        for index, row in input_df.iterrows():
            sample_accession = row[identifier]
            prep_type = row['library_strategy']
            if prep_type not in sample_count_dict.keys() :
                sample_count_dict[prep_type]= {sample_accession:0}
            elif sample_accession not in sample_count_dict[prep_type].keys():
                sample_count_dict[prep_type][sample_accession]= 0
            else:
                sample_count_dict[prep_type][sample_accession] = sample_count_dict[prep_type][sample_accession] + 1

            sampleUrl = "http://www.ebi.ac.uk/ena/data/view/" + sample_accession \
                    + "&display=xml"
            if DEBUG: logger.info(sampleUrl)

            response = requests.get(sampleUrl)
            xml_dict=parse(response.content)
            if 'SAMPLE_SET' in xml_dict.keys():

                input_df.at[index,'sample_title_specific']=xml_dict['SAMPLE_SET']['SAMPLE']['TITLE']

                sn=xml_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_NAME']
                for s in sn.keys():
                    col = scrub_special_chars(s).lower()
                    #print(col)
                    input_df.at[index,col]=sn[s]

                sa=xml_dict['SAMPLE_SET']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE']
                for s in sa:
                    col = scrub_special_chars(s['TAG']).lower()
                    #print(col)
                    try:
                        input_df.at[index,col]=s['VALUE']
                    except:
                        logger.warning('No value found for sample attribute: ' + col + '. Setting to "not provided".')
                        input_df.at[index,col]='not provided'
                input_df.at[index,'prep_file']=prep_type + '_' + str(sample_count_dict[prep_type][sample_accession])
            else:
                logger.warning('No metadata found for sample named: ' + sample_accession + ' omitting.')
        #set sample_name based on identifier column
        input_df['sample_name']=input_df[identifier]

        #need to catch common issue where identifier is identical for unique samples
        #for now assume that in this case, libarary_name will be unique, and if it isn't combined sample and run names

        if len(input_df) > 1 and input_df[identifier].nunique() == 1:
            if input_df['library_name'].nunique() != 1:
                #print(str(len(_input_df)) + " is length and input_df[identifier].nunique() = " + str(input_df[identifier].nunique()))
                input_df['sample_name']=input_df['library_name'].apply(lambda x: scrub_special_chars(x))
            else:
                input_df['sample_name']=input_df[identifier]+ '.' + input_df[run_accession]

        input_df['run_prefix']=input_df[run_accession]

        #the loop above takes the most time so write out a placeholder file for faster re-running if interrupted
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

    return output_df

def normalize_types(md,mapping):

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
    
def create_details_file(study_details_df, study_accession,mode='ebi',prefix='',file_suffix="_detail"):
    """Returns the details of the EBI/SRA study

    If the accession ID is valid, generate a .details.txt, and return the
    detail file name of this EBI/SRA study. Else return None

    Parameters
    ----------
    study_accession : string
        The accession ID of the EBI/SRA study

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
    
def write_info_files(final_df,max_prep,prefix=''):
    if max_prep > len(final_df):
        max_prep = len(final_df)
        
    prep_info_columns = ['run_prefix','experiment_accession','platform','instrument_model','library_strategy',
                         'library_source','library_layout','library_selection','fastq_ftp','ena_checklist',
                         'ena_spot_count','ena_base_count','ena_first_public','ena_last_update','instrument_platform',
                         'submitted_format','sequencing_method','target_gene','target_subfragment','primer']
    amplicon_min_prep_list=['target_gene','target_subfragment','primer']
    amplicon_type_preps = ['AMPLICON','OTHER']
    final_df.columns =[scrub_special_chars(col).lower() for col in final_df.columns]
    if DEBUG: logger.info(final_df.columns)
    
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
        prep_df = final_df[final_df['prep_file']==prep_file]
        if prep_file in amplicon_type_preps: #check to see if the prep is amplicon-style, specified by list above
            for min_prep in amplicon_min_prep_list: #if amplicon-style, enforce presence or null values for minimum prep info information
                if min_prep not in prep_df.columns():
                    prep_df[min_prep]='XXEBIXX' #will throw warning, but okay with current pandas
        prep_df= prep_df[prep_df.columns[prep_df.columns.isin(prep_info_columns)]].set_index('sample_name')
        prep_df=prep_df.dropna(axis=1,how='all')
        prep_df_list = [prep_df[i:i+max_prep] for i in range(0,prep_df.shape[0],max_prep)]
        prep_count=0
        for prep in prep_df_list:
            prep.to_csv(prefix+'_prep_info_'+ prep_file + '_part' + str(prep_count) +'.tsv',sep='\t',index=True,index_label='sample_name')
            prep_count += 1


def run_host_depletion(fastq_file_list,filter_db='',cpus=4,output_dir='./',method='bowtie2'):
    db_dict = {'bowtie2':'/databases/bowtie/Human_phiX174/Human_phix174',
                'minimap2':'/databases/minimap2/human-phix-db.mmi'
              }
    if filter_db == '':
        filter_db = db_dict[method]

    read1_fastq=str(fastq_file_list[0])
    filtered_fastq_1=read1_fastq.replace('.fastq.gz','.R1.filtered')

    if len(fastq_file_list) == 2:
        read2_fastq=str(fastq_file_list[1])
    else:
        read2_fastq = '' #to reduce code complexity for checking for paired vs unpaired, dummy file string here

    filtered_fastq_2=read2_fastq.replace('.fastq.gz','.R2.filtered')

    if method == 'bowtie2':
        bowtie2_args=['bowtie2', '-p', str(cpus), '-x', filter_db]
        stv_args_1 =['samtools', 'view', '-f', '12', '-F', '256']
        sts_args =['samtools', 'sort', '-@', str(cpus),'-n']
        stv_args_2 = ['samtools', 'view','-bS']
        btb_args = ['bedtools', 'bamtofastq', '-i', '-', '-fq', filtered_fastq_1]

        if len(fastq_file_list) == 2:
            bowtie2_args=bowtie2_args + ['-1',read1_fastq,'-2',read2_fastq]
            btb_args = btb_args + ['-fq2',filtered_fastq_2]
        else:
            bowtie2_args=bowtie2_args + ['-U',read1_fastq]

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
            subprocess.run(['rm', read1_fastq])
        if path.isfile(filtered_fastq_2):
            subprocess.run(['gzip', filtered_fastq_2])
            subprocess.run(['rm', read2_fastq])

    elif method == 'minimap2':

        output_fastq1=filtered_fastq_1 + '.gz'
        fastp_args=['fastp','-l','100','-i', read1_fastq]
        minimap2_args = ['minimap2','-ax', 'sr', '-t', str(cpus), filter_db,"-","-a"]
        stf_args=['samtools','fastq', '-@', str(cpus), '-f', '12', '-F', '256', '-1',output_fastq1]

        if len(fastq_file_list)==2:
            fastp_args = fastp_args + ['-I',read2_fastq]
            output_fastq2 =filtered_fastq_2 + '.gz'
            stf_args= stf_args+ ['-2',output_fastq2]

        #add final flags for fastp after determining input file count
        fastp_args = fastp_args + ['-w',str(cpus),'--stdout']

        fastp_ps = subprocess.Popen(fastp_args, stdout=subprocess.PIPE)
        minimap2_ps = subprocess.Popen(minimap2_args,stdin=fastp_ps.stdout, stdout=subprocess.PIPE)
        stf_ps = subprocess.Popen(stf_args,stdin=minimap2_ps.stdout, stdout=subprocess.PIPE)
        stf_ps.wait()

        #logger.warning("minimap2 not yet supported. Skipping host depletion.")
    else:
        logger.warning("Selected depletion method '" + method + "' not currently supported. Please select either bowtie2 or minimap2.")

def fetch_sequencing_data(download_df,output_dir="./",mode='ebi',host_deplete=False,host_db='/databases/bowtie/Human_phiX174/Human_phiX174',cpus=4):
    """Fetch all the meta file(s) for EBI study

    Parameters
    ----------
    study_accession : string
        The accession ID of the EBI study

    study_details : string
        study detail file name

    """

    logger.info("Downloading the fastqs")
    # Download the fastqs
    
    for index, row in download_df.iterrows():
        if mode == 'ebi':
            try:
                files = row['fastq_ftp'].split(';')
            except:
                try:
                    files = row['download_path'].split(';')
                except:
                    logger.warning("Skipping sample:" + row['sample_name']
                                    + ", run: " + row['run_accession'] + "No fastq ftp found.")
            if host_deplete:
                if len(files) > 2 and host_deplete:
                    logger.warning("More than 2 files in the fastq download path. Skipping host depletion for "
                     + row['sample_name'])
            hd_file_list =[]
            for f in files:
                fq_path = output_dir + "/" + f.split('/')[-1]
                hd_file_list.append(fq_path)
                if DEBUG: logger.info(f)
                if type(f) != str:
                    logger.warning("Skipping sample:" + row['sample_name']
                                   + ", run: " + row['run_accession'] + "fastq ftp path is not string.")      
                elif path.isfile(fq_path):
                    logger.warning("Skipping " + fq_path)
                    logger.warning("File exists")
                else:
                    urlretrieve("ftp://" +f, fq_path)

            valid_hd_prep_types = ['WGS','WGA','WXS','ChIP-Seq','RNA-Seq']
            if host_deplete:
                if len(hd_file_list) > 2:
                    logger.warning("More than 2 files in the fastq download path. Skipping host depletion for "
                         + row['sample_name'])
                else:
                    if row['library_strategy'] in valid_hd_prep_types:
                        logger.info("Running host depletion on " + row['sample_name'])
                        run_host_depletion(hd_file_list,host_db,cpus,output_dir)
                    elif force_hd:
                        logger.warning("Forcing host depletion for " + + " with library_strategy " + + ". This may break...")
                        run_host_depletion(hd_file_list,host_db,cpus,output_dir)
                    else:
                        logger.warning("Skipping host depletion for " + row['sample_name'] + " Library Strategy is "
                                     + row['library_strategy'] + ". Valid formats are: "+ str(valid_hd_prep_types))
        elif mode =='sra':
            subprocess.run(['fastq-dump', '-I', '--split-files', '--gzip', '--outdir', output_dir, row['run_accession']])
            #TODO implement host depletion for sra
        else:
            raise Exception(mode + " is not a valid repository")


### END of methods ###

if __name__ == '__main__':
    # parse the flags and initialize output file names
    # TODO:Reword help info

    #set up logging
    handler = logging.StreamHandler()
    fmt_str = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handler.setFormatter(logging.Formatter(fmt_str))
    logger = logging.getLogger(__name__)
    logger.addHandler(handler)

    parser = ArgumentParser(description='Please note that the following ' +
                            'packages have to be installed for running this ' +
                            'script: 1)lxml 2)pandas 3)glob 4)csv 5)sys ' +
                            '6)urllib 7)argparse 8)requests 9)xmltodict ' +
                            '10)subprocess 11)bioconda 12)sra-tools 13)os ' +
                            '14)entrez-direct 15)pyyaml')
    parser.add_argument("-project","--project", nargs='*',help="EBI/ENA project or study accession(s) " +
                        "to retrieve")
    parser.add_argument("-o","--output", default='./',help='directory for output files. Default is working directory.')
    parser.add_argument("-mode", "--mode", default='ebi', help="sra accession " +
                        "repository to be queried.", choices=['ebi','sra'])
    parser.add_argument("-prefix", "--prefix", nargs='*', help="prefix(es) to prepend to output info files")
    parser.add_argument("-src","--sources",nargs='*', help="list of one or more sources for restricting sample selection.")
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
    parser.add_argument("-name","--scientific_names",nargs='*',help="List of scientific_names to restrict for selection.")
    parser.add_argument("-yaml","--validators",nargs='*', help="one or more yaml files in QIIMP format for validation.")
    parser.add_argument("-yaml-dir","--yaml_dir", default ='./', help="One or more yaml files in QIIMP format for validation. Loads yml files in ./ by default.")
    parser.add_argument("-no-seqs", "--no_seqs", default=False,action='store_true', help="Omit download of fastq files.")
    parser.add_argument("-v", "--verbose", default=False, action='store_true', help="Output additional messages.")
    parser.add_argument("-log", "--log",default='./output.log',help="filename for logger. Defaults to [output_dir]/[ProjectID]_output.log")
    parser.add_argument("-prep-max","--prep_max",type=int, default=10000,help="Max number of samples per prep info file.")
    parser.add_argument("-f","--force_yaml",default=False, action='store_true', help="Advanced: force use of specified yaml for validation.")
    parser.add_argument("-hd","--host_deplete",default=False, action='store_true', help="Advanced: host deplete using bowtie2. Uses human_PhiX db on barnacle by default.")
    parser.add_argument("-db","--host_db",default='/databases/bowtie/Human_phiX174/Human_phix174', help="Advanced: specify the path to the host database for depletion.")
    parser.add_argument("-p","--cpus",type=int, default=4,help="Number of processors to use during host depletion. Default is 4.")
    parser.add_argument("-fhd","--force_host_depletion",default=False, action='store_true', help="Advanced: force host depletion for non-supported libary strategies. May break.")
    parser.add_argument("-map","--empo_mapping",help="Advanced: .tsv with sample_type to EMPO mappings to use for normalization. Should contain the following columns: [sample_type, simple_sample_type, empo_1, empo_2,empo_3]")
    parser.add_argument("-method","--depletion_method",default='bowtie2',help="Advanced: set host depletion method. bowtie2 by default. TODO implement minimap2 as second option.")

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
                    --prep_max [Max number of samples per prep info file: https://qiita.ucsd.edu/static/doc/html/faq.html?highlight=size#how-should-i-split-my-samples-within-preparations]
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
        host_db = args.host_db
        proc= args.cpus
        force_hd = args.force_host_depletion
        hd_method=args.depletion_method

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

            #tidy metadata
            md=get_sample_info(study,mode,platforms,strategies,yaml_validator_dict,file_prefix,names,sources,empo_mapping_dict)

            #write out files
            write_info_files(md,max_prep,file_prefix)
            logger.info(host_db)
            if not omit_seqs:
                fetch_sequencing_data(md,output,mode,host_deplete,host_db,proc)

            #kludge to handle output log for now
            if path.isfile('./output.log'):
                shutil.move('./output.log',file_prefix + '_output.log')
