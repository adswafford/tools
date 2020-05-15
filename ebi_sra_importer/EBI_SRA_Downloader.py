#!/usr/bin/env python3
# conda command to install all dependencies:
#   conda create -n ebi_sra_importer pandas requests entrez-direct sra-tools xmltodict lxml -c bioconda -c conda-forge -y
#
# pip command to install all dependencies:
#   pip install csv glob requests subprocess xmltodict sys lxml os urllib
#   pip install argparse pandas bioconda sra-tools entrez-direct
#
# Instruction:
#       -ebi, -sra flags allow the user to choose download source
#           python3 EBI_SRA_Downloader.py -ebi {ebiaccession_number}
#           python3 EBI_SRA_Downloader.py -sra {sraaccession_number}
#       -sample {name} flag allows the user to specify the sample file name
#       -prep {name} flag allows the user to specify the prep file name
#       -study {name} flag allows the user to specify the study info file name
#       -all-seqs allows the script to accept all sample types
#       -all-platforms allows the script to accept samples from all platforms
#       -debug true flag to enter debug mode (not download fastq files)
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
from os import path, makedirs
from urllib.request import urlretrieve
from argparse import ArgumentParser
from pandas import read_csv, DataFrame

DEBUG = False
ALL_SEQS = False
ALL_PLATFORMS = False
handler = logging.StreamHandler()
fmt_str = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
handler.setFormatter(logging.Formatter(fmt_str))
logger = logging.getLogger(__name__)
logger.addHandler(handler)


def ebi_create_details_file(study_accession, file_suffix="_detail"):
    """Returns the details of the EBI study

    If the accession ID is valid, generate a .details.txt, and return the
    detail file name of this EBI study. Else return None

    Parameters
    ----------
    study_accession : string
        The accession ID of the EBI study

    file_suffix : string
        The suffix for the output study detail file

    Returns
    -------
    string
        study details file name
    """
    # Grab the details related to the given accession
    study_details = study_accession + file_suffix + ".txt"
    host = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
    read_type = "&result=read_run&"
    fields = "library_name,secondary_sample_accession,run_accession," + \
             "experiment_accession,fastq_ftp,library_source," + \
             "instrument_platform,submitted_format,library_strategy," +\
             "library_layout"
    url = ''.join([host, study_accession, read_type, "fields=", fields])
    response = requests.get(url)

    # Check for valid accessions
    if not response:
        raise Exception(study_accession + " is not a valid EBI study ID")

    logger.info("generating study detail file")
    detailWoHeader = response.content.partition(b'\n')[2]  # remove header
    detailWoHeader = detailWoHeader.split(b'\n')
    finalDetail = ""
    for row in detailWoHeader:
        if len(row) == 0:
            continue
        row_list = row.decode("utf-8").split('\t')
        if not ALL_SEQS and row_list[5].upper() != "METAGENOMIC":
            logger.warning("Library source is " + row_list[5] +
                            " not Metagenomic for " +
                           row_list[1] + ". Omitting " + row_list[1])
            continue
            # skip row
        elif not ALL_PLATFORMS and row_list[6].lower() != "illumina":
            logger.warning("Instrument platform is " + row_list[6] + 
                            " not Illumina for " +
                           row_list[1] + ". Omitting " + row_list[1])
            continue
            # skip row
        else:
            for i in range(len(row_list)):
                if len(row_list[i]) == 0:
                    row_list[i] = "unspecified"
            if ALL_PLATFORMS:
                row_string = '\t'.join(row_list)
            else:
                row_string = '\t'.join(row_list[:6]) + "\tIllumina\t" + \
                    '\t'.join(row_list[7:])

        if len(finalDetail) == 0:
            finalDetail = row_string
        else:
            finalDetail = finalDetail + '\n' + row_string
    if len(finalDetail) == 0:
        if ALL_SEQS:
            raise Exception(study_accession + " has no sample or run that" +
                            " is from Illumina")
        elif ALL_PLATFORMS:
            raise Exception(study_accession + " has no sample or run that" +
                            " is METAGENOMIC")
        else:
            raise Exception(study_accession + " has no sample or run that" +
                            " is METAGENOMIC or from Illumina")
    with open(study_details, 'w') as file:
        file.write(finalDetail)
    return study_details


def sra_create_details_file(study_accession, file_suffix="_detail"):
    """Returns the details of the SRA study

    If the accession ID is valid, generate a .details.txt, and return the
    detail file name of this SRA study. Else return None

    Parameters
    ----------
    study_accession : string
        The accession ID of the SRA study

    file_suffix : string
        The suffix for the output study detail file

    Returns
    -------
    string
        study details file name
    """
    p1 = subprocess.Popen(['esearch', '-db', 'sra', '-query', study_accession],
                          stdout=subprocess.PIPE)
    for i in p1.stdout:
        if "<Count>0</Count>" in i.decode("utf-8"):
            p1.stdout.close()
            raise Exception(study_accession + " is not a valid SRA study ID")

    p1 = subprocess.Popen(['esearch', '-db', 'sra', '-query', study_accession],
                          stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['efetch', '-format', 'runinfo'], stdin=p1.stdout,
                          stdout=subprocess.PIPE)
    p1.stdout.close()
    logger.info("generating study detail file")
    # since the index of each filed efetch returns may differ
    # we use a dictionary to get the indices of the fileds we need in the
    # return content
    indices = {
        'library_name': -1,
        'sample_accession': -1,
        'run_accession': -1,
        'experiment_accession': -1,
        'download_path': -1,
        'library_source': -1,
        'instrument_platform': -1,
        'submitted_format': -1,
        'library_strategy': -1,
        'library_layout': -1
    }
    isHeader = True
    output = ""
    for i in p2.stdout:
        if isHeader:
            # the header (the first line) returned by efetch has all the
            # fields' names. we map the fields we need to the corresponding
            # locations in the header, and get the indices we need for fields.
            header = i.decode("utf-8").strip().split(',')
            indices['library_name'] = header.index('LibraryName')
            indices['sample_accession'] = header.index('Sample')
            indices['run_accession'] = header.index('Run')
            indices['experiment_accession'] = header.index('Experiment')
            indices['download_path'] = header.index('download_path')
            indices['library_source'] = header.index('LibrarySource')
            indices['instrument_platform'] = header.index('Platform')
            indices['library_strategy'] = header.index('LibraryStrategy')
            indices['library_layout'] = header.index('LibraryLayout')
            isHeader = False
        else:
            line = i.decode("utf-8").strip().split(',')
            if len(line) < indices[max(indices, key=indices.get)]:
                continue
            if not ALL_SEQS:
                if line[indices['library_source']].upper() != "METAGENOMIC":
                    logger.warning(line[indices['library_source']])
                    logger.warning("Library source is " +
                                    line[indices['library_source']] +
                                    " not Metagenomic for " +
                                   line[indices['run_accession']] +
                                   ". Omitting " +
                                   line[indices['run_accession']])
                    continue
            elif not ALL_PLATFORMS:
                if line[indices['instrument_platform']].lower() != "illumina":
                    logger.warning("Instrument platform is " + 
                                    line[indices['instrument_platform']] +
                                    " not Illumina for " +
                                    line[indices['run_accession']] + ". Omitting "
                                    + line[indices['run_accession']])
                    continue

            for key in indices:
                if indices[key] == -1 or len(line[indices[key]]) < 1:
                    output = output + "unspecified\t"
                else:
                    output = output + line[indices[key]] + "\t"
            output += "\n"
    p2.stdout.close()
    if len(output) == 0:
        if ALL_SEQS:
            raise Exception(study_accession + " has no sample or run that" +
                            " is from Illumina")
        elif ALL_PLATFORMS:
            raise Exception(study_accession + " has no sample or run that" +
                            " is METAGENOMIC")
        else:
            raise Exception(study_accession + " has no sample or run that is "
                            " is METAGENOMIC or from Illumina")
    with open(study_accession + file_suffix + ".txt", "w") as f:
        f.write(output)
    return study_accession + file_suffix + ".txt"


def ebi_create_info_file(study_file, study_accession, xmlFile='feed.xml'):
    """fetch the study information of the EBI study based on accession ID
        and generate a study file for the information

    Parameters
    ----------
    study_accession : string
        The accession ID of the EBI study

    study_file : string
        The study_info file name
    """
    # pull xml from website: https://www.ebi.ac.uk/ena/data/view/ERP012804
    URL = "http://www.ebi.ac.uk/ena/data/view/{0}&display=xml".\
        format(study_accession)
    response = requests.get(URL)

    if len(etree.fromstring(response.content).getchildren()) == 0:
        logger.warning(study_accession + " is a valid EBI study ID, but its "
                       + "information (http://www.ebi.ac.uk/ena/data/view/" +
                       study_accession + ") page is not supported.\nskipping" +
                       " creating study_info.txt")
        return
    is_study = etree.fromstring(response.content).getchildren()[0]
    if is_study.tag != 'STUDY':
        if is_study.find('IDENTIFIERS').find('SECONDARY_ID') is None:
            raise Exception(study_accession + " is not a valid EBI study ID")
        else:
            study_accession = is_study.find('IDENTIFIERS').\
                find('SECONDARY_ID').text
        URL = "http://www.ebi.ac.uk/ena/data/view/{0}&display=xml".\
            format(study_accession)
        response = requests.get(URL)

    logger.info("generating feed.xml file")
    # Creating a dummy xml to dump the output
    with open(xmlFile, 'wb') as file:
        file.write(response.content)
    with open(xmlFile) as fd:
        doc = parse(fd.read())

    logger.info("generating study info file")
    study_title = doc['ROOT']['STUDY']['DESCRIPTOR']['STUDY_TITLE']
    alias = doc['ROOT']['STUDY']['@alias']
    study_abstract = doc['ROOT']['STUDY']['DESCRIPTOR']['STUDY_ABSTRACT'] \
        if 'STUDY_ABSTRACT' in doc['ROOT']['STUDY']['DESCRIPTOR'] else "None"
    description = doc['ROOT']['STUDY']['DESCRIPTOR']['STUDY_DESCRIPTION'] \
        if 'STUDY_DESCRIPTION' in doc['ROOT']['STUDY']['DESCRIPTOR'] \
        else "None"
    # Creating default values for PI and env
    # TODO: Come up with a better way to handle these fields
    PI = "EBI import"  # default principal investigator
    env = "miscellaneous natural or artificial environment"
    # default environmental package

    # Write to the study_file
    file = open(study_file, "w")
    file.write("STUDY_TITLE" + "\t" + "ALIAS" + "\t" + "STUDY_ABSTRACT" +
               "\t" + "STUDY_DESCRIPTION" + "\t" +
               "PRINCIPAL_INVESTIGATOR" + "\t" +
               "ENVIRONMENTAL_PACKAGES" + "\n")
    file.write(study_title + "\t" + alias + "\t" + study_abstract + "\t"
               + description + "\t" + PI + "\t" + env + "\n")
    file.close()


def sra_create_info_file(study_file, study_accession):
    """fetch the study information of the SRA study based on accession ID
        and generate a study file for the information

    Parameters
    ----------
    study_accession : string
        The accession ID of the SRA study

    study_file : string
        The study_info file name
    """
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

    logger.info("generating study info file")
    study_title = doc['STUDY']['DESCRIPTOR']['STUDY_TITLE']
    alias = doc['STUDY']['@alias']
    study_abstract = doc['STUDY']['DESCRIPTOR']['STUDY_ABSTRACT'] \
        if 'STUDY_ABSTRACT' in doc['STUDY']['DESCRIPTOR'] else "None"
    description = doc['STUDY']['DESCRIPTOR']['STUDY_DESCRIPTION'] \
        if 'STUDY_DESCRIPTION' in doc['STUDY']['DESCRIPTOR'] \
        else "None"
    # Creating default values for PI and env
    # TODO: Come up with a better way to handle these fields
    PI = "SRA import"  # default principal investigator
    env = "miscellaneous natural or artificial environment"
    # default environmental package

    # Write to the study_file
    file = open(study_file, "w")
    file.write("STUDY_TITLE" + "\t" + "ALIAS" + "\t" + "STUDY_ABSTRACT" +
               "\t" + "STUDY_DESCRIPTION" + "\t" +
               "PRINCIPAL_INVESTIGATOR" + "\t" +
               "ENVIRONMENTAL_PACKAGES" + "\n")
    file.write(study_title + "\t" + alias + "\t" + study_abstract + "\t"
               + description + "\t" + PI + "\t" + env + "\n")
    file.close()


def create_prep_file(prep_file, study_details):
    """Generate the prep file(s) for EBI study

    Parameters
    ----------
    prep_file : string
        prep file name

    study_details : string
        study detail file name

    """
    logger.info("generating prep file(s)")
    file = open(study_details, 'r')
    details = [line.strip().split('\t') for line in file]
    file.close()
    details = sorted(details, key=lambda l: l[8])

    temp_library_strategy = details[0][8]
    prep_file_name = prep_file[:prep_file.rfind(".")] + "_" +\
        temp_library_strategy + prep_file[prep_file.rfind("."):]
    file = open(prep_file_name, "w")
    file.write("sample_name" + "\t" + "run_prefix" + "\t"
               + "experiment_accession" + "\t" + "instrument_platform" + "\t"
               + "library_strategy" + "\t" + "library_source" + "\t"
               + "library_layout" + "\n")
    file.close()

    sample_accession = []
    for row in details:
        # find which library_strategy it is
        library_strategy = row[8]
        if library_strategy != temp_library_strategy:
            temp_library_strategy = library_strategy
            sample_accession = []
            prep_file_name = prep_file[:prep_file.rfind(".")] + "_" +\
                temp_library_strategy + prep_file[prep_file.rfind("."):]
            file = open(prep_file_name, "w")
            file.write("sample_name" + "\t" + "run_prefix" + "\t"
                       + "experiment_accession" + "\t" + "instrument_platform"
                       + "\t" + "library_strategy" + "\t" + "library_source"
                       + "\t" + "library_layout" + "\n")
            file.close()

        sample = row[1]
        run_prefix = sample + "." + row[2]
        # find which prep file to write
        if sample in sample_accession:
            # sample ID appeared before
            sample_accession.append(row[1])
            sample_count = sample_accession.count(row[1])
            append2 = prep_file_name.rfind(".")
            next_prep_file = prep_file_name[:append2] + str(sample_count) +\
                prep_file_name[append2:]
            file_path = "./" + next_prep_file
            if not path.exists(file_path):
                f1 = open(next_prep_file, "w")
                f1.write("sample_name" + "\t" + "run_prefix" + "\t"
                         + "experiment_accession" + "\t"
                         + "instrument_platform" + "\t"
                         + "library_strategy" + "\t"
                         + "library_source" + "\t"
                         + "library_layout" + "\n")
                f1.close()
            write_prep = next_prep_file
        else:
            # sample ID not appeared before
            sample_accession.append(row[1])
            write_prep = prep_file_name
        # write to the selected prep file
        write_string = (row[1]) + "\t" + run_prefix + "\t" + (row[3]) + "\t" \
            + (row[6]) + "\t" + (row[8]) + "\t" + (row[5])\
            + "\t" + (row[9]) + "\n"
        with open(write_prep) as f1:
            if write_string in f1.read():
                continue
        with open(write_prep, 'a') as f1:
            f1.write(write_string)


def ebi_create_sample_file(sample_file, study_accession, study_details):
    """fetch the study information of each sample in EBI study
        and generate a sample file for these information

    Parameters
    ----------
    study_accession : string
        The accession ID of the EBI study

    sample_file : string
        The sample file name

    study_details : string
        The study detail file name
    """
    def xml_to_dict(xml_fp):
        # Converts xml string to a dictionary
        root = etree.parse(xml_fp).getroot()
        sample = root.getchildren()[0]
        metadata = {}

        attributes = sample.find('SAMPLE_ATTRIBUTES')
        for node in attributes.iterfind('SAMPLE_ATTRIBUTE'):
            tag = node.getchildren()[0]
            value = node.getchildren()[1]
            if value.text is None:
                metadata[tag.text.strip('" ').upper()] = 'not provided'
            else:
                metadata[tag.text.strip('" ').upper()] \
                    = value.text.strip('" ')

        #adding loops to look for additional data
        title= sample.find('TITLE')
        if title.text is None:
            metadata['title'] = 'not provided'
        else:
            metadata['title'] = title.text.strip('" ')
        description = sample.find('DESCRIPTION')
        try:
            if description.text is None:
                metadata['description'] = 'not provided'
            else:
                if sep in description.text:
                    split_desc=description.text.strip('" ').split(sep)
                    counter=0
                    for i in split_desc:
                        metadata['description_field_' + str(counter)] = i
                        counter += 1
                else:
                    metadata['description'] = description.text.strip('" ')
        except:
            metadata['description'] = 'not provided'

        nameInfo = sample.find('SAMPLE_NAME')
        for node in nameInfo:
            tag = node.tag
            value = node.text
            if value is None:
                metadata[tag.text.strip('" ').upper()] = 'not provided'
            else:
                metadata[tag.strip('" ').upper()] = value.strip('" ')

        idInfo = sample.find('IDENTIFIERS')
        for node in idInfo:
            value = node.text
            d = node.attrib
            if len(d) > 0:
                for k in d.keys():
                    tag = node.tag + "_" + d[k]
                    if value is None:
                        metadata[tag.text.strip('" ').upper()] = 'not provided'
                    else:
                        metadata[tag.strip('" ').upper()] = value.strip('" ')
            else:
                tag = node.tag

                if value is None:
                    metadata[tag.text.strip('" ').upper()] = 'not provided'
                else:
                    metadata[tag.strip('" ').upper()] = value.strip('" ')

        return metadata

    logger.info("downloading sample.txt file for each sample")
    details_df = read_csv(study_details, sep='\t', header=None)
    for row in details_df.iterrows():
        library_name = row[1][0]
        current_path = "./" + study_accession + "/" + library_name

        sample_accession = row[1][1]
        if sample_accession == 'unspecified': # and not DEBUG:
            raise Exception(sample_accession + " does not contain metadata")
        if path.exists(current_path + "/" + sample_accession + ".txt"):
            continue

        if not path.exists(current_path):
            makedirs(current_path)

        tempUrl = "http://www.ebi.ac.uk/ena/data/view/" + sample_accession \
            + "&display=xml"
        temp_response = requests.get(tempUrl)
        if not temp_response:
            raise Exception(sample_accession + " is not a valid EBI sample ID")
        with open(current_path + "/" + sample_accession + ".txt", 'wb') as f:
            f.write(temp_response.content)

    logger.info("generating sample file")
    # generates content for the sample file
    all_samp = []
    all_samp_exp = []
    # For each sample file containing metadata in xml format,
    # convert to dictionary format
    samples = []
    for fp in glob.glob('%s/*/*.txt' % study_accession):
        id_, _ = path.splitext(fp)
        a = id_.rfind("/")
        newid = id_[a+1:]
        if newid in samples:
            continue
        samples.append(newid)
        parsed = xml_to_dict(fp)
        parsed["sample_name"] = newid
        if '.' in id_:
            all_samp_exp.append(parsed)
        else:
            all_samp.append(parsed)

    all_samp = DataFrame(all_samp, dtype=str).fillna('Not provided')
    all_samp_exp = DataFrame(all_samp_exp, dtype=str).fillna('Not provided')

    result = (all_samp, all_samp_exp)

    newlist = []
    for item in result:
        header = item.columns[:, ].values.astype(str).tolist()
        values = item.values.tolist()
        values = [val+["EBI"] for val in values]
        listitem = [header] + values
        newlist.append(listitem)
    # C ontains header names
    header = newlist[0][0]
    final = []
    # The column names should not involve these characters
    # TODO: Construct a more exclusive list
    for i in header:
        temp = str(i).lower().replace(" ", "_").replace("-", "_")\
            .replace("(", "").replace(")", "").replace("/", "")
        final.append(temp)
    final.append("public_import_source")
    newlist[0][0] = final
    # Write to the sample file
    with open(study_accession + "_sample_info.txt", "w") as f:
        for items in newlist:
            wr = csv.writer(f, delimiter="\t")
            print(items)
            wr.writerows(items)


def sra_create_sample_file(sample_file, study_accession, study_details):
    """fetch the study information of each sample in SRA study
        and generate a sample file for these information

    Parameters
    ----------
    study_accession : string
        The accession ID of the SRA study

    sample_file : string
        The sample file name

    study_details : string
        The study detail file name
    """
    def xml_to_dict(xml_fp):
        # Converts xml string to a dictionary
        root = etree.parse(xml_fp).getroot()
        sample = root.getchildren()[0].find('SAMPLE')
        metadata = {}
        attributes = sample.find('SAMPLE_ATTRIBUTES')
        for node in attributes.iterfind('SAMPLE_ATTRIBUTE'):
            tag = node.getchildren()[0]
            value = node.getchildren()[1]
            if value.text is None:
                metadata[tag.text.strip('" ').upper()] = 'not provided'
            else:
                metadata[tag.text.strip('" ').upper()] \
                    = value.text.strip('" ')
					
		#adding loops to look for additional data
        title= sample.find('TITLE')
        if title.text is None:
            metadata['title'] = 'not provided'
        else:
            metadata['title'] = title.text.strip('" ')
        description = sample.find('DESCRIPTION')
        try:
            if description.text is None:
                metadata['description'] = 'not provided'
            else:
                if sep in description.text:
                    split_desc=description.text.strip('" ').split(sep)
                    counter=0
                    for i in split_desc:
                        metadata['description_field_' + str(counter)] = i
                        counter += 1
                else:
                    metadata['description'] = description.text.strip('" ')
        except:
            metadata['description'] = 'not provided'

        nameInfo = sample.find('SAMPLE_NAME')
        for node in nameInfo:
            tag = node.tag
            value = node.text
            if value is None:
                metadata[tag.text.strip('" ').upper()] = 'not provided'
            else:
                metadata[tag.strip('" ').upper()] = value.strip('" ')

        idInfo = sample.find('IDENTIFIERS')
        for node in idInfo:
            value = node.text
            d = node.attrib
            if len(d) > 0:
                for k in d.keys():
                    tag = node.tag + "_" + d[k]
                    if value is None:
                        metadata[tag.text.strip('" ').upper()] = 'not provided'
                    else:
                        metadata[tag.strip('" ').upper()] = value.strip('" ')
            else:
                tag = node.tag

                if value is None:
                    metadata[tag.text.strip('" ').upper()] = 'not provided'
                else:
                    metadata[tag.strip('" ').upper()] = value.strip('" ')

        return metadata

    logger.info("downloading sample.txt file for each sample")
    details_df = read_csv(study_details, sep='\t', header=None)
    for row in details_df.iterrows():
        library_name = row[1][0]
        current_path = "./" + study_accession + "/" + library_name
        run_accession = row[1][2]

        p1 = subprocess.Popen(['esearch', '-db', 'sra', '-query',
                               run_accession], stdout=subprocess.PIPE)
        for i in p1.stdout:
            if "<Count>0</Count>" in i.decode("utf-8"):
                p1.stdout.close()
                raise Exception(run_accession + " is not a valid SRA run ID")

        if path.exists(current_path + "/" + run_accession + ".txt"):
            continue
        if not path.exists(current_path):
            makedirs(current_path)

        p1 = subprocess.Popen(['esearch', '-db', 'sra', '-query',
                               run_accession], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['efetch', '-format', 'native'], stdin=p1.stdout,
                              stdout=subprocess.PIPE)
        p1.stdout.close()

        with open(current_path + "/" + run_accession + ".txt", 'wb') as f:
            for i in p2.stdout:
                f.write(i)
        p2.stdout.close()

    logger.info("generating sample file")
    # generates content for the sample file
    all_samp = []
    all_samp_exp = []
    # For each sample file containing metadata in xml format,
    # convert to dictionary format
    samples = []
    for fp in glob.glob('%s/*/*.txt' % study_accession):
        id_, _ = path.splitext(fp)
        a = id_.rfind("/")
        newid = id_[a+1:]
        if newid in samples:
            continue
        samples.append(newid)
        parsed = xml_to_dict(fp)
        parsed["sample_name"] = newid
        if '.' in id_:
            all_samp_exp.append(parsed)
        else:
            all_samp.append(parsed)

    all_samp = DataFrame(all_samp, dtype=str).fillna('Not provided')
    all_samp_exp = DataFrame(all_samp_exp, dtype=str).fillna('Not provided')

    result = (all_samp, all_samp_exp)

    newlist = []
    for item in result:
        header = item.columns[:, ].values.astype(str).tolist()
        values = item.values.tolist()
        values = [val+["NCBI"] for val in values]
        listitem = [header] + values
        newlist.append(listitem)

    # C ontains header names
    header = newlist[0][0]
    final = []
    # The column names should not involve these characters
    # TODO: Construct a more exclusive list
    for i in header:
        temp = str(i).lower().replace(" ", "_").replace("-", "_")\
            .replace("(", "").replace(")", "").replace("/", "")
        final.append(temp)
    final.append("public_import_source")
    newlist[0][0] = final
    # Write to the sample file
    with open(study_accession + "_sample_info.txt", "w") as f:
        for items in newlist:
            wr = csv.writer(f, delimiter="\t")
            wr.writerows(items)


def ebi_fetch_data_file(study_accession, study_details):
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
    # import the details as a dataframe
    details_df = read_csv(study_details, sep='\t', header=None)
    for row in details_df.iterrows():
        submitted_format = row[1][7]
        library_name = row[1][0]
        current_path = "./" + study_accession + "/" + library_name
        sample_accession = row[1][1]
        run_accession = row[1][2]
        fastq_ftp = row[1][4]
        if type(row[1][4]) is not str:
            logger.warning("No fastq ftp found for sample: " + sample_accession
                           + ", run: " + run_accession)
            logger.warning("Skipping sample: " + sample_accession + ", run: "
                           + run_accession)
            continue

        if not path.exists(current_path):
            makedirs(current_path)

        fastq_ftp = fastq_ftp.split(';')
        if isinstance(submitted_format, str):
            submitted_format = submitted_format.split(';')
        else:
            submitted_format = None
        for i in range(len(fastq_ftp)):
            # Check for the format(sff or fastq)
            fq_path = current_path + "/" + sample_accession + "." + \
                run_accession + "_R"
            if len(fastq_ftp) > 1:
                fq_path = fq_path + str(i + 1)
            if submitted_format is None and ".sff" in fastq_ftp[i]:
                fq_path = fq_path + ".sff"
            elif submitted_format is None and ".fastq.gz" in fastq_ftp[i]:
                fq_path = fq_path + ".fastq.gz"
            elif len(submitted_format) < i and submitted_format[i] == "SFF":
                fq_path = fq_path + ".sff"
            else:
                fq_path = fq_path + ".fastq.gz"
            if path.isfile(fq_path):
                logger.warning("Skipping " + fq_path)
                logger.warning("File exists")
                continue
            urlretrieve("ftp://" + fastq_ftp[i], fq_path)


def sra_fetch_data_file(study_details):
    """Fetch all the fastq file(s) for SRA study

    Parameters
    ----------
    study_details : string
        study detail file name

    """

    logger.info("Downloading the fastqs")
    # Download the fastqs
    # import the details as a dataframe
    details_df = read_csv(study_details, sep='\t', header=None)
    for row in details_df.iterrows():
        subprocess.run(['fastq-dump', '-I', '--split-files', row[1][2]])


if __name__ == '__main__':

    # parse the flags and initialize output file names
    # TODO:Reword help info
    parser = ArgumentParser(description='Please note that the following ' +
                            'packages has to be installed for running this ' +
                            'script: 1)lxml 2)pandas 3)glob 4)csv 5)sys ' +
                            '6)urllib 7)argparse 8)requests 9)xmltodict ' +
                            '10)subprocess 11)bioconda 12)sra-tools 13)os' +
                            '14)entrez-direct')
    parser.add_argument("-ebi", "--ebiaccession", help="ebi accession " +
                        "number whose info is to be processed")
    parser.add_argument("-sra", "--sraaccession", help="sra accession " +
                        "number whose info is to be processed")
    parser.add_argument("-sample", "--sample_fileName", help="sample_file" +
                        "Name which will contain per sample information")
    parser.add_argument("-prep", "--prep_fileName", help="prep_fileName" +
                        " that contains info like: sample-name,prefix of" +
                        " fastqs/sffs,library-source(metagenomic/genomic/" +
                        "transcriptional),instrument_platform")
    parser.add_argument("-study", "--study_fileName", help="Study_file" +
                        " that contains study information")
    parser.add_argument("-debug", "--debug", action='store_true', help="Debug mode: don't " +
                        "download fastq files")
    parser.add_argument("-all-seqs", "--all_seqs", action='store_true', help="Accept " +
                        "all type of sequence samples")
    parser.add_argument("-all-platforms", "--all_platforms", action='store_true', help="Accept " +
                        "all platform samples")
    parser.add_argument("-sep","--sep",help="separator for description, default is ';' ")
    args = parser.parse_args()

    if args.ebiaccession is None and args.sraaccession is None:
        print("""
                python EBI_SRA_Downloader.py -ebi [accession]
                    Generate the study info, study detail, prep, and  sample
                    files for the entered EBI accession, and download the
                    FASTQ files.
                python EBI_SRA_Downloader.py -sra [accession]
                    Generate the study info, study detail, prep, and  sample
                    files for the entered SRA accession, and download the
                    FASTQ files.
                Optional flags:
                    -sample_info [sample_info_file_name]
                    -prep_info [prep_info_file_name]
                    -study_info [study_info_file_name]
                    -debug
                    -all-seqs
                    -all-platforms
                    -sep
               """)
        sys.exit(2)

    DEBUG = args.debug
    ALL_SEQS = args.all_seqs
    ALL_PLATFORMS = args.all_platforms

    if args.ebiaccession is not None:
        # Output file names
        sample_file_name = args.ebiaccession + "_sample_info.txt" \
            if args.sample_fileName is None else args.sample_fileName
        prep_file_name = args.ebiaccession + "_prep_info.txt" \
            if args.prep_fileName is None else args.prep_fileName
        study_file_name = args.ebiaccession + "_study_info.txt" \
            if args.study_fileName is None else args.study_fileName

        #set parser settings
        sep= ';' if args.sep is None else args.sep
        # Call create_details_file to generate .details.txt
        study_details = ebi_create_details_file(args.ebiaccession)

        # Call create_info_file to generate feed.xml file and study_info file
        ebi_create_info_file(study_file_name, args.ebiaccession)

        # Call create_prep_file to generate prep files (based on .detail.txt)
        create_prep_file(prep_file_name, study_details)

        # Call create_sample_file to generate sample file
        # (based on .detail.txt)
        ebi_create_sample_file(sample_file_name, args.ebiaccession,
                               study_details)
        if DEBUG:
            sys.exit()

        # Call fetch_data_file to download all the fastqs files
        # (based on .details)
        # Get metadata info(internally converts xml to tuple format)
        # ebi_fetch_data_file(args.ebiaccession, 'PRJEB26419_detail.txt')
        ebi_fetch_data_file(args.ebiaccession, study_details)

    if args.sraaccession is not None:
        # Output file names
        sample_file_name = args.sraaccession + "_sample_info.txt" \
            if args.sample_fileName is None else args.sample_fileName
        prep_file_name = args.sraaccession + "_prep_info.txt" \
            if args.prep_fileName is None else args.prep_fileName
        study_file_name = args.sraaccession + "_study_info.txt" \
            if args.study_fileName is None else args.study_fileName
        # Call create_details_file to generate .details.txt
        study_details = sra_create_details_file(args.sraaccession)
        # Call create_info_file to generate feed.xml file and study_info file
        sra_create_info_file(study_file_name, args.sraaccession)
        # Call create_prep_file to generate prep files (based on .detail.txt)
        create_prep_file(prep_file_name, study_details)
        # Call create_sample_file to generate sample file
        # (based on .detail.txt)
        sra_create_sample_file(sample_file_name, args.sraaccession,
                               study_details)
        if DEBUG:
            sys.exit()

        # Call fetch_data_file to download all the fastqs files
        # (based on .details)
        # Get metadata info(internally converts xml to tuple format)
        sra_fetch_data_file(study_details)
