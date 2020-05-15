#!/usr/bin/env python3
#
# To run the test cases:
#   python3 test_EBI_SRA_Downloader.py
#
# conda command to install all dependencies:
#   conda create -n ebi_sra_importer pandas requests entrez-direct sra-tools xmltodict lxml -c bioconda -c conda-forge -y
# pip command to install all dependencies:
#   pip install csv glob requests subprocess xmltodict sys lxml os urllib
#   pip install argparse pandas bioconda sra-tools entrez-direct
#
# We also need EBI_SRA_Downloader.py to be under the same directory
#
# libraries used
import glob
import unittest
from shutil import rmtree
from pandas import read_csv
from os import path, remove
from EBI_SRA_Downloader import (ebi_create_details_file,
                                sra_create_details_file, ebi_create_info_file,
                                sra_create_info_file, create_prep_file,
                                ebi_create_sample_file, sra_create_sample_file,
                                ebi_fetch_data_file, sra_fetch_data_file)


class EBI_SRA_DownloaderTests(unittest.TestCase):
    def setUp(self):
        pass

    def test_ebi_create_details_file(self):
        # invalid accession
        with self.assertRaises(Exception):
            ebi_create_details_file("0")
        # valid accession but no entry
        with self.assertRaises(Exception):
            ebi_create_details_file("ERP0")
        # valid accession but not from Illumia
        with self.assertRaises(Exception):
            ebi_create_details_file("ERP012804")
        # valid accession
        detail_file = ebi_create_details_file("PRJNA390646")
        self.assertTrue(path.exists(detail_file)
                        and path.getsize(detail_file) > 0)
        if path.exists(detail_file):
            remove(detail_file)

    def test_sra_create_details_file(self):
        # invalid accession
        with self.assertRaises(Exception):
            sra_create_details_file("0")
        # valid accession but not metagenomic
        with self.assertRaises(Exception):
            sra_create_details_file("PRJNA312948")
        # valid accession
        detail_file = sra_create_details_file("PRJNA545409")
        self.assertTrue(path.exists(detail_file)
                        and path.getsize(detail_file) > 0)
        if path.exists(detail_file):
            remove(detail_file)

    def test_ebi_create_info_file(self):
        # invalid accession
        info_file = "0_study_info.txt"
        ebi_create_info_file(info_file, "0")
        self.assertFalse(path.exists(info_file))
        # valid accession but no entry
        info_file = "PRJNA295847_study_info.txt"
        ebi_create_info_file(info_file, "PRJNA295847")
        self.assertFalse(path.exists(info_file))
        self.assertFalse(path.exists('feed.xml'))
        # valid accession
        info_file = "PRJNA390646_study_info.txt"
        ebi_create_info_file(info_file, "PRJNA390646")
        self.assertTrue(path.exists(info_file)
                        and path.getsize(info_file) > 0)
        if path.exists(info_file):
            remove(info_file)
        self.assertTrue(path.exists('feed.xml')
                        and path.getsize('feed.xml') > 0)
        if path.exists('feed.xml'):
            remove('feed.xml')

    def test_sra_create_info_file(self):
        # invalid accession
        info_file = "0_study_info.txt"
        with self.assertRaises(Exception):
            sra_create_info_file(info_file, "0")
        # valid accession
        info_file = "PRJNA545409_study_info.txt"
        sra_create_info_file(info_file, "PRJNA545409")
        self.assertTrue(path.exists(info_file)
                        and path.getsize(info_file) > 0)
        if path.exists(info_file):
            remove(info_file)

    def test_create_prep_file(self):
        # invalid study detail file
        prep_file = "PRJNA390646_prep_info.txt"
        with self.assertRaises(FileNotFoundError):
            create_prep_file(prep_file, "PRJNA390646_detail.txt")
        # valid study detail file
        detail_file = ebi_create_details_file("PRJNA390646")
        create_prep_file(prep_file, detail_file)
        prep_fname = [fp for fp in glob.glob('PRJNA390646_prep_info*.txt')]
        self.assertTrue(len(prep_fname) > 0)
        for fp in prep_fname:
            self.assertTrue(path.getsize(fp) > 0)
            if path.exists(fp):
                remove(fp)
        remove(detail_file)

        # invalid study detail filfale
        prep_file = "0_prep_info.txt"
        with self.assertRaises(FileNotFoundError):
            create_prep_file(prep_file, "0_detail.txt")
        # valid study detail file
        detail_file = sra_create_details_file("PRJNA545409")
        prep_file = "PRJNA545409_prep_info.txt"
        create_prep_file(prep_file, detail_file)
        prep_fname = [fp for fp in glob.glob('PRJNA545409_prep_info*.txt')]
        self.assertTrue(len(prep_fname) > 0)
        for fp in prep_fname:
            self.assertTrue(path.getsize(fp) > 0)
            if path.exists(fp):
                remove(fp)
        remove(detail_file)

    def test_ebi_create_sample_file(self):
        # invalid study detail file
        sample_file = "0_sample_info.txt"
        with self.assertRaises(FileNotFoundError):
            ebi_create_sample_file(sample_file, "0",
                                   "0_detail.txt")
        # valid study detail file
        sample_file = "ERP106503_sample_info.txt"
        detail_file = ebi_create_details_file("ERP106503")
        ebi_create_sample_file(sample_file, "ERP106503", detail_file)
        remove(detail_file)
        self.assertTrue(path.exists(sample_file)
                        and path.getsize(sample_file) > 0)
        if path.exists(sample_file):
            remove(sample_file)
        sample_fname = [fp for fp in glob.glob('ERP106503/*/*')]
        self.assertTrue(len(sample_fname) > 0)
        for fp in sample_fname:
            self.assertTrue(path.getsize(fp) > 0)
            if path.exists(fp):
                remove(fp)
        rmtree("ERP106503")

    def test_sra_create_sample_file(self):
        # invalid study detail file
        sample_file = "0_sample_info.txt"
        with self.assertRaises(FileNotFoundError):
            sra_create_sample_file(sample_file, "0",
                                   "0_detail.txt")
        # valid study detail file
        sample_file = "PRJNA545409_sample_info.txt"
        detail_file = sra_create_details_file("PRJNA545409")
        sra_create_sample_file(sample_file, "PRJNA545409", detail_file)
        remove(detail_file)
        self.assertTrue(path.exists(sample_file)
                        and path.getsize(sample_file) > 0)
        if path.exists(sample_file):
            remove(sample_file)
        sample_fname = [fp for fp in glob.glob('PRJNA545409/*/*')]
        self.assertTrue(len(sample_fname) > 0)
        for fp in sample_fname:
            self.assertTrue(path.getsize(fp) > 0)
            if path.exists(fp):
                remove(fp)
        rmtree("PRJNA545409")

    def test_ebi_fetch_data_file(self):
        # invalid study detail file
        with self.assertRaises(FileNotFoundError):
            ebi_fetch_data_file("ERP106503", "ERP106503_detail.txt")
        # valid study detail file
        detail_file = ebi_create_details_file("ERP106503")
        ebi_fetch_data_file("ERP106503", detail_file)
        remove(detail_file)
        fastq_fname = [fp for fp in glob.glob('ERP106503/*/*')]
        self.assertTrue(len(fastq_fname) > 0)
        for fp in fastq_fname:
            self.assertTrue(path.getsize(fp) > 0)
            if path.exists(fp):
                remove(fp)
        rmtree("ERP106503")

    def test_sra_fetch_data_file(self):
        # invalid study detail file
        with self.assertRaises(FileNotFoundError):
            sra_fetch_data_file("0_detail.txt")
        # valid study detail file
        detail_file = ebi_create_details_file("PRJNA533871")
        sra_fetch_data_file(detail_file)
        details_df = read_csv(detail_file, sep='\t', header=None)
        run_accession = []
        for row in details_df.iterrows():
            run_accession.append(row[1][2])
        remove(detail_file)
        fastq_fname = [fp for fp in glob.glob('*')
                       if fp[:fp.rfind("_")] in run_accession]
        self.assertTrue(len(fastq_fname) > 0)
        for fp in fastq_fname:
            self.assertTrue(path.getsize(fp) > 0)
            if path.exists(fp):
                remove(fp)


if __name__ == '__main__':
    unittest.main()
