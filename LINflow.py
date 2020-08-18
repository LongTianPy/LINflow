#!/usr/bin/python
"""
Usage:
python LINflow.py initiate ./test
python LINflow.py add_scheme ./test
python ./LINflow.py add_genomes $HOME/idea_projects/linbase/LINflow/test/ -i \
$HOME/idea_projects/linbase/LINflow/fasta/ \
-m $HOME/idea_projects/linbase/LINflow/fasta/meta_data2.meta -s 1
"""

# IMPORT
import os
import sys
import uuid
import argparse
import shutil
import pandas as pd
from os.path import isfile, isdir, join
from distutils.spawn import find_executable
import re
import numpy as np
import decimal
from pyani import pyani_files, anib
from pyani.pyani_config import FRAGSIZE, ALIGNDIR
from pyani import run_multiprocessing as run_mp

# used to create and manage fasta signatures (hashes)
from sourmash import MinHash, SourmashSignature, save_signatures
from sourmash_custom import load_dbs_and_sigs
from sourmash.search import search_databases

# used to read fasta files
import screed
import glob

from LINdb import connect_to_db
# used to sleep process
import time

# used for tree building
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

# used to create lock files
from filelock import Timeout, FileLock


# CLASSES
class LINbase(object):
    """
    Object to interact with the MySQL database
    """

    def __init__(self, database_name='LINflow'):
        self.conn, self.cursor = connect_to_db()
        self.db_name = database_name
        if database_exists(database_name):
            self.cursor.execute("USE `{}`".format(self.db_name))

    def get_cursor(self):
        return self.cursor

    def commit(self):
        self.conn.commit()

    def initialize(self):
        if database_exists(self.db_name):
            print("Database with the name {} exists please delete the database"
                  " or select a new name".format(self.db_name))

        self.cursor.execute("CREATE DATABASE `{}`".format(self.db_name))
        #self.cursor.execute(
        #    "CREATE TABLESPACE `{0}` ADD DATAFILE './{0}.ibd' Engine='InnoDB';".format(self.db_name))
        self.conn.commit()
        self.cursor.execute("USE `{}`".format(self.db_name))

        self.cursor.execute("CREATE TABLE Genome (Genome_ID INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,"
                            "FilePath VARCHAR(3000) NOT NULL)")
        self.cursor.execute("CREATE TABLE Scheme (Scheme_ID INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,"
                            "Cutoff MEDIUMTEXT NOT NULL,"
                            "LabelNum INT UNSIGNED NOT NULL,"
                            "Description TEXT)")
        self.cursor.execute("CREATE TABLE ANI (ANI_ID INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,"
                            "Genome_ID INT UNSIGNED NOT NULL,"
                            "SubjectGenome INT UNSIGNED NOT NULL,"
                            "ANI DOUBLE NOT NULL)")
        self.cursor.execute("CREATE TABLE LIN (LIN_ID INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,"
                            "Genome_ID INT UNSIGNED NOT NULL,"
                            "Scheme_ID INT UNSIGNED NOT NULL,"
                            "LIN MEDIUMTEXT NOT NULL)")
        self.cursor.execute("CREATE TABLE Taxonomy (Taxonomy_ID INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,"
                            "Genome_ID INT UNSIGNED NOT NULL,"
                            "Genus VARCHAR(255) NOT NULL,"
                            "Species VARCHAR(255) NOT NULL,"
                            "Strain VARCHAR(255) NOT NULL)")
        self.add_scheme(
            "70,75,80,85,90,95,96,97,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999",
            20, "The default scheme of LINbase.org")

        fine_scheme = ','.join([str(i / float(10)) for i in range(700, 1000, 1)])
        self.add_scheme(fine_scheme, 300, "A scheme used to approximate the tree")
        # self.conn.commit()
        fine_scheme = ','.join([str(i / float(100)) for i in range(7000, 10000, 1)])
        self.add_scheme(fine_scheme, 3000, "A scheme used fo a precise tree")
        self.conn.commit()

    # obj is of type getLIN
    def add_lin_async(self, obj, genome_id):
        # genome_id, scheme_id, num_lins, conserved_lin
        tail_length = obj.label_num - len(obj.conserved_LIN.split(','))

        # if there is a tail
        if obj.idx_to_change == -1:
            self.add_lin(genome_id, obj.Scheme_ID, obj.conserved_LIN)
            return
        elif obj.idx_to_change == 0:
            head = ''
            tail = "".join([',0'] * (obj.label_num - 1))
            wildcard = "%" + "".join([',0'] * (obj.label_num - 1))
        else:
            tail = "".join([',0'] * (tail_length - 1))
            wildcard = '%'
            head = obj.conserved_LIN + ","

        query = "INSERT INTO LIN (Genome_ID,Scheme_ID,LIN) " \
                "SELECT {0}, {1}, CONCAT(\"{2}\", " \
                "MAX(SUBSTRING_INDEX(SUBSTRING_INDEX(LIN,',',{3}),',',1)) + 1, \"{4}\") " \
                "FROM LIN WHERE Scheme_ID={1} AND LIN LIKE \"{2}{5}\"".format(
            genome_id, obj.Scheme_ID, head, -1 * tail_length, tail, wildcard)

        try:
            self.cursor.execute(query)
        except Exception as e:
            print("SQL ERROR: " + query + "\n")
            print(e)

    def add_genome(self, filepath: str) -> int:
        self.cursor.execute("INSERT INTO Genome (FilePath) VALUES ('{0}')".format(filepath))
        return self.cursor.lastrowid

    def add_scheme(self, LIN_percents, number_of_labels, description) -> int:
        self.cursor.execute(
            "INSERT INTO Scheme (Cutoff, LabelNum, Description) VALUES "
            "('{0}', {1}, '{2}')".format(LIN_percents, number_of_labels, description)
        )
        return self.cursor.lastrowid

    def add_lin(self, genome_id, scheme_id, lin):
        self.cursor.execute(
            "INSERT INTO LIN (Genome_ID,Scheme_ID,LIN) VALUES ({0},{1},'{2}')".format(
                genome_id, scheme_id, lin)
        )

    def add_taxonomy(self, genome_id, strain, genus='', species=''):
        self.cursor.execute(
            "INSERT INTO Taxonomy (Genome_ID,Genus,Species,Strain) VALUES "
            "({0},'{1}', '{2}','{3}')".format(genome_id, genus, species, strain)
        )

    def add_ani(self, genome_id, subject_genome, ani):
        self.cursor.execute(
            "INSERT INTO ANI (Genome_ID, SubjectGenome, ANI) VALUES ({0},{1},{2})".format(
                genome_id, subject_genome, ani)
        )

    def get_genome_path(self, genome_id) -> str:
        self.cursor.execute(
            "SELECT FilePath FROM Genome WHERE Genome_ID = {0}".format(genome_id)
        )
        return self.cursor.fetchone()[0]

    def get_genomeID(self, filename):
        self.cursor.execute(
            "SELECT Genome_ID FROM Genome WHERE FilePath like '{0}'".format(filename)
        )
        results = self.cursor.fetchone()
        if results is not None:
            return int(results[0])
        return None

    def get_LIN(self, genome_id, scheme_id):
        self.cursor.execute(
            "SELECT LIN FROM LIN WHERE Genome_ID={0} AND Scheme_ID={1}".format(genome_id, scheme_id)
        )
        results = self.cursor.fetchone()
        if results is not None and len(results) > 0:
            return results[0]
        return None

    def len_scheme(self, scheme_id) -> int:
        self.cursor.execute(
            "SELECT LabelNum FROM Scheme WHERE SCHEME_ID={0}".format(scheme_id)
        )
        return self.cursor.fetchone()[0]

    def query_first(self, query: str):
        self.cursor.execute(query)
        return self.cursor.fetchone()

    def query_all(self, query: str):
        self.cursor.execute(query)
        return self.cursor.fetchall()

    def close(self):
        self.conn.close()


class getLIN(object):
    """
    Read the LIN of the top similar genome.
    """

    def __init__(self, Genome_ID, Scheme_ID, similarity, db: LINbase):
        """
        Initialize a LIN object from the database by parses the returned DB content

        Parameters
        ----------
        Genome_ID : str
            ID of the Genome
        Scheme_ID : str
            Represents the used LIN scheme which describes the LIN's similarity rating
        similarity : float
            Similarity measure for the LIN (between 0 and 1) 0 being dissimilar and 1 being identical
        c : SQLlite database connection handle
            Established database connection handle used to find the LIN in
        """
        self.Genome_ID = Genome_ID
        self.Scheme_ID = Scheme_ID
        self.db = db
        tmp = db.query_first("SELECT LabelNum from Scheme WHERE Scheme_ID={0}".format(Scheme_ID))
        self.label_num = int(tmp[0])
        self.similarity = float(similarity) * 100
        self.idx_to_change = None
        self.conserved_LIN = None
        self.parse()

    def parse(self, Genome_ID=None, Scheme_ID=None, similarity=None, db=None):
        """
        Determine highest LIN to change based on the given

        Parameters
        ----------
        Genome_ID : str, default (uses object field)
            ID of the referenced genome
        Scheme_ID : str, default (uses object field)
            Represents the used LIN scheme which describes the LIN's similarity rating
        similarity : float, default (uses object field)
            Similarity measure for the LIN (between 0 and 1) 0 being dissimilar and 1 being identical
        db : LINbase
            Linbase object for handling database interaction
        """
        if not Genome_ID:
            Genome_ID = self.Genome_ID
        if not Scheme_ID:
            Scheme_ID = self.Scheme_ID
        if not similarity:
            similarity = self.similarity
        if not db:
            db = self.db
        # Read the LIN of this genome
        lin = db.get_LIN(int(Genome_ID), int(Scheme_ID))
        self.LIN = lin.split(",")
        # Read the cutoff (LIN representation) of this scheme
        result = db.query_first('SELECT Cutoff from Scheme where Scheme_ID={0}'.format(Scheme_ID))
        cutoff = [float(i) for i in result[0].split(',')]

        if similarity < cutoff[0]:
            self.idx_to_change = 0
            self.conserved_LIN = ''
        elif similarity >= cutoff[-1]:
            self.idx_to_change = -1
            self.conserved_LIN = lin
        else:
            for i in range(len(cutoff) - 1):
                if cutoff[i] <= similarity < cutoff[i + 1]:
                    self.idx_to_change = i + 1
                    self.conserved_LIN = ",".join(self.LIN[:i + 1])
                    break

        '''self.idx_to_change = idx_to_change
        if idx_to_change == 0:
            self.conserved_LIN = ''
        elif idx_to_change == 'n/a':
            self.conserved_LIN = lin
        else:
            self.conserved_LIN = lin[:idx_to_change]'''


class Assign_LIN(object):
    """
    Get the biggest number assigned to the idx_to_change with the same conserved part of LIN
    """

    def __init__(self, getLIN_object, db):
        self.idx_to_change = getLIN_object.idx_to_change
        self.conserved_LIN = getLIN_object.conserved_LIN
        self.label_num = getLIN_object.label_num
        self.scheme_id = getLIN_object.Scheme_ID
        self.db = db
        self.new_LIN = None
        self.assign()

    def assign(self, idx_to_change=None, conserved_LIN=None, label_num=None, db=None, scheme_id=None):
        """
            Generate a new LIN number based on the current LINs

            Parameters
            ----------
            idx_to_change : int, default (uses object field)
                LIN index where the query is different from te reference
            conserved_LIN : str, default (uses object field)
                LIN value where the query is different from the reference
            label_num : int, default (uses object field)
                Number of LINs in the schema
            db : LINbase
                Linbase object for handling database interaction
            scheme_id : str, default (uses object field)
                LIN scheme ID
        """
        if not idx_to_change:
            idx_to_change = self.idx_to_change
        if not conserved_LIN:
            conserved_LIN = self.conserved_LIN
        if not label_num:
            label_num = self.label_num
        if not db:
            db = self.db
        if not scheme_id:
            scheme_id = self.scheme_id

        # Find representative LINs to compare with based on LIN index
        query = ""
        if not idx_to_change:
            # Get all genomes when similarity < first scheme base
            query = "SELECT LIN.LIN FROM LIN where Scheme_ID={0}".format(scheme_id)
        else:
            # Get genomes in the LINgroup
            query = "SELECT LIN.LIN from LIN WHERE LIN.LIN LIKE '{0},%' and Scheme_ID={1}".format(conserved_LIN,
                                                                                                  scheme_id)
        representatives = db.query_all(query)

        if idx_to_change >= 0:
            # If the index is already determined generate a new LIN number in the group
            # (last assigned LIN in th group + 1)
            LINs = [int(i[0].split(',')[idx_to_change]) for i in representatives]
            num_to_assign = str(max(LINs) + 1)

        if idx_to_change >= 0 and idx_to_change != label_num - 1:
            # If the index is not the last LIN position (last similarity metric) aka normal case
            tail = ['0'] * (label_num - 1 - idx_to_change)
            tail = ','.join(tail)
            new_LIN = conserved_LIN + ',%s,' % num_to_assign + tail
        elif idx_to_change < 0:
            # Too high for LINs to differentiate then same as closest LIN
            new_LIN = conserved_LIN
        else:
            # When the last LIN can differentiate the genomes just change last LIN number
            new_LIN = conserved_LIN + ',%s,' % num_to_assign

        self.new_LIN = new_LIN.strip(',')


def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="LINflow"
    )
    parser.add_argument(
        "function",
        type=str,
        choices=[
            'initiate',
            'show_schemes',
            'add_scheme',
            'add_genomes',
            'LIN_distance',
            'build_trees',
            'cleanup'])
    parser.add_argument(
        "workspace",
        type=str,
        help="The location of the workspace")
    parser.add_argument(
        "-s",
        dest="Scheme_ID",
        type=int,
        default=0,
        help="The Scheme based on which LINs are going to be assigned.", )
    parser.add_argument(
        "-i",
        dest="input_dir",
        default='',
        help="The directory of genomes going to be added.")
    parser.add_argument(
        "-m",
        dest="metadata",
        default='',
        help="The metadata corresponding to the genomes. "
             "Download the sample in https://bit.ly/2Y6Pw3R, and save as CSV "
             "(comma separated values) format file.")
    parser.add_argument(
        "-d",
        dest="df",
        default='',
        help="The file name of the ANI matrix.")
    parser.add_argument(
        "-l",
        dest='lingroup',
        type=str,
        default='',
        help="The LINgroup of genomes to show distance. Default: '', to see all.")
    parser.add_argument(
        "-t",
        dest='taxonomy',
        type=str,
        default='',
        help="The taxon of selection. Default: ''.")
    parser.add_argument(
        "-g",
        dest='genome_labels',
        type=int,
        default=7,
        help="Genome labels included in distance matrix 1(Genome) 2(Species)"
             "4(Strain). Add numbers to combine. Default: 7(All) 1+2+4=7")

    # parser.add_argument("-p", dest="privacy", help="Is it private information")
    args = parser.parse_args()
    return args


def database_exists(db_name) -> bool:
    conn, cursor = connect_to_db()
    cursor.execute("SELECT SCHEMA_NAME FROM information_schema.schemata "
                   "where SCHEMA_NAME =  '{}'".format(db_name))
    db_exists = cursor.fetchall()
    if len(db_exists) > 0:
        return True
    return False


def initiate(workspace) -> LINbase:
    """
    Initializes an empty directory with the database and default schema

    Parameters
    ----------
    workspace : str
        location of the workspace new or update old workspace
    """
    if os.path.isdir(workspace):
        print("Provided workspace is not empty, please use a clean environment")
        exit(-1)
    else:
        os.mkdir(workspace)
        os.chdir(workspace)
        # mysql_uid = getpwnam('mysql')[2]
        # os.chown('./', mysql_uid, mysql_uid)
        # os.chmod('./data', 775)
        db = LINbase(os.path.basename(workspace))
        db.initialize()
        for name in ['Genomes', "ANI", 'Signatures']:
            os.mkdir(name)
        os.chdir('Signatures')
        for name in ['All', "tmp_sig", 'tmp_result', 'rep_bac']:
            os.mkdir(name)
        return db


def create_sketch(filepath, dest):
    """
    Creates a sourmash signature

    Parameters
    ----------
    filepath : str
        path to the FASTA file
    dest : str
        path to save sketch file
    """
    signatures = []
    for k in [21, 51]:
        minhash = MinHash(n=2000, ksize=k)
        for read in screed.open(filepath):
            minhash.add_sequence(read.sequence, True)
        signatures.append(SourmashSignature(minhash, filename=os.path.abspath(filepath)))
    with open(dest, "w+") as signature_file:
        save_signatures(signatures, signature_file)
    # cmd = "sourmash compute -o {0} {1} -k 21,51 -n 2000 -q".format(dest, filepath)
    return signatures


def compare_sketch(query_sigs, lingroup, k) -> pd.DataFrame:
    """
    Compare signatures using sourmash

    Parameters
    ----------
    query_sigs : List of SourmashSignature
        list of sourmash signatures
    lingroup : str
        determined LINgroup membership to save under the result under the appropriate path
    k : int
        window size for comparing signatures

    Returns
    -------
     df_list  : pandas DataFrame
        Sourmash result format
        File header [similarity,name,filename,md5] e.g. 0.xxx,/path/to/query/zzz.fasta,
        path/to/workspace/Signatures/rep_back/1.sig,md5 of query file
    """
    if lingroup == "rep_bac":
        dest = rep_bac_dir
    else:
        dest = join(sourmash_dir, lingroup)
    folder_size = len([file for file in os.listdir(dest) if isfile(join(dest, file))])

    # get the appropriate k-mer signature
    query = None
    for sig in query_sigs:
        if sig.minhash.ksize == int(k):
            query = sig

    all_sigs = glob.glob(dest + "/*.sig")
    registered_sigs = [file for file in all_sigs if re.match(r'.*\/[0-9]+.sig$', file)]
    new_sigs = np.setdiff1d(np.array(all_sigs), np.array(registered_sigs)).tolist()
    while True:
        try:
            ref_signatures = load_dbs_and_sigs(registered_sigs + new_sigs, query,
                                               is_similarity_query=False, traverse=False)
        # this is thrown when a file in the list is removed from the directory during loading
        except (EnvironmentError, ValueError) as e:
            continue
        break

    results = search_databases(query, ref_signatures, threshold=0.0001,
                               do_containment=False, best_only=False, ignore_abundance=False)
    # cmd = "sourmash search {0} {1} -n {2} -k {4} -q --threshold 0.0001 -o {3}"
    # cmd = cmd.format(query, join(dest, '*.sig'), folder_size, join(sourmash_result, "tmp_result.txt"), k)

    # Organize results into a dataframe
    df_list = []
    index = []
    for result in results:
        record = dict()

        record['similarity'] = result.similarity
        record['name'] = result.name
        record['filename'] = result.filename
        record['md5'] = result.md5
        # should an int be appended?
        index.append(result.filename.split("/")[-1].split(".")[0])
        df_list.append(record)
    return pd.DataFrame(df_list, index)


def parse_result(result_file) -> pd.DataFrame:
    """
    Parse Sourmash result file

    Parameters
    ----------
    result_file : str
        path to the Sourmash result file

    Returns
    -------
     df : Pandas Dataframe
        Sourmash results in Dataframe format
    """

    df = pd.read_csv(result_file, sep=",", header=0)
    if df.empty:
        return df
    else:
        ids = []
        for each in df['filename']:
            id = int(each.split('/')[-1].split('.')[0])
            ids.append(id)
        df.index = ids
        return df


def add_first_genome(filename, taxonomy, target_filename, scheme_id, db):
    """
    Add first genome to the database

    Parameters
    ----------
    filename : str
        path to the query FASTA file
    taxonomy : dict
        dictionary of taxonomy data e.g. genome name, strain name
    target_filename : str
        path to the query FASTA file's signature
    scheme_id : str
        LIN scheme id
    db : LINbase
        LIN database object
    """
    db.add_genome(target_filename)
    shutil.copy(filename, target_filename)
    db.add_taxonomy(
        genome_id=1, genus=taxonomy['genus'], species=taxonomy['species'], strain=taxonomy['strain'])
    db.add_ani(1, 1, 1)
    num = db.query_first("select count(LabelNum) FROM Scheme")[0]
    # num = int(c.fetchone()[0])
    for scheme_num in range(1, num + 1):
        zeros = db.len_scheme(scheme_num)
        db.add_lin(1, scheme_num, ",".join(["0"] * zeros))

    LIN_dir = ",".join(["0"] * 6)
    os.mkdir(join(sourmash_dir, LIN_dir))
    sig_path = join(sourmash_dir, LIN_dir, "1.sig")
    create_sketch(filename, sig_path)
    os.symlink(sig_path, join(rep_bac_dir, "1.sig"), False)
    # shutil.copy(filename, target_filename)


def add_genome(filename, taxonomy, target_filename, scheme_id, db):
    """
    Add a genome to the database

    Parameters
    ----------
    filename : str
        path to the query FASTA file
    target_filename : str
        path to the query FASTA file
    taxonomy : dict
        dictionary of taxonomy data e.g. genome name, strain name
    scheme_id : str
        LIN scheme id
    db : LINbase
        LIN database object
    """
    # Need to remove the file if it is a duplicate (100% ANI)
    shutil.copy(filename, target_filename)
    tmp_filename = os.path.splitext(os.path.basename(target_filename))[0]
    tmp_sig = create_sketch(filename, join(rep_bac_dir, tmp_filename + ".sig"))
    tmp_sig_path = join(rep_bac_dir, tmp_filename + ".sig")

    while True:
        df: pd.DataFrame = compare_sketch(tmp_sig, "rep_bac", 21)
        ANIb_result = 0.6
        SubjectGenome = 1

        if tmp_filename in df.index:
            df.drop(tmp_filename, inplace=True)

        if not df.empty:
            rep_bac_Genome_ID = None
            # if Subject genome is not already registered (has LIN and genomeID)
            # removes similar parallel x <=> y deadlock where each waits for the other
            if not str(df.index[0]).isnumeric():
                query_fasta_file_age = os.path.getctime(target_filename)
                # find dependant
                for idx in range(0, len(df)):
                    match_sig = df['filename'][idx]
                    # best match has a genomeID use it
                    if df.index[idx].isnumeric():
                        rep_bac_Genome_ID = int(df.index[idx])
                        break
                    subject_file_path = join(
                        genome_dir, os.path.basename(match_sig).replace('.sig', '.fasta'))
                    subject_file_age = os.path.getctime(subject_file_path)
                    if subject_file_age < query_fasta_file_age:
                        while rep_bac_Genome_ID is None:
                            rep_bac_Genome_ID = db.get_genomeID(subject_file_path)
                            time.sleep(10)
                            print("stuck in genomeID " + subject_file_path)

            if rep_bac_Genome_ID is None:
                rep_bac_Genome_ID = int(df.index[0])

            while True:
                rep_bac_LIN = db.get_LIN(rep_bac_Genome_ID, 1)
                # if the SubjectGenome has a LIN continue onward
                if rep_bac_LIN is not None:
                    break
                time.sleep(10)
                print("stuck in LIN")

            rep_bac_LINgroup = ",".join(rep_bac_LIN.split(",")[:6])
            jaccard_similarity = df.loc[str(rep_bac_Genome_ID), "similarity"]

            # increase sensitivity (k= 21 -> 51) since the similarity is high
            if jaccard_similarity > 0.2475:
                df = compare_sketch(tmp_sig, rep_bac_LINgroup, 51)
            elif 0.2475 >= jaccard_similarity > 0.0025:
                df = compare_sketch(tmp_sig, rep_bac_LINgroup, 21)

            ANIb_result = -1
            SubjectGenome = ""
            for each_subject_genome_ID in df.index[:3]:
                subject_genome_file = db.get_genome_path(int(each_subject_genome_ID))
                this_ANIb_result = pyani_pairwise_adapter(subject_genome_file, filename)

                #  find the genome which has the highest ANI value
                if this_ANIb_result > ANIb_result:
                    ANIb_result = this_ANIb_result
                    SubjectGenome = int(each_subject_genome_ID)

        genome_id = db.add_genome(target_filename)
        new_LIN_object = getLIN(Genome_ID=SubjectGenome, Scheme_ID=scheme_id, similarity=ANIb_result, db=db)
        db.add_lin_async(new_LIN_object, genome_id)
        db.add_ani(genome_id, SubjectGenome, ANIb_result)
        db.add_taxonomy(genome_id, genus=taxonomy['genus'], species=taxonomy['species'], strain=taxonomy['strain'])

        # Technically schemes to update
        schemes_to_update = np.setdiff1d(np.array([1, 2, 3]), np.array([scheme_id]))
        for scheme_num in schemes_to_update:
            new_LIN_object_other = getLIN(
                Genome_ID=SubjectGenome, Scheme_ID=scheme_num, similarity=ANIb_result, db=db)
            db.add_lin_async(new_LIN_object_other, genome_id)
            # db.add_lin(genome_id, scheme_num, new_LIN_other)
        # Update signature location
        while True:
            LIN_1 = db.get_LIN(genome_id, 1)
            # if the SubjectGenome has a LIN continue onward
            if LIN_1 is not None:
                break
            time.sleep(10)
            print("stuck in LIN2")

        lingroup = ",".join(LIN_1.split(",")[:6])

        if not isdir(join(sourmash_dir, lingroup)):
            waited_redo = False
            while True:
                try:
                    rep_lock.acquire()
                    break
                except Timeout:
                    time.sleep(10)
                    waited_redo = True
            if waited_redo:
                continue
            os.mkdir(join(sourmash_dir, lingroup))
            shutil.copyfile(tmp_sig_path, join(sourmash_dir, lingroup, "{0}.sig".format(str(genome_id))))
            # change this to symbolic link
            os.rename(tmp_sig_path, join(rep_bac_dir, "{0}.sig".format(str(genome_id))))
            rep_lock.release()
        else:
            shutil.copyfile(tmp_sig_path, join(sourmash_dir, lingroup, "{0}.sig".format(str(genome_id))))
            os.remove(tmp_sig_path)
        break


def show_schemes(workspace):
    os.chdir(workspace)
    project_name = os.path.basename(workspace.strip('/'))
    db = LINbase(project_name)
    tmp = db.query_all('SELECT Scheme_ID,Cutoff,Description FROM Scheme')
    print('Scheme_ID\tDescription\tScheme')
    for i in tmp:
        print('{0}\t{1}\t{2}'.format(i[0], i[2], i[1]))
    db.close()


def add_scheme(workspace):
    os.chdir(workspace)
    project_name = os.path.basename(workspace.strip('/'))
    db = LINbase(project_name)
    scheme = input(
        "Please enter a new scheme delimited by comma (,)"
        + " without the percentage mark and do not include 100 at the end, e.g. 70,80,90,99"
        + "\nOR\nin the format MIN,STEP,MAX:")
    scheme_list = scheme.split(",")
    scheme_list = [float(num.strip()) for num in scheme_list]

    if len(scheme_list) == 3 and scheme_list[1] < scheme_list[0]:
        scheme_arr = np.arange(scheme_list[0], scheme_list[2] + scheme_list[1], scheme_list[1])
        decimal_places = abs(decimal.Decimal(str(scheme_list[1])).as_tuple().exponent)
        scheme_arr = [round(num, decimal_places) for num in scheme_arr.tolist()]
        scheme = re.sub(r'\n?\s+', '', str(scheme_arr))[1:-1]
    elif len(scheme_list) == 3 and scheme_list[1] >= scheme_list[0]:
        print("Numbers should be ordered or be in the format \"Min,Step,Max\"")
        sys.exit()
    else:
        scheme = ",".join(str(num) for num in scheme_list)

    length = len(scheme.split(","))
    description = input("Please enter a short description (optional):")
    scheme_id = db.add_scheme(scheme, length, description)
    # c.execute("SELECT Scheme_ID,Description,Cutoff FROM Scheme WHERE Scheme_ID=(SELECT max(Scheme_ID) from Scheme)")
    # tmp = c.fetchone()
    print('Scheme_ID\tDescription\tCutoff')
    print('{0}\t{1}\t{2}'.format(scheme_id, description, length))


def check_arguments(args):
    if args.input_dir == '':
        print("Please provide a directory with genomes to be added with -i")
        exit(1)
    if args.metadata == '':
        print("Please provide the metadata file with -m")
        exit(2)
    if args.Scheme_ID <= 0:
        print("Please select one Scheme_ID with -s")
        exit(3)
    if not find_executable("sourmash"):
        print("Please install sourmash first")
        exit(4)
    if not find_executable("average_nucleotide_identity.py"):
        print("Please install pyani first")
        exit(5)

    return args.input_dir, args.metadata, args.Scheme_ID


def pyani_pairwise_adapter(ref_file_path, query_file_path):
    # ani_path = find_executable("average_nucleotide_identity.py")
    # sys.path.append(os.path.dirname(ani_path))
    # from average_nucleotide_identity import unified_anib
    results = -1.0
    infiles = [ref_file_path, query_file_path]
    infile_lengths = pyani_files.get_sequence_lengths(infiles)

    task_working_dir = join(os.getcwd(), str(uuid.uuid4()))
    blastdir = os.path.join(task_working_dir, ALIGNDIR['ANIb'])
    if not isdir(task_working_dir):
        os.mkdir(task_working_dir)
        os.mkdir(blastdir)
    # create fragmented files
    fragfiles, fraglengths = anib.fragment_fasta_files(infiles, blastdir, FRAGSIZE)
    # create BLAST jobs
    jobgraph = anib.make_job_graph(infiles, fragfiles,
                                   anib.make_blastcmd_builder('ANIb', blastdir))
    # run BLAST jobs in multiprocessing mode and return the cumulative exit status
    cumval = run_mp.run_dependency_graph(jobgraph)
    if cumval:
        raise Exception("At least one BLAST run failed.")

    # process BLAST data
    try:
        data = anib.process_blast(blastdir, infile_lengths, fraglengths=FRAGSIZE, mode="ANIb")
        for dfr, filestem in data.data:
            if "ANIb_percentage_identity" in filestem:
                results = dfr.iloc[0, 1]
                break
    except ZeroDivisionError:
        raise Exception("ANI calculation failed between reference: \"{}\" and query \"{}\"".format(
            ref_file_path, query_file_path))
    shutil.rmtree(task_working_dir)
    return results


# delete generated files + folders + database
def cleanup(workspace):
    project_folder = os.path.abspath(workspace)
    dbname = os.path.basename(workspace.strip('/'))
    if os.path.isdir(project_folder):
        shutil.rmtree(project_folder)
    conn, c = connect_to_db()
    c.execute('DROP DATABASE IF EXISTS {}'.format(dbname))
    try:
        c.execute('DROP TABLESPACE {}'.format(dbname))
    except Exception as e:
        pass


def rename_clades(clade):
    if not clade.clades:
        return
    else:
        clade.name = None
        for c in clade.clades:
            rename_clades(c)


def save_tree(tree, output_file):
    # You can change the Tree figure sizes here figsize=(width,height) default 20,30
    fig, ax = plt.subplots(figsize=(80, 90))
    Phylo.draw(tree, axes=ax, branch_labels=None, do_show=False)
    xlim = ax.get_xlim()
    xmax = 5 * (xlim[1] - xlim[0]) / 4
    ax.set_xlim(right=xmax)

    plt.rc('font', size=18)
    plt.rc('axes', titlesize=23)
    plt.rc('axes', labelsize=23)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.rc('legend', fontsize=23)
    plt.rc('figure', titlesize=23)
    plt.savefig(output_file)


# MAIN
if __name__ == '__main__':
    argv = sys.argv
    args = get_parsed_args()
    method = args.function
    workspace = os.path.abspath(args.workspace)
    project_name = os.path.basename(workspace.strip('/'))
    if method == 'initiate':

        if database_exists(project_name):
            print("ERROR: Project folder names should be unique.\n"
                  "A Database/Project with the name {} exists".format(project_name))
            exit(-3)
        initiate(project_name)
        exit(0)
    if method == 'cleanup':
        cleanup(workspace)
        exit(0)

    if not database_exists(project_name):
        print("Please initiate first")
        exit(6)

    genome_dir = join(workspace, "Genomes")
    ani_dir = join(workspace, "ANI")
    rep_bac_dir = join(workspace, "Signatures", "rep_bac")
    sourmash_dir = join(workspace, "Signatures", "All")
    tmp_sig_dir = join(workspace, "Signatures", "tmp_sig")
    sourmash_result = join(workspace, "Signatures", "tmp_result")
    rep_lock_path = join(workspace, "rep_lock.lock")
    rep_lock = FileLock(rep_lock_path)

    if method == 'show_schemes':
        show_schemes(workspace)

    elif method == 'add_scheme':
        add_scheme(workspace)

    elif method == 'add_genomes':
        input_dir, meta = [os.path.abspath(p) for p in check_arguments(args)[:2]]
        scheme_id = check_arguments(args)[2]
        df = pd.read_csv(meta, sep=",", header=0, index_col=0)
        os.chdir(workspace)
        db = LINbase(project_name)
        for i in df.index:
            filename = join(input_dir, str(i))
            uuid_filename = str(uuid.uuid4()) + ".fasta"
            target_filename = join(genome_dir, uuid_filename)
            genus = df.loc[i, "genus"]
            species = df.loc[i, "species"]
            strain = df.loc[i, "strain"]
            taxonomy = {"genus": genus, "species": species, "strain": strain}
            first_flag = False
            # size = c.fetchone()[0]
            size = db.query_first("select COUNT(Genome_ID) from Genome")[0]
            while size == 0:
                size = db.query_first("select COUNT(Genome_ID) from Genome")[0]
                try:
                    if size == 0 and rep_lock.acquire(timeout=1):
                        add_first_genome(
                            filename, taxonomy, target_filename, scheme_id, db)
                        rep_lock.release()
                        first_flag = True
                    else:
                        time.sleep(30)
                except Timeout:
                    time.sleep(30)
                break

            size = db.query_first("select COUNT(Genome_ID) from Genome")[0]
            if size > 0 and not first_flag:
                add_genome(filename, taxonomy, target_filename, scheme_id, db)
                print("Added " + filename)
        db.close()

    elif method == 'LIN_distance':
        if args.df == '':
            print("Please provide the name of the output distance matrix.")
            exit(7)

        df = os.path.abspath(args.df)
        lingroup = args.lingroup
        taxonomy = args.taxonomy
        os.chdir(workspace)
        db = LINbase(project_name)
        scheme_id = 3
        if args.Scheme_ID:
            scheme_id = args.Scheme_ID

        result = db.query_first("SELECT Cutoff FROM Scheme WHERE Scheme_ID={}".format(scheme_id))[0]
        scheme = result.split(',')
        scheme = [float(i) / 100 for i in scheme]
        sql = "SELECT LIN.Genome_ID,Taxonomy.Genus,Taxonomy.Species, Taxonomy.Strain,LIN.LIN " \
              "FROM LIN,Taxonomy,Genome " \
              "WHERE LIN.Genome_ID=Taxonomy.Genome_ID AND Genome.Genome_ID=LIN.Genome_ID AND " \
              "LIN.Scheme_ID={}".format(scheme_id)

        sql_extra = ""
        if lingroup != '':
            sql_extra = "AND LIN.Genome_ID IN (SELECT Genome_ID FROM LIN WHERE LIN LIKE '{},%')".format(lingroup)
        if taxonomy != '':
            taxonomy_list = taxonomy.split(' ')
            genus = taxonomy_list[0]
            sql_extra += "AND Taxonomy.Genus='{}' ".format(genus)
            if len(taxonomy.split(' ')) > 1:
                species = taxonomy_list[1]
                sql_extra += "AND Taxonomy.Species='{}'".format(species)
        sql += sql_extra

        tmp = db.query_all(sql)

        # select labels to be shown
        requested_labels = [1, 2, 3]
        if args.genome_labels:
            binary_num = bin(args.genome_labels)
            requested_labels = []
            for index in range(1, len(binary_num) - 1):
                if int(binary_num[-index]) > 0:
                    requested_labels += [index]

        names = [" ".join([row[col] for col in requested_labels]) for row in tmp]

        dm = pd.DataFrame(0, columns=names, index=names)
        for i in range(len(names)):
            for j in range(len(names)):
                idx = names[i]
                col = names[j]
                if i == j:
                    dm.loc[idx, col] = 1
                else:
                    lin_idx = tmp[i][4].split(',')
                    lin_col = tmp[j][4].split(',')
                    if lin_idx[0] != lin_col[0]:
                        dm.loc[idx, col] = 0.6
                    else:
                        scheme_limit = 0
                        for scheme_limit in range(0, len(scheme)):
                            if lin_idx[scheme_limit] == lin_col[scheme_limit]:
                                continue
                            else:
                                break
                        threshold = scheme[scheme_limit - 1]
                        dm.loc[idx, col] = threshold
        dm.to_csv(df, sep='\t')

    elif method == 'build_trees':
        if args.df == '':
            print("Please provide the name of the input distance matrix.")
            exit(8)

        distance_file_path = os.path.abspath(args.df)
        df = pd.read_csv(distance_file_path, sep="\t", header=0)
        df = df.drop(['Unnamed: 0'], axis=1)

        os.chdir(workspace)
        if not os.path.isdir('trees'):
            os.mkdir("trees")
        os.chdir("trees")

        # Distance matrix to similarity matirx
        similarity_matirx = 1 - df.values
        strains = list(df.columns)
        constructor = DistanceTreeConstructor()

        # create a list of list with only th upper triangle values
        # The inner lists have decreasing lengths)
        tril = np.tril(similarity_matirx).tolist()
        for i in range(len(tril)):
            for j in range(len(tril[i])):
                if j == i:
                    tril[i] = tril[i][0:j + 1]
                    break

        dm = DistanceMatrix(strains, tril)
        nj_tree = constructor.nj(dm)
        upgma_tree = constructor.upgma(dm)
        rename_clades(nj_tree.clade)
        rename_clades(upgma_tree.clade)

        save_tree(nj_tree, "nj.png")
        save_tree(upgma_tree, "upgma.png")
        Phylo.write(nj_tree, "nj.tree", 'newick')
        Phylo.write(upgma_tree, "upgma.tree", 'newick')

        nj_tree.ladderize()
        upgma_tree.ladderize()

        save_tree(nj_tree, "nj_ladderized.png")
        save_tree(upgma_tree, "upgma_ladderized.png")
        Phylo.write(nj_tree, "nj_ladderized.tree", 'newick')
        Phylo.write(upgma_tree, "upgma_ladderized.tree", 'newick')

    exit(0)
