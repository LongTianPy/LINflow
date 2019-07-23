#!/usr/bin/python
"""
"""

# IMPORT
import os
import sys
import uuid
import sqlite3
import argparse
import shutil
import pandas as pd
from os.path import isfile, isdir, join
from distutils.spawn import find_executable

# OBJECTS
class getLIN(object):
    """
    Read the LIN of the top similar genome.
    """
    def __init__(self, Genome_ID, Scheme_ID, similarity,c):
        self.Genome_ID = Genome_ID
        self.Scheme_ID = Scheme_ID
        self.c=c
        c.execute("SELECT LabelNum from Scheme WHERE Scheme_ID={0}".format(Scheme_ID))
        self.label_num = int(c.fetchone()[0])
        self.similarity = float(similarity)*100
        self.parse()
    def parse(self, Genome_ID = None, Scheme_ID = None, similarity = None, c = None):
        if not Genome_ID:
            Genome_ID = self.Genome_ID
        if not Scheme_ID:
            Scheme_ID = self.Scheme_ID
        if not similarity:
            similarity = self.similarity
        if not c:
            c = self.c
        # Read the LIN of this genome
        c.execute('SELECT LIN from LIN where Genome_ID = {0} and LIN.Scheme_ID={1}'.format(int(Genome_ID),int(Scheme_ID)))
        lin = c.fetchone()[0].split(',')
        self.LIN = lin
        # Read the cutoff of this scheme
        c.execute('SELECT Cutoff from Scheme where Scheme_ID={0}'.format(Scheme_ID))
        cutoff = c.fetchone()[0].split(',')
        cutoff = [float(i) for i in cutoff]
        if similarity < cutoff[0]:
            idx_to_change = 0
        elif similarity >= cutoff[-1]:
            idx_to_change = "n/a"
        else:
            for i in range(len(cutoff)-1):
                if cutoff[i] <= similarity < cutoff[i+1]:
                    idx_to_change = i+1
                else:
                    continue
        self.idx_to_change = idx_to_change
        if idx_to_change == 0:
            self.conserved_LIN = ''
        elif idx_to_change == 'n/a':
            self.conserved_LIN = lin
        else:
            self.conserved_LIN = lin[:idx_to_change]

class Assign_LIN(object):
    """ Get the biggest number assigned to the idx_to_change with the same conserved part of LIN
    """
    def __init__(self, getLIN_object,c):
        self.idx_to_change = getLIN_object.idx_to_change
        self.conserved_LIN = ','.join(getLIN_object.conserved_LIN)
        self.label_num = getLIN_object.label_num
        self.scheme_id = getLIN_object.Scheme_ID
        self.c=c
        self.assign()
    def assign(self, idx_to_change=None, conserved_LIN=None, label_num=None,c=None,scheme_id=None):
        if not idx_to_change:
            idx_to_change = self.idx_to_change
        if not conserved_LIN:
            conserved_LIN = self.conserved_LIN
        if not label_num:
            label_num = self.label_num
        if not c:
            c=self.c
        if not scheme_id:
            scheme_id = self.scheme_id
        if conserved_LIN == '':
            c.execute("SELECT LIN.LIN FROM LIN where Scheme_ID=4")
            tmp = c.fetchall()
        else:
            sql="SELECT LIN.LIN from LIN WHERE LIN.LIN LIKE '{0}%' and Scheme_ID={1}".format(conserved_LIN,scheme_id)
            # print sql
            c.execute(sql)
            tmp = c.fetchall()
        if type(idx_to_change) == int:
            LINs = [int(i[0].split(',')[idx_to_change]) for i in tmp]
            num_to_assign = str(max(LINs)+1)
        if type(idx_to_change) == int and idx_to_change != label_num - 1:
            tail = ['0'] * (label_num - 1 - idx_to_change)
            tail = ','.join(tail)
            new_LIN = conserved_LIN + ',%s,'%num_to_assign + tail
        elif type(idx_to_change) != int:
            new_LIN = conserved_LIN
        else:
            new_LIN = conserved_LIN + ',%s'%num_to_assign
        if new_LIN.startswith(','):
            new_LIN = new_LIN[1:]
        else:
            new_LIN = new_LIN
        self.new_LIN = new_LIN


# FUNCTIONS
def connect_to_db():
    conn = sqlite3.connect('LINbase.db')
    c = conn.cursor()
    return conn, c

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="LINflow"
    )
    parser.add_argument("function", type=str, choices=['initiate','show_schemes','add_scheme','add_genomes'])
    parser.add_argument("workspace", type=str, help="The location of the workspace")
    parser.add_argument("-s", dest="Scheme_ID", help="The Scheme based on which LINs are going to be assigned.", type=int,default=0)
    parser.add_argument("-i", dest="input_dir", help="The directory of genomes going to be added.",default='')
    parser.add_argument("-m", dest="metadata", default='', help="The metadata corresponding to the genomes. Download the sample in https://bit.ly/2Y6Pw3R, and save as CSV (comma separated values) format file.")
    # parser.add_argument("-p", dest="privacy", help="Is it private information")
    args = parser.parse_args()
    return args

def initiate(workspace):
    if os.path.isdir(workspace):
        print("Provided workspace is not empty, please use a clean environment")
    else:
        os.mkdir(workspace)
        os.chdir(workspace)
        conn, c = connect_to_db()
        c.execute('CREATE TABLE Genome (Genome_ID INTEGER PRIMARY KEY AUTOINCREMENT,'
                  'FilePath TEXT NOT NULL)')
        c.execute('CREATE TABLE Scheme (Scheme_ID INTEGER PRIMARY KEY AUTOINCREMENT,'
                  'Cutoff text(255) NOT NULL,'
                  'LabelNum int NOT NULL,'
                  'Description TEXT)')
        c.execute('CREATE TABLE ANI (ANI_ID INTEGER PRIMARY KEY AUTOINCREMENT,'
                  'Genome_ID INT NOT NULL,'
                  'SubjectGenome INT NOT NULL,'
                  'ANI DOUBLE NOT NULL)')
        c.execute('CREATE TABLE LIN (LIN_ID INTEGER PRIMARY KEY AUTOINCREMENT,'
                  'Genome_ID INT NOT NULL,'
                  'Scheme_ID INT NOT NULL,'
                  'LIN TEXT NOT NULL)')
        c.execute('CREATE TABLE Taxonomy (Taxonomy_ID INTEGER PRIMARY KEY AUTOINCREMENT,'
                  'Genome_ID INT NOT NULL,'
                  'Genus TEXT NOT NULL,'
                  'Species TEXT NOT NULL,'
                  'Strain TEXT NOT NULL)')
        c.execute(
                "INSERT INTO Scheme (Cutoff, LabelNum, Description) values ('70,75,80,85,90,95,96,97,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999', 20, 'The default scheme of LINbase.org')")
        conn.commit()
        fine_scheme = ','.join([str(i/float(10)) for i in range(700,1000,1)])
        c.execute("INSERT INTO Scheme (Cutoff, LabelNum, Description) VALUES ('{}',300,'A scheme used to approximate the tree')".format(fine_scheme))
        conn.commit()
        conn.close()
        os.mkdir('Genomes')
        os.mkdir("ANI")
        os.mkdir('Signatures')
        os.chdir('Signatures')
        os.mkdir('All')
        os.mkdir("tmp_sig")
        os.mkdir("tmp_result")
        os.mkdir('rep_bac')

def create_sketch(filepath,dest):
    cmd = "sourmash compute -o {0} {1} -k 21,51 -n 2000 -q".format(dest, filepath)
    os.system(cmd)
    return dest

def compare_sketch(query,LINgroup,k):
    if LINgroup == "rep_bac":
        dest = rep_bac_dir
    else:
        dest = join(sourmash_dir, LINgroup)
    folder_size = len([file for file in os.listdir(dest) if isfile(join(dest,file))])
    cmd = "sourmash search {0} {1} -n {2} -k 21 -q --threshold 0.0001 -o {3}"
    cmd = cmd.format(query, join(dest,'*.sig'), folder_size, join(sourmash_result,"tmp_result.txt"))
    os.system(cmd)
    return join(sourmash_result,"tmp_result.txt")

def parse_result(result_file):
    df = pd.read_csv(result_file,sep=",",header=0)
    if df.empty:
        return df
    else:
        ids = []
        for each in df['filename']:
            id = int(each.split('/')[-1].split('.')[0])
            ids.append(id)
        df.index = ids
        return df

def add_genome(filename, taxonomy, target_filename,scheme_id):
    tmp_sig = create_sketch(filename,join(tmp_sig_dir,"tmp.sig"))
    result_file = compare_sketch(tmp_sig, "rep_bac",'21')
    df = parse_result(result_file)
    if df.empty:
        new_LIN_object = getLIN(Genome_ID=1, Scheme_ID=scheme_id,similarity=0.6,c=c)
        new_LIN = Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
        ANIb_result = 0.6
        SubjectGenome=1
    else:
        rep_bac_Genome_ID = int(df.index[0])
        c.execute("select LIN from LIN where Genome_ID={0} and Scheme_ID=1".format(rep_bac_Genome_ID))
        rep_bac_LIN = c.fetchone()[0]
        rep_bac_LINgroup = ",".join(rep_bac_LIN.split(",")[:6])
        jaccard_similarity = df.loc[rep_bac_Genome_ID,"similarity"]
        if jaccard_similarity > 0.2475:
            result_file = compare_sketch(tmp_sig,rep_bac_LINgroup,'51')
            df = parse_result(result_file)
            SubjectGenome = int(df.index[0])
            c.execute("select FilePath from Genome where Genome_ID={0}".format(SubjectGenome))
            subject_genome_file = c.fetchone()[0]
            sub_working_dir = join(ani_dir, str(uuid.uuid4()))
            if not isdir(sub_working_dir):
                os.mkdir(sub_working_dir)
            shutil.copyfile(filename,join(sub_working_dir,"tmp.fasta"))
            shutil.copyfile(subject_genome_file,join(sub_working_dir,"{0}.fasta".format(str(SubjectGenome))))
            pyani_cmd = "average_nucleotide_identity.py " \
                        "-i {0} -o {1} -m ANIb --nocompress -f".format(sub_working_dir, join(sub_working_dir, 'output'))
            os.system(pyani_cmd)
            ANIb_result = pd.read_csv(join(sub_working_dir, "output", "ANIb_percentage_identity.tab"), sep="\t",
                                             header=0,
                                             index_col=0).loc['tmp', str(SubjectGenome)]
            shutil.rmtree(sub_working_dir)
            new_LIN_object = getLIN(Genome_ID=SubjectGenome, Scheme_ID=scheme_id, similarity=ANIb_result,c=c)
            new_LIN = Assign_LIN(getLIN_object=new_LIN_object,c=c).new_LIN
        elif jaccard_similarity <= 0.2475 and jaccard_similarity > 0.0025:
            result_file = compare_sketch(tmp_sig,rep_bac_LINgroup,'21')
            df = parse_result(result_file)
            ANIb_result = 0
            SubjectGenome = 0
            for each_subject_genome_ID in df.index[:3]:
                sub_working_dir = ani_dir + str(uuid.uuid4()) + "/"
                if not isdir(sub_working_dir):
                    os.mkdir(sub_working_dir)
                c.execute("select FilePath from Genome where Genome_ID={0}".format(each_subject_genome_ID))
                subject_genome_file = c.fetchone()[0]
                shutil.copyfile(filename, join(sub_working_dir, "tmp.fasta"))
                shutil.copyfile(subject_genome_file,
                                join(sub_working_dir, "{0}.fasta".format(each_subject_genome_ID)))
                pyani_cmd = "average_nucleotide_identity.py " \
                            "-i {0} -o {1} -m ANIb --nocompress -f".format(sub_working_dir,
                                                                           join(sub_working_dir, 'output'))
                os.system(pyani_cmd)
                this_ANIb_result = pd.read_csv(join(sub_working_dir, "output", "ANIb_percentage_identity.tab"),
                                                 sep="\t",
                                                 header=0,
                                                 index_col=0).loc['tmp', str(each_subject_genome_ID)]
                shutil.rmtree(sub_working_dir)
                if this_ANIb_result > ANIb_result:
                    ANIb_result = this_ANIb_result
                    SubjectGenome = each_subject_genome_ID
            new_LIN_object = getLIN(Genome_ID=SubjectGenome, Scheme_ID=scheme_id, similarity=ANIb_result, c=c)
            new_LIN = Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
        else:
            SubjectGenome = int(df.index[0])
            c.execute("select FilePath from Genome where Genome_ID={0}".format(SubjectGenome))
            subject_genome_file = c.fetchone()[0]
            sub_working_dir = join(ani_dir, str(uuid.uuid4()))
            if not isdir(sub_working_dir):
                os.mkdir(sub_working_dir)
            shutil.copyfile(filename, join(sub_working_dir, "tmp.fasta"))
            shutil.copyfile(subject_genome_file, join(sub_working_dir, "{0}.fasta".format(str(SubjectGenome))))
            pyani_cmd = "average_nucleotide_identity.py " \
                        "-i {0} -o {1} -m ANIb --nocompress -f".format(sub_working_dir, join(sub_working_dir, 'output'))
            os.system(pyani_cmd)
            ANIb_result = pd.read_csv(join(sub_working_dir, "output", "ANIb_percentage_identity.tab"), sep="\t",
                                        header=0,
                                        index_col=0).loc['tmp', str(SubjectGenome)]
            shutil.rmtree(sub_working_dir)
            new_LIN_object = getLIN(Genome_ID=SubjectGenome, Scheme_ID=scheme_id, similarity=ANIb_result, c=c)
            new_LIN = Assign_LIN(getLIN_object=new_LIN_object, c=c).new_LIN
    c.execute("INSERT INTO Genome (FilePath) VALUES ('{0}')".format(target_filename))
    conn.commit()
    c.execute("SELECT Genome_ID FROM Genome WHERE FilePath='{0}'".format(target_filename))
    genome_id = c.fetchone()[0]
    c.execute("INSERT INTO ANI (Genome_ID,SubjectGenome,ANI) VALUES ({0},{1},{2})".format(genome_id,SubjectGenome,ANIb_result))
    conn.commit()
    c.execute("INSERT INTO Taxonomy (Genome_ID, Genus, Species, Strain) VALUES ({0},'{1}','{2}','{3}')".format(genome_id,taxonomy['genus'],taxonomy['species'],taxonomy['strain']))
    conn.commit()
    c.execute("INSERT INTO LIN (Genome_ID,Scheme_ID,LIN) VALUES ({0},{1},'{2}')".format(genome_id,scheme_id,new_LIN))
    conn.commit()
    if scheme_id not in [1,2]:
        new_LIN_object_1 = getLIN(Genome_ID=SubjectGenome,Scheme_ID=1,similarity=ANIb_result,c=c)
        new_LIN_1 = Assign_LIN(getLIN_object=new_LIN_object_1,c=c).new_LIN
        c.execute(
            "INSERT INTO LIN (Genome_ID,Scheme_ID,LIN) VALUES ({0},{1},'{2}')".format(genome_id, 1, new_LIN_1))
        conn.commit()
        new_LIN_object_2 = getLIN(Genome_ID=SubjectGenome, Scheme_ID=2, similarity=ANIb_result, c=c)
        new_LIN_2 = Assign_LIN(getLIN_object=new_LIN_object_2, c=c).new_LIN
        c.execute(
                "INSERT INTO LIN (Genome_ID,Scheme_ID,LIN) VALUES ({0},{1},'{2}')".format(genome_id, 2, new_LIN_2))
        conn.commit()
    else:
        the_other_scheme_id = 3-scheme_id
        new_LIN_object_other = getLIN(Genome_ID=SubjectGenome, Scheme_ID=the_other_scheme_id, similarity=ANIb_result, c=c)
        new_LIN_other = Assign_LIN(getLIN_object=new_LIN_object_other, c=c).new_LIN
        c.execute(
                "INSERT INTO LIN (Genome_ID,Scheme_ID,LIN) VALUES ({0},{1},'{2}')".format(genome_id, the_other_scheme_id, new_LIN_other))
        conn.commit()
    # Update signature
    c.execute("SELECT LIN FROM LIN WHERE Genome_ID={0} and Scheme_ID=1".format(genome_id))
    LIN_1 = c.fetchone()[0]
    lingroup = ",".join(LIN_1.split(",")[:6])
    if not isdir(join(sourmash_dir,lingroup)):
        os.mkdir(join(sourmash_dir,lingroup))
        create_sketch(filename, join(rep_bac_dir,"{0}.sig".format(str(genome_id))))
        shutil.copyfile(join(rep_bac_dir,"{0}.sig".format(str(genome_id))),join(sourmash_dir,lingroup,"{0}.sig".format(str(genome_id))))
    else:
        create_sketch(filename, join(sourmash_dir,lingroup,"{0}.sig".format(str(genome_id))))
    shutil.copy(filename, target_filename)


# MAIN
if __name__ == '__main__':
    argv = sys.argv
    args = get_parsed_args()
    method = args.function
    workspace = args.workspace
    if method == 'initiate':
        initiate(workspace)
    else:
        if not isfile(join(workspace,"LINbase.db")):
            print("Please initiate first")
            sys.exit()
        else:
            genome_dir = join(workspace,"Genomes")
            ani_dir = join(workspace,"ANI")
            rep_bac_dir = join(workspace,"Signatures","rep_bac")
            sourmash_dir = join(workspace, "Signatures", "All")
            tmp_sig_dir = join(workspace, "Signatures", "tmp_sig")
            sourmash_result = join(workspace, "Signatures", "tmp_result")
            if method == 'show_schemes':
                os.chdir(workspace)
                conn,c=connect_to_db()
                c.execute('SELECT Scheme_ID,Cutoff,Description FROM Scheme')
                tmp = c.fetchall()
                print('Scheme_ID\tDescription\tScheme')
                for i in tmp:
                 print('{0}\t{1}\t{2}'.format(i[0],i[2],i[1]))
                conn.close()
            elif method == 'add_scheme':
                os.chdir(workspace)
                conn,c=connect_to_db()
                scheme = input("Please enter a new scheme delimited by comma (,) without the percentage mark and do not include 100 at the end, e.g. 70,80,90,99: ")
                num = len(scheme.split(","))
                description = input("Please enter a short description (optional):")
                c.execute("INSERT INTO Scheme (Cutoff,LabelNum,Description) VALUES ('{0}',{1},'{2}')".format(scheme,num,description))
                conn.commit()
                c.execute("SELECT Scheme_ID,Description,Cutoff FROM Scheme WHERE Scheme_ID=(SELECT max(Scheme_ID) from Scheme)")
                tmp = c.fetchone()
                print('{0}\t{1}\t{2}'.format(tmp[0],tmp[1],tmp[2]))
            elif method == 'add_genomes':
                if args.input_dir == '':
                    print("Please provide a directory with genomes to be added with -i")
                    sys.exit()
                else:
                    input_dir = args.input_dir
                if args.metadata == '':
                    print("Please provide the metadata file with -m")
                    sys.exit()
                else:
                    meta = args.metadata
                if args.Scheme_ID == 0:
                    print("Please select one Scheme_ID with -s")
                    sys.exit()
                else:
                    scheme_id = args.Scheme_ID
                df = pd.read_csv(meta,sep=",",header=0,index_col=0)
                if not find_executable("sourmash"):
                    print("Please install sourmash first")
                    sys.exit()
                if not find_executable("average_nucleotide_identity.py"):
                    print("Please install pyani first")
                    sys.exit()
                os.chdir(workspace)
                conn, c = connect_to_db()
                for i in df.index:
                    filename = join(input_dir,str(i))
                    uuid_filename = str(uuid.uuid4()) + ".fasta"
                    target_filename = join(genome_dir,uuid_filename)
                    genus = df.loc[i,"genus"]
                    species = df.loc[i,"species"]
                    strain = df.loc[i,"strain"]
                    taxonomy = {"genus":genus,"species":species,"strain":strain}
                    c.execute("select count(Genome_ID) from Genome")
                    size = c.fetchone()[0]
                    if size == 0:
                        genome_id = 1
                        uuid_filename = str(uuid.uuid4()) + ".fasta"
                        c.execute("insert into Genome (FilePath) values ('{0}')".format(target_filename))
                        conn.commit()
                        c.execute("insert into Taxonomy (Genome_ID,Genus,Species,Strain) values (1,'{0}','{1}', '{2}')".format(taxonomy['genus'],taxonomy['species'],taxonomy['strain']))
                        conn.commit()
                        c.execute("insert into ANI (Genome_ID, SubjectGenome, ANI) values (1,1,1)")
                        c.execute("insert into LIN (Genome_ID, Scheme_ID, LIN) values (1,1,'{0}')".format(",".join(["0"] * 20)))
                        c.execute("insert into LIN (Genome_ID,Scheme_ID,LIN) values (1,2,'{0}')".format(",".join(["0"] * 300)))
                        os.mkdir(join(sourmash_dir,"0,0,0,0,0,0"))
                        create_sketch(filename,join(sourmash_dir,"0,0,0,0,0,0","1.sig"))
                        shutil.copyfile(join(sourmash_dir,"0,0,0,0,0,0","1.sig"),join(rep_bac_dir, "1.sig"))
                        shutil.copy(filename, target_filename)
                    else:
                        new_LIN, SubjectGenome, ANIb_result =  add_genome(filename,taxonomy,target_filename,scheme_id)
                conn.close()




