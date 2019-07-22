#!/usr/bin/python
"""
"""

# IMPORT
import os
import sys
import sqlite3
import argparse

# FUNCTIONS
def connect_to_db('LINbase.db'):
    conn = sqlite3.connect('LINbase.db')
    c = conn.cursor()
    return conn, c

def get_parsed_args():
    parser = argparse.ArgumentParser(
        description="LINflow"
    )
    parser.add_argument("function", type=str, choices=['initiate','show_schemes','add_scheme','add_genomes'], required=True, dest='function')
    parser.add_argument("workspace",dest="workspace", help="The location of the workspace",required=True)
    parser.add_argument("--scheme_id", dest="Scheme_ID", help="The Scheme based on which LINs are going to be assigned.")
    parser.add_argument("-s", dest="Interest_ID", help="Interest ID")
    parser.add_argument("-t", dest="Taxonomy", help="Taxonomy")
    parser.add_argument("-a", dest="Attributes",help="Attributes")
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
        c.execute('CREATE TABLE Genome (Genome_ID INT NOT NULL AUTO_INCREMENT,'
                  'Interest_ID INT NOT NULL,'
                  'FilePath TEXT NOT NULL,'
                  'PRIMARY KEY (Genome_ID))')
        c.execute('CREATE TABLE Scheme (Scheme_ID int NOT NULL AUTO_INCREMENT,'
                  'Cutoff text(255) NOT NULL,'
                  'LabelNum int NOT NULL,'
                  'Description TEXT,'
                  'PRIMARY KEY (Scheme_ID))')
        c.execute('CREATE TABLE ANI (ANI_ID INT NOT NULL AUTO_INCREMENT,'
                  'Genome_ID INT NOT NULL,'
                  'SubjectGenome_ID INT NOT NULL,'
                  'ANI DOUBLE NOT NULL,'
                  'PRIMARY KEY (ANI_ID))')
        c.execute('CREATE TABLE LIN (LIN_ID INT NOT NULL AUTO_INCREMENT,'
                  'Genome_ID INT NOT NULL,'
                  'Scheme_ID INT NOT NULL,'
                  'LIN TEXT NOT NULL,'
                  'PRIMARY KEY (LIN_ID))')
        c.execute('CREATE TABLE Taxonomy (Taxonomy_ID INT NOT NULL AUTO_INCREMENT,'
                  'Genome_ID INT NOT NULL,'
                  'Genus TEXT NOT NULL,'
                  'Species TEXT NOT NULL,'
                  'Strain TEXT NOT NULL,'
                  'PRIMARY KEY (Taxonomy_ID))')
        c.execute(
                "INSERT INTO Scheme (Cutoff, LabelNum, Description) values ('70,75,80,85,90,95,96,97,98,98.5,99,99.25,99.5,99.75,99.9,99.925,99.95,99.975,99.99,99.999', 20, 'The default scheme of LINbase.org')")
        conn.commit()
        fine_scheme = ','.join([str(i/float(10)) for i in range(700,1000,1)])
        c.execute("INSERT INTO Scheme (Cutoff, LabelNum) VALUES (?,300,'A scheme used to approximate the tree')",fine_scheme)
        conn.commit()
        conn.close()
        os.mkdir('Genomes')
        os.mkdir('Signatures')
        os.chdir('Signatures')
        os.mkdir('All')
        os.mkdir('rep_bac')

# MAIN
if __name__ == '__main__':
    argv = sys.argv
    args = get_parsed_args()
    method = args.function
    workspace = args.workspace
    if method == 'initiate':
        initiate(workspace)
    elif method == 'show_schemes':
         os.chdir(workspace)
         conn,c=connect_to_db()
         c.execute('SELECT Genome_ID,Cutoff,Description FROM Scheme')
         tmp = c.fetchall()
         print('Genome_ID\tDescription\tScheme')
         for i in tmp:
             print('{0}\t{1}\t{2}'.format(i[0],i[2],i[1]))
         conn.close()
    elif method == 'add_scheme':
        os.chdir(workspace)
        conn,c=connect_to_db()
        scheme = input("Please enter a new scheme delimited by comma (,) without the percentage mark, e.g. 70,80,90,100: ")
        num = len(scheme.split(","))
        description = input("Please enter a short description (optional):")
        c.execute("INSERT INTO Scheme (Cutoff,LabelNum,Description) VALUES ('{0}',{1},'{2}')".format(scheme,num,description))
        conn.commit()
        c.execute("SELECT Genome_ID,Description,Cutoff FROM Scheme WHERE Scheme_ID=(SELECT max(Scheme_ID) from Scheme)")
        tmp = c.fetchone()
        print('{0}\t{1}\t{2}'.format(tmp[0],tmp[1],tmp[2]))
    elif method == 'add_genomes':