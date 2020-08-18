# LINflow
LINflow is a standalone Life Identification Number (LIN) assignment 
workflow used to discover the genomic relatedness of bacteria.  

## Dependencies
* Bash Shell
* MySQL 5.6.14 
* BLAST
* MUMmer
* Conda with Python 3.6+
* Sourmash
* MySQL python connector
* Filelock Python3 module

**Important**: PyANI requires BLAST and MUMmer executables on PATH

### BLAST
The latest version of BLAST can be found [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 
and installed using your package manager (apt, yum, brew) or 
extracted from the tar.gz files as precompiled binaries.

### MUMmer
The latest version of MUMmer can be found [here](https://sourceforge.net/projects/mummer/files/) 
and installed using the following commands:
```shell script
tar xzf MUMmer{version}.tar.gz
cd MuMer{version}
./configure --prefix=/path/to/MUMer/installation/directory
make
make install
```
You need to use ```sudo make install``` if you are installing at a system location.
You can find MUMmer's full installation instructions [here](https://github.com/mummer4/mummer/blob/master/INSTALL.md).

### Configuring PATH
You need to add the BLAST and MUMmer binaries to your PATH. 
On Linux you can add the following line to your "~/.bashrc"
(~/.bash_profile in OSX) to make this change permanent, 
or run it to make the two available within the current terminal.
The bin directory can be found within the installation directory:
```shell script
export PATH=$PATH:/path/to/BLAST/bin/directory:/path/to/MUMer/bin/directory
```
To ensure your PATH setup was successful, the two commands below 
should show you the usage options of BLAST and MUMmer without any 
error messages

`blastn -h`
```
USAGE
  blastn [-h] [-help] [-import_search_strategy filename]
    [-export_search_strategy filename] [-task task_name] [-db database_name]
    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
    [-negative_gilist filename] [-negative_seqidlist filename]
    [-taxids taxids] [-negative_taxids taxids] [-taxidlist filename]
    [-negative_taxidlist filename] [-entrez_query entrez_query]
    ...

DESCRIPTION
   Nucleotide-Nucleotide BLAST 2.9.0+

Use '-help' to print detailed descriptions of command line arguments
```


`mummer -h`
```
Usage: mummer [options] <reference-file> <query-files>

Find and output (to stdout) the positions and length of all
sufficiently long maximal matches of a substring in
<query-file> and <reference-file>

Options:
-mum           compute maximal matches that are unique in both sequences
-mumcand       same as -mumreference
```

### Conda with Python 3.6+
To get access to the latest Sourmash release you need to install Conda. 
A minimalist approach is to install Miniconda for python3 available [here](https://docs.conda.io/en/latest/miniconda.html#linux-installers).
You can then create and activate a virtual environment and install the required modules:
```shell script
conda create -n linflow
conda activate linflow
conda install -c conda-forge -c bioconda sourmash pyani
conda install mysqlclient filelock
```

### MySQL Connection
Create a file at `~/.my.cnf` with the syntax below. You should be able to log into MySQL
by typing `mysql` without a username or password after creating the file. Make sure you have permission to create
databases and tables when logging in using this method.
```
[client]
user=yourusername
password=yourpassowrd
socket=/var/run/mysqld/mysqld.sock
```
You can find your socket file at /etc/mysql/mysql.conf.d/mysqld.cnf (Ubuntu)
or /var/lib/mysql/mysql.sock (CentOS)

## Usage
```
Usage: LINflow.py [-h] [-s SCHEME_ID] [-i INPUT_DIR] [-m METADATA]
                  {initiate,show_schemes,add_scheme,add_genomes,infer_distance_by_LIN,cleanup} workspace
LINflow
positional arguments:
  {initiate,show_schemes,add_scheme,add_genomes}
  workspace             The location of the workspace
optional arguments:
  -h, --help            show this help message and exit
  -s SCHEME_ID          The Scheme based on which LINs are going to be
                        assigned.
  -i INPUT_DIR          The directory of genomes going to be added.
  -m METADATA           The metadata corresponding to the genomes. Download
                        the sample in https://bit.ly/m2Y6Pw3R, and save as CSV
                        (comma separated values) format file.
  -d SIMILARITY_MATRIX  Path to the similarity matrix output. The accuracy of the
                        matrix will be based on the the number of positions in
                        the selected scheme.
  -g GENOME_LABEL       Genome label included in distance matrix 1(Genome) 2(Species)
                        4(Strain). Add numbers to combine. Default: 7(All) 1+2+4=7
```
## Functions
### Initiate the database
Initiate a workspace for LINflow with a database and file structure.  
```shell script
python LINflow.py initiate path/to/workspace
```  
This will create a SQLite database and folders to store genome files and signatures.  

### Show all schemes
Show current schemes in the database, based on which LINs are assigned.
```shell script
python LINflow.py show_schemes path/to/workspace
```
The first column is the Scheme ID as integers.

### Add a new scheme
Add a new scheme to the database.
```shell script
python LINflow.py add_scheme
```
This will ask users to add a string of numbers delimited by comma (,) 
and a descriptive sentence of this new scheme or if evenly spaced, 
a MIN, MAX range with a step size could be provided as well.  
When adding manually, the added scheme should be a list of percentages 
without the percent mark (%), for example, if you wish to add a scheme 
of 70%, 80%, 90%, 95%, 99%, 99.5%, 99.9%, enter "70,80,90,95,99,99.5,99.9". 
Use the `show_schemes` function to see examples.

### Add genomes
Add genomes to the database and get LINs assigned.
```shell script
python LINflow.py add_genomes /path/to/workspace \
-i /path/to/folder/of/genomes \
-m /path/to/metadata/file \
-s scheme_id
```
Use `-i` for the directory with only the genome files going to be added 
to the database. Use `-m` for the taxonomic information of the genomes 
going to be added, including genus, species, and strain names, use [https://bit.ly/2Y6Pw3R] 
as the template and save it as a CSV file. Use `-s` for indicating the Scheme_ID 
by entering the integer corresponding to the scheme, see result from the `show_schemes` function.

### Export similarity matrix
Exports a similarity matrix (tab separated) based on the selected Scheme_ID (denoted by `-s`).
The Genome labels (denoted by `-g`) that will be used in the matrix based on the numbers 
1(Genome) 2(Species) 4(Strain). Two or more labels can be used together by adding their numbers.  
e.g. 6 (=2+4) will include the *species* and *strain* name but not the *genome* name. 
 7 (=1+2+4) includes Taxon data for each label. 
```shell script
python LINflow.py infer_distance_by_LIN /path/to/workspace \
-s schema scheme_id \
-d /path/to/output/matrix_file \
-g genome_labels_to_include
```

### Build tree
Creates two NJ and UPGMA trees each, one ladderized and one not. The trees will be saved 
under the workspace directory in the "trees" directory.
```shell script
python LINflow.py build_trees /path/to/workspace
-d /path/to/input/matrix_file
```


### Clean Workspace
Deletes all files and folders in workspace and its corresponding MySQL database. 
(**ALL DATA WILL BE ERASED**)
```shell script
python LINflow.py cleanup /path/to/workspace
```

