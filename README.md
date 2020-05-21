# LINflow
LINflow is a standalone Life Identification Number (LIN) assignment workflow used to compute the similarity among a set of prokaryotic genomes.

## Dependencies
sourmash  
pyani

## Usage
```
usage: LINflow.py [-h] [-s SCHEME_ID] [-i INPUT_DIR] [-m METADATA]
                  {initiate,show_schemes,add_scheme,add_genomes} workspace

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
                        the sample in https://bit.ly/2Y6Pw3R, and save as CSV
                        (comma separated values) format file.

```
## Functions
### Initiate the database
Initiate a workspace for LINflow with a database and file structure.  
```shell
python ./LINflow.py initiate path/to/workspace
```  
This will create a SQLite database and folders to store genome files and signatures.  

### Show schemes
Show current schemes in the database, based on which LINs are assigned.
```
python ./LINflow.py show_schemes
```
The first column is the Scheme ID as integers.

### Add a new scheme
Add a new scheme to the database.
```
python ./LINflow.py add_scheme
```
This will ask users to add a string of numbers delimited by comma (,) and a descriptive sentence of this new scheme.  
The added scheme should be a list of percentages without the percent mark (%), for example, if you wish to add a scheme 
of 70%, 80%, 90%, 95%, 99%, 99.5%, 99.9%, enter "70,80,90,95,99,99.5,99.9". Use `python ./LINflow.py show_schemes` to 
see examples.  

### Add genomes
Add genomes to the database and get LINs assigned.
```
python ./LINflow.py add_genomes /path/to/workspace -i /path/to/folder/of/genomes -m /path/to/metadata/file -s scheme_id
```
Use `-i` for the directory with only the genome files going to be added to the database.  Use `-m` for the taxonomic 
information of the genomes going to be added, including genus, species and strain names, use [https://bit.ly/2Y6Pw3R] 
as the template and save it as a CSV file. Use `-s` for indicating the Scheme ID by entering the integer corresponding 
to the scheme, see result from the `show_schemes` function.
