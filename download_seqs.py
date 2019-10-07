# Note: Don't download more than 100 sequences during peak usage times
# Make sure you have installed biopython: https://biopython.org/wiki/Download
# Make sure you are using Python 3

# import modules
from Bio import Entrez
from Bio import SeqIO
import os

####### Gather user input #######
# ask user for email (tells NCBI who they are)
Entrez.email = input("Enter your email: \n")
# ask user for a csv-formatted list of genbank IDs    
seq_list = input("Paste accession numbers in csv format:\n").split(',')
# ask user to choose what format they want
choice = input("What format do you want? \n   \"1\" = nucleotide, genbank \n   \"2\" = nucleotide, fasta \n   \"3\" = protein, fasta \n")

####### MAKE DIRECTORY ############
# define the name of the directory to be created
current_path = os.getcwd()
new_directory = "seq-files/"
new_directory_full_path = str(current_path) + "/" + str(new_directory)
# make new sub-directory
try:
    access_rights = 0o755 # define the access rights
    os.mkdir(new_directory, access_rights)
except OSError:
    print ("Creation of the directory %s failed" % new_directory)
else:
    print ("Successfully created the directory %s " % new_directory)

############ DOWNLOAD FILES #############
# you can instead choose fasta format by changing the script as indicated below
# reference for "db": https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
# reference for "rettype": https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
i = 0
for seq_acc in seq_list:
    if choice == '1': # genbank format
        handle = Entrez.efetch(db="nuccore", id=seq_acc, rettype="gbwithparts", retmode="text")
        myfilename = seq_acc + '.gb'
        myfile = open(os.path.join(new_directory, myfilename), 'w')
        myfile.write(handle.read())
        myfile.close()
        i += 1
    if choice == '2': # fasta format
        handle = Entrez.efetch(db="nucleotide", id=seq_acc, rettype="fasta", retmode="text")
        myfilename = seq_acc + '.fasta'
        myfile = open(os.path.join(new_directory, myfilename), 'w')
        myfile.write(handle.read())
        myfile.close()
        i += 1
    if choice == '3': # proteins
        handle = Entrez.efetch(db="protein", id=seq_acc, rettype="fasta", retmode="text")
        myfilename = seq_acc + '.fasta'
        myfile = open(os.path.join(new_directory, myfilename), 'w')
        myfile.write(handle.read())
        myfile.close()
        i += 1

# print the number of files downloaded and the directory they're downloaded to
print ("Successfully downloaded", i, "sequence files in", new_directory_full_path)

