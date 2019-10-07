# Note: Don't download more than 100 sequences during peak usage times
# Make sure you have installed biopython: https://biopython.org/wiki/Download
# Make sure you are using Python 3

# import modules
from Bio import Entrez
from Bio import SeqIO
import os

# tell NCBI who you are (put in your email here)
Entrez.email = "jvandal@vandals.uidaho.edu" 

# ask user for input in csv-formatted list of genbank IDs    
myinput = input("Please specify the GenBankIDs you want to save (in csv format):\n").split(',')

# define the name of the directory to be created
current_path = os.getcwd()
new_directory = "seq-files/"
new_directory_full_path = str(current_path) + "/" + str(new_directory)

# define the access rights
access_rights = 0o755

# make new sub-directory
try:
    os.mkdir(new_directory, access_rights)
except OSError:
    print ("Creation of the directory %s failed" % new_directory)
else:
    print ("Successfully created the directory %s " % new_directory)

# downloads genbank format files into folder running
i = 0
for myfilename in myinput:
    handle = Entrez.efetch(db="nucleotide", id=myfilename, rettype="gb", retmode="text")
    myfilename = myfilename + '.gb'
    #myfile = open(str(myfilename)+'.txt', 'w')
    myfile = open(os.path.join(new_directory, myfilename), 'w')
    myfile.write(handle.read())
    myfile.close()
    i += 1

# print the number of files downloaded and the directory they're downloaded to
print ("Successfully downloaded", i, "sequence files in", new_directory_full_path)
