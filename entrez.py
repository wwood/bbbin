#!/usr/bin/python
#
# Script to read FASTA sequences from GenBank
#
# This script is in the public domain.
#
# Usage: entrez.py   database  search_term output_file
#
# eg  entrez.py protein nematode nematode.fna
#
# Note: user needs to place email address in the script - name@domain.edu.au, line 30
#
# Ross Hall 2009
#


from Bio import Entrez
import os
import sys

argc = len(sys.argv)
if argc != 4: 
	print 'Usage: entrez.py database search_term ouputfile'
	print 'Eg: entrez.py protein "Mollusca [ORGN]" mollusca_protein.fasta' 
	exit(1);
	

database = sys.argv[1]
search_term = sys.argv[2]
output_file = sys.argv[3]


# Need your email address when querying Entrez
email_address = "name@domain.edu.au"


# Query the Entrez database
Entrez.email = email_address
#search_handle = Entrez.esearch(db=database,term=search_term,usehistory="y", retmax=50000)
search_handle = Entrez.esearch(db=database,term=search_term,usehistory="y", retmax=5)
search_results = Entrez.read(search_handle)
search_handle.close()

# Write number of sequences found
gi_list = search_results["IdList"]
count = int(search_results["Count"])
print count,
print ' sequences found.'

# Use the history of the above search
# to provide queries for download
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"] 


# Download batches of 5000 sequences at atime
#batch_size = 5000
batch_size = 5

# File where the sequences are written to
out_handle = open(output_file, "w")

# Download the sequences
for start in range(0,count,batch_size) :
    end = min(count, start+batch_size)
    print "Going to download record %i to %i" % (start+1, end)
    fetch_handle = Entrez.efetch(db=database, rettype="fasta",
                                 retstart=start, retmax=batch_size,
                                 webenv=webenv, query_key=query_key)
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()
