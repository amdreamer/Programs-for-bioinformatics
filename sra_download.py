#!/usr/bin/env python
import subprocess
import argparse
import csv
import sys

# add arguments
download=argparse.ArgumentParser(description="downlaoding sra files and transfer into fastq files, input sra ID, output fq files.")
download.add_argument('-sraID_file',type=argparse.FileType('r'),help='passing a file with sra IDs, one ID for each line.')
download.add_argument('-ID',nargs='+',help='Instead of passing a file with sraID for each line, user can also input sraID directly by using this parameter. ex: -ID SRR2121685 SRR2121686 SRR2121687')
options=download.parse_args()

# Get sra IDs from parameters.
# if input as sraID_file
if options.sraID_file is not None:
    sraID_file = csv.reader(options.sraID_file,quotechar=None)
    sra_numbers=[]
    for row in sraID_file:
        sra_numbers.append(row[0])
    print("The files ready for download:" + sra_numbers)

# if input by directly asigning sraIDs
if options.ID is not None:
    sra_numbers=options.ID
    print("The files ready for download:" + sra_numbers)

# example: sra_numbers=["SRR2121685", "SRR2121686", "SRR2121687", "SRR2121688"]
if all(v is None for v in [options.sraID_file, options.ID]):
    print("Please input sra IDs!")
    sys.exit()


# Start downloading
print("Start downloading sra files:")
# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = " ~/install/sratoolkit.2.11.0-mac64/bin/prefetch " + sra_id
    # replace the correct path of sratoolkit in your system!
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

print("Start generating fastq files:")
# this will extract the .sra files from above into the currently folder
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = " ~/install/sratoolkit.2.11.0-mac64/bin/fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ~/Downloads/GEO/sra/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)
