#!/usr/bin/env python3
import os
import sys

from Bio import SeqIO
from Bio.Blast import NCBIXML


def parseVSGblast(inputfile): 
    v = inputfile+'.xml'
    n = inputfile+'_nonVSG.xml' 
    fasta_file = inputfile+'.fasta'
    cdhit_fasta_file = inputfile+'_cdhit.fasta'

    # Check which file exists
    if os.path.exists(cdhit_fasta_file):
        s = cdhit_fasta_file
    elif os.path.exists(fasta_file):
        s = fasta_file
    else:
        raise FileNotFoundError("Neither .fasta nor _cdhit.fasta file found.")

    record_dict = SeqIO.index(s, "fasta")
    hit_list = []
    exclude_list = []

    with open(n) as nonVSGresult_handle:
        for blast_record_nonVSG in NCBIXML.parse(nonVSGresult_handle):
            for alignment in blast_record_nonVSG.alignments:
                for hsp in alignment.hsps:
                    percent_identity = (100.0 * hsp.identities) / alignment.length
                    percent_query_identity = (100.0 * hsp.identities) / blast_record_nonVSG.query_letters
                    if ((percent_query_identity > 30 and hsp.identities > 300) or (percent_identity > 90)) \
                            and blast_record_nonVSG.query not in exclude_list:
                        exclude_list.append(str(blast_record_nonVSG.query))

    with open(v) as result_handle, open(inputfile.split('.')[0]+'_VSGs.fasta', 'w') as outfile:
        for blast_record in NCBIXML.parse(result_handle):
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 1.0e-10 \
                            and blast_record.query not in hit_list \
                            and blast_record.query not in exclude_list:
                        hit_list.append(str(blast_record.query))
                        SeqIO.write(record_dict[blast_record.query], outfile, "fasta")

if __name__ == "__main__":
    parseVSGblast(sys.argv[1])
