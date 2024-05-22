import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO

def parseVSGblast(inputfile): 
    v = inputfile+'.xml'
    n = inputfile+'_nonVSG.xml' 
    s = inputfile+"_cdhit.fasta"
    
    result_handle = open(v)
    blast_records = NCBIXML.parse(result_handle) 
    record_dict = SeqIO.index(s,"fasta")
    
    outfile = open(inputfile.split('.')[0]+'_VSGs.fasta', 'w')
    hit_list = []

    nonVSGresult_handle = open(n)
    blast_records_nonVSG = NCBIXML.parse(nonVSGresult_handle) 
    exclude_list = []
                
    for blast_record_nonVSG in blast_records_nonVSG:
        for alignment in blast_record_nonVSG.alignments:
            for hsp in alignment.hsps: 
                percent_identity = (100.0 * hsp.identities) / alignment.length 
                percent_query_identity = (100.0 * hsp.identities) / blast_record_nonVSG.query_letters
                if (percent_query_identity > 30 and hsp.identities > 300) or (percent_identity > 90):
                    if not blast_record_nonVSG.query in exclude_list:
                        exclude_list.append(str(blast_record_nonVSG.query))

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 1.0e-10: 
                    if not blast_record.query in hit_list: 
                        if not blast_record.query in exclude_list: 
                            hit_list.append(str(blast_record.query))
                            SeqIO.write(record_dict[blast_record.query], outfile, "fasta")
                            
if __name__ == "__main__":
    parseVSGblast(sys.argv[1])
