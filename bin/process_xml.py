import argparse
from Bio.Blast import NCBIXML

def makescoredict(xmlName): 
    result_handle = open(xmlName, 'r')
    blast_records = NCBIXML.parse(result_handle)
    hit_list = []
    scoredict = {}
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 1.0e-10: 
                    if not blast_record.query in hit_list: 
                        hit_list.append(str(blast_record.query)) 
                        scoredict[str(blast_record.query)] = ('\t'+str(alignment.title)+'\t'+str((100.0 * hsp.identities) / blast_record.query_letters)+'\t'+str((100.0 * hsp.identities) / alignment.length)+'\t'+str(alignment.length)+'\t'+str(blast_record.query_letters)+'\n')
    return scoredict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process BLAST XML output.')
    parser.add_argument('xmlName', type=str, help='The name of the XML file to process.')
    args = parser.parse_args()

    result = makescoredict(args.xmlName)
    print(result)
