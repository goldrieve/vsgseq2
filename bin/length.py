from Bio import SeqIO
import csv
import argparse

def count_vsg_lengths(fasta_file, csv_file):
    """
    Counts the length of each VSG sequence in a multifasta file and
    prints the output as a CSV with the fasta header name and its sequence length.

    Args:
        fasta_file (str): Path to the input multifasta file.
        csv_file (str): Path to the output CSV file.
    """

    with open(csv_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['vsg_id', 'length'])  # Write header row

        for record in SeqIO.parse(fasta_file, "fasta"):
            fasta_header = record.id.replace("cluster_", "")  # Extract the fasta header (sequence ID)
            sequence_length = len(record.seq)  # Calculate the sequence length
            csvwriter.writerow([fasta_header, sequence_length])  # Write to CSV

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Count VSG sequence lengths from a FASTA file and output to a CSV.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file.")
    parser.add_argument("csv_file", help="Path to the output CSV file.")

    # Parse command line arguments
    args = parser.parse_args()

    # Call the function with the provided arguments
    count_vsg_lengths(args.fasta_file, args.csv_file)

    print(f"VSG lengths written to {args.csv_file}")
