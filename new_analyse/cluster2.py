import pandas as pd
import re

def parse_cluster_file(cluster_file):
    """
    Parse the cluster file to create a mapping of transcript IDs to cluster numbers.
    """
    cluster_dict = {}
    with open(cluster_file, 'r') as f:
        current_cluster = None
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                # Extract the cluster number
                current_cluster = line.split()[1]
            elif line:
                # Extract the transcript ID
                match = re.search(r'>([^\s\.]+)', line)  # Match transcript ID after '>'
                if match and current_cluster is not None:
                    transcript_id = match.group(1)  # Extract the transcript ID
                    cluster_dict[transcript_id] = current_cluster  # Map to the current cluster
    return cluster_dict

def update_tpm_with_clusters(tpm_file, cluster_file, output_file):
    """
    Update the TPM file with cluster numbers from the cluster file.
    """
    # Parse the cluster file to get the mapping of transcript IDs to cluster numbers
    cluster_dict = parse_cluster_file(cluster_file)
    
    # Read the TPM file
    tpm_df = pd.read_csv(tpm_file, sep=',', header=0)
    
    # Add a new 'cluster' column by mapping the first column (transcript IDs) to the cluster numbers
    tpm_df['cluster'] = tpm_df.iloc[:, 0].map(lambda x: cluster_dict.get(x.split('.')[0], 'NA'))
    
    # Save the updated TPM file
    tpm_df.to_csv(output_file, sep=',', index=False)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python update_tpm_with_clusters.py <tpm_file> <cluster_file> <output_file>")
        sys.exit(1)
    
    tpm_file = sys.argv[1]
    cluster_file = sys.argv[2]
    output_file = sys.argv[3]
    
    update_tpm_with_clusters(tpm_file, cluster_file, output_file)
