import re

def add_cluster_to_lines(cluster_file, output_file):
    with open(cluster_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_cluster = None
        for line in infile:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = line.split()[1]  # Extract cluster number
                outfile.write(line + '\n')  # Write the cluster header
            elif line:
                match = re.search(r'>([^\s]+)', line)
                if match:
                    transcript_id = match.group(1)
                    outfile.write(f"{line}\t{current_cluster}\n")  # Add cluster to the line
                else:
                    outfile.write(line + '\n')  # Write other lines as is

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python add_cluster_to_lines.py <cluster_file> <output_file>")
        sys.exit(1)

    cluster_file = sys.argv[1]
    output_file = sys.argv[2]
    add_cluster_to_lines(cluster_file, output_file)
