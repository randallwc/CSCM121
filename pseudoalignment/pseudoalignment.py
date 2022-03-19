import networkx as nx
import csv

# CONSTANTS
KMER_LENGTH = 30
TRANSCRIPTOME_FASTA = 'chr11_transcriptome.fasta'
READS_FASTA = 'reads.fasta'


def get_num_lines(file_name):
    with open(file_name, 'r') as f:
        return sum(1 for _ in f)


current_num_lines = get_num_lines(TRANSCRIPTOME_FASTA)
current_line_num = 0


def read_line(f):
    global current_num_lines
    global current_line_num
    print(f'\r{current_line_num/current_num_lines * 100:.2f} %', end='')
    current_line_num += 1
    return f.readline()


print('CREATE DEBRUIJN GRAPH')
graph = nx.DiGraph()  # directed graph
kmer_dict = {}
file = open(TRANSCRIPTOME_FASTA, 'r')
current_line = read_line(file)
while True:
    if current_line == '':
        break
    elif current_line[0] == '>':
        current_isoform = current_line[1:-1]
        current_line = read_line(file)
        current_read = ''
        while current_line != '' and current_line[0] != '>':
            current_read += current_line[:-1]
            current_line = read_line(file)
        previous_node = None
        for i in range(len(current_read) - KMER_LENGTH + 1):
            current_kmer = current_read[i:i + KMER_LENGTH]
            if current_kmer in graph:
                kmer_dict[current_kmer].add(current_isoform)
            else:
                graph.add_node(current_kmer)
                kmer_dict[current_kmer] = {current_isoform}
            if previous_node is not None:
                graph.add_edge(previous_node, current_kmer)
            previous_node = current_kmer
    else:
        exit(1)
file.close()

# GET EQUIVALENCE CLASSES
print()
print('GETTING EQUIVALENCE CLASSES')

current_num_lines = get_num_lines(READS_FASTA)
current_line_num = 0
file = open(READS_FASTA, 'r')
no_equivalence_class = 'NA'
equivalence_dict = {no_equivalence_class: [0, 0]}
for line in file:
    print(f'\r{current_line_num/current_num_lines * 100:.2f} %', end='')
    current_line_num += 1
    if line[0] == '>':
        continue
    current_read = line[:-1]
    current_kmer = current_read[:KMER_LENGTH]
    if current_kmer not in graph:
        equivalence_dict[no_equivalence_class][0] += 1
        continue
    current_isoforms = kmer_dict[current_kmer]
    previous_kmer = current_kmer
    for i in range(1, len(current_read) - KMER_LENGTH + 1):
        current_kmer = current_read[i:i + KMER_LENGTH]
        if current_kmer in graph.successors(previous_kmer):
            current_isoforms = current_isoforms & kmer_dict[current_kmer]
            previous_kmer = current_kmer
        else:
            current_isoforms = set()
            break
    if len(current_isoforms) > 0:
        current_isoforms_s = ','.join(sorted(list(current_isoforms)))
        if current_isoforms_s in equivalence_dict:
            equivalence_dict[current_isoforms_s][0] += 1
        else:
            equivalence_dict[current_isoforms_s] = [1, len(current_isoforms)]
    else:
        equivalence_dict[no_equivalence_class][0] += 1
file.close()

print()
print('SORTING THE EQUIVALENCE CLASSES')

output = sorted(
    [(count[0], count[1], equivalence_class) for equivalence_class, count in equivalence_dict.items()],
    key=lambda a: (a[1], a[2])
)

print('WRITING OUTPUT TO FILE')
with open(f'pseudoalignment-output-{KMER_LENGTH}.csv', 'w') as file:
    fhandle = csv.writer(file)
    heading = [('counts', 'number of items in equivalence class', 'isoforms in equivalence class')]
    fhandle.writerows(heading)
    fhandle.writerows(output)

print('DONE')
