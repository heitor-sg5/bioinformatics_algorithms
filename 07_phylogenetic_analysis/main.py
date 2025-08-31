import tkinter as tk
from tkinter import filedialog
import time
from phylogenetic_analysis import DistanceMatrixBase, UPGMA, NeighborJoining, SmallParsimonyAndNNI

def get_user_inputs():
    matrix_choice = input("Enter matrix type (0=Hamming (default), 1=Kimura, 2=Indel-Free): ").strip()
    matrix_type = int(matrix_choice) if matrix_choice and matrix_choice.isdigit() and int(matrix_choice) in [0, 1, 2] else 0
    matrix_name = {0: "Hamming Distance Matrix", 1: "Kimura Distance Matrix", 2: "Indel-Free Sequences"}[matrix_type]
    return matrix_type, matrix_name

def get_fasta_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")]
    )
    if not file_path:
        print("No file selected. Using default 'test.txt'.")
        file_path = 'test.txt'
        try:
            sequences = load_fasta(file_path)
        except FileNotFoundError:
            print("test.txt not found.")
            return
    else:
        sequences = load_fasta(file_path)
    if not sequences:
        print("No valid sequences found in the file.")
        return
    print(f"Loaded {len(sequences)} sequences.")
    return sequences

def load_fasta(file_path):
    sequences = []
    current_seq = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append(''.join(current_seq))
    return sequences

def main():
    print("Select a list of genomes FASTA for tree construction")
    sequences = get_fasta_file()
    matrix_type, matrix_name = get_user_inputs()
    results = []
    overall_start = time.time()
    distance_base = DistanceMatrixBase()
    if matrix_type == 2:
        print("Running Small Parsimony with NNI...")
        algo_name = "Small Parsimony with NNI"
        start_time = time.time()
        aligned_seqs = distance_base.remove_indel_columns(sequences)
        matrix_str = "Indel-free Matrix:\n"
        for seq in aligned_seqs:
            matrix_str += "".join(seq) + "\n"
        results.append(matrix_str + "\n")
        algo = SmallParsimonyAndNNI()
        tree = algo.run(aligned_seqs)
        runtime = time.time() - start_time
        results.append(f"--- {algo_name} ---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append("Tree:\n")
        tree_str = algo.format_tree(tree)
        results.append(tree_str + "\n")
        results.append("Newick Tree:\n")
        newick_str = algo.format_newick(tree)
        results.append(newick_str + "\n")
        results.append("\n")
    else:
        start_time = time.time()
        print("Building distance matix...")
        dist_matrix = distance_base.build_distance_matrix(sequences, matrix_type)
        matrix_str = f"{matrix_name}:\n"
        for row in dist_matrix:
            matrix_str += f"{' '.join(map(str, row))}\n"
        results.append(matrix_str + "\n")
        algorithms = [
            ("UPGMA", UPGMA()),
            ("Neighbor Joining", NeighborJoining())
        ]
        for algo_name, algo in algorithms:
            print(f"Running {algo_name}...")
            start_time = time.time()
            tree = algo.run(dist_matrix, list(range(len(sequences))))
            runtime = time.time() - start_time
            results.append(f"--- {algo_name} ---\n")
            results.append(f"Runtime: {runtime:.1f} seconds\n")
            results.append("Tree:\n")
            tree_str = algo.format_tree(tree)
            results.append(tree_str + "\n")
            results.append("Newick Tree:\n")
            newick_str = algo.format_newick(tree)
            results.append(newick_str + "\n")
            results.append("\n")
    total_runtime = time.time() - overall_start
    results.append(f"Total Runtime: {total_runtime:.1f} seconds\n")
    with open('results.txt', 'w') as f:
        f.writelines(results)
    print("Results written to results.txt.")
    print(f"Total Runtime: {total_runtime:.1f} seconds")

if __name__ == "__main__":
    main()