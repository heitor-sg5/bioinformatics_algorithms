import tkinter as tk
from tkinter import filedialog
import time
from motif_search import GreedyMotifSearch, RandomMotifSearch, GibbsSamplerMotifSearch

def get_user_inputs():
    k_input = input("Enter k (default 15): ").strip()
    n_input = input("Enter N (default 50): ").strip()
    k = int(k_input) if k_input else 15
    n = int(n_input) if n_input else 50
    return k, n

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
    print("Select a list of genomes FASTA for search")
    sequences = get_fasta_file()
    k, n = get_user_inputs()
    algorithms = [
        ('Greedy Motif Search', GreedyMotifSearch()),
        ('Random Motif Search', RandomMotifSearch()),
        ('Gibbs Sampler Motif Search', GibbsSamplerMotifSearch())
    ]
    results = []
    overall_start = time.time()
    for algo_name, algo in algorithms:
        print(f"Running {algo_name}...")
        start_time = time.time()
        motifs = algo.run(sequences, k, len(sequences), n)
        runtime = time.time() - start_time
        consensus_seq = algo.consensus(motifs)
        score_val = algo.score(motifs)
        results.append(f"---{algo_name}---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append(f"Consensus Sequence: {consensus_seq}\n")
        results.append(f"Score: {score_val}\n")
        results.append("Motifs:\n")
        for i, motif in enumerate(motifs, 1):
            results.append(f"  Sequence {i}: {motif}\n")
        results.append("\n")
    total_runtime = time.time() - overall_start
    results.append(f"Total Runtime: {total_runtime:.1f} seconds\n")
    with open('results.txt', 'w') as f:
        f.writelines(results)
    print("Results written to results.txt.")
    print(f"Total Runtime: {total_runtime:.1f} seconds")

if __name__ == "__main__":
    main()