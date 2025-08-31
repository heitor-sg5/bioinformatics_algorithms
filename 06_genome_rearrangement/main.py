import tkinter as tk
from tkinter import filedialog
import time
from genome_rearrangement import SyntenyBlockConstruction, BreakpointSort, TwoBreakSort

def get_user_inputs():
    k_input = input("Enter k (default 12): ").strip()
    max_distance_input = input("Enter max distance (default 100): ").strip()
    min_size_input = input("Enter minimum size (default 5): ").strip()
    scale_input = input("Enter scale factor for dot plot (10, 100, 1000, 10000, 100000, 1000000, default 1000): ").strip()
    k = int(k_input) if k_input and k_input.isdigit() else 12
    max_distance = int(max_distance_input) if max_distance_input and max_distance_input.isdigit() else 100
    min_size = int(min_size_input) if min_size_input and min_size_input.isdigit() else 5
    scale_options = [10, 100, 1000, 10000, 100000, 1000000]
    scale = int(scale_input) if scale_input and scale_input.isdigit() and int(scale_input) in scale_options else 1000
    scale_unit = {10: 'b', 100: 'hb', 1000: 'kb', 10000: '10kb', 100000: '100kb', 1000000: 'Mb'}[scale]
    return k, max_distance, min_size, scale, scale_unit

def get_fasta_files():
    root = tk.Tk()
    root.withdraw()
    print("Select the first genome FASTA file:")
    file_path1 = filedialog.askopenfilename(
        title="Select the FIRST FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")]
    )
    if not file_path1:
        print("First file not selected. Using default 'bp_test.txt'.")
        file_path1 = 'bp_test1.txt'
    
    print("Select the second genome FASTA file:")
    file_path2 = filedialog.askopenfilename(
        title="Select the SECOND FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")]
    )
    if not file_path2:
        print("Second file not selected. Using default 'bp_test.txt'.")
        file_path2 = 'bp_test2.txt'
    try:
        seqs1 = load_fasta(file_path1)
        seqs2 = load_fasta(file_path2)
    except FileNotFoundError:
        print(f"{file_path1} or {file_path2} not found.")
        return
    return seqs1, seqs2

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
    print("Select two genomes FASTA for rearrangement")
    seqs1, seqs2 = get_fasta_files()
    k, max_distance, min_size, scale, scale_unit = get_user_inputs()
    is_pairwise = len(seqs1) == 1 and len(seqs2) == 1
    print(f"{'Pairwise' if is_pairwise else 'Multi-chromosome'} rearrangement selected.")
    results = []
    overall_start = time.time()
    print("Building synteny blocks...")
    synteny = SyntenyBlockConstruction(k, max_distance, min_size)
    if is_pairwise:
        print("Running Breakpoint Sort...")
        algo_name = "Breakpoint Sort with Synteny Blocks"
        start_time = time.time()
        shared_kmers = synteny.find_shared_kmers(seqs1[0], seqs2[0])
        blocks = synteny.synteny_blocks(shared_kmers)
        perm1, perm2 = synteny.signed_permutations(blocks)
        sorter = BreakpointSort()
        distance, steps = sorter.run(perm2)
        runtime = time.time() - start_time
        results.append(f"---{algo_name}---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append(f"Initial Genomes:\n")
        results.append(f"P: {sorter.format_perm(perm1)}\n")
        results.append(f"Q: {sorter.format_perm(perm2)}\n")
        for i, step in enumerate(steps, 1):
            results.append(f"Step {i}: {sorter.format_perm(step)}\n")
        results.append(f"Reversal Distance: {distance}\n\n")
        synteny.plot_dotplot(shared_kmers, len(seqs1[0]), len(seqs2[0]), scale, scale_unit)
    else:
        print("Running Two-Break Sort...")
        algo_name = "Two-Break Sort with Synteny Blocks"
        start_time = time.time()
        shared_kmers = synteny.find_shared_kmers_multichr(seqs1, seqs2)
        P_lists, Q_lists = synteny.permutations_grouped_by_chromosomes(seqs1, seqs2)
        sorter = TwoBreakSort()
        distance, steps = sorter.run(P_lists, Q_lists)
        runtime = time.time() - start_time
        results.append(f"---{algo_name}---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append(f"Initial Genomes:\n")
        results.append(f"P: {sorter.format_genome(P_lists)}\n")
        results.append(f"Q: {sorter.format_genome(Q_lists)}\n")
        for i, step in enumerate(steps, 1):
            results.append(f"Step {i}: {sorter.format_genome(step)}\n")
        results.append(f"2-Break Distance: {distance}\n\n")
        synteny.plot_dotplot_multichr(shared_kmers, seqs1, seqs2, scale, scale_unit)
    total_runtime = time.time() - overall_start
    results.append(f"Total Runtime: {total_runtime:.1f} seconds\n")
    with open('results.txt', 'w') as f:
        f.writelines(results)
    print("Results written to results.txt.")
    print(f"Total Runtime: {total_runtime:.1f} seconds")

if __name__ == "__main__":
    main()