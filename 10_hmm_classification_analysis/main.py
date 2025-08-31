import tkinter as tk
from tkinter import filedialog
import time
import numpy as np
from hmm_classification import ProfileHMM, FindPath, BaumWelch, Viterbi, MostDivergent

def get_user_inputs():
    iterations_input = input("Enter iterations (default 5): ").strip()
    theta_input = input("Enter theta (default 0.5): ").strip()
    pseudocount_input = input("Enter pseudocount (default 0.01): ").strip()
    iterations = int(iterations_input) if iterations_input and iterations_input.replace('.', '', 1).isdigit() else 5
    theta = float(theta_input) if theta_input and theta_input.replace('.', '', 1).isdigit() else 0.5
    pseudocount = float(pseudocount_input) if pseudocount_input and pseudocount_input.replace('.', '', 1).isdigit() else 0.01
    return theta, pseudocount, iterations

def get_fasta_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select FASTA file with single sequence (X)",
        filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
    )
    if not file_path:
        print("No single-sequence file selected. Using default 'test1.txt'.")
        file_path = "test2.txt"
    try:
        text = load_fasta(file_path)
    except FileNotFoundError:
        print(f"{file_path} not found.")
        return
    print(f"Genome loaded.")
    return text

def load_fasta(file_path):
    with open(file_path, 'r') as file:
        text = ''.join(line.strip() for line in file if not line.startswith(">"))
    return text

def get_fasta_files():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select FASTA file with multiple sequences",
        filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
    )
    if not file_path:
        print("No multi-sequence file selected. Using default 'test1.txt'.")
        file_path = "test1.txt"
    try:
        sequences = load_fastas(file_path)
    except FileNotFoundError:
        print(f"{file_path} not found.")
        return
    print(f"Loaded {len(sequences)} sequences.")
    return sequences

def load_fastas(file_path):
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
    print("Select a list of genomes FASTA for alignment")
    sequences = get_fasta_files()
    print("Select a genome FASTA for classification")
    x = get_fasta_file()
    theta, pseudocount, iterations = get_user_inputs()
    phmm = ProfileHMM()
    print("Building profile HMM from training alignment...")
    hmm, states, transition, emission = phmm.build_profile(x, theta, pseudocount, sequences)
    with open("hmm_matrices.txt", "w") as f:
        f.write("---Initialized HMM---\n")
        f.write(hmm + "\n")
    results = []
    overall_start = time.time()
    start_time = time.time()
    print("Finding hidden path...")
    path, prob = FindPath(states, transition, emission, iterations).run(x)
    runtime = time.time() - start_time
    results.append("---Initialized HMM---\n")
    results.append(f"Runtime: {runtime:.1f} seconds\n")
    results.append(f"Hidden Path: {path}\n")
    results.append(f"Probability: {prob:.3g}\n")
    results.append("\n")
    algorithms = [
        ("Baum-Welch HMM", BaumWelch(states, transition, emission, iterations)),
        ("Viterbi HMM", Viterbi(states, transition, emission, iterations)),
        ("Most Divergent", MostDivergent(states, transition, emission, iterations))
    ]
    for algo_name, algo in algorithms:
        print(f"Running {algo_name}...")
        if algo_name == "Most Divergent":
            start_time = time.time()
            prob, index = algo.run(sequences)
            runtime = time.time() - start_time
            results.append(f"---{algo_name}---\n")
            results.append(f"Runtime: {runtime:.1f} seconds\n")
            results.append(f"Sequence: {index}\n")
            results.append(f"Probability: {prob:.3g}\n\n")
            continue
        start_time = time.time()
        hmm, prob, path = algo.run(x)
        runtime = time.time() - start_time
        with open("hmm_matrices.txt", "a") as f:
            f.write(f"---{algo_name} Improved HMM---\n")
            f.write(hmm + "\n")
        results.append(f"---{algo_name}---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append(f"Hidden Path: {path}\n")
        results.append(f"Probability: {prob:.3g}\n\n")
    total_runtime = time.time() - overall_start
    results.append(f"Total Runtime: {total_runtime:.1f} seconds\n")
    with open("results.txt", "w") as f:
        f.writelines(results)
    print("Results written to results.txt.")
    print(f"Total Runtime: {total_runtime:.1f} seconds")

if __name__ == "__main__":
    main()