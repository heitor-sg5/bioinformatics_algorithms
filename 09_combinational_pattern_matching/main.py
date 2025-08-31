import tkinter as tk
from tkinter import filedialog
import time
from pattern_matching import PatternMatchingBase, BWTMatching, SuffixArrayMatching, PrefixTrieMatching, ApproximatePatternMatching, SuffixTree, GeneralizedSuffixTree

def get_fasta_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")]
    )
    if not file_path:
        print("No file selected. Using default 'test1.txt'.")
        file_path = 'test1.txt'
        try:
            text = load_fasta(file_path)
        except FileNotFoundError:
            print("test1.txt not found.")
            return
    else:
        text = load_fasta(file_path)
    print("Genome loaded.")
    return text

def load_fasta(file_path):
    with open(file_path, 'r') as file:
        text = ''.join(line.strip() for line in file if not line.startswith(">"))
    return text

def get_reads_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select reads file",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    if not file_path:
        print("No file selected. Using default 'reads.txt'.")
        file_path = 'reads.txt'
        try:
            reads = load_reads(file_path)
        except FileNotFoundError:
            print("reads.txt not found.")
            return
    else:
        reads = load_reads(file_path)
    print("Reads loaded.")
    return reads

def load_reads(file_path):
    with open(file_path, 'r') as f:
        return f.read().strip().split()

def main():
    analysis_type = 2
    base = PatternMatchingBase()
    results = []
    if analysis_type == 1:
        print("Select a genome FASTA for search")
        text = get_fasta_file()
        pattern = input("Enter pattern to find: ").strip()
        overall_start = time.time()
        bwt = base.burrows_wheeler_transform(text)
        results.append(f"BWT: {bwt}\n")
        results.append("\n")
        algorithms = [
            ("BWT Matching", BWTMatching()),
            ("Suffix Array Matching", SuffixArrayMatching())
        ]
        for algo_name, algo in algorithms:
            print(f"Running {algo_name}...")
            start_time = time.time()
            positions, count = algo.run(text, pattern)
            runtime = time.time() - start_time
            results.append(f"--- {algo_name} ---\n")
            results.append(f"Runtime: {runtime:.1f} seconds\n")
            results.append(f"Pattern: {pattern}\n")
            results.append(f"Total occurrences: {count}\n")
            results.append(f"Positions: {positions}\n\n")
    elif analysis_type == 2:
        print("Select a genome FASTA for search")
        text = get_fasta_file()
        print("Select a text file with reads for search")
        reads = get_reads_file()
        d_input = input("Enter d (default 1): ").strip()
        d = int(d_input) if d_input.isdigit() else 1
        overall_start = time.time()
        bwt = base.burrows_wheeler_transform(text)
        results.append(f"BWT: {bwt}\n")
        results.append("\n")
        algorithms = [
            ("Prefix Trie Matching", PrefixTrieMatching()),
            ("Approximate Pattern Matching", ApproximatePatternMatching(d))
        ]
        for algo_name, algo in algorithms:
            print(f"Running {algo_name}...")
            start_time = time.time()
            result = algo.run(text, reads)
            runtime = time.time() - start_time
            results.append(f"--- {algo_name} ---\n")
            results.append(f"Runtime: {runtime:.1f} seconds\n")
            results.append("Patterns and Positions:\n")
            for pattern, (positions, count) in result.items():
                if count != 0:
                    results.append(f"Pattern: {pattern}\n")
                    results.append(f"Total occurrences: {count}\n")
                    results.append(f"Positions: {positions}\n")
            results.append("\n")
    else:
        print("Select a genome FASTA for comparison")
        text1 = get_fasta_file()
        print("Select a second genome FASTA for comparison")
        text2 = get_fasta_file()
        overall_start = time.time()
        start_time = time.time()
        suffix_tree = SuffixTree()
        longest_substr, length = suffix_tree.run(text1)
        runtime = time.time() - start_time
        results.append("--- Suffix Tree (Text 1) ---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append(f"Longest repeated substring: {longest_substr}\n")
        results.append(f"Length: {length}\n\n")
        start_time = time.time()
        longest_substr, length = suffix_tree.run(text2)
        runtime = time.time() - start_time
        results.append("--- Suffix Tree (Text 2) ---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append(f"Longest repeated substring: {longest_substr}\n")
        results.append(f"Length: {length}\n\n")
        start_time = time.time()
        gen_suffix_tree = GeneralizedSuffixTree()
        longest_shared, shortest_non_shared = gen_suffix_tree.run(text1, text2)
        runtime = time.time() - start_time
        results.append("--- Generalized Suffix Tree ---\n")
        results.append(f"Runtime: {runtime:.1f} seconds\n")
        results.append(f"Longest shared substring: {longest_shared}\n")
        results.append(f"Shortest non-shared substring: {shortest_non_shared}\n\n")
    total_runtime = time.time() - overall_start
    results.append(f"Total Runtime: {total_runtime:.1f} seconds\n")
    with open('results.txt', 'w') as f:
        f.writelines(results)
    print("Results written to results.txt.")
    print(f"Total Runtime: {total_runtime:.1f} seconds")

if __name__ == "__main__":
    main()