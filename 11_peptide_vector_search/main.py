import tkinter as tk
from tkinter import filedialog
import ast
import time
from peptide_search import PeptideBase, PeptideSearch, StatisticalSignificance, SpectralAlignment

def get_user_inputs():
    k_input = input("Enter k (default 3): ").strip()
    d_input = input("Enter d (default 1): ").strip()
    k = int(k_input) if k_input else 3
    d = int(d_input) if d_input else 1
    return k, d

def get_spectra_file():
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select spectra file",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    if not file_path:
        print("No file selected. Using default 'test1.txt'.")
        file_path = 'test1.txt'
        try:
            spectra = load_spectra(file_path)
        except FileNotFoundError:
            print("test2.txt not found.")
            return
    else:
        spectra = load_spectra(file_path)
    if not spectra:
        print("No valid spectrum found in the file.")
        return
    print(f"Loaded {len(spectra)} spectra.")
    return spectra

def load_spectra(file_path):
    spectra = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            spectrum = ast.literal_eval(line)
            spectra.append(spectrum)
    return spectra

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
            peptides = load_fasta(file_path)
        except FileNotFoundError:
            print("test2.txt not found.")
            return
    else:
        peptides = load_fasta(file_path)
    if not peptides:
        print("No valid sequences found in the file.")
        return
    print(f"Loaded {len(peptides)} peptides.")
    return peptides

def load_fasta(file_path):
    peptides = []
    current_seq = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    peptides.append(''.join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            peptides.append(''.join(current_seq))
    return peptides

def main():
    print("Select a text file with your spectra for search")
    spectra = get_spectra_file()
    print("Select a peptide database FASTA for search")
    peptides = get_fasta_file()
    k, d = get_user_inputs()
    print("Building vectors...")
    spectral_vectors = [PeptideBase.spectrum_to_vector(s) for s in spectra]
    results = []
    overall_start = time.time()
    print("Searching for peptides...")
    ps = PeptideSearch(k, peptides, spectral_vectors)
    pms = ps.run()
    for i, match in enumerate(pms, 1):
        start_time = time.time()
        results.append(f"---PMS {i}---\n")
        results.append(f"Spectrum vector: {match[1]}\n")
        results.append(f"Best matching peptide: {match[0]}\n")
        results.append(f"Score: {match[2]}\n")
        if match[2] != 0:
            stat = StatisticalSignificance(k, match[1])
            align = SpectralAlignment(d, match[0], match[1])
            p = stat.run()
            mod = align.run()
            results.append(f"Probability: {p:.4f}\n")
            results.append(f"Spectral alignment: {mod}\n")
        runtime = time.time() - start_time
        results.append(f"Runtime: {runtime:.1f} seconds\n\n")
    total_runtime = time.time() - overall_start
    results.append(f"Total runtime: {total_runtime:.1f} seconds\n")
    with open("results.txt", "w") as f:
        f.writelines(results)
    print("Results written to results.txt.")
    print(f"Total runtime: {total_runtime:.1f} seconds")

if __name__ == "__main__":
    main()
