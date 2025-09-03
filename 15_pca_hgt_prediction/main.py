import tkinter as tk
from tkinter import filedialog
import time
from Bio import SeqIO
import sys
from hgt_prediction import ParseFiles, GeneAnalysis, WindowAnalysis, Clusters

def get_user_inputs():
    root = tk.Tk()
    root.withdraw()
    print("Select Reference GFF3 file")
    gff3 = filedialog.askopenfilename(title="Select Reference GFF3 file")
    print("Select Reference FASTA file (same chromosome as GFF3)")
    ref_fa = filedialog.askopenfilename(title="Select Reference FASTA file")
    return gff3, ref_fa

def ask_analysis_mode():
    while True:
        select_mode = input("Select mode (gene (1) or window (0), default 1): ").strip()
        mode = int(select_mode) if select_mode else 1
        return mode

def main():
    gff3, ref_fa = get_user_inputs()
    ref_records = list(SeqIO.parse(ref_fa, 'fasta'))
    if len(ref_records) == 0:
        print('Reference FASTA appears empty.')
        sys.exit(1)
    ref_record = ref_records[0]
    parse = ParseFiles()
    start_time = time.time()
    genes_exons, used_seqid = parse.parse_gff3(gff3, seqid_filter=ref_record.id)
    if not genes_exons:
        genes_exons, used_seqid = parse.parse_gff3(gff3, seqid_filter=None)
        if not genes_exons:
            print('No exon features parsed from GFF3.')
            sys.exit(1)
    print(f'Parsed {len(genes_exons)} genes from GFF3.')
    mode = ask_analysis_mode()
    window = WindowAnalysis()
    gene = GeneAnalysis()
    genes = []
    if mode == 1:
        genes = gene.run_gene_based_hgt(ref_record, genes_exons)
    else:
        try:
            window_input = int(input("Enter sliding window size (default 1000): "))
            window_size = window_input if window_input else 1000
        except ValueError:
            print("Invalid window size. Exiting.")
            sys.exit(1)
        genes = window.run_window_based_hgt(ref_record, window_size)
    cluster = Clusters()
    clustered_genes, outliers = cluster.cluster_genes_by_scores(genes, eps=0.6, min_samples=5)
    cluster.plot_hgt_results(clustered_genes, outliers, genome_length=len(ref_record.seq))
    is_gene_analysis = (mode == 1)
    results = cluster.write_outlier_results(outliers, ref_record, is_gene=is_gene_analysis)
    total_runtime = time.time() - start_time
    results.append(f"Total Runtime: {total_runtime:.1f} seconds\n")
    with open('results.txt', 'w') as f:
        f.writelines(results)
    print("Results written to results.txt.")
    print(f"Total Runtime: {total_runtime:.1f} seconds")

if __name__ == "__main__": 
    main()