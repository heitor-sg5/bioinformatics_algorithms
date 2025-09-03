from collections import defaultdict
import numpy as np
import itertools
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN

class ParseFiles:
    def parse_gff3(self, gff3_path, seqid_filter=None):
        genes_exons = defaultdict(list)
        used_seqid = None
        with open(gff3_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                seqid, source, ftype, start, end, score, strand, phase, attributes = parts[:9]
                if seqid_filter and seqid != seqid_filter:
                    continue
                if used_seqid is None:
                    used_seqid = seqid
                if ftype.lower() in ('exon','cds'):
                    start_i = int(start)
                    end_i = int(end)
                    attr = {k:v for k,v in (item.split('=') for item in attributes.split(';') if '=' in item)}
                    gene_id = None
                    if 'Parent' in attr:
                        gene_id = attr['Parent']
                    elif 'gene_id' in attr:
                        gene_id = attr['gene_id']
                    elif 'gene' in attr:
                        gene_id = attr['gene']
                    elif 'ID' in attr:
                        gene_id = attr['ID']
                    else:
                        gene_id = f"unknown_{seqid}_{start}_{end}"
                    genes_exons[gene_id].append((start_i, end_i))
        for gid, exons in genes_exons.items():
            exons_sorted = sorted(exons, key=lambda x: x[0])
            merged = []
            for s,e in exons_sorted:
                if not merged or s > merged[-1][1] + 1:
                    merged.append([s,e])
                else:
                    merged[-1][1] = max(merged[-1][1], e)
            genes_exons[gid] = [(a,b) for a,b in merged]
        return genes_exons, used_seqid
    
class Scores:
    def gc_content(self, seq):
        seq = seq.upper()
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq) if len(seq) > 0 else 0.0

    def extract_codon_usage(self, seq):
        seq = seq.upper()
        codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
        usage = defaultdict(int)
        for codon in codons:
            if set(codon) <= {'A', 'C', 'G', 'T'}:
                usage[codon] += 1
        return usage

    def compute_pca_features(self, codon_usage_list):
        codons = [''.join(p) for p in itertools.product('ACGT', repeat=3)]
        codon_index = {codon: i for i, codon in enumerate(codons)}
        mat = np.zeros((len(codon_usage_list), len(codons)), dtype=float)
        for i, usage in enumerate(codon_usage_list):
            for codon, count in usage.items():
                mat[i, codon_index[codon]] = count
        row_sums = mat.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        mat = mat / row_sums
        pca = PCA(n_components=2)
        coords = pca.fit_transform(mat)
        return coords

    def build_2nd_order_markov(self, sequence):
        sequence = sequence.upper().replace("U", "T")
        valid_bases = {'A', 'C', 'G', 'T'}
        sequence = ''.join([b for b in sequence if b in valid_bases])
        trinucs = [''.join(p) for p in itertools.product('ACGT', repeat=3)]
        trinuc_to_idx = {tri: idx for idx, tri in enumerate(trinucs)}
        transition_counts = np.zeros((64, 64), dtype=np.float64)
        for i in range(len(sequence) - 5):
            tri1 = sequence[i:i+3]
            tri2 = sequence[i+3:i+6]
            if tri1 in trinuc_to_idx and tri2 in trinuc_to_idx:
                idx1 = trinuc_to_idx[tri1]
                idx2 = trinuc_to_idx[tri2]
                transition_counts[idx1, idx2] += 1
        alpha = 1.0
        transition_counts += alpha
        row_sums = transition_counts.sum(axis=1, keepdims=True)
        transition_probs = transition_counts / row_sums
        return transition_probs, trinuc_to_idx

    def log_likelihood_markov(self, seq, transition_probs, trinuc_to_idx):
        seq = seq.upper().replace("U", "T")
        ll = 0.0
        for i in range(len(seq) - 5):
            tri1 = seq[i:i+3]; tri2 = seq[i+3:i+6]
            if tri1 in trinuc_to_idx and tri2 in trinuc_to_idx:
                idx1 = trinuc_to_idx[tri1]; idx2 = trinuc_to_idx[tri2]
                p = transition_probs[idx1, idx2]
                if p > 0:
                    ll += float(np.log(p))
        return ll
    
class WindowAnalysis(Scores):
    def run_window_based_hgt(self, ref_record, window_size):
        print(f"Running window-based HGT detection with window size {window_size}...")
        genome_seq = str(ref_record.seq).upper()
        genome_length = len(genome_seq)
        transition_probs, trinuc_to_idx = self.build_2nd_order_markov(genome_seq)
        genome_gc = self.gc_content(genome_seq)
        window_results = []
        codon_usage_list = []
        for start in range(0, genome_length, window_size):
            end = min(start + window_size, genome_length)
            window_seq = genome_seq[start:end]
            if len(window_seq) < 6:
                continue
            gc = self.gc_content(window_seq)
            gc_dev = abs(gc - genome_gc)
            markov_ll = self.log_likelihood_markov(window_seq, transition_probs, trinuc_to_idx)
            codon_usage = self.extract_codon_usage(window_seq)
            codon_usage_list.append(codon_usage)
            window_results.append({
                "gene_id": f"window_{start+1}_{end}",
                "start": start + 1,
                "end": end,
                "gc_content": gc,
                "gc_deviation": gc_dev,
                "markov_log_likelihood": markov_ll,
                "pca1": None,
                "pca2": None
            })
        if window_results:
            pca_coords = self.compute_pca_features(codon_usage_list)
            for i, coords in enumerate(pca_coords):
                window_results[i]["pca1"] = coords[0]
                window_results[i]["pca2"] = coords[1]
        print(f"Computed scores for {len(window_results)} windows.")
        return window_results
    
class GeneAnalysis(Scores):
    def run_gene_based_hgt(self, ref_record, genes_exons):
        print("Running gene-based HGT detection...")
        genome_seq = str(ref_record.seq)
        genome_gc = self.gc_content(genome_seq)
        transition_probs, trinuc_to_idx = self.build_2nd_order_markov(genome_seq)
        gene_results = []
        codon_usage_list = []
        for gene_id, exons in genes_exons.items():
            gene_seq = ''.join([genome_seq[start-1:end] for start, end in exons])
            if len(gene_seq) < 6:
                continue
            gc = self.gc_content(gene_seq)
            gc_dev = abs(gc - genome_gc)
            ll = self.log_likelihood_markov(gene_seq, transition_probs, trinuc_to_idx)
            codon_usage = self.extract_codon_usage(gene_seq)
            codon_usage_list.append(codon_usage)
            gene_results.append({
                "gene_id": gene_id,
                "start": exons[0][0],
                "end": exons[-1][1],
                "gc_content": gc,
                "gc_deviation": gc_dev,
                "markov_log_likelihood": ll,
                "pca1": None,
                "pca2": None
            })
        pca_coords = self.compute_pca_features(codon_usage_list)
        for i, coords in enumerate(pca_coords):
            gene_results[i]["pca1"] = coords[0]
            gene_results[i]["pca2"] = coords[1]
        print(f"Computed scores for {len(gene_results)} genes.")
        return gene_results
    
class Clusters:
    def cluster_genes_by_scores(self, genes, eps, min_samples):
        if not genes:
            return [], []
        X = np.asarray([
            [g.get("gc_content", 0.0),
            g.get("markov_log_likelihood", 0.0),
            g.get("pca1", 0.0)]
            for g in genes
        ], dtype=float)
        Xs = StandardScaler().fit_transform(X)
        db = DBSCAN(eps=eps, min_samples=min_samples, metric="euclidean", n_jobs=-1)
        labels = db.fit_predict(Xs)
        outliers = []
        for i, g in enumerate(genes):
            g["cluster_id"] = int(labels[i])
            g["outlier"] = (labels[i] == -1)
            if g["outlier"]:
                outliers.append(g)
        return genes, outliers

    def plot_hgt_results(self, genes, outliers, genome_length, bin_size=100000):
        if not genes:
            print("No genes/windows to plot.")
            return
        genes_p = [g for g in genes if g.get("pca1") is not None and g.get("pca2") is not None]
        outs_p  = [g for g in outliers if g.get("pca1") is not None and g.get("pca2") is not None]
        gc_list    = [g["gc_content"] for g in genes]
        markov_list= [g["markov_log_likelihood"] for g in genes]
        pca1_list  = [g["pca1"] for g in genes_p]
        pca2_list  = [g["pca2"] for g in genes_p]
        out_pca1   = [g["pca1"] for g in outs_p]
        out_pca2   = [g["pca2"] for g in outs_p]
        out_gc     = [g["gc_content"] for g in outliers]
        out_markov = [g["markov_log_likelihood"] for g in outliers]
        fig = plt.figure(figsize=(14, 10))
        gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])
        ax0 = fig.add_subplot(gs[0, 0])
        if pca1_list:
            ax0.scatter(pca1_list, pca2_list, alpha=0.6, label='Genes/Windows')
        if outliers and out_pca1:
            ax0.scatter(out_pca1, out_pca2, label='Outliers')
        ax0.set_xlabel("PCA1"); ax0.set_ylabel("PCA2"); ax0.set_title("PCA Scatter"); ax0.legend()
        ax1 = fig.add_subplot(gs[0, 1])
        ax1.scatter(gc_list, markov_list, alpha=0.6, label='Genes/Windows')
        if outliers:
            ax1.scatter(out_gc, out_markov, label='Outliers')
        ax1.set_xlabel("GC Content"); ax1.set_ylabel("Markov Log-Likelihood")
        ax1.set_title("GC% vs Markov Score"); ax1.legend()
        ax2 = fig.add_subplot(gs[1, :])
        if genome_length is None:
            genome_length = max(g["end"] for g in genes)
        bins = np.arange(0, genome_length + bin_size, bin_size)
        hgt_counts = np.zeros(len(bins)-1, dtype=int)
        for g in outliers:
            s, e = g["start"], g["end"]
            bi = min((s-1)//bin_size, len(hgt_counts)-1)
            bj = min((e-1)//bin_size, len(hgt_counts)-1)
            for k in range(bi, bj+1):
                hgt_counts[k] += 1
        centers = bins[:-1] + bin_size/2
        ax2.bar(centers, hgt_counts, width=bin_size*0.8)
        ax2.set_xlabel("Genome Position (bp)"); ax2.set_ylabel("HGT Count")
        ax2.set_title(f"HGT Count per {bin_size} bases")
        plt.suptitle("HGT Analysis")
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()

    def write_outlier_results(self, outliers, ref_record, is_gene=True):
        results = []
        results.append(f"Total outliers detected: {len(outliers)}\n\n")
        for idx, g in enumerate(sorted(outliers, key=lambda x: x["start"])):
            start = g["start"]
            end = g["end"]
            length = end - start + 1
            gc = g["gc_content"]
            markov = g.get("markov_log_likelihood", 0.0)
            pca1 = g.get("pca1", 0.0)
            pca2 = g.get("pca2", 0.0)
            seq = str(ref_record.seq[start-1:end])
            if is_gene:
                results.append(
                    f"Gene: {g['gene_id']} | Position: {start}-{end} | Length: {length} | "f"GC: {gc:.3f} | Score: {markov:.3f} | PCA: ({pca1:.3f}, {pca2:.3f})\n"
                    f"Sequence:\n{seq}\n\n"
                )
            else:
                results.append(
                    f"Window: {idx+1} | Position: {start}-{end} | Length: {length} | "
                    f"GC: {gc:.3f} | Score: {markov:.3f} | PCA: ({pca1:.3f}, {pca2:.3f})\n"
                    f"Sequence:\n{seq}\n\n"
                )
        return results