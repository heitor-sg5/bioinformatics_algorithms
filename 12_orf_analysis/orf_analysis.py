import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import heapq
import math

class ORFBase:
    base_to_int_map = {'A':0,'C':1,'G':2,'T':3}

    codons = [
        'AAA','AAC','AAG','AAT','ACA','ACC','ACG','ACT',
        'AGA','AGC','AGG','AGT','ATA','ATC','ATG','ATT', 
        'CAA','CAC','CAG','CAT','CCA','CCC','CCG','CCT',
        'CGA','CGC','CGG','CGT','CTA','CTC','CTG','CTT',
        'GAA','GAC','GAG','GAT','GCA','GCC','GCG','GCT',
        'GGA','GGC','GGG','GGT','GTA','GTC','GTG','GTT',
        'TAA','TAC','TAG','TAT','TCA','TCC','TCG','TCT',
        'TGA','TGC','TGG','TGT','TTA','TTC','TTG','TTT'
    ]

    def reverse_complement(self, seq):
        complement = {'A':'T','T':'A','C':'G','G':'C'}
        return ''.join(complement[b] for b in reversed(seq))
    
    def build_2nd_markov(self, leaderboard, codon_to_vector):
        transition_counts = np.zeros((64, 64, 64))
        for entry in leaderboard:
            seq = entry["seq"]
            for i in range(0, len(seq) - 9, 3):
                c1 = codon_to_vector[seq[i:i+3]]
                c2 = codon_to_vector[seq[i+3:i+6]]
                c3 = codon_to_vector[seq[i+6:i+9]]
                transition_counts[c1, c2, c3] += 1
        transition_probs = transition_counts + 1e-8
        transition_probs /= transition_probs.sum(axis=2, keepdims=True)
        return transition_probs

    def markov_score(self, orf, codon_to_vector, transition_probs):
        log_prob = 0.0
        transitions = 0
        for i in range(0, len(orf) - 9, 3):
            c1 = codon_to_vector[orf[i:i+3]]
            c2 = codon_to_vector[orf[i+3:i+6]]
            c3 = codon_to_vector[orf[i+6:i+9]]
            log_prob += math.log(transition_probs[c1, c2, c3])
            transitions += 1
        return log_prob / max(transitions, 1)

class FirstPassORF(ORFBase):
    def genome_to_vectors(self, genome, rev_genome, codon_to_vector):
        frames = []
        for seq in [genome, rev_genome]:
            for offset in range(3):
                codons = np.array([codon_to_vector[seq[i:i+3]] for i in range(offset, len(seq)-2, 3)])
                frames.append(codons)
        return frames
    
    def find_start_stop_codons(self, frames, codon_to_vector):
        start_codons = {codon_to_vector[c] for c in ['GTG', 'TTG', 'ATG']}
        stop_codons = {codon_to_vector[c] for c in ['TGA', 'TAG', 'TAA']}
        orfs_per_frame = [[] for _ in range(6)]
        for idx, frame in enumerate(frames):
            start_mask = np.isin(frame, list(start_codons))
            stop_mask = np.isin(frame, list(stop_codons))
            start_positions = np.flatnonzero(start_mask)
            stop_positions = np.flatnonzero(stop_mask)
            stop_idx = 0
            for start in start_positions:
                while stop_idx < len(stop_positions) and stop_positions[stop_idx] <= start:
                    stop_idx += 1
                if stop_idx >= len(stop_positions):
                    continue
                stop = stop_positions[stop_idx]
                orfs_per_frame[idx].append((start, stop))
        return orfs_per_frame
    
    def filter_orfs(self, orfs, genome, min_size):
        new_orfs = []
        for orf_list in orfs:
            for start_idx, stop_idx in orf_list:
                nt_start = start_idx * 3
                nt_stop = stop_idx * 3 + 3
                orf_seq = genome[nt_start:nt_stop]
                orf_len = len(orf_seq) // 3
                if orf_len < min_size:
                    continue
                score = 0
                new_orfs.append({"seq": orf_seq, "score": score})
        return new_orfs

class SecondPassORF(ORFBase):
    def build_leaderboard(self, orfs, genome, min_size, L, codon_to_vector, transition_probs):
        leaderboard = []
        for orf_list in orfs:
            for start_idx, stop_idx in orf_list:
                nt_start = start_idx * 3
                nt_stop = stop_idx * 3 + 3
                orf_seq = genome[nt_start:nt_stop]
                orf_len = len(orf_seq) // 3
                if orf_len < min_size:
                    continue
                score = self.markov_score(orf_seq, codon_to_vector, transition_probs)
                if len(leaderboard) < L:
                    heapq.heappush(leaderboard, (score, {"seq": orf_seq, "score": score}))
                else:
                    if score > leaderboard[0][0]:
                        heapq.heappushpop(leaderboard, (score, {"seq": orf_seq, "score": score}))
        return [item[1] for item in sorted(leaderboard, key=lambda x: x[0], reverse=True)]

    def find_orfs(self, orfs, genome, rev_genome, codon_to_vector, min_size, transition_probs):
        n = len(genome)
        new_orfs = []
        for frame_idx, orf_list in enumerate(orfs):
            if frame_idx < 3:
                strand = "+"
                frame_number = frame_idx + 1
                for start_idx, stop_idx in orf_list:
                    nt_start = start_idx * 3 + frame_idx
                    nt_stop = stop_idx * 3 + frame_idx + 3
                    orf_seq = genome[nt_start:nt_stop]
                    orf_len = len(orf_seq) // 3
                    if orf_len < min_size:
                        continue
                    score = self.markov_score(orf_seq, codon_to_vector, transition_probs)
                    new_orfs.append({
                        "seq": orf_seq,
                        "len": orf_len,
                        "start": nt_start,
                        "end": nt_stop,
                        "strand": strand,
                        "frame": frame_number,
                        "score": score
                    })
            else:
                strand = "-"
                frame_number = frame_idx - 2
                for start_idx, stop_idx in orf_list:
                    nt_stop = n - start_idx * 3 - frame_number
                    nt_start = n - stop_idx * 3 - frame_number - 3
                    orf_seq = rev_genome[start_idx*3 + (frame_idx-3): stop_idx*3 + (frame_idx-3) + 3]
                    orf_len = len(orf_seq) // 3
                    if orf_len < min_size:
                        continue
                    score = self.markov_score(orf_seq, codon_to_vector, transition_probs)
                    new_orfs.append({
                        "seq": orf_seq,
                        "len": orf_len,
                        "start": nt_start,
                        "end": nt_stop,
                        "strand": strand,
                        "frame": frame_number,
                        "score": score
                    })
        return new_orfs

    def orf_overlap(self, orfs, max_overlap, t_type):
        scores = np.array([orf['score'] for orf in orfs])
        if t_type == 0:
            mean = np.mean(scores)
            std = np.std(scores)
            t = mean + 1 * std
        else:
            t = np.percentile(scores, t_type)
        if max_overlap is None:
            return orfs
        orfs_sorted = sorted(orfs, key=lambda x: x["score"], reverse=True)
        filtered_orfs = []
        used_regions = []
        for orf in orfs_sorted:
            overlaps = 0
            if orf["score"] < t:
                continue
            for region in used_regions:
                if not (orf["end"] <= region[0] or orf["start"] >= region[1]):
                    overlaps += 1
                    if overlaps >= max_overlap:
                        break
            if overlaps < max_overlap:
                filtered_orfs.append(orf)
                used_regions.append((orf["start"], orf["end"]))
        return sorted(filtered_orfs, key=lambda x: x["start"])

class TwoPassORF(FirstPassORF, SecondPassORF):
    def first_scan(self, genome, rev_genome, codon_to_vector, min_size):
        frames = self.genome_to_vectors(genome, rev_genome, codon_to_vector)
        orfs = self.find_start_stop_codons(frames, codon_to_vector)
        filtered_orfs = self.filter_orfs(orfs, genome, min_size)
        transitions = self.build_2nd_markov(filtered_orfs, codon_to_vector)
        return transitions, orfs
    
    def second_scan(self, genome, rev_genome, codon_to_vector, min_size, max_overlap, t, L, orfs, transitions):
        leaderboard = self.build_leaderboard(orfs, genome, min_size, L, codon_to_vector, transitions)
        new_transitions = self.build_2nd_markov(leaderboard, codon_to_vector)
        new_orfs = self.find_orfs(orfs, genome, rev_genome, codon_to_vector, min_size, new_transitions)
        filtered = self.orf_overlap(new_orfs, max_overlap, t)
        return filtered
    
    def two_pass(self, genome, min_size, max_overlap, t, L):
        codon_to_vector = {c: self.base_to_int_map[c[0]]*16 + self.base_to_int_map[c[1]]*4 + self.base_to_int_map[c[2]] for c in self.codons}
        rev_genome = self.reverse_complement(genome)
        transitions, orfs = self.first_scan(genome, rev_genome, codon_to_vector, min_size)
        new_orfs = self.second_scan(genome, rev_genome, codon_to_vector, min_size, max_overlap, t, L, orfs, transitions)
        return new_orfs

class Charts(ORFBase):
    def score_dist(self, scores):
        plt.hist(scores, bins=100)
        plt.xlabel("Score")
        plt.ylabel("ORFs")
        plt.show()
    
    def orf_density(self, genome, orfs):
        window_size = 10000
        genome_len = len(genome)
        density = [0] * (genome_len // window_size + 1)
        for orf in orfs:
            start_bin = orf['start'] // window_size
            end_bin = orf['end'] // window_size
            for b in range(start_bin, end_bin + 1):
                density[b] += 1
        plt.bar(range(len(density)), density, width=1.0)
        plt.xlabel(f"Position (per {window_size} bp)")
        plt.ylabel("ORFs")
        plt.title("ORF Density")
        plt.show()

    def codon_freq(self, orfs):
        codon_counts = Counter()
        for orf in orfs:
            seq = orf['seq']
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                codon_counts[codon] += 1
        total_codons = sum(codon_counts.values())
        return {c: count/total_codons for c, count in codon_counts.items()}
    
    def summary_info(self, orfs, scores, lengths, frequencies):
        results = []
        results.append(f"---Summary---\n")
        results.append(f"Mean score: {np.mean(scores):.1f} | Min: {np.min(scores):.1f} | Max: {np.max(scores):.1f}\n")
        results.append(f"Mean length: {np.mean(lengths) * 3:.0f} ({np.mean(lengths):.0f}) | Min: {np.min(lengths) * 3:.0f} ({np.min(lengths)}) | Max: {np.max(lengths) * 3:.0f} ({np.max(lengths)})\n")
        pos = sum(1 for orf in orfs if orf['strand'] == '+')
        neg = sum(1 for orf in orfs if orf['strand'] == '-')
        f1 = sum(1 for orf in orfs if orf['frame'] == 1)
        f2 = sum(1 for orf in orfs if orf['frame'] == 2)
        f3 = sum(1 for orf in orfs if orf['frame'] == 3)
        results.append(f"Strand distribution (+/-): {pos}/{neg} | Frame ratio (1:2:3): {f1}:{f2}:{f3}\n")
        results.append(f"Codon frequencies:\n")
        line = []
        for i, codon in enumerate(self.codons, start=1):
            freq = frequencies.get(codon, 0)
            line.append(f"{codon}: {freq:.2f}")
            if i % 8 == 0:
                results.append(" | ".join(line)+ "\n")
                line = []
        if line:
            results.append(" | ".join(line)+ "\n")
        results.append("\n")
        return results
    
    def display(self, orfs, genome):
        scores = [orf["score"] for orf in orfs]
        lengths = [orf["len"] for orf in orfs]
        self.score_dist(scores)
        self.orf_density(genome, orfs)
        freq = self.codon_freq(orfs)
        summary = self.summary_info(orfs, scores, lengths, freq)
        return summary