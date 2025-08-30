import numpy as np
import matplotlib.pyplot as plt

class CpGBase:
    base_to_int_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    def gc_content(self, seq):
        seq_array = np.frombuffer(seq.encode(), dtype=np.uint8)
        gc_count = np.sum((seq_array == ord('G')) | (seq_array == ord('C')))
        return gc_count / len(seq) if len(seq) > 0 else 0.0

    def obs_exp_cpg(self, seq):
        seq_array = np.frombuffer(seq.encode(), dtype=np.uint8)
        L = len(seq)
        if L < 2:
            return 0.0
        c_count = np.sum(seq_array == ord('C'))
        g_count = np.sum(seq_array == ord('G'))
        cpg_count = np.sum((seq_array[:-1] == ord('C')) & (seq_array[1:] == ord('G')))
        if c_count == 0 or g_count == 0 or L == 0:
            return 0.0
        return (cpg_count * L) / (c_count * g_count)

    def build_1st_markov(self, seq):
        seq_array = np.frombuffer(seq.encode(), dtype=np.uint8)
        trans_counts = np.zeros((4, 4), dtype=np.float64)
        valid_bases = np.isin(seq_array, [ord(b) for b in self.base_to_int_map])
        base_indices = np.array([self.base_to_int_map.get(chr(b), -1) for b in seq_array[valid_bases]], dtype=np.int32)
        if len(base_indices) < 2:
            trans_counts += 1.0
            return trans_counts / trans_counts.sum(axis=1, keepdims=True)
        pairs = np.vstack((base_indices[:-1], base_indices[1:])).T
        np.add.at(trans_counts, (pairs[:, 0], pairs[:, 1]), 1)
        alpha = 1.0
        trans_counts += alpha
        return trans_counts / trans_counts.sum(axis=1, keepdims=True)

    def llr_score(self, region_seq, p_region, p_bg):
        seq_array = np.frombuffer(region_seq.encode(), dtype=np.uint8)
        valid_bases = np.isin(seq_array, [ord(b) for b in self.base_to_int_map])
        base_indices = np.array([self.base_to_int_map.get(chr(b), -1) for b in seq_array[valid_bases]], dtype=np.int32)
        if len(base_indices) < 2:
            return 0.0
        pairs = np.vstack((base_indices[:-1], base_indices[1:])).T
        valid_pairs = (pairs[:, 0] >= 0) & (pairs[:, 1] >= 0)
        pairs = pairs[valid_pairs]
        llr = np.sum(np.log(p_region[pairs[:, 0], pairs[:, 1]] / np.maximum(p_bg[pairs[:, 0], pairs[:, 1]], 1e-10)))
        return llr

    def get_threshold(self, genome, p_bg, L, s, t):
        n_windows = (len(genome) - L + 1) // s + 1
        scores = np.zeros(n_windows, dtype=np.float64)
        for idx, i in enumerate(range(0, len(genome) - L + 1, s)):
            window = genome[i:i + L]
            p_region = self.build_1st_markov(window)
            scores[idx] = self.llr_score(window, p_region, p_bg)
        return np.percentile(scores, t)

class FindCpGIslands(CpGBase):
    def extend_island(self, genome, start, end, p_bg, w, p, min_gc, min_obsexp, min_llr, min_length, max_length):
        start = max(0, start)
        end = min(len(genome), end)
        region = genome[start:end]
        if len(region) < min_length:
            return None
        if len(region) > max_length:
            end = start + max_length
            region = genome[start:end]
        p_region = self.build_1st_markov(region)
        gc = self.gc_content(region)
        obsexp = self.obs_exp_cpg(region)
        llr = self.llr_score(region, p_region, p_bg)
        if gc < min_gc or obsexp < min_obsexp or llr < min_llr:
            return None
        if p <= 0:
            return {"seq": region, 
                    "len": len(region), 
                    "start": start, 
                    "end": end, 
                    "gc": gc, 
                    "obsexp": obsexp, 
                    "llr": llr
                }
        new_start = max(0, start - w)
        if new_start < start:
            left_extended = self.extend_island(genome, new_start, end, p_bg, w, p - 1, min_gc, min_obsexp, min_llr, min_length, max_length)
            if left_extended:
                start, end, region, gc, obsexp, llr = (left_extended["start"], left_extended["end"], left_extended["seq"], left_extended["gc"], left_extended["obsexp"], left_extended["llr"])
        new_end = min(len(genome), end + w)
        if new_end > end:
            right_extended = self.extend_island(genome, start, new_end, p_bg, w, p - 1, min_gc, min_obsexp, min_llr, min_length, max_length)
            if right_extended:
                start, end, region, gc, obsexp, llr = (right_extended["start"], right_extended["end"], right_extended["seq"], right_extended["gc"], right_extended["obsexp"], right_extended["llr"])
        return {"seq": region, 
                "len": len(region), 
                "start": start, 
                "end": end, 
                "gc": gc, 
                "obsexp": obsexp, 
                "llr": llr}

    def find_islands(self, genome, p_bg, L, s, w, p, min_gc, min_obsexp, min_llr, min_length, max_length):
        cpg_islands = []
        for i in range(0, len(genome) - L + 1, s):
            window = genome[i:i + L]
            p_region = self.build_1st_markov(window)
            gc = self.gc_content(window)
            obsexp = self.obs_exp_cpg(window)
            llr = self.llr_score(window, p_region, p_bg)
            if gc >= min_gc and obsexp >= min_obsexp and llr >= min_llr:
                island = self.extend_island(genome, i, i + L, p_bg, w, p, min_gc, min_obsexp, min_llr, min_length, max_length)
                if island:
                    cpg_islands.append(island)
        return cpg_islands

    def filter_overlap(self, cpg_islands, merge_gap, max_length, p_bg):
        if not cpg_islands:
            return []
        sorted_islands = sorted(cpg_islands, key=lambda x: x["start"])
        merged_islands = []
        current = sorted_islands[0].copy()
        for island in sorted_islands[1:]:
            if island["start"] <= current["end"] + merge_gap:
                new_start = current["start"]
                new_end = min(max(current["end"], island["end"]), new_start + max_length)
                seq_end = new_end - new_start
                new_seq = (current["seq"] + island["seq"])[0:seq_end]
                gc = self.gc_content(new_seq)
                obsexp = self.obs_exp_cpg(new_seq)
                p_region = self.build_1st_markov(new_seq)
                llr = self.llr_score(new_seq, p_region, p_bg)
                if gc >= min(current["gc"], island["gc"]) and obsexp >= min(current["obsexp"], island["obsexp"]) and llr >= min(current["llr"], island["llr"]):
                    current.update({
                        "seq": new_seq,
                        "len": len(new_seq),
                        "start": new_start,
                        "end": new_end,
                        "gc": gc,
                        "obsexp": obsexp,
                        "llr": llr
                    })
                else:
                    merged_islands.append(current)
                    current = island.copy()
            else:
                merged_islands.append(current)
                current = island.copy()
        merged_islands.append(current)
        return merged_islands

    def refine_island(self, genome, island, p_bg, min_gc, min_obsexp, min_llr, min_length, max_length):
        orig_start, orig_end = island["start"], island["end"]
        orig_gc, orig_obsexp, orig_llr = island["gc"], island["obsexp"], island["llr"]
        best_island = island.copy()
        best_score = (orig_gc / max(orig_gc, 1e-10)) + (orig_obsexp / max(orig_obsexp, 1e-10)) + (orig_llr / max(orig_llr, 1e-10))
        max_step = max(1, (orig_end - orig_start) // 10)
        min_step = 1
        step = max_step
        while step >= min_step:
            improved = False
            for delta_start, delta_end in [(step, 0), (0, step), (step, step), (0, 0)]:
                new_start = max(0, orig_start + delta_start)
                new_end = min(len(genome), orig_end - delta_end)
                if new_end - new_start < min_length or new_end - new_start > max_length:
                    continue
                new_seq = genome[new_start:new_end]
                if not new_seq:
                    continue
                p_region = self.build_1st_markov(new_seq)
                gc = self.gc_content(new_seq)
                obsexp = self.obs_exp_cpg(new_seq)
                llr = self.llr_score(new_seq, p_region, p_bg)
                if gc < min_gc or obsexp < min_obsexp or llr < min_llr:
                    continue
                score = (gc / max(orig_gc, 1e-10)) + (obsexp / max(orig_obsexp, 1e-10)) + (llr / max(orig_llr, 1e-10))
                if score > best_score:
                    best_score = score
                    best_island = {
                        "seq": new_seq,
                        "len": len(new_seq),
                        "start": new_start,
                        "end": new_end,
                        "gc": gc,
                        "obsexp": obsexp,
                        "llr": llr
                    }
                    improved = True
            if not improved:
                step //= 2
            else:
                orig_start, orig_end = best_island["start"], best_island["end"]
                orig_gc, orig_obsexp, orig_llr = best_island["gc"], best_island["obsexp"], best_island["llr"]
                best_score = (orig_gc / max(orig_gc, 1e-10)) + (orig_obsexp / max(orig_obsexp, 1e-10)) + (orig_llr / max(orig_llr, 1e-10))
        return best_island
    
    def run(self, genome, L, s, w, p, min_gc, min_obsexp, min_length, max_length, t, merge_gap):
        p_bg = self.build_1st_markov(genome)
        min_llr = self.get_threshold(genome, p_bg, L, s, t)
        cpg_islands = self.find_islands(genome, p_bg, L, s, w, p, min_gc, min_obsexp, min_llr, min_length, max_length)
        merged_cpg_islands = self.filter_overlap(cpg_islands, merge_gap, max_length, p_bg)
        refined_cpg_islands = []
        for island in merged_cpg_islands:
            refined_cpg_islands.append(self.refine_island(genome, island, p_bg, min_gc, min_obsexp, min_llr, min_length, max_length))
        return refined_cpg_islands
    
class Charts:
    def cpg_distribution(self, genome, cpg_islands):
        window_size = 500000
        genome_len = len(genome)
        density = [0] * (genome_len // window_size + 1)
        for island in cpg_islands:
            start_bin = island['start'] // window_size
            end_bin = island['end'] // window_size
            for b in range(start_bin, end_bin + 1):
                density[b] += 1
        plt.figure(figsize=(14, 5))
        plt.bar(range(len(density)), density, width=1.0, color="lightgreen", edgecolor="black")
        plt.xlabel(f"Position (per {window_size} bp)")
        plt.ylabel("CpG Islands")
        plt.title("CpG Island Density")
        plt.tight_layout()
        plt.show()

    def summary(self, cpg_islands):
        gcs = [island['gc'] for island in cpg_islands]
        oes = [island['obsexp'] for island in cpg_islands]
        lengths = [island['len'] for island in cpg_islands]
        scores = [island['llr'] for island in cpg_islands]
        results = []
        results.append("---Summary---\n")
        results.append(f"Mean length: {np.mean(lengths):.0f} | Min: {np.min(lengths)} | Max: {np.max(lengths)}\n")
        results.append(f"Mean GC%: {np.mean(gcs):.1f} | Min: {np.min(gcs):.1} | Max: {np.max(gcs):.1f}\n")
        results.append(f"Mean expected ratio: {np.mean(oes):.1f} | Min: {np.min(oes):.1} | Max: {np.max(oes):.1f}\n")
        results.append(f"Mean score: {np.mean(scores):.1f} | Min: {np.min(scores):.1f} | Max: {np.max(scores):.1f}\n")
        return results
    
    def display(self, cpg_islands, genome):
        self.cpg_distribution(genome, cpg_islands)
        summary = self.summary(cpg_islands)
        return summary