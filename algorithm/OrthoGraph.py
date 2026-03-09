import pandas as pd
from collections import defaultdict
from algorithm.SequenceSpace import SequenceSpace

class OrthoGraph:
    def __init__(self, ss: SequenceSpace, k: int):
        self.ss = ss
        self.k = k
        self.t = k+1 # forbidden k-mer length
        self.init_neighbors()

        print(f"Initialized graph for {self.ss.N} sequences with k={self.k}.")

    def _tmers(self, s:str) -> set[str]:
        """Return the set of all t-mers in the sequence s.
        Rk: t-mer is a substring of length t, where t = k+1 (the forbidden k-mer length)"""
        return {s[i:i+self.t] for i in range(len(s) - self.t + 1)}

    def init_neighbors(self):
        """
        Compute the neighbors for each sequence. Suffices to compute all
        t=k+1-mers for each sequence, then declare as neighbors any sequences 
        that share any t-mer.

        Fills in self.neighbors and self.self_conflicting in place.
        """
        s1 = self.ss.s1
        s2 = self.ss.s2
        N = self.ss.N

        self.neighbors = [set() for _ in range(N)]
        self.self_conflicting = set() # track sequences that could bind to themselves

        # precompute t-mers for every sequence
        t_mers_s1 = [self._tmers(s) for s in s1]
        t_mers_s2 = [self._tmers(s) for s in s2]

        # create a map from t-mer to the set of sequences which contain it
        tmer_to_s1 = defaultdict(set) # {t_mer: set of indices in s1 that contain t_mer}
        tmer_to_s2 = defaultdict(set) # ... for s2
        for s in range(N):
            for t_mer in t_mers_s1[s]:
                tmer_to_s1[t_mer].add(s)
            for t_mer in t_mers_s2[s]:
                tmer_to_s2[t_mer].add(s)

        # find neighbors based on shared t-mers
        # achieves the logical-OR of conflicts from (s1 x s1) and (s1 x s2)
        for s in range(N): # for string s in s1
            bad = set()
            for t_mer in t_mers_s1[s]:
                bad |= tmer_to_s1[t_mer] # AA conflicts
                bad |= tmer_to_s2[t_mer] # AB conflicts

            # > detect self-complementarity
            if t_mers_s1[s] & t_mers_s2[s]:
                self.self_conflicting.add(s)
                bad.add(s)

            self.neighbors[s] = bad


    def is_orthogonal(self, s1: str, s2: str):
        return s1 not in self.neighbors[s2]
    
    def get_orthogonal_sequences(self) -> tuple[list[str], list[str]]:
        """
        Greedy algorithm to find a large set of orthogonal sequences.
        
        Starts from the sequence with the fewest neighbors, adds it to the orthogonal
        set, removes its neighboring nodes and their incident edges, and repeats until no 
        sequences remain.
        """
        remaining = set(range(self.ss.N)) - self.self_conflicting
        orth1, orth2 = [], []

        while remaining:
            best = min(remaining, key=lambda s: len(self.neighbors[s] & remaining))
            orth1.append(self.ss.s1[best])
            orth2.append(self.ss.s2[best])

            # remove chosen vertex and all its conflicts
            remaining -= self.neighbors[best]
            remaining.discard(best)

        return orth1, orth2
    
    def export_orthogonal_set(self, filename: str):
        """Export the orthogonal set to a CSV file with columns "Sequence" and "Complement"."""
        orth1, orth2 = self.get_orthogonal_sequences()
        df = pd.DataFrame({"Sequence": orth1, "Complement": orth2})
        df.to_csv(filename, index=False)
        print(f"Exported {len(orth1)} orthogonal sequences to {filename}.")