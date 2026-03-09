from itertools import product
import pandas as pd

class SequenceSpace:
    def __init__(self, seq_length: int=8):
        """Initialize our binding library"""
        self.BASES = "ATGC"
        self.COMPLEMENT = str.maketrans("ATGC", "TACG")
        self.SEQ_LENGTH = seq_length
        self.init_temp_space()
        self.init_space(seq_length)
        self.N = len(self.s1)

    def init_temp_space(self):
        """Initialize the melting temperature map from the excel file
        with all melting temps for sequences of length SEQ_LENGTH
        
        Final output:
            -> self.tm[sequence] = melting temp
        """
        # read in the excel file with melting temps
        self.tm_df = pd.read_excel(f"data/melting_temps/{self.SEQ_LENGTH}bp.xlsx")
        self.tm = {}
        for _, row in self.tm_df.iterrows():
            seq = row['Sequence']
            comp = row['Complement']
            temp = row['Melting Temp']
            self.tm[seq] = temp
            self.tm[comp] = temp

    def init_space(self, seq_length: int):
        """Initialize the binding space. Creates s1 and s2, which are lists of the
        color sequences and their complements respectively."""
        self.seq_length = seq_length
        sequences = [''.join(p) for p in product(self.BASES, repeat=seq_length)]
        self.s1, self.s2 = self.sort_sequences(sequences)

    def complement(self, seq: str, reverse=True) -> str:
        """Return the complement of a sequence"""
        # we reverse cus the excel reverses it
        comp = seq.translate(self.COMPLEMENT)
        return comp if not reverse else comp[::-1]

    def sort_sequences(self, sequences: list[str]) -> tuple[list[str], list[str]]:
        """Sort sequences into color and complements"""
        seen = set()
        s1, s2 = [], []
        for s in sequences:
            if s in seen:
                continue
            c = self.complement(s, reverse=True)
            # choose a canonical representative
            color, comp = (s, c) if s < c else (c, s)

            # add to our lists and seen set
            s1.append(color)
            s2.append(comp)
            seen.add(s)
            seen.add(c)

        return s1, s2
    
    def filter_sequences(self, gc_range: tuple[float, float]=(0.3, 0.7), 
                         temp_range: tuple[float, float]=(5.0, 50.0), 
                         exclude_ggg=True):
        """Filter the sequence space based on the provided criteria.
        
        Args:
            gc_range: Min/max GC proportion per sequence
            temp_range: Min/max melting temperature per sequence
            exclude_ggg: Whether to exclude sequences containing "GGG" or "CCC"
        """
        s1, s2 = [], []
        for s, c in zip(self.s1, self.s2):
            
            good_ggg = self._excludes_ggg(s, exclude_ggg)
            good_gc = self._filter_gc(s, gc_range)
            good_temp = self._filter_temp(s, temp_range)
            good_palindrome = not self._is_palindrome(s)

            if good_gc and good_temp and good_ggg and good_palindrome:
                s1.append(s)
                s2.append(c)
        
        print(f"Filtered from {self.N} to {len(s1)} sequences.")
        self.N = len(s1)
        self.s1, self.s2 = s1, s2

    def _filter_gc(self, s: str, gc_range: tuple[float, float]=(0.3, 0.7)):
        """Returns True if the sequence 's' passes the gc_range filter check."""
        gc_count = s.count("G") + s.count("C")
        gc_frac = gc_count / self.SEQ_LENGTH

        if gc_frac >= gc_range[0] and gc_frac <= gc_range[1]:
            return True
        return False

    def _excludes_ggg(self, s: str, exclude_ggg=True):
        """Returns True if the sequence 's' passes the excludes "GGG" check."""
        if not exclude_ggg:
            return True
        
        if "GGG" in s or "CCC" in s:
            return False 
        return True

    def _filter_temp(self, s: str, temp_range: tuple[float, float]=(5, 50)):
        """Returns True if the sequence 's' passes the melting temperature range check."""
        temp = self.tm.get(s, None)
        if not temp:
            return False
        if temp_range[0] <= temp <= temp_range[1]:
            return True
        return False
    
    def _is_palindrome(self, s: str) -> bool:
        """Check if a sequence is a palindrome"""
        return s == self.complement(s)