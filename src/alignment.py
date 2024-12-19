# src/alignment.py

from Bio import Align
from Bio.Seq import Seq
from typing import Dict, Tuple, List, Optional
import pandas as pd
import numpy as np

class SequenceAligner:
    """
    A class to handle sequence alignment operations
    """
    def __init__(self):
        self.aligner = Align.PairwiseAligner()
        self._configure_aligner()
        self.alignments = {}

    def _configure_aligner(self):
        """Configure default alignment parameters"""
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -0.5
        self.aligner.extend_gap_score = -0.1

    def align_sequences(self, seq1: str, seq2: str, seq1_id: str = "seq1", seq2_id: str = "seq2") -> Dict:
        """
        Perform pairwise sequence alignment
        
        Args:
            seq1: First sequence string
            seq2: Second sequence string
            seq1_id: Identifier for first sequence
            seq2_id: Identifier for second sequence
            
        Returns:
            Dictionary containing alignment results
        """
        try:
            # Perform the alignment
            alignments = self.aligner.align(seq1, seq2)
            best_alignment = alignments[0]  # Get the best scoring alignment

            # Create alignment strings
            aligned_seq1, aligned_seq2 = best_alignment.aligned

            # Generate formatted alignment representation
            alignment_result = self._format_alignment(seq1, seq2, best_alignment)

            # Store the results
            self.alignments[f"{seq1_id}_{seq2_id}"] = {
                "score": best_alignment.score,
                "aligned_seq1": aligned_seq1,
                "aligned_seq2": aligned_seq2,
                "identity": self._calculate_identity(best_alignment),
                "alignment_result": alignment_result
            }

            return self.alignments[f"{seq1_id}_{seq2_id}"]

        except Exception as e:
            print(f"Error during alignment: {e}")
            return {}

    def _calculate_identity(self, alignment) -> float:
        """Calculate the percentage identity between aligned sequences"""
        matches = sum(a == b for a, b in zip(str(alignment[0]), str(alignment[1])))
        total = len(str(alignment[0]))
        return (matches / total) * 100 if total > 0 else 0

    def _format_alignment(self, seq1: str, seq2: str, alignment) -> str:
        """Format the alignment for display"""
        aligned_seq1, aligned_seq2 = alignment.aligned
        
        # Convert alignment coordinates to strings
        s1 = list(seq1)
        s2 = list(seq2)
        
        # Create the alignment string
        alignment_str = ""
        match_str = ""
        
        for idx in range(len(aligned_seq1[0])):
            if idx < len(s1) and idx < len(s2):
                if s1[idx] == s2[idx]:
                    match_str += "|"
                else:
                    match_str += " "
            else:
                match_str += " "
                
        return f"{seq1}\n{match_str}\n{seq2}"

    def get_alignment_summary(self) -> pd.DataFrame:
        """
        Get a summary of all alignments as a DataFrame
        
        Returns:
            DataFrame containing alignment statistics
        """
        data = []
        for pair_id, result in self.alignments.items():
            data.append({
                "pair_id": pair_id,
                "score": result["score"],
                "identity_percentage": result["identity"],
            })
        return pd.DataFrame(data)