# src/sequence_processor.py

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction  # Changed from GC to gc_fraction
import pandas as pd
from pathlib import Path
from typing import Dict, List, Union, Optional

class SequenceProcessor:
    """
    A class to handle basic DNA sequence processing operations.
    """
    def __init__(self):
        self.sequences = {}  # Store sequences with their IDs
        self.stats = {}      # Store sequence statistics

    def load_sequence(self, file_path: Union[str, Path], file_format: str = "fasta") -> bool:
        """
        Load sequences from a file (FASTA/FASTQ)
        
        Args:
            file_path: Path to the sequence file
            file_format: Format of the file ('fasta' or 'fastq')
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            for record in SeqIO.parse(str(file_path), file_format):
                self.sequences[record.id] = record
                # Calculate basic statistics for the sequence
                self._calculate_sequence_stats(record)
            return True
        except Exception as e:
            print(f"Error loading sequence file: {e}")
            return False

    def _calculate_sequence_stats(self, record) -> None:
        """
        Calculate basic statistics for a sequence
        """
        sequence = str(record.seq)
        self.stats[record.id] = {
            "length": len(sequence),
            "gc_content": gc_fraction(sequence) * 100,  # Convert to percentage
            "base_counts": {
                "A": sequence.count('A'),
                "T": sequence.count('T'),
                "G": sequence.count('G'),
                "C": sequence.count('C')
            }
        }

    def get_sequence_stats(self, sequence_id: Optional[str] = None) -> Dict:
        """
        Get statistics for a specific sequence or all sequences
        
        Args:
            sequence_id: ID of the sequence (optional)
            
        Returns:
            Dict containing sequence statistics
        """
        if sequence_id:
            return self.stats.get(sequence_id, {})
        return self.stats

    def get_sequences_as_dataframe(self) -> pd.DataFrame:
        """
        Convert sequence statistics to a pandas DataFrame
        
        Returns:
            DataFrame containing sequence statistics
        """
        data = []
        for seq_id, stats in self.stats.items():
            row = {
                "sequence_id": seq_id,
                "length": stats["length"],
                "gc_content": stats["gc_content"],
                **{f"count_{base}": count for base, count in stats["base_counts"].items()}
            }
            data.append(row)
        return pd.DataFrame(data)