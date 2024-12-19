# src/pattern_finder.py

from Bio import motifs
from Bio.Seq import Seq
import re
from collections import Counter
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np

class PatternFinder:
    """
    A class to handle DNA sequence pattern analysis and motif finding
    """
    def __init__(self):
        self.patterns = {}
        self.common_motifs = {
            'TATA_box': 'TATAAA',
            'CAP_signal': 'CATTTT',
            'Kozak_sequence': 'ACCATGG',
            'GC_box': 'GGGCGG'
        }

    def find_pattern(self, sequence: str, pattern: str, min_length: int = 4) -> List[Dict]:
        """
        Find all occurrences of a specific pattern in the sequence
        
        Args:
            sequence: Input DNA sequence
            pattern: Pattern to search for
            min_length: Minimum length of pattern to consider
            
        Returns:
            List of dictionaries containing pattern information
        """
        if len(pattern) < min_length:
            return []

        results = []
        for match in re.finditer(pattern, sequence):
            results.append({
                'pattern': pattern,
                'start': match.start(),
                'end': match.end(),
                'length': len(pattern)
            })
        return results

    def find_repeats(self, sequence: str, min_length: int = 3, max_length: int = 8) -> Dict[str, List[Dict]]:
        """
        Find repetitive elements in the sequence
        
        Args:
            sequence: Input DNA sequence
            min_length: Minimum length of repeat
            max_length: Maximum length of repeat
            
        Returns:
            Dictionary of repeat patterns and their locations
        """
        repeats = {}
        
        for length in range(min_length, max_length + 1):
            kmers = [sequence[i:i+length] for i in range(len(sequence) - length + 1)]
            kmer_counts = Counter(kmers)
            
            for kmer, count in kmer_counts.items():
                if count > 1:  # Only consider repeated patterns
                    repeat_info = self.find_pattern(sequence, kmer)
                    if repeat_info:
                        repeats[kmer] = {
                            'count': count,
                            'locations': repeat_info
                        }
        
        return repeats

    def find_regulatory_elements(self, sequence: str) -> Dict[str, List[Dict]]:
        """
        Find known regulatory elements in the sequence
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            Dictionary of regulatory elements and their locations
        """
        regulatory_elements = {}
        
        for motif_name, motif_seq in self.common_motifs.items():
            locations = self.find_pattern(sequence, motif_seq)
            if locations:
                regulatory_elements[motif_name] = locations
                
        return regulatory_elements

    def find_custom_motif(self, sequences: List[str], motif_length: int = 6) -> Dict:
        """
        Find potential common motifs in a list of sequences
        
        Args:
            sequences: List of DNA sequences
            motif_length: Length of motif to look for
            
        Returns:
            Dictionary containing motif information
        """
        if not sequences:
            return {}

        # Convert sequences to Seq objects
        instances = [Seq(seq) for seq in sequences]
        
        # Create a motif object
        m = motifs.create(instances)
        
        return {
            'consensus': str(m.consensus),
            'pwm': m.counts.normalize(),
            'sequences': sequences,
            'motif_length': motif_length
        }

    def get_pattern_summary(self, sequence: str) -> pd.DataFrame:
        """
        Generate a summary of all patterns found in the sequence
        
        Args:
            sequence: Input DNA sequence
            
        Returns:
            DataFrame containing pattern statistics
        """
        # Find all types of patterns
        repeats = self.find_repeats(sequence)
        regulatory = self.find_regulatory_elements(sequence)
        
        # Combine results into a DataFrame
        data = []
        
        # Add repeat patterns
        for pattern, info in repeats.items():
            for loc in info['locations']:
                data.append({
                    'pattern_type': 'repeat',
                    'pattern': pattern,
                    'start': loc['start'],
                    'end': loc['end'],
                    'length': loc['length'],
                    'count': info['count']
                })
        
        # Add regulatory elements
        for element, locations in regulatory.items():
            for loc in locations:
                data.append({
                    'pattern_type': 'regulatory',
                    'pattern': element,
                    'start': loc['start'],
                    'end': loc['end'],
                    'length': loc['length'],
                    'count': 1
                })
        
        return pd.DataFrame(data)