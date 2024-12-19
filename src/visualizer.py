# src/visualizer.py

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq
from pathlib import Path

class SequenceVisualizer:
    """
    A class to handle visualization of sequence data and analysis results
    """
    # src/visualizer.py
    def __init__(self):
        # Set style using a built-in style
        plt.style.use('default')  # Using default style instead of seaborn
        # Set color palette manually
        colors = ['#2ecc71', '#e74c3c', '#3498db', '#f1c40f', '#9b59b6']
        sns.set_palette(colors)
        
    def plot_sequence_composition(self, sequences: Dict[str, Dict], 
                                save_path: Optional[Path] = None) -> None:
        """
        Create a bar plot of base composition for each sequence
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        data = []
        for seq_id, stats in sequences.items():
            for base, count in stats['base_counts'].items():
                data.append({
                    'Sequence': seq_id,
                    'Base': base,
                    'Count': count
                })
        
        df = pd.DataFrame(data)
        sns.barplot(data=df, x='Sequence', y='Count', hue='Base', ax=ax)
        
        plt.title('Sequence Base Composition')
        plt.xlabel('Sequence ID')
        plt.ylabel('Base Count')
        
        if save_path:
            plt.savefig(save_path)
        plt.close()

    def plot_gc_content(self, sequences: Dict[str, Dict], 
                       window_size: int = 10,
                       save_path: Optional[Path] = None) -> None:
        """
        Create a line plot of GC content along the sequence
        """
        fig, ax = plt.subplots(figsize=(12, 6))
        
        for seq_id, stats in sequences.items():
            sequence = ''.join([base for base, count in stats['base_counts'].items() 
                              for _ in range(count)])
            
            # Calculate GC content in windows
            gc_content = []
            positions = []
            
            for i in range(0, len(sequence) - window_size + 1):
                window = sequence[i:i + window_size]
                gc = (window.count('G') + window.count('C')) / window_size * 100
                gc_content.append(gc)
                positions.append(i + window_size // 2)
            
            ax.plot(positions, gc_content, label=seq_id)
        
        plt.title(f'GC Content Distribution (Window Size: {window_size})')
        plt.xlabel('Position')
        plt.ylabel('GC Content (%)')
        plt.legend()
        
        if save_path:
            plt.savefig(save_path)
        plt.close()

    def visualize_alignment(self, alignment_result: Dict,
                          save_path: Optional[Path] = None) -> None:
        """
        Create a visualization of sequence alignment
        """
        fig, ax = plt.subplots(figsize=(15, 4))
        
        alignment_str = alignment_result['alignment_result'].split('\n')
        
        # Plot sequences as colored blocks
        y_positions = [2, 0]  # Position for each sequence
        
        for idx, seq in enumerate([alignment_str[0], alignment_str[2]]):
            x_position = 0
            for base in seq:
                color = self._get_base_color(base)
                ax.add_patch(plt.Rectangle((x_position, y_positions[idx]), 
                                         1, 0.8, 
                                         facecolor=color))
                x_position += 1
        
        # Add sequence labels
        ax.text(-1, 2, 'Seq 1', va='center')
        ax.text(-1, 0, 'Seq 2', va='center')
        
        # Customize plot
        ax.set_ylim(-0.5, 3)
        ax.set_xlim(-1.5, len(alignment_str[0]) + 0.5)
        ax.set_title(f'Sequence Alignment (Identity: {alignment_result["identity"]:.1f}%)')
        ax.axis('off')
        
        if save_path:
            plt.savefig(save_path)
        plt.close()

    def plot_pattern_distribution(self, pattern_summary: pd.DataFrame,
                                save_path: Optional[Path] = None) -> None:
        """
        Create a visualization of pattern distributions along the sequence
        """
        if pattern_summary.empty:
            return
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Plot patterns as colored blocks
        patterns = pattern_summary['pattern'].unique()
        colors = sns.color_palette("husl", len(patterns))
        pattern_colors = dict(zip(patterns, colors))
        
        for idx, row in pattern_summary.iterrows():
            ax.add_patch(plt.Rectangle((row['start'], 0),
                                     row['length'],
                                     0.8,
                                     facecolor=pattern_colors[row['pattern']],
                                     alpha=0.6,
                                     label=row['pattern']))
        
        # Remove duplicate labels
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(),
                  title="Patterns", bbox_to_anchor=(1.05, 1))
        
        plt.title('Pattern Distribution Along Sequence')
        plt.xlabel('Position')
        plt.ylim(-0.2, 1)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
        plt.close()

    @staticmethod
    def _get_base_color(base: str) -> str:
        """Return color for DNA base"""
        colors = {
            'A': '#2ecc71',  # Green
            'T': '#e74c3c',  # Red
            'G': '#f1c40f',  # Yellow
            'C': '#3498db',  # Blue
            '-': '#95a5a6'   # Gray
        }
        return colors.get(base, '#ffffff')  # White for unknown bases