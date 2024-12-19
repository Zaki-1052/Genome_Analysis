# main.py

from src.sequence_processor import SequenceProcessor
from src.alignment import SequenceAligner
from src.pattern_finder import PatternFinder
from src.visualizer import SequenceVisualizer
from pathlib import Path

def main():
    # Initialize all processors
    processor = SequenceProcessor()
    aligner = SequenceAligner()
    pattern_finder = PatternFinder()
    visualizer = SequenceVisualizer()
    
    # Create output directory for plots
    output_dir = Path("data/output")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load sequences
    data_path = Path("data/example_sequences/test.fasta")
    
    if processor.load_sequence(data_path):
        print("Sequences loaded successfully!")
        
        # Get statistics and create visualizations
        stats = processor.get_sequence_stats()
        print("\nSequence Statistics:")
        print(stats)
        
        # Plot sequence composition
        visualizer.plot_sequence_composition(
            stats,
            save_path=output_dir / "sequence_composition.png"
        )
        print("\nGenerated sequence composition plot")
        
        # Plot GC content
        visualizer.plot_gc_content(
            stats,
            save_path=output_dir / "gc_content.png"
        )
        print("Generated GC content plot")
        
        # Perform alignment
        sequences = list(processor.sequences.values())
        if len(sequences) >= 2:
            print("\nPerforming sequence alignment...")
            alignment_result = aligner.align_sequences(
                str(sequences[0].seq),
                str(sequences[1].seq),
                sequences[0].id,
                sequences[1].id
            )
            
            # Visualize alignment
            visualizer.visualize_alignment(
                alignment_result,
                save_path=output_dir / "alignment.png"
            )
            print("Generated alignment visualization")
        
        # Pattern analysis and visualization
        print("\nAnalyzing sequence patterns...")
        first_seq = str(sequences[0].seq)
        pattern_summary = pattern_finder.get_pattern_summary(first_seq)
        
        if not pattern_summary.empty:
            visualizer.plot_pattern_distribution(
                pattern_summary,
                save_path=output_dir / "pattern_distribution.png"
            )
            print("Generated pattern distribution plot")

        print("\nAll visualizations have been saved to the 'data/output' directory")

if __name__ == "__main__":
    main()