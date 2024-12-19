# main.py

from src.sequence_processor import SequenceProcessor
from src.alignment import SequenceAligner
from src.pattern_finder import PatternFinder
from pathlib import Path

def main():
    # Initialize the processors
    processor = SequenceProcessor()
    aligner = SequenceAligner()
    pattern_finder = PatternFinder()
    
    # Load sequences
    data_path = Path("data/example_sequences/test.fasta")
    
    if processor.load_sequence(data_path):
        print("Sequences loaded successfully!")
        
        # Get statistics
        stats = processor.get_sequence_stats()
        print("\nSequence Statistics:")
        print(stats)
        
        # Get DataFrame representation
        df = processor.get_sequences_as_dataframe()
        print("\nSequence DataFrame:")
        print(df)
        
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
            
            print("\nAlignment Result:")
            print(f"Score: {alignment_result['score']}")
            print(f"Identity: {alignment_result['identity']:.2f}%")
            print("\nAlignment:")
            print(alignment_result['alignment_result'])

        # Pattern analysis
        print("\nAnalyzing sequence patterns...")
        first_seq = str(sequences[0].seq)
        
        # Find repeats
        repeats = pattern_finder.find_repeats(first_seq)
        print("\nRepetitive Elements:")
        for pattern, info in repeats.items():
            print(f"Pattern '{pattern}' found {info['count']} times")
            
        # Find regulatory elements
        regulatory = pattern_finder.find_regulatory_elements(first_seq)
        if regulatory:
            print("\nRegulatory Elements Found:")
            for element, locations in regulatory.items():
                print(f"{element}: {len(locations)} occurrences")
        
        # Get pattern summary
        pattern_summary = pattern_finder.get_pattern_summary(first_seq)
        if not pattern_summary.empty:
            print("\nPattern Summary:")
            print(pattern_summary)

if __name__ == "__main__":
    main()