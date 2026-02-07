#!/usr/bin/env python3
"""
Homology Analysis Script
========================
This script performs homology-based analysis using NCBI BLAST to:
- Find similar proteins in public databases
- Identify evolutionary relationships
- Provide functional clues based on homologs
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import os
import sys

# File paths
INPUT_FILE = "E:/Functional_Sequence_Characterization/Data/input/input.fasta"
OUTPUT_FILE = "E:/Functional_Sequence_Characterization/Data/Results/blast_results.txt"
XML_OUTPUT = "E:/Functional_Sequence_Characterization/Data/Results/blast_results.xml"

# BLAST parameters
BLAST_DATABASE = "nr"  # Non-redundant protein database
BLAST_PROGRAM = "blastp"  # Protein-protein BLAST
EXPECT_THRESHOLD = 0.001  # E-value cutoff
MAX_HITS = 10  # Maximum number of hits to retrieve

# Ensure results directory exists
os.makedirs("E:/Functional_Sequence_Characterization/Data/Results", exist_ok=True)


def run_blast_search(sequence):
    """
    Perform BLAST search against NCBI database.
    
    """
    print("Submitting BLAST query to NCBI...")
    print("This may take 2-5 minutes depending on server load...")
    print("Please be patient...")
    print()
    
    try:
        # Submit BLAST query
        result_handle = NCBIWWW.qblast(
            program=BLAST_PROGRAM,
            database=BLAST_DATABASE,
            sequence=sequence,
            expect=EXPECT_THRESHOLD,
            hitlist_size=MAX_HITS,
            format_type="XML"
        )
        
        # Save XML results
        with open(XML_OUTPUT, 'w') as save_file:
            blast_results = result_handle.read()
            save_file.write(blast_results)
        
        result_handle.close()
        
        print("✓ BLAST search completed successfully!")
        return blast_results
        
    except Exception as e:
        print(f"ERROR during BLAST search: {e}")
        print("This could be due to:")
        print("  - No internet connection")
        print("  - NCBI server is down")
        print("  - Network firewall blocking the request")
        sys.exit(1)


def parse_blast_results(xml_file):
    """
    Parse BLAST XML results and extract meaningful information.
    """
    print("Parsing BLAST results...")
    
    hits_data = []
    
    try:
        with open(xml_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        hit_info = {
                            'title': alignment.title,
                            'length': alignment.length,
                            'e_value': hsp.expect,
                            'score': hsp.score,
                            'identities': hsp.identities,
                            'positives': hsp.positives,
                            'gaps': hsp.gaps,
                            'query_length': blast_record.query_length,
                            'query_seq': hsp.query,
                            'subject_seq': hsp.sbjct,
                            'match_seq': hsp.match,

                        }
                        
                        # Calculate percentage identity
                        hit_info['percent_identity'] = (hsp.identities / len(hsp.query)) * 100
                        hit_info['percent_positives'] = (hsp.positives / len(hsp.query)) * 100
                        
                        hits_data.append(hit_info)
        
        print(f"✓ Found {len(hits_data)} significant hits")
        return hits_data
        
    except Exception as e:
        print(f"ERROR parsing BLAST results: {e}")
        sys.exit(1)


def extract_species_from_title(title):
    """
    Extract species information from BLAST hit title.  """
    if '[' in title and ']' in title:
        start = title.rfind('[')
        end = title.rfind(']')
        return title[start+1:end]
    return "Unknown"


def analyze_homologs(hits_data):
    """
    Analyze homologous sequences to extract functional insights.
    """
    if not hits_data:
        return {
            'species_distribution': {},
            'functional_keywords': {},
            'avg_identity': 0,
            'conservation_level': 'No hits found'
        }
    
    # Species distribution
    species_count = {}
    for hit in hits_data:
        species = extract_species_from_title(hit['title'])
        species_count[species] = species_count.get(species, 0) + 1
    
    
    
    # Average identity
    avg_identity = sum(hit['percent_identity'] for hit in hits_data) / len(hits_data)
    
    # Conservation level
    if avg_identity > 80:
        conservation = "Highly conserved protein"
    elif avg_identity > 50:
        conservation = "Moderately conserved protein"
    elif avg_identity > 30:
        conservation = "Poorly conserved protein"
    else:
        conservation = "Highly divergent protein"
    
    return {
        'species_distribution': species_count,
        'avg_identity': avg_identity,
        'conservation_level': conservation
    }

# MAIN ANALYSIS FUNCTION


def perform_homology_analysis():
    """
    Main function to perform homology analysis.
    """
    print("=" * 80)
    print("HOMOLOGY ANALYSIS USING BLAST")
    print("=" * 80)
    print()
    
    # Check if input file exists
    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: Input file not found: {INPUT_FILE}")
        sys.exit(1)
    
    # Read the sequence
    print(f"Reading sequence from: {INPUT_FILE}")
    try:
        record = SeqIO.read(INPUT_FILE, "fasta")
        sequence = str(record.seq)
        description = record.description
    except Exception as e:
        print(f"ERROR reading FASTA file: {e}")
        sys.exit(1)
    
    print(f"Sequence ID: {record.id}")
    print(f"Length: {len(sequence)} amino acids")
    print()
    
    # Run BLAST search
    print("Step 1: Running BLAST search...")
    print(f"Database: {BLAST_DATABASE}")
    print(f"E-value threshold: {EXPECT_THRESHOLD}")
    print(f"Maximum hits: {MAX_HITS}")
    print()
    
    blast_results = run_blast_search(sequence)
    print()
    
    # Parse results
    print("Step 2: Parsing BLAST results...")
    hits_data = parse_blast_results(XML_OUTPUT)
    print()
    
    # Analyze homologs
    print("Step 3: Analyzing homologous sequences...")
    analysis = analyze_homologs(hits_data)
    print("✓ Analysis complete")
    print()
    
   
    # GENERATE OUTPUT REPORT
    
    
    print(f"Step 4: Writing results to: {OUTPUT_FILE}")
    
    with open(OUTPUT_FILE, 'w') as f:
        # Header
        f.write("=" * 80 + "\n")
        f.write("BLAST HOMOLOGY ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        # Query Information
        f.write("QUERY SEQUENCE INFORMATION\n")
        f.write("-" * 80 + "\n")
        f.write(f"Sequence ID:        {record.id}\n")
        f.write(f"Description:        {description}\n")
        f.write(f"Length:             {len(sequence)} amino acids\n")
        f.write(f"BLAST Program:      {BLAST_PROGRAM}\n")
        f.write(f"Database:           {BLAST_DATABASE}\n")
        f.write(f"E-value threshold:  {EXPECT_THRESHOLD}\n")
        f.write("\n")
        
        # Summary Statistics
        f.write("SUMMARY STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total hits found:   {len(hits_data)}\n")
        f.write(f"Average identity:   {analysis['avg_identity']:.2f}%\n")
        f.write(f"Conservation:       {analysis['conservation_level']}\n")
        f.write("\n")
        
        # Top BLAST Hits
        f.write("TOP BLAST HITS\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Rank':<6}{'E-value':<12}{'Identity':<12}{'Description'}\n")
        f.write("-" * 80 + "\n")
        
        for i, hit in enumerate(hits_data[:10], 1):  # Show top 10
            title_short = hit['title'][:60] + "..." if len(hit['title']) > 60 else hit['title']
            f.write(f"{i:<6}{hit['e_value']:<12.2e}{hit['percent_identity']:<11.1f}% {title_short}\n")
        
        f.write("\n")
        
        # Species Distribution
        f.write("SPECIES DISTRIBUTION (Top 10)\n")
        f.write("-" * 80 + "\n")
        sorted_species = sorted(analysis['species_distribution'].items(), 
                               key=lambda x: x[1], reverse=True)
        
        for species, count in sorted_species[:10]:
            f.write(f"{species:<50} {count:>3} hits\n")
        f.write("\n")
        
        # Detailed Alignments
        f.write("DETAILED ALIGNMENT INFORMATION (Top 5 Hits)\n")
        f.write("=" * 80 + "\n\n")
        
        for i, hit in enumerate(hits_data[:5], 1):
            f.write(f"Hit #{i}\n")
            f.write("-" * 80 + "\n")
            f.write(f"Description: {hit['title']}\n")
            f.write(f"E-value:     {hit['e_value']:.2e}\n")
            f.write(f"Score:       {hit['score']:.1f}\n")
            f.write(f"Identities:  {hit['identities']}/{len(hit['query_seq'])} " +
                   f"({hit['percent_identity']:.1f}%)\n")
            f.write(f"Positives:   {hit['positives']}/{len(hit['query_seq'])} " +
                   f"({hit['percent_positives']:.1f}%)\n")
            f.write(f"Gaps:        {hit['gaps']}/{len(hit['query_seq'])}\n")
           
            
            # Format alignment in blocks of 60
            query = hit['query_seq']
            match = hit['match_seq']
            subject = hit['subject_seq']
            
            for j in range(0, len(query), 60):
                f.write(f"\nQuery  {j+1:>4}  {query[j:j+60]}\n")
                f.write(f"             {match[j:j+60]}\n")
                f.write(f"Sbjct  {j+1:>4}  {subject[j:j+60]}\n")
            
            f.write("\n" + "=" * 80 + "\n\n")
        
        # Functional Prediction
        f.write("FUNCTIONAL PREDICTION BASED ON HOMOLOGY\n")
        f.write("-" * 80 + "\n")
        
        if hits_data:
            top_hit = hits_data[0]
            f.write("Based on the top BLAST hit:\n")
            f.write(f"  Description: {top_hit['title']}\n")
            f.write(f"  Identity:    {top_hit['percent_identity']:.1f}%\n")
            f.write(f"  E-value:     {top_hit['e_value']:.2e}\n\n")
            
            if top_hit['percent_identity'] > 90:
                confidence = "Very High"
                note = "Likely same or very similar function"
            elif top_hit['percent_identity'] > 70:
                confidence = "High"
                note = "Likely similar function with possible variations"
            elif top_hit['percent_identity'] > 50:
                confidence = "Moderate"
                note = "Related function, but may have diverged"
            elif top_hit['percent_identity'] > 30:
                confidence = "Low"
                note = "Distantly related, function may differ significantly"
            else:
                confidence = "Very Low"
                note = "Weak homology, function uncertain"
            
            f.write(f"  Confidence:  {confidence}\n")
            f.write(f"  Note:        {note}\n")
        else:
            f.write("No significant homologs found.\n")
            f.write("This could indicate:\n")
            f.write("  - Novel protein with no characterized homologs\n")
            f.write("  - Highly divergent sequence\n")
            f.write("  - Sequence errors\n")
        
        f.write("\n")
    
        
        # Footer
        f.write("=" * 80 + "\n")
        f.write("END OF BLAST REPORT\n")
        f.write("=" * 80 + "\n")
    
    print("✓ Results saved successfully!")
    print()
    print("=" * 80)
    print("HOMOLOGY ANALYSIS COMPLETE")
    print("=" * 80)
    print()
    
    if hits_data:
        print(f"Key Findings:")
        print(f"  - {len(hits_data)} significant homologs found")
        print(f"  - Average identity: {analysis['avg_identity']:.1f}%")
        print(f"  - {analysis['conservation_level']}")
        print(f"  - Top hit: {hits_data[0]['title'][:60]}...")
    else:
        print("No significant homologs found")
    
    print()
    print(f"Full report saved to: {OUTPUT_FILE}")
    print(f"XML results saved to: {XML_OUTPUT}")
    print()


# ====================================================================================
# SCRIPT EXECUTION
# ====================================================================================

if __name__ == "__main__":
    try:
        perform_homology_analysis()
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)