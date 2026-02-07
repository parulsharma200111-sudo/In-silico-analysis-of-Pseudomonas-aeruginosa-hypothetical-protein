"""
Sequence Quality Control Script
================================
This script performs comprehensive quality control on protein sequences including:
- Sequence validation
- Physicochemical properties calculation
- Amino acid composition analysis
- Basic sequence statistics
"""

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import molecular_weight
import os
import sys

# File paths
INPUT_FILE = "E:/Functional_Sequence_Characterization/Data/input/input.fasta"
OUTPUT_FILE = "E:/Functional_Sequence_Characterization/Data/Results/qc_summary.txt"

# Functions

def validate_protein_sequence(sequence):
    """ Validate that the sequence contains only valid amino acid charaters"""
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    sequence_chars = set(sequence.upper())
    invalid_chars = sequence_chars - valid_amino_acids

    return len(invalid_chars) == 0, invalid_chars

def calculate_basic_properties(sequence):
    """ Caluculate basic physiochemical properties of the protein sequence"""

    analysed_seq = ProteinAnalysis(sequence)
    properties = {
        'length': len(sequence),
        'molecular_weight': analysed_seq.molecular_weight(),
        'aromaticity': analysed_seq.aromaticity(),
        'instability_index': analysed_seq.instability_index(),
        'isoelectric_point': analysed_seq.isoelectric_point(),
        'gravy': analysed_seq.gravy(),  # Grand average of hydropathicity
        'secondary_structure': analysed_seq.secondary_structure_fraction(),
        'molar_extinction_coefficient': analysed_seq.molar_extinction_coefficient(),
        'flexibility' : analysed_seq.flexibility()

    }

    return properties


def calculate_amino_acid_composition(sequence):
    """Calculate the percentage composition of each amino acid."""

    analysed_seq = ProteinAnalysis(sequence)
    # Some versions expose `amino_acids_percent` as a dict attribute,
    # others as a callable method. Handle both.
    aa_percent = analysed_seq.amino_acids_percent
    if callable(aa_percent):
        return aa_percent()
    return aa_percent

def interpret_properties(properties):
    """ Provide biological interpretation of calculated properties"""

    interpretations = {}
    
    # molecular weight interpretation 
    if properties['molecular_weight'] < 10000:
        interpretations['size'] = "Small_peptide"
    elif properties['molecular_weight'] < 50000:
        interpretations['size'] = "Medium-size protein"
    else:
        interpretations['size'] = "Large protein"

    # Instability index interpretation
    if properties['gravy'] > 0:
        interpretations['hydropathy'] = "Hydrophobic (positive GRAVY)"
    else:
        interpretations['hydropathy'] = "Hydrophilic (negative GRAVY)"

    # Stability interpretation using instability index (ProtParam rule of thumb)
    # instability_index < 40 -> stable; >= 40 -> unstable
    instability = properties.get('instability_index')
    try:
        if instability is not None and float(instability) < 40.0:
            interpretations['stability'] = "Stable (instability index < 40)"
        else:
            interpretations['stability'] = "Unstable (instability index >= 40)"
    except Exception:
        interpretations['stability'] = "Unknown"

    # pI interpretation 

    if properties['isoelectric_point'] < 7:
        interpretations['charge'] = "Acidic protein (pI < 7)"
    elif properties['isoelectric_point'] > 7:
        interpretations['charge'] = "Basic protein (pI > 7)"
    else:
        interpretations['charge'] = "Neutral protein (pI = 7)"

    return interpretations

# MAIN ANALYSIS FUNCTION 

def perform_qc_analysis():
    """ Main function to perform quality control analysis on the input sequence """
    print("*" * 80)
    print("PROTEIN SEQUENCE QULAITY CONTROL ANALYSIS")
    print("*" * 80)

    # Read the sequence

    record  = SeqIO.read(INPUT_FILE, "fasta")
    sequence = str(record.seq)
    description = record.description

    print(f"Sequence ID: {record.id}")
    print(f"Description: {description}")
    print()

    # Validate sequence
    print("Step 1: Validating sequence...")
    is_valid, invalid_chars = validate_protein_sequence(sequence)
    
    if not is_valid:
        print(f"WARNING: Invalid amino acid characters found: {invalid_chars}")
        print("Proceeding with analysis, but results may be unreliable.")
    else:
        print("✓ Sequence validated successfully - all characters are valid amino acids")
    print()
    
    # Calculate properties
    print("Step 2: Calculating physicochemical properties...")

    # Remove any non-standard amino acids (e.g., 'X') for analysis to avoid
    # ProteinAnalysis errors; keep uppercase for consistency
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    analysis_sequence = ''.join([ch for ch in sequence.upper() if ch in valid_amino_acids])
    if not analysis_sequence:
        print("ERROR: No valid amino acids found after cleaning sequence. Aborting.")
        sys.exit(1)

    if not is_valid:
        print(f"Note: invalid characters removed; using cleaned sequence length {len(analysis_sequence)} for analysis.")

    properties = calculate_basic_properties(analysis_sequence)
    print()
    
    # Calculate amino acid composition
    print("Step 3: Analyzing amino acid composition...")
    aa_composition = calculate_amino_acid_composition(analysis_sequence)
    print()
  
    # Get interpretations
    interpretations = interpret_properties(properties)

    # GENERATE OUTPUT REPORT

    print(f"Step 4: Writing results to: {OUTPUT_FILE}")
    
    with open (OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write("PROTEIN SEQUENCE QUALITY CONTROL REPORT\n\n")
        f.write("SEQUENCE INFORMATION\n\n")
        f.write(f"Sequence ID:          {record.id}\n")
        f.write(f"Description:          {description}\n")
        f.write(f"Sequence Length:      {properties['length']} amino acids\n")
        f.write(f"Validation Status:    {'PASSED' if is_valid else 'FAILED'}\n")
        if not is_valid:
            f.write(f"Invalid Characters:   {invalid_chars}\n")
        f.write("\n")

        # PHYSIOCHEMICAL PROPERTIES
        f.write("PHYSIOCHEMICAL PROPERTIES\n\n")
        f.write(f"Molecular Weight:  {properties['molecular_weight']:.2f} Da\n")
        f.write(f"                      ({properties['molecular_weight']/1000:.2f} kDa)\n")
        f.write(f"Isoelectric Point:    {properties['isoelectric_point']:.2f}\n")
        f.write(f"Aromaticity:          {properties['aromaticity']:.4f}\n")
        f.write(f"Instability Index:    {properties['instability_index']:.2f}\n")
        f.write(f"GRAVY:                {properties['gravy']:.4f}\n")
        f.write("\n")

        # Secondary Structure Prediction
        helix, turn, sheet = properties['secondary_structure']
        f.write("SECONDARY STRUCTURE PREDICTION (Fraction)\n\n")
       
        f.write(f"Helix:                {helix:.4f} ({helix*100:.2f}%)\n")
        f.write(f"Turn:                 {turn:.4f} ({turn*100:.2f}%)\n")
        f.write(f"Sheet:                {sheet:.4f} ({sheet*100:.2f}%)\n")
        f.write("\n")

        # Biological Interpretation
        f.write("BIOLOGICAL INTERPRETATION\n\n")
        f.write(f"Size Category:        {interpretations['size']}\n")
        f.write(f"Stability:            {interpretations['stability']}\n")
        f.write(f"Hydropathy:           {interpretations['hydropathy']}\n")
        f.write(f"Charge:               {interpretations['charge']}\n")
        f.write("\n")

        # Amino Acid Composition
        f.write("AMINO ACID COMPOSITION\n\n")
        f.write("AA\tCount\tPercentage\t""\n")

        # Sort amino acids by frequency
        sorted_aa = sorted(aa_composition.items(), key=lambda x: x[1], reverse=True)
        
        for aa, percentage in sorted_aa:
            count = int(percentage * properties['length'])
            f.write(f"{aa}\t{count}\t{percentage*100:.2f}%\n")
        
        f.write("\n")

         # Summary Statistics
        f.write("SUMMARY STATISTICS\n\n")
        
        
        # Count charged, polar, and hydrophobic residues
        charged_aa = sum([aa_composition.get(aa, 0) for aa in 'DEKR'])
        polar_aa = sum([aa_composition.get(aa, 0) for aa in 'STNQ'])
        hydrophobic_aa = sum([aa_composition.get(aa, 0) for aa in 'AVILMFW'])
        
        f.write(f"Charged residues:     {charged_aa*100:.2f}% (D,E,K,R)\n")
        f.write(f"Polar residues:       {polar_aa*100:.2f}% (S,T,N,Q)\n")
        f.write(f"Hydrophobic residues: {hydrophobic_aa*100:.2f}% (A,V,I,L,M,F,W)\n")
        f.write("\n\n\n")
        
        # Footer
       
        f.write("END OF REPORT\n")
        
print("QUALITY CONTROL ANALYSIS COMPLETE")
    
if __name__ == "__main__":
    
    perform_qc_analysis()

































