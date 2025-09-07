#!/usr/bin/env python3
"""
PyMOL script to analyze two protein structures:
1. Load two complex structures from different directories with the same filename
2. Align structures and calculate RMSD
3. Calculate SASA (Solvent Accessible Surface Area)
4. Find interacting residues within 4.0Å distance

Usage: pymol -c pymol_analysis.py -- <pdb_filename> <dir1> <dir2>
Example: pymol -c pymol_analysis.py -- protein.pdb /path/to/dir1 /path/to/dir2
"""

import sys
import os
import csv
from datetime import datetime
from pymol import cmd, stored

class AnalysisResults:
    """Class to store all analysis results"""
    def __init__(self):
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.filename = ""
        self.dir1 = ""
        self.dir2 = ""
        self.rmsd = 0.0
        self.backbone_rmsd = 0.0
        self.aligned_atoms = 0
        self.backbone_aligned_atoms = 0
        self.sasa1 = 0.0
        self.sasa2 = 0.0
        self.sasa_difference = 0.0
        self.interacting_pairs = []
        self.interface1_residues = []
        self.interface2_residues = []
        self.interface1_count = 0
        self.interface2_count = 0

def load_structures(filename, dir1, dir2, results):
    """Load two PDB structures from different directories"""
    
    # Construct full paths
    pdb1_path = os.path.join(dir1, filename)
    pdb2_path = os.path.join(dir2, filename)
    
    # Check if files exist
    if not os.path.exists(pdb1_path):
        print(f"Error: File {pdb1_path} does not exist!")
        return False
    if not os.path.exists(pdb2_path):
        print(f"Error: File {pdb2_path} does not exist!")
        return False
    
    # Load structures with distinct names
    cmd.load(pdb1_path, "struct1")
    cmd.load(pdb2_path, "struct2")
    
    # Store paths in results
    results.filename = filename
    results.dir1 = dir1
    results.dir2 = dir2
    
    print(f"Loaded structure 1 from: {pdb1_path}")
    print(f"Loaded structure 2 from: {pdb2_path}")
    
    return True

def align_and_calculate_rmsd(results):
    """Align structures and calculate RMSD"""
    
    # Align struct2 to struct1
    alignment_result = cmd.align("struct2", "struct1")
    results.rmsd = alignment_result[0]
    results.aligned_atoms = alignment_result[1]
    
    # Also calculate RMSD for backbone atoms only
    backbone_result = cmd.align("struct2 and name CA+C+N+O", "struct1 and name CA+C+N+O")
    results.backbone_rmsd = backbone_result[0]
    results.backbone_aligned_atoms = backbone_result[1]
    
    print(f"RMSD calculation completed")
    
    return results.rmsd, results.backbone_rmsd

def calculate_sasa(results):
    """Calculate SASA for both structures"""
    
    # Calculate SASA for structure 1
    results.sasa1 = cmd.get_area("struct1", state=1)
    
    # Calculate SASA for structure 2
    results.sasa2 = cmd.get_area("struct2", state=1)
    
    # Calculate SASA difference
    results.sasa_difference = abs(results.sasa1 - results.sasa2)
    
    print(f"SASA calculation completed")
    
    return results.sasa1, results.sasa2

def find_interacting_residues(results, cutoff=4.0):
    """Find interacting residues within cutoff distance"""
    
    print(f"Finding interacting residues within {cutoff}Å...")
    
    # Get all residues from both structures
    stored.residues1 = []
    stored.residues2 = []
    
    cmd.iterate("struct1 and name CA", "stored.residues1.append((resi, resn, chain))")
    cmd.iterate("struct2 and name CA", "stored.residues2.append((resi, resn, chain))")
    
    # Find interacting residue pairs
    for res1 in stored.residues1:
        resi1, resn1, chain1 = res1
        selection1 = f"struct1 and resi {resi1} and chain {chain1}"
        
        for res2 in stored.residues2:
            resi2, resn2, chain2 = res2
            selection2 = f"struct2 and resi {resi2} and chain {chain2}"
            
            # Calculate minimum distance between residues
            distance = cmd.distance("temp_dist", selection1, selection2)
            cmd.delete("temp_dist")
            
            if distance <= cutoff:
                results.interacting_pairs.append({
                    'struct1_res': f"{resn1}{resi1}:{chain1}",
                    'struct2_res': f"{resn2}{resi2}:{chain2}",
                    'distance': distance
                })
    
    # Sort by distance
    results.interacting_pairs.sort(key=lambda x: x['distance'])
    
    print(f"Found {len(results.interacting_pairs)} interacting residue pairs")
    
    return results.interacting_pairs

def analyze_interface_residues(results):
    """Analyze interface residues in more detail"""
    
    print(f"Analyzing interface residues...")
    
    # Find residues at the interface (within 5Å of the other structure)
    cmd.select("interface1", "struct1 within 5.0 of struct2")
    cmd.select("interface2", "struct2 within 5.0 of struct1")
    
    # Count interface residues
    results.interface1_count = cmd.count_atoms("interface1 and name CA")
    results.interface2_count = cmd.count_atoms("interface2 and name CA")
    
    # Get interface residue details
    stored.interface1_residues = []
    stored.interface2_residues = []
    
    cmd.iterate("interface1 and name CA", "stored.interface1_residues.append(f'{resn}{resi}:{chain}')")
    cmd.iterate("interface2 and name CA", "stored.interface2_residues.append(f'{resn}{resi}:{chain}')")
    
    results.interface1_residues = sorted(set(stored.interface1_residues))
    results.interface2_residues = sorted(set(stored.interface2_residues))
    
    # Clean up selections
    cmd.delete("interface1")
    cmd.delete("interface2")

def save_aligned_structures(output_dir="./"):
    """Save aligned structures"""
    
    output_path1 = os.path.join(output_dir, "aligned_struct1.pdb")
    output_path2 = os.path.join(output_dir, "aligned_struct2.pdb")
    
    cmd.save(output_path1, "struct1")
    cmd.save(output_path2, "struct2")
    
    print(f"Aligned structures saved to {output_path1} and {output_path2}")

def write_summary_csv(results, output_file="pymol_analysis_summary.csv"):
    """Write summary results to CSV"""
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header and basic info
        writer.writerow(['Analysis Summary'])
        writer.writerow(['Timestamp', results.timestamp])
        writer.writerow(['PDB Filename', results.filename])
        writer.writerow(['Directory 1', results.dir1])
        writer.writerow(['Directory 2', results.dir2])
        writer.writerow([])
        
        # Write RMSD results
        writer.writerow(['RMSD Analysis'])
        writer.writerow(['Metric', 'Value', 'Unit'])
        writer.writerow(['Total RMSD', f'{results.rmsd:.3f}', 'Å'])
        writer.writerow(['Backbone RMSD', f'{results.backbone_rmsd:.3f}', 'Å'])
        writer.writerow(['Aligned Atoms', results.aligned_atoms, 'atoms'])
        writer.writerow(['Backbone Aligned Atoms', results.backbone_aligned_atoms, 'atoms'])
        writer.writerow([])
        
        # Write SASA results
        writer.writerow(['SASA Analysis'])
        writer.writerow(['Structure', 'SASA', 'Unit'])
        writer.writerow(['Structure 1', f'{results.sasa1:.2f}', 'Ų'])
        writer.writerow(['Structure 2', f'{results.sasa2:.2f}', 'Ų'])
        writer.writerow(['SASA Difference', f'{results.sasa_difference:.2f}', 'Ų'])
        writer.writerow([])
        
        # Write interface analysis
        writer.writerow(['Interface Analysis'])
        writer.writerow(['Structure', 'Interface Residues Count'])
        writer.writerow(['Structure 1', results.interface1_count])
        writer.writerow(['Structure 2', results.interface2_count])
        writer.writerow([])
        
        # Write interface residues
        writer.writerow(['Structure 1 Interface Residues'])
        writer.writerow([', '.join(results.interface1_residues)])
        writer.writerow([])
        writer.writerow(['Structure 2 Interface Residues'])
        writer.writerow([', '.join(results.interface2_residues)])
        writer.writerow([])
        
        # Write interacting pairs summary
        writer.writerow(['Interacting Residue Pairs Summary'])
        writer.writerow(['Total Pairs within 4.0Å', len(results.interacting_pairs)])
        writer.writerow([])

def write_interactions_csv(results, output_file="pymol_interacting_residues.csv"):
    """Write detailed interacting residues to CSV"""
    
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        writer.writerow(['Structure1_Residue', 'Structure2_Residue', 'Distance_Angstrom'])
        
        # Write all interacting pairs
        for pair in results.interacting_pairs:
            writer.writerow([
                pair['struct1_res'],
                pair['struct2_res'],
                f"{pair['distance']:.2f}"
            ])

def main():
    """Main analysis function"""
    
    # Parse command line arguments
    if len(sys.argv) < 4:
        print("Usage: pymol -c pymol_analysis.py -- <pdb_filename> <dir1> <dir2>")
        print("Example: pymol -c pymol_analysis.py -- protein.pdb /path/to/dir1 /path/to/dir2")
        sys.exit(1)
    
    filename = sys.argv[1]
    dir1 = sys.argv[2]
    dir2 = sys.argv[3]
    
    print(f"PyMOL Structure Analysis")
    print(f"========================")
    print(f"PDB filename: {filename}")
    print(f"Directory 1: {dir1}")
    print(f"Directory 2: {dir2}")
    
    # Initialize results storage
    results = AnalysisResults()
    
    # Initialize PyMOL
    cmd.reinitialize()
    
    # Load structures
    if not load_structures(filename, dir1, dir2, results):
        sys.exit(1)
    
    # Perform analysis
    print("Performing structural alignment and RMSD calculation...")
    align_and_calculate_rmsd(results)
    
    print("Calculating SASA...")
    calculate_sasa(results)
    
    print("Finding interacting residues...")
    find_interacting_residues(results, cutoff=4.0)
    
    print("Analyzing interface residues...")
    analyze_interface_residues(results)
    
    # Save aligned structures
    print("Saving aligned structures...")
    save_aligned_structures()
    
    # Write results to CSV files
    print("Writing results to CSV files...")
    summary_file = f"pymol_analysis_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    interactions_file = f"pymol_interacting_residues_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
    
    write_summary_csv(results, summary_file)
    write_interactions_csv(results, interactions_file)
    
    print(f"\n=== ANALYSIS COMPLETED ===")
    print(f"Summary results saved to: {summary_file}")
    print(f"Interacting residues saved to: {interactions_file}")
    print(f"Aligned structures saved to: aligned_struct1.pdb and aligned_struct2.pdb")
    
    # Print brief summary to console
    print(f"\n=== BRIEF SUMMARY ===")
    print(f"Total RMSD: {results.rmsd:.3f} Å")
    print(f"Backbone RMSD: {results.backbone_rmsd:.3f} Å")
    print(f"Structure 1 SASA: {results.sasa1:.2f} Ų")
    print(f"Structure 2 SASA: {results.sasa2:.2f} Ų")
    print(f"Interacting residue pairs (≤4.0Å): {len(results.interacting_pairs)}")
    print(f"Interface residues: Struct1={results.interface1_count}, Struct2={results.interface2_count}")

if __name__ == "__main__":
    main()