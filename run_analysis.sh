#!/bin/bash

# Batch script to run PyMOL analysis
# Usage: ./run_analysis.sh or submit with sbatch run_analysis.sh

#SBATCH -N 1
#SBATCH --gpus=1
#SBATCH -J pymol_analysis
#SBATCH -o pymol_analysis_%j.out
#SBATCH -e pymol_analysis_%j.err

module purge
module load miniforge
source activate pymol

# Define directories and script path
DIR1=/data/run01/scz0sfc/PPI_eval/relax/data/original
DIR2=/data/run01/scz0sfc/PPI_eval/relax/data/af3
SCRIPT_PATH=/data/run01/scz0sfc/PPI_eval/relax/pymol_analysis.py

echo "Running PyMOL analysis..."
echo "Directory 1: $DIR1"
echo "Directory 2: $DIR2"
echo "Script path: $SCRIPT_PATH"

# Check if script exists
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: Script not found at $SCRIPT_PATH"
    exit 1
fi

# Check if directories exist
if [ ! -d "$DIR1" ]; then
    echo "Error: Directory 1 not found: $DIR1"
    exit 1
fi

if [ ! -d "$DIR2" ]; then
    echo "Error: Directory 2 not found: $DIR2"
    exit 1
fi

# Run PyMOL with the analysis script for each PDB file
echo "Processing 1OTR.pdb..."
pymol -c ${SCRIPT_PATH} -- 1OTR.pdb "$DIR1" "$DIR2"

echo "Processing 1UEL.pdb..."
pymol -c ${SCRIPT_PATH} -- 1UEL.pdb "$DIR1" "$DIR2"

echo "Processing 2G3Q.pdb..."
pymol -c ${SCRIPT_PATH} -- 2G3Q.pdb "$DIR1" "$DIR2"

echo "Processing 1buh.pdb..."
pymol -c ${SCRIPT_PATH} -- 1buh.pdb "$DIR1" "$DIR2"

echo "All analyses completed!"

# Optional: List generated files
echo "Generated files:"
ls -la pymol_*_*.csv aligned_*.pdb 2>/dev/null || echo "No output files found in current directory"