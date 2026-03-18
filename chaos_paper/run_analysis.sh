#!/bin/bash

# Arguments passed from the command line
DIM=$1
GEOM=$2

# Check if both arguments are provided
if [ -z "$DIM" ] || [ -z "$GEOM" ]; then
    echo "Error: Missing arguments."
    echo "Usage: ./run_analysis.sh <dim> <geom>"
    echo "Example: ./run_analysis.sh 2D full"
    exit 1
fi

echo "========================================"
echo " Initializing Analysis for $DIM - $GEOM"
echo "========================================"

OUTPUT_DIR_NAME="post_processing"

# 1. Execute Python script (Create CSV, Pointplot, Curves and AUC Bars)
echo "[1/2] Running Python..."
python3 post_processing_script.py --dim "$DIM" --geom "$GEOM" --output_dir "$OUTPUT_DIR_NAME"

# Check if the Python script executed successfully
if [ $? -ne 0 ]; then
    echo "Error in execution of Python script. Aborting."
    exit 1
fi

# Directory where the R script will look for the CSV and save results
OUTPUT_DIR="./$DIM/$GEOM/$OUTPUT_DIR_NAME"

# 2. Execute R script (Inferential Statistics and Forest Plot)
echo "[2/2] Running R..."
Rscript statistical_analysis.R "$OUTPUT_DIR"

echo "========================================"
echo " Analysis Completed for $DIM - $GEOM"
echo " Results saved in: $OUTPUT_DIR"
echo "========================================"
