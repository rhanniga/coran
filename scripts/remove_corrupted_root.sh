#!/bin/bash

# This script removes corrupted ROOT files from a specified directory.
# A ROOT file is considered corrupted if it cannot be opened correctly,
# which in ROOT terminology is called a "zombie" file.

# Check if a directory is provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory_path>"
    exit 1
fi

TARGET_DIR=$1

# Check if the provided directory exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' not found."
    exit 1
fi

echo "Scanning for corrupted ROOT files in: $TARGET_DIR"

# Find all files ending in .root in the target directory and its subdirectories
find "$TARGET_DIR" -type f -name "*.root" | while read -r filepath; do
    # Use ROOT to check if the file is a "zombie".
    # The -l flag prevents the splash screen, -b runs in batch mode, and -q exits after processing.
    # The -e flag executes the provided ROOT C++ code.
    # TFile f(...) attempts to open the file.
    # If f.IsZombie() is true, the file is corrupted, and we exit ROOT with a non-zero status code (1).
    root -l -b -q -e "TFile f(\"$filepath\"); if (f.IsZombie()) { gSystem->Exit(1); }"
    
    # Check the exit status of the last command
    if [ $? -ne 0 ]; then
        echo "Found corrupted file: $filepath. Removing it."
        rm "$filepath"
    fi
done

echo "Scan complete."
