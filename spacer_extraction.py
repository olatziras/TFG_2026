#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CRISPR Array Extractor

Description:
    This script parses FASTQ files to extract CRISPR-Cas 
    arrays based on specific repeat sequences. It identifies full and partial 
    repeats, extracts the intervening spacers, filters them by length, and 
    deduplicates the arrays. Finally, it generates FASTA, FASTQ, and CSV 
    summary files for downstream analysis.

Author: Andrea Fernández
Date: April 2026

Usage:
    python spacers_analysis_complete.py input_file.fastq [additional_files.fastq ...]
"""

import csv
import os
import argparse
from collections import Counter

# --- Script Configuration ---

# Repeat sequence definitions
# The full repeat is considered the 29 bp "stem", as the last nucleotide
# is treated as part of the spacer.
FULL_REPEAT_SEQ = "GAGTTCCCCGCGCCAGCGGGGATAAACC" # Length: 29
FULL_REPEAT_LEN = len(FULL_REPEAT_SEQ)

PARTIAL_REPEAT_SEQ = "GAGTTCCCCGCGCCAGCGG"      # Length: 19
PARTIAL_REPEAT_LEN = len(PARTIAL_REPEAT_SEQ)

# Search and filtering parameters
MAX_MISMATCHES = 1
# Limit for the "real" spacer length (excluding the 1st nucleotide of the PAM)
MAX_REAL_SPACER_LENGTH = 38

# --- Helper Functions ---

def reverse_complement(dna_seq):
    """Computes the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(dna_seq.upper()))

def count_mismatches(s1, s2):
    """Counts mismatches between two sequences of equal length."""
    return sum(1 for char1, char2 in zip(s1, s2) if char1 != char2)

def find_next_repeat(sequence, search_start_pos, pattern, max_mismatches):
    """
    Generic function to find the next repeat (full or partial).
    Returns (absolute_start_index, matched_sequence) or (-1, None) if not found.
    """
    pattern_len = len(pattern)
    for i in range(search_start_pos, len(sequence) - pattern_len + 1):
        sub_sequence = sequence[i : i + pattern_len]
        if count_mismatches(sub_sequence, pattern) <= max_mismatches:
            return i, sub_sequence
    return -1, None

# --- Core Logic Functions ---

def extract_arrays_from_read(read_id, sequence, orientation):
    """
    Extracts a single CRISPR array from a sequence, assuming a maximum of one per read.
    An array is represented as a list of dictionaries, where each dictionary contains spacer data.
    Returns a list containing the found array, or an empty list if none is found.
    """
    # 1. Find the start of the array (the first full repeat)
    first_repeat_idx, first_repeat_seq = find_next_repeat(sequence, 0, FULL_REPEAT_SEQ, MAX_MISMATCHES)

    if first_repeat_idx == -1:
        return [] # No array found in this read

    current_array_spacers = []
    prev_repeat_idx = first_repeat_idx
    prev_repeat_seq = first_repeat_seq
    
    # The spacer sequence begins immediately AFTER the 29 bp repeat
    search_pos_for_next = prev_repeat_idx + len(prev_repeat_seq)

    while True:
        # Priority 1: Search for the next full repeat to continue the array
        next_full_idx, next_full_seq = find_next_repeat(sequence, search_pos_for_next, FULL_REPEAT_SEQ, MAX_MISMATCHES)

        is_terminal = False
        chosen_next_repeat_idx = -1
        chosen_next_repeat_seq = None

        if next_full_idx != -1:
            # Full repeat found, the array continues
            chosen_next_repeat_idx = next_full_idx
            chosen_next_repeat_seq = next_full_seq
        else:
            # No more full repeats found, search for a partial repeat to terminate the array
            next_partial_idx, next_partial_seq = find_next_repeat(sequence, search_pos_for_next, PARTIAL_REPEAT_SEQ, MAX_MISMATCHES)
            if next_partial_idx != -1:
                chosen_next_repeat_idx = next_partial_idx
                chosen_next_repeat_seq = next_partial_seq
                is_terminal = True
            else:
                # No repeat found, the array ends abruptly
                break

        # Extract the spacer sequence
        spacer_seq = sequence[search_pos_for_next : chosen_next_repeat_idx]

        # Avoid appending empty spacers if repeats are adjacent
        if not spacer_seq: 
            if is_terminal: break
            prev_repeat_idx = chosen_next_repeat_idx
            prev_repeat_seq = chosen_next_repeat_seq
            search_pos_for_next = prev_repeat_idx + len(prev_repeat_seq)
            continue

        # Calculate the "real" length for downstream filtering
        real_spacer_len = len(spacer_seq) - 1 if len(spacer_seq) > 0 else 0

        current_array_spacers.append({
            "Read_ID": f"{read_id}_{orientation}",
            "Previous_Repeat": prev_repeat_seq,
            "Spacer": spacer_seq,
            "Next_Repeat": chosen_next_repeat_seq,
            "Real_Spacer_Length": real_spacer_len
        })

        if is_terminal:
            break

        # Move indices forward for the next iteration
        prev_repeat_idx = chosen_next_repeat_idx
        prev_repeat_seq = chosen_next_repeat_seq
        search_pos_for_next = prev_repeat_idx + len(prev_repeat_seq)

    # Return the array wrapped in a list if spacers were successfully extracted
    if current_array_spacers:
        return [current_array_spacers]
    else:
        return []

def read_fastq(filepath):
    """Generator to read FASTQ records (yields ID and sequence only)."""
    try:
        with open(filepath, 'r') as f:
            while True:
                header = f.readline()
                if not header: break
                if not header.startswith('@'): continue
                seq = f.readline().strip()
                f.readline() # Skip '+' line
                f.readline() # Skip quality scores
                yield header.strip().split()[0][1:], seq
    except FileNotFoundError:
        print(f"Error: FASTQ file not found at {filepath}")
        return

def process_file(fastq_filename):
    """
    Processes a single FASTQ file, extracting CRISPR arrays and generating output files.
    """
    print(f"\n--- Processing file: {fastq_filename} ---")
    base_name = os.path.basename(fastq_filename).removesuffix('.fastq').removesuffix('.fq')

    # --- 1. Initial Array Extraction ---
    initial_arrays = []
    reads_processed = 0
    for read_id, sequence in read_fastq(fastq_filename):
        reads_processed += 1
        if reads_processed % 10000 == 0:
            print(f"  ... {reads_processed} reads processed")

        # Check both forward and reverse complement orientations
        initial_arrays.extend(extract_arrays_from_read(read_id, sequence, "FWD"))
        initial_arrays.extend(extract_arrays_from_read(read_id, reverse_complement(sequence), "REV"))

    initial_spacer_count = sum(len(arr) for arr in initial_arrays)
    print(f"  Found {len(initial_arrays)} initial arrays containing a total of {initial_spacer_count} spacers.")

    # --- 2. Array Filtering and Deduplication ---
    final_arrays = []
    seen_array_signatures = set()

    for array in initial_arrays:
        # Filter by spacer length (discard arrays with excessively long spacers)
        is_too_long = any(spacer['Real_Spacer_Length'] > MAX_REAL_SPACER_LENGTH for spacer in array)
        if is_too_long:
            continue

        # Deduplicate based on the sequence of spacers within the array
        array_signature = tuple(spacer['Spacer'] for spacer in array)
        if array_signature in seen_array_signatures:
            continue

        seen_array_signatures.add(array_signature)
        final_arrays.append(array)

    final_spacer_count = sum(len(arr) for arr in final_arrays)
    print(f"  After filtering and deduplication, {len(final_arrays)} unique arrays remain with {final_spacer_count} spacers.")

    if not final_arrays:
        print("  No arrays remained after filtering. Output files will not be generated.")
        return

    # --- 3. Output Generation ---

    # Detailed spacers CSV file
    csv_spacers_filename = f"spacers_{base_name}.csv"
    print(f"  Generating detailed CSV file: {csv_spacers_filename}")
    with open(csv_spacers_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Array_ID", "Spacer_Position", "Previous_Repeat", "Spacer_Sequence", "Next_Repeat"])
        array_id_counter = 1
        for array in final_arrays:
            for i, spacer_info in enumerate(array, 1):
                writer.writerow([
                    array_id_counter,
                    i,
                    spacer_info['Previous_Repeat'],
                    spacer_info['Spacer'],
                    spacer_info['Next_Repeat']
                ])
            array_id_counter += 1

    # Summary statistics CSV file
    summary_filename = f"summary_{base_name}.csv"
    print(f"  Generating summary file: {summary_filename}")
    with open(summary_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # a. General Statistics
        writer.writerow(["Statistic", "Value"])
        writer.writerow(["Reads Processed", reads_processed])
        writer.writerow(["Initial Arrays Found", len(initial_arrays)])
        writer.writerow(["Initial Spacers Found", initial_spacer_count])
        writer.writerow(["Arrays after Filtering/Deduplication", len(final_arrays)])
        writer.writerow(["Spacers in Final Arrays", final_spacer_count])
        writer.writerow([]) # Blank line

        # b. Array Length Distribution
        writer.writerow(["Array Length (spacers)", "Count"])
        array_len_counts = Counter(len(arr) for arr in final_arrays)
        for length, count in sorted(array_len_counts.items()):
            writer.writerow([length, count])
        writer.writerow([])

        # c. Spacer Length Distribution (Binned)
        writer.writerow(["Spacer Length (real)", "Count"])
        spacer_len_bins = Counter()
        for array in final_arrays:
            for spacer in array:
                real_len = spacer['Real_Spacer_Length']
                if real_len < 30:
                    spacer_len_bins['<30'] += 1
                elif real_len > 35:
                    spacer_len_bins['>35'] += 1
                else:
                    spacer_len_bins[str(real_len)] += 1
        
        # Sort bins for output
        order = ['<30', '30', '31', '32', '33', '34', '35', '>35']
        for bin_name in order:
            if bin_name in spacer_len_bins:
                 writer.writerow([bin_name, spacer_len_bins[bin_name]])


    # FASTA and FASTQ files for downstream alignment (e.g., BWA)
    fastq_filename_out = f"unique_spacers_{base_name}.fastq"
    fasta_filename_out = f"unique_spacers_{base_name}.fasta"
    print(f"  Generating FASTA/FASTQ files for alignment: {fasta_filename_out}, {fastq_filename_out}")
    with open(fastq_filename_out, 'w') as fq_out, open(fasta_filename_out, 'w') as fa_out:
        array_id_counter = 1
        for array in final_arrays:
            for i, spacer_info in enumerate(array, 1):
                spacer_seq = spacer_info['Spacer']
                header_id = f"{array_id_counter}_{i}"
                
                # FASTQ format (using dummy 'I' for high quality scores)
                fq_out.write(f"@{header_id}\n")
                fq_out.write(f"{spacer_seq}\n")
                fq_out.write("+\n")
                fq_out.write("I" * len(spacer_seq) + "\n")

                # FASTA format
                fa_out.write(f">{header_id}\n")
                fa_out.write(f"{spacer_seq}\n")
            array_id_counter += 1

    print(f"--- Processing of {fastq_filename} complete. ---")


# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extracts, filters, and analyzes CRISPR-Cas arrays from FASTQ files.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'fastq_files',
        metavar='FILE.fastq',
        nargs='+',  # Accepts one or more files
        help="One or more FASTQ files to analyze."
    )
    args = parser.parse_args()

    for filename in args.fastq_files:
        if not os.path.exists(filename):
            print(f"Error: Input file '{filename}' not found. Skipping.")
            continue
        process_file(filename)

    print("\nAll files have been processed successfully.")