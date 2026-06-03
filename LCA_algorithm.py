#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LCA (Lowest Common Ancestor) Analysis Tool

Description:
    This script calculates the Lowest Common Ancestor (LCA) for BLAST hits. 
    It groups alignments by query, applies a "Best-Hit" strategy based on 
    the maximum bitscore, and determines the common taxonomic ancestor. 
    It also features a fallback mechanism to fetch 
    missing TaxIDs directly from the NCBI Entrez API.

Author: Andrea Fernández
Date: April 2026

Usage:
    python LCA_algorithm.py <blast_output.tsv> [output_filename.tsv]

BLAST Format Requirement:
    qseqid sseqid pident length mismatch gapopen evalue bitscore staxids stitle
"""

import sys
import csv
import time
import json
import urllib.request
import urllib.error
from collections import defaultdict
from ete3 import NCBITaxa

# Initialize NCBI taxonomy database
ncbi = NCBITaxa()

# --- COLUMN CONFIGURATION ---
COL_QSEQID   = 0
COL_SSEQID   = 1
COL_PIDENT   = 2
COL_BITSCORE = 7  
COL_TAXID    = 8  

def fetch_taxids_from_ncbi(accessions):
    """
    Retrieves TaxIDs from NCBI for accessions that returned '0' or missing values.
    """
    url_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    mapping = {}
    unique_accs = list(set([a for a in accessions if a]))
    
    if not unique_accs:
        return mapping

    print(f"   > API: Correcting {len(unique_accs)} TaxIDs (value 0) using NCBI...")
    batch_size = 200
    
    for i in range(0, len(unique_accs), batch_size):
        batch = unique_accs[i:i+batch_size]
        ids_str = ",".join(batch)
        url = f"{url_base}?db=nucleotide&id={ids_str}&retmode=json"
        
        try:
            with urllib.request.urlopen(url) as response:
                data = json.loads(response.read().decode())
                if 'result' in data:
                    for uid in data['result']:
                        if uid == 'uids': continue
                        item = data['result'][uid]
                        taxid = int(item.get('taxid', 0))
                        caption = item.get('caption', '')
                        acc_ver = item.get('accessionversion', '')
                        if taxid > 0:
                            if caption: mapping[caption] = taxid
                            if acc_ver: mapping[acc_ver] = taxid
            time.sleep(0.35) # Delay to respect NCBI API limits
            sys.stdout.write(f"\r     Processed: {min(i+batch_size, len(unique_accs))}/{len(unique_accs)}")
            sys.stdout.flush()
        except Exception as e:
            print(f"\n     [Warning] API Error: {e}")
            
    print("\n   > Recovery finished.")
    return mapping

def perform_lca_analysis(blast_file, output_file, identity_threshold=95.0):
    """
    Main function to parse BLAST results and compute the LCA.
    """
    print(f"1. Reading BLAST file: {blast_file}")
    
    file_buffer = []
    missing_accessions = set()
    
    # Pass 1: Read file and detect missing TaxIDs
    with open(blast_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) <= COL_TAXID: continue
            
            file_buffer.append(row)
            
            # Check if staxids is "0", empty, or N/A
            raw_taxid = row[COL_TAXID].strip()
            if raw_taxid == "0" or not raw_taxid or raw_taxid == "N/A":
                missing_accessions.add(row[COL_SSEQID])

    # Pass 2: Correct TaxIDs via API if needed
    recovered_map = {}
    if missing_accessions:
        print(f"2. Detected {len(missing_accessions)} hits with missing TaxIDs.")
        recovered_map = fetch_taxids_from_ncbi(list(missing_accessions))
    else:
        print("2. No missing TaxIDs detected.")

    # Pass 3: Group hits and apply Best-Hit strategy
    print("3. Applying Best-Hit strategy (maximum bitscore)...")
    blast_hits = defaultdict(list)
    
    for row in file_buffer:
        try:
            query_id = row[COL_QSEQID]
            subject_acc = row[COL_SSEQID]
            identity = float(row[COL_PIDENT])
            bitscore = float(row[COL_BITSCORE])
            
            # Extract TaxID (handling potential multiple IDs separated by ';')
            raw_val = row[COL_TAXID].strip()
            valid_taxids = []
            
            if raw_val != "0" and raw_val != "" and raw_val != "N/A":
                valid_taxids = [int(tid) for tid in raw_val.split(';') if tid.isdigit()]
            
            # If missing, use the recovery map
            if not valid_taxids:
                taxid_api = recovered_map.get(subject_acc, 0)
                if taxid_api == 0 and '.' in subject_acc:
                    taxid_api = recovered_map.get(subject_acc.split('.')[0], 0)
                if taxid_api > 0:
                    valid_taxids = [taxid_api]
            
            # Filter by identity threshold and store the TaxIDs associated with this bitscore
            if identity >= identity_threshold and valid_taxids:
                blast_hits[query_id].append((bitscore, valid_taxids))
                
        except (ValueError, IndexError):
            continue

    # Pass 4: LCA Calculation
    print(f"4. Calculating LCA for {len(blast_hits)} queries...")
    with open(output_file, 'w') as out:
        out.write("QueryID\tLCA_TaxID\tScientific_Name\tRank\tLineage\n")
        
        for qid, hits in blast_hits.items():
            if not hits: continue
            
            # Strategy: Keep only hits with the maximum bitscore
            max_score = max(hits, key=lambda x: x[0])[0]
            
            # Collect all TaxIDs with the maximum bitscore
            best_taxids = []
            for score, tids in hits:
                if score >= max_score:
                    best_taxids.extend(tids)
            
            # Remove duplicates for LCA calculation
            unique_taxids = list(set(best_taxids))
            
            if not unique_taxids:
                continue

            try:
                # If only 1 TaxID exists, it is the LCA
                if len(unique_taxids) == 1:
                    lca_taxid = unique_taxids[0]
                    
                    # Get scientific name
                    name_dict = ncbi.get_taxid_translator([lca_taxid])
                    lca_name = name_dict.get(lca_taxid, "Unknown")
                    
                    # Get rank
                    rank_dict = ncbi.get_rank([lca_taxid])
                    lca_rank = rank_dict.get(lca_taxid, "no rank")
                
                # If multiple TaxIDs exist, build topology to find the common ancestor
                else:
                    tree = ncbi.get_topology(unique_taxids)
                    lca_node = tree.get_tree_root()
                    lca_taxid = lca_node.taxid
                    lca_name = lca_node.name
                    lca_rank = lca_node.rank
                
                # Retrieve full lineage names
                lineage = ncbi.get_lineage(lca_taxid)
                names = ncbi.get_taxid_translator(lineage)
                lineage_names = ";".join([names[t] for t in lineage])
                
                out.write(f"{qid}\t{lca_taxid}\t{lca_name}\t{lca_rank}\t{lineage_names}\n")
                
            except Exception as e:
                print(f"\n[Warning] Error calculating LCA for {qid} (TaxIDs: {unique_taxids}): {e}")
                continue

    print(f"\n✅ Analysis complete. Results saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 LCA_algorithm.py <blast_output.tsv> [output_file.tsv]")
        sys.exit(1)

    input_f = sys.argv[1]
    output_f = sys.argv[2] if len(sys.argv) >= 3 else input_f.rsplit('.', 1)[0] + "_LCA_BestHit.tsv"
    perform_lca_analysis(input_f, output_f)