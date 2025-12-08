# GeneAffinity
This is a program developed by Cole Kuznitz as part of the McHugh Lab at UCSD.  It was developed as a tool to further look at binding between our genes of interest, and expand upon our search.  

# FASTA Binding Pipeline

## Overview
This is a Streamlit web app to find potential RNA-RNA binding regions using a sliding-window approach.

## Features
- Upload FASTA sequences for Gene A and Gene B
- Set window size, step, flank, energy cutoff, and threads
- Runs fast candidate search (RiSearch2) and rescoring (IntaRNA)
- Outputs aggregated results in CSV format

## Design
Input sequences
   ↓
Encode / hash windows (fast)
   ↓
Candidate binding regions
   ↓
Deep nucleotide-level analysis (slow but rare)
   ↓
Scored binding sites

## Installation
1. Install Python 3.9+  
2. Install dependencies:


```bash
pip install -r requirements.txt
