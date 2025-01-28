# BioSeqAnalyzer: Biological Sequence Analysis and Visualization with Python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from fpdf import FPDF


def selectFile():
    # Get all FASTA files in the current directory
    fastaFiles = [file for file in os.listdir() if file.endswith(".fasta")]

    # Check if no FASTA files are found
    if not fastaFiles:
        print("No FASTA files found in the current directory.")
        return
     
    # Display available FASTA files
    print("\nAvailable FASTA files:")
    for index, fileName in enumerate(fastaFiles, 1):
        print(f"{index}. {fileName}")
    
      # Prompt user to select a file by number
    userChoice = input("\nEnter the number of the FASTA file you want to use: ")

    # Validate user input
    if userChoice.isdigit():
        fileIndex = int(userChoice) - 1
        if 0 <= fileIndex and fileIndex < len(userChoice):
            selectedFile = fastaFiles[fileIndex]
            print(f"\nYou selected: {selectedFile}")
            return selectedFile
        else:
            print("\nInvalid choice. Please select a valid number from the list.")
            return
    else:
        print("\nInvalid input. Please enter a number corresponding to the file.")
        return

# Function to load a FASTA file and return sequences as a dictionary
def getSequences(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

# Function to calculate GC content using BioPython
def calculateGCContent(sequence):
    return gc_fraction(sequence) * 100

# Function to get the reverse complement of a DNA sequence using BioPython
def reverseComplement(sequence):
    seqObj = Seq(sequence)
    return str(seqObj.reverse_complement())

# Function to transcribe DNA to RNA (T -> U)
def transcribeToRNA(sequence):
    seqObj = Seq(sequence)
    return str(seqObj.transcribe())

# Function to calculate nucleotide composition
def nucleotideComposition(sequence):
    total = len(sequence)
    return {
        'A (%)': (sequence.count('A') / total) * 100,
        'T (%)': (sequence.count('T') / total) * 100,
        'G (%)': (sequence.count('G') / total) * 100,
        'C (%)': (sequence.count('C') / total) * 100,
    }

def main():
    
    selectedFile = selectFile()
    
    # Load sequences from the FASTA file
    sequences = getSequences(selectedFile)
    if not sequences:
        print("No sequences found in the FASTA file.")
        return

    print(sequences)

main()