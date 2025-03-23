from Bio import Entrez, SeqIO
import streamlit as st
from io import BytesIO

def fetch_fasta_sequence(accession):
    """
    Retrieves a FASTA sequence from the NCBI database.

    Parameters:
        accession (str): Accession number of the sequence.

    Returns:
        str: FASTA sequence.
    """
    Entrez.email = "abhinavrana18july@gmail.com"  # Replace with your email

    # Search for the accession in the NCBI database
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()
        return record
    except Exception as e:
        return f"Failed to retrieve {accession}: {str(e)}"


def generate_fasta_file(input_file):
    """
    Generates a FASTA file from a list of accession numbers provided in an input file.

    Parameters:
        input_file (file-like object): A file-like object containing accession numbers (one per line). The file should be in plain text format.

    Returns:
        BytesIO: In-memory FASTA file containing the retrieved sequences.
    """
    if input_file is not None:
        try:
            # Read accession numbers from the input file
            accessions = [line.strip() for line in input_file.read().decode("utf-8").splitlines()]
        except Exception as e:
            st.error("Error reading the file. Please ensure it is a text file containing accession numbers.")
            return None

        # Retrieve the FASTA sequences
        results = {}
        total_accessions = len(accessions)
        progress_bar = st.progress(0)  # Initialize the progress bar
        percentage_text = st.empty()  # Placeholder for percentage text

        for i, accession in enumerate(accessions):
            results[accession] = fetch_fasta_sequence(accession)
            progress_percentage = (i + 1) / total_accessions  # Calculate percentage
            progress_bar.progress(progress_percentage)  # Update the progress bar
            percentage_text.text(f"Progress: {int(progress_percentage * 100)}%")  # Update percentage text

        # Save the retrieved sequences into a BytesIO object
        fasta_io = io.BytesIO()
        with fasta_io as file:
            for accession, sequence in results.items():
                file.write(f">{accession}\n{sequence}\n".encode('utf-8'))

        fasta_io.seek(0)  # Reset the file pointer to the beginning
        return fasta_io, total_accessions
