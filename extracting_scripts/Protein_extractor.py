from Bio import SeqIO

# Input GenBank file
genbank_file = "input_file.gbk"  # Replace with your GenBank file path

# Output FASTA file for amino acid sequences
output_fasta = "amino_acid_sequences.fasta"

# Open the output file for writing
with open(output_fasta, "w") as fasta_out:
    # Parse the GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            # Check if the feature is a CDS (coding sequence)
            if feature.type == "CDS":
                # Extract the protein sequence if available
                if "translation" in feature.qualifiers:
                    protein_seq = feature.qualifiers["translation"][0]
                    # Extract the protein ID or locus tag for naming
                    protein_id = feature.qualifiers.get("protein_id", ["unknown_protein"])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", ["unknown_locus"])[0]
                    # Write the protein sequence to the FASTA file
                    fasta_out.write(f">{protein_id} {locus_tag}\n{protein_seq}\n")

print(f"Amino acid sequences have been extracted to {output_fasta}")