# !pip install biopython

import warnings
from Bio import SeqIO, Entrez, Align
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.PDB import PDBList, PDBParser
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# ------------------------------------------------------------
# Suppress harmless warnings
# ------------------------------------------------------------
warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', PDBConstructionWarning)

# ------------------------------------------------------------
# Set email for NCBI
# ------------------------------------------------------------
Entrez.email = "saijaswanth561@gmail.com"

# ------------------------------------------------------------
# 1. Load Nucleotide and Protein Sequences
# ------------------------------------------------------------
nucleotide_record = SeqIO.read("NM_000633.fasta", "fasta")
protein_record = SeqIO.read("2F5N.fasta", "fasta")

nucleotide_seq = nucleotide_record.seq
protein_seq = protein_record.seq

print("=== Sequence Information ===")
print("Nucleotide ID:", nucleotide_record.id)
print("Nucleotide Length:", len(nucleotide_seq), "bases")
print("Protein ID:", protein_record.id)
print("Protein Length:", len(protein_seq), "amino acids\n")

# ------------------------------------------------------------
# 2. Nucleotide Composition & GC Content
# ------------------------------------------------------------
print("=== Nucleotide Composition ===")
print("A:", nucleotide_seq.count("A"))
print("T:", nucleotide_seq.count("T"))
print("G:", nucleotide_seq.count("G"))
print("C:", nucleotide_seq.count("C"))
gc_content = gc_fraction(nucleotide_seq) * 100
print(f"GC Content: {gc_content:.2f}%\n")

# ------------------------------------------------------------
# 3. Translation of DNA to Protein (handles partial codon issue)
# ------------------------------------------------------------
print("=== Translation of DNA to Protein ===")
trimmed_seq = nucleotide_seq[:len(nucleotide_seq)//3 * 3]  # trim to multiple of 3
translated_protein = trimmed_seq.translate(to_stop=True)

print("Translated Protein Sequence (first 100 aa):")
print(translated_protein[:100], "...")
print("Protein Length:", len(translated_protein), "amino acids\n")

# ------------------------------------------------------------
# 4. Genome Annotation (from GenBank)
# ------------------------------------------------------------
print("=== Genome Annotation (from GenBank) ===")
try:
    handle = Entrez.efetch(db="nucleotide", id=nucleotide_record.id,
                           rettype="gb", retmode="text")
    gb_record = SeqIO.read(handle, "genbank")
    handle.close()

    for feature in gb_record.features:
        if feature.type in ["CDS", "5'UTR", "3'UTR", "exon"]:
            print(f"Feature Type: {feature.type}")
            if "gene" in feature.qualifiers:
                print("  Gene:", feature.qualifiers["gene"][0])
            if "product" in feature.qualifiers:
                print("  Product:", feature.qualifiers["product"][0])
            if feature.type == "CDS":
                print("  CDS Location:", feature.location)
                print("  Protein ID:", feature.qualifiers.get("protein_id", ['N/A'])[0])
            print("-" * 50)
    print()
except Exception as e:
    print("Genome annotation not available. Error:", e, "\n")

# ------------------------------------------------------------
# 5. Pairwise Alignment
# ------------------------------------------------------------
print("=== Pairwise Alignment ===")
aligner = Align.PairwiseAligner()
aligner.mode = 'global'
alignment = aligner.align(translated_protein, protein_seq)
print("Alignment Score:", alignment.score)
print(alignment[0], "\n")

# ------------------------------------------------------------
# 6. BLASTN (Nucleotide vs nt)
# ------------------------------------------------------------
print("=== Running BLASTN === (this may take a few minutes)")
result_handle_n = NCBIWWW.qblast("blastn", "nt", nucleotide_seq)

with open("blastn_results.xml", "w") as out_n:
    out_n.write(result_handle_n.read())
result_handle_n.close()

print("BLASTN search completed. Results saved as 'blastn_results.xml'.\n")

# Parse BLASTN results
with open("blastn_results.xml") as blastn_out:
    blastn_record = NCBIXML.read(blastn_out)

print("Top 3 BLASTN Hits:\n")
for alignment in blastn_record.alignments[:3]:
    print("Title:", alignment.title)
    print("Length:", alignment.length)
    for hsp in alignment.hsps:
        print(f" Score: {hsp.score},  E-value: {hsp.expect}")
    print("-" * 50)

# ------------------------------------------------------------
# 7. BLASTP (Protein vs SwissProt)
# ------------------------------------------------------------
print("\n=== Running BLASTP === (this may take a few minutes)")
result_handle_p = NCBIWWW.qblast("blastp", "swissprot", protein_seq)

with open("blastp_results.xml", "w") as out_p:
    out_p.write(result_handle_p.read())
result_handle_p.close()

print("BLASTP search completed. Results saved as 'blastp_results.xml'.\n")

# Parse BLASTP results
with open("blastp_results.xml") as blastp_out:
    blastp_record = NCBIXML.read(blastp_out)

print("Top 3 BLASTP Hits:\n")
for alignment in blastp_record.alignments[:3]:
    print("Title:", alignment.title)
    print("Length:", alignment.length)
    for hsp in alignment.hsps:
        print(f" Score: {hsp.score},  E-value: {hsp.expect}")
    print("-" * 50)

# ------------------------------------------------------------
# 8. Protein Structure Retrieval and Parsing
# ------------------------------------------------------------
print("\n=== Protein Structure Parsing (PDB) ===")
pdbl = PDBList()
pdb_id = protein_record.id.split('|')[1]  # Correctly extract PDB ID
pdb_file = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb')
print("Downloaded PDB file:", pdb_file)

parser = PDBParser()
structure = parser.get_structure("Protein", pdb_file)
print("Model IDs:", [m.id for m in structure])
print("Chains:", [c.id for c in structure[0]])
if 'A' in structure[0]:
    print("Residues in Chain A:", len([r for r in structure[0]['A']]))