#!/usr/bin/env python3
import csv

# Standard Genetic Code Dictionary
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def parse_fasta(filepath):
    """Parses FASTA files into header-sequence pairs."""
    records = []
    header, seq = None, []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(seq)))
                    seq = []
                header = line[1:].split()[0]
            else:
                seq.append(line.upper())
        if header:
            records.append((header, "".join(seq)))
    return records

def calc_gc(seq):
    if not seq:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    return round(((g + c) / len(seq)) * 100, 2)

def translate(seq):
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        protein.append(CODON_TABLE.get(codon, "X"))
    return "".join(protein)

def main():
    sample_fasta = "sample.fasta"
    
    # 1. Create a demo FASTA if none exists
    with open(sample_fasta, "w") as f:
        f.write(">Seq1_Synthetic gene target\nATGCGTACGTTAGCCTAGCTAGCTAATCGATCG\n")
        f.write(">Seq2_High_GC control\nATGCCCGGGCCCGGGATGCGCGCGTAA\n")

    # 2. Parse & Analyze
    records = parse_fasta(sample_fasta)
    print(f"==> Parsed {len(records)} records from {sample_fasta}\n")

    # 3. Export metrics to CSV (for analysis.R)
    output_csv = "sequence_metrics.csv"
    with open(output_csv, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["ID", "Length", "GC_Content", "Protein"])
        
        for name, sequence in records:
            length = len(sequence)
            gc = calc_gc(sequence)
            prot = translate(sequence)
            writer.writerow([name, length, gc, prot])
            print(f"[{name}] Length: {length}bp | GC: {gc}% | Protein: {prot}")

    print(f"\n==> Saved sequence metrics to '{output_csv}'")

if __name__ == "__main__":
    main()
