import argparse
from Bio import SeqIO, SearchIO
from Bio.Blast.Applications import NcbiblastxCommandline

def perform_blastx(query_file, ref_protein, output_file):
    blastx_cline = NcbiblastxCommandline(query=query_file, subject=ref_protein, outfmt=5, out=output_file)
    blastx_cline()

def find_amino_acid_substitutions(blast_output, ref_protein):
    substitutions = []
    blast_records = SearchIO.parse(blast_output, "blast-xml")
    
    for record in blast_records:
        for hit in record.hits:
            for hsp in hit.hsps:
                ref_seq = hsp.hit.seq
                query_seq = hsp.query.seq

                # Find amino acid substitutions
                for i in range(len(ref_seq)):
                    if ref_seq[i] != query_seq[i] and ref_seq[i] != '-':
                        substitutions.append({
                            'position': hsp.hit_start // 3 + (i // 3) + 1,
                            'ref_aa': ref_seq[i],
                            'query_aa': query_seq[i]
                        })
    return substitutions

def main():
    parser = argparse.ArgumentParser(description='Perform BLASTX and identify amino acid substitutions.')
    parser.add_argument('-a', '--assembly', required=True, help='Input assembly file')
    parser.add_argument('-r', '--reference', required=True, help='Reference protein sequence file')
    parser.add_argument('-o', '--output', required=True, help='BLAST output file')

    args = parser.parse_args()

    # Perform BLASTX
    perform_blastx(args.assembly, args.reference, args.output)

    # Find amino acid substitutions
    substitutions = find_amino_acid_substitutions(args.output, args.reference)

    print("Amino acid substitutions found:")
    for sub in substitutions:
        print(f"Position {sub['position']}: {sub['ref_aa']} -> {sub['query_aa']}")

if __name__ == "__main__":
    main()
