#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process blastx {
    input:
    path assembly
    path reference
    path output

    output:
    path output

    script:
    """
    blastx -query ${assembly} -subject ${reference} -outfmt 5 -out ${output}
    """
}

process find_amino_acid_substitutions {
    input:
    path blast_output
    path reference

    output:
    file 'substitutions.txt'

    script:
    """
    python -c '
import Bio.SearchIO as SearchIO

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

substitutions = find_amino_acid_substitutions("${blast_output}", "${reference}")

with open("substitutions.txt", "w") as f:
    for sub in substitutions:
        f.write(f"Position {sub['position']}: {sub['ref_aa']} -> {sub['query_aa']}\n")
    '
    """
}

workflow {
    assembly = file(params.assembly)
    reference = file(params.reference)
    output = file('blast_output.xml')

    blastx(assembly, reference, output)

    find_amino_acid_substitutions(output, reference)
}
