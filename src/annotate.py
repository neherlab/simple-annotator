#!/bin/env python3
import os
import copy
from Bio import SeqIO, AlignIO
import numpy as np

from defaults import reference_sequences

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Annotate sequences using a genbank reference')
    parser.add_argument('--reference', help='Genbank accession of reference sequence (will be fetched from genbank)')
    parser.add_argument('--virus', help='Virus name, default reference for the virus will be fetched from genbank')
    parser.add_argument('--sequences', required=True, help='Fasta file of sequences to annotate')
    parser.add_argument('--output-format', choices=['genbank', 'gff3', 'tbl'], default='genbank', type=str, help='Output format')
    parser.add_argument('--output-dir', required=True, type=str, help='Output directory')
    return parser.parse_args()

def get_reference_sequence(accession):
    from Bio import Entrez
    Entrez.email = "hello@nextstrain.org"
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    print(f"Fetching reference sequence {accession} from genbank")
    return SeqIO.read(handle, "genbank")

def get_coordinate_map(reference, record):
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmp_dir:
        SeqIO.write([reference, record], f"{tmp_dir}/sequences.fasta", "fasta")

        os.system(f"mafft --quiet {tmp_dir}/sequences.fasta > {tmp_dir}/aligned.fasta")
        ref_aligned, qry_aligned = AlignIO.read(f"{tmp_dir}/aligned.fasta", "fasta")

    aln_to_qry = np.cumsum(np.array(qry_aligned.seq) !='-')
    ref_to_qry = aln_to_qry[np.array(ref_aligned.seq) !='-'] - 1

    positions_mapping_to_end_of_query = np.where(ref_to_qry == len(record.seq)-1)[0]
    if len(positions_mapping_to_end_of_query) > 1:
        ref_to_qry[positions_mapping_to_end_of_query[1]:] = len(record.seq)

    return ref_to_qry


def annotate_sequence(reference, record):
    from Bio import SeqFeature

    coordinate_map = get_coordinate_map(reference, record)

    new_features = []
    # annotate features
    for feature in reference.features:
        new_feat = copy.deepcopy(feature)
        if new_feat.type == 'source':
            new_feat.location = SeqFeature.FeatureLocation(0, len(record.seq))
            new_feat.qualifiers = {}
            for v in ['organism', 'mol_type']:
                if v in feature.qualifiers:
                    new_feat.qualifiers[v] = feature.qualifiers[v]
        else:
            new_feat.location = SeqFeature.FeatureLocation(
                    max(0,int(coordinate_map[feature.location.start])),
                    min(int(coordinate_map[feature.location.end-1]+1), len(record.seq)))

            # if coordinate_map[feature.location.start]==-1:
            #     continue
            # if coordinate_map[feature.location.end]==len(record.seq):
            #     continue

            if new_feat.type == 'CDS':
                new_feat.qualifiers['translation'] = str(new_feat.extract(record.seq).translate())

            if 'locus_tag' in new_feat.qualifiers:
                new_feat.qualifiers.pop('locus_tag')

        if new_feat.location.end - new_feat.location.start > 0.3*(feature.location.end - feature.location.start):
            new_features.append(new_feat)

    record.features = new_features
    record.annotations = {k:v for k,v in reference.annotations.items() if k in ['organism', 'taxonomy', 'molecule_type']}

    return record

def tbl_entry(start, end, key, qualifiers):
    return f"{start}\t{end}\t{key}\t\t\n" + "\n".join([f"\t\t\t{q[0]}\t{q[1]}" for q in qualifiers]) + "\n"

def feature_to_tbl(feat):
    if feat.type == 'source':
        return tbl_entry(feat.location.start, feat.location.end, 'source',
                         [('mol_type', feat.qualifiers['mol_type'][0])])
    elif feat.type == 'CDS':
        return tbl_entry(feat.location.start, feat.location.end, 'CDS',
                        [('product', feat.qualifiers['product'][0]),
                         ('protein_id', feat.qualifiers['protein_id'][0])])
    elif feat.type == 'gene':
        return tbl_entry(feat.location.start, feat.location.end, 'gene',
                        [('gene', feat.qualifiers['gene'][0])])
    else:
        return ''

def write_tbl(record, filename):
    with open(filename, 'w') as f:
        f.write(f">Feature {record.id}\n")
        for feat in record.features:
            f.write(feature_to_tbl(feat))

if __name__=="__main__":

    args = parse_args()
    if not os.path.isdir(args.output_dir):
        print("Output directory does not exist!")
        exit()

    if args.virus is None and args.reference is None:
        print("Please specify a reference sequence accession or a virus name.")
        exit()

    if args.virus is not None and args.virus not in reference_sequences:
        print("Virus name not found in reference sequences.")
        exit()

    ref_id = reference_sequences[args.virus] if args.virus else args.reference

    reference_sequence = get_reference_sequence(ref_id)

    annotated_sequences = []

    for record in SeqIO.parse(args.sequences, "fasta"):
        annotated_sequences.append(annotate_sequence(reference_sequence, record))
        print(f"Annotated {record.id}")
        if args.output_format == 'genbank':
            SeqIO.write(record, f"{args.output_dir}/{record.id}.gb", "genbank")
        elif args.output_format == 'gff3':
            from BCBio.GFF import GFF3Writer
            gffwriter = GFF3Writer()
            with open(f"{args.output_dir}/{record.id}.gff3", "w") as handle:
                gffwriter.write([record], handle)
        elif args.output_format == 'tbl':
            write_tbl(record, f"{args.output_dir}/{record.id}.tbl")

