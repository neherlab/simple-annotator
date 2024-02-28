import pandas as pd
import argparse
from Bio import SeqIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequences", help="input sequences")
    parser.add_argument("--fields", nargs='+', help="fields encoded in fasta header")
    parser.add_argument("--organism",  help="organism")
    parser.add_argument("--header-separator", default='|', help="metadata")
    parser.add_argument("--strain-separator", default='/', help="metadata")
    parser.add_argument("--output-sequences", help="output sequences")
    parser.add_argument("--output-metadata", help="output sequences")
    args = parser.parse_args()

    metadata = []
    out_seqs = []
    for seq in SeqIO.parse(args.sequences, "fasta"):
        entries = seq.id.split(args.header_separator)
        isolate = entries[args.fields.index("isolate")]
        accession = entries[args.fields.index("accession")]
        collection_date = entries[args.fields.index("collection-date")]
        isolate_details = isolate.split("/")
        country = isolate_details[2]
        division = isolate_details[3].split("-")[0]
        seq_id = isolate_details[3]
        metadata.append({"Sequence_ID":seq_id, "isolate": isolate, "note": accession, "collection-date": collection_date, 
                         "country": f"{country}: {division}","organism": args.organism})


        seq.id = seq_id
        seq.description = " ".join(f"[{k}={v}]" for k,v in metadata[-1].items() if k != "Sequence_ID") + ' ' + isolate
        out_seqs.append(seq)

    metadata = pd.DataFrame(metadata)
    metadata.to_csv(args.output_metadata, index=False, sep="\t") 
    SeqIO.write(out_seqs, args.output_sequences, "fasta")
