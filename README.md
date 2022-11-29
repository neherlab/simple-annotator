## Simple viral annotator

This simple script aligns viral sequences to a reference genome that is fetched from genbank and lifts the annotation from the reference to the query. This is only expected to work for sequences that are very similar to the reference sequence.

Usage:
```
./src/annotate.py --sequences example_data/example_rsv-b.fasta  --virus rsv-b --output-dir test_output --output-format tbl
```
### Available output formats

 * genbank
 * gff3
 * tbl


### Current Limitations

 * no support for partial features
 * no support for complex features (slippage, splicing, etc)

