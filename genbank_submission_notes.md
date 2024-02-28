# Submission to genbank

## Create a template file with author and submitter information

Navigate to [submit.ncbi.nlm.nih.gov/genbank/template/submission/](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) to create a file that contains the author and contact information associated with the submission. 
This file will have the file ending `.sbt`. Download this file and place it into this directory.

## Workflow

To facilitate preparing a submission, this directory contains a workflow that runs the necessary steps.
The workflow is currently written around RSV. 

Open the `Snakefile` in this directory and adjust the fields in the top section. The workflow expects a fasta file with metadata coded into the fasta header with fields `strain|accession|collection_date` where geographic location are derived from the strain name expected to be of the form
`virus/type/country/state-uniqueID/year`. 

The workflow will 
 - split the input into different virus types (RSV-A vs RSV-B)
 - parse the metadata and produce a fasta file with ending `fsa` and metadata coded into the headline in a gb compliant manner
 - annotated the genome
 - assembles the submission files using `table2asn`

 

## Create metadata
Metadata can be provided as a tab separated file with file ending `.src`.
The first column has to be named `Sequence_ID` and contain the fasta identifier of the corresponding sequence.
See [here for an example](https://www.ncbi.nlm.nih.gov/WebSub/html/help/sample_files/source-table-sample.txt).
A list of common field names are available on the web: https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html

Common ones include 
 - Collection_date
 - Country
 - Isolation_source
 - Lat_Long
 - Host
 - Isolate
 - Strain


`Country` can include additional information on location after a colon as in "USA:CA".

These fields can also be included in the fasta header of each sequence with [fields described here](https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/).


## run the annotation

```
python simple-annotator/src/annotate.py --sequences consensus_sequences/all_seqs.fasta --reference NC_001802.1 --output-format tbl --output-dir annotations
```




