input_data = "data/sequences.fasta"
template_file = "ncbi_author_template.sbt"
pathogens = ['rsv/a', 'rsv/b']
platform = 'linux64'  # or 'mac'

rule all:
    input:
        expand("{pathogen}/submission.sqn", pathogen=pathogens)


rule get_table2asn:
    output:
        table2asn = 'table2asn'
    params:
        bin = "linux64.table2asn.gz" if platform == 'linux64' else "mac.table2asn.gz"
    shell:
        '''
        curl https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn/{params.bin} -o {params.bin}
        mv {params.bin} {output.table2asn}.gz
        gunzip {output.table2asn}.gz
        chmod +x {output}
        '''

def get_organisms(w):
    return {'rsv/a': 'Human respiratory syncytial virus A', 'rsv/b': 'Human respiratory syncytial virus B'}.get(w.pathogen,'')

rule split:
    input:
        input_data
    output:
        expand("split_by_pathogen/nextstrain/{pathogen}/sequences.fasta", pathogen=pathogens)
    shell:
        '''
        nextclade3 sort --output-dir split_by_pathogen  {input}
        '''


rule parse:
    input:
        "split_by_pathogen/nextstrain/{pathogen}/sequences.fasta"
    output:
        sequences = "split_by_pathogen/nextstrain/{pathogen}/submission.fsa",
        metadata = "split_by_pathogen/nextstrain/{pathogen}/metadata.tsv"
    params:
        fields = "isolate accession collection-date",
        organism = get_organisms
    shell:
        '''
        python src/parse_metadata.py --sequences {input} --output-sequences {output.sequences} \
            --output-metadata {output.metadata} --organism {params.organism:q} --fields {params.fields}
        '''


rule annotate:
    input:
        sequences = "split_by_pathogen/nextstrain/{pathogen}/submission.fsa"
    output:
        annotation = "split_by_pathogen/nextstrain/{pathogen}/submission.tbl"
    params:
        out_dir =  "split_by_pathogen/nextstrain/{pathogen}/annotations",
        virus = lambda w: w.pathogen.replace('/', '-'),
    shell:
        '''
        mkdir -p {params.out_dir}
        python src/annotate.py \
            --sequences {input.sequences} --virus {params.virus} \
            --output-format tbl --output-dir {params.out_dir} 
        cat {params.out_dir}/*.tbl > {output.annotation}
        '''


rule assemble_submission:
    input:
        sequences = "split_by_pathogen/nextstrain/{pathogen}/submission.fsa",
        annotation = "split_by_pathogen/nextstrain/{pathogen}/submission.tbl",
        template = template_file,
        table2asn = "table2asn"
    output:
        submission_file = "{pathogen}/submission.sqn"
    params:
        indir = "split_by_pathogen/nextstrain/{pathogen}",
        prefix = lambda w: w.pathogen.split('/')[-1]
    shell:
        '''
        ./table2asn -indir {params.indir} -t {input.template}  -M n -Z -split-dr || true
        mkdir -p {wildcards.pathogen}
        mv {params.indir}/submission.* {wildcards.pathogen}/
        mv {params.prefix}.stats {wildcards.pathogen}/submission.stats
        '''



