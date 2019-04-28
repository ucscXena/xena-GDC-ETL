# Map xena_dtype to corresponding metadata template.
METADATA_TEMPLATE = {
    'htseq_counts': 'template.rna.meta.json',
    'htseq_fpkm': 'template.rna.meta.json',
    'htseq_fpkm-uq': 'template.rna.meta.json',
    'mirna': 'template.mirna.meta.json',
    'mirna_isoform': 'template.mirna_isoform.meta.json',
    'cnv': 'template.cnv.meta.json',
    'masked_cnv': 'template.cnv.meta.json',
    'muse_snv': 'template.snv.meta.json',
    'mutect2_snv': 'template.snv.meta.json',
    'somaticsniper_snv': 'template.snv.meta.json',
    'varscan2_snv': 'template.snv.meta.json',
    'GDC_phenotype': 'template.phenotype.meta.json',
    'survival': 'template.survival.meta.json',
    'methylation27': 'template.methylation.meta.json',
    'methylation450': 'template.methylation.meta.json'
}
# Jinja2 template variables for corresponding "xena_dtype".
METADATA_VARIABLES = {
    'htseq_counts': {'gdc_type': 'HTSeq - Counts'},
    'htseq_fpkm': {'gdc_type': 'HTSeq - FPKM',
                   'unit': 'fpkm'},
    'htseq_fpkm-uq': {'gdc_type': 'HTSeq - FPKM-UQ',
                      'unit': 'fpkm-uq'},
    'mirna': {'gdc_type': 'miRNA Expression Quantification'},
    'mirna_isoform': {'gdc_type': 'Isoform Expression Quantification'},
    'cnv': {'gdc_type': 'Copy Number Segment'},
    'masked_cnv': {'gdc_type': 'Masked Copy Number Segment'},
    'muse_snv': {'gdc_type': 'MuSE Variant Aggregation and Masking'},
    'mutect2_snv': {
            'gdc_type': 'MuTect2 Variant Aggregation and Masking'
        },
    'somaticsniper_snv': {
            'gdc_type': 'SomaticSniper Variant Aggregation and Masking'
        },
    'varscan2_snv': {
            'gdc_type': 'VarScan2 Variant Aggregation and Masking'
        },
    'methylation27': {'platform_num': '27'},
    'methylation450': {'platform_num': '450'}
}
valid_dtype = [
    'htseq_counts', 'htseq_fpkm', 'htseq_fpkm-uq', 'mirna',
    'masked_cnv', 'muse_snv', 'mutect2_snv',
    'somaticsniper_snv', 'varscan2_snv', 'GDC_phenotype',
    'survival', 'methylation27', 'methylation450'
]
