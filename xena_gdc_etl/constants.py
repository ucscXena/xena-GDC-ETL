GDC_RELEASE_URL = (
    'https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/'
)

# Map GDC project_id to Xena specific cohort name.
GDC_XENA_COHORT = {
    'TCGA-BRCA': 'GDC TCGA Breast Cancer (BRCA)',
    'TCGA-LUAD': 'GDC TCGA Lung Adenocarcinoma (LUAD)',
    'TCGA-UCEC': 'GDC TCGA Endometrioid Cancer (UCEC)',
    'TCGA-LGG': 'GDC TCGA Lower Grade Glioma (LGG)',
    'TCGA-HNSC': 'GDC TCGA Head and Neck Cancer (HNSC)',
    'TCGA-PRAD': 'GDC TCGA Prostate Cancer (PRAD)',
    'TCGA-LUSC': 'GDC TCGA Lung Squamous Cell Carcinoma (LUSC)',
    'TCGA-THCA': 'GDC TCGA Thyroid Cancer (THCA)',
    'TCGA-SKCM': 'GDC TCGA Melanoma (SKCM)',
    'TCGA-OV': 'GDC TCGA Ovarian Cancer (OV)',
    'TCGA-STAD': 'GDC TCGA Stomach Cancer (STAD)',
    'TCGA-COAD': 'GDC TCGA Colon Cancer (COAD)',
    'TCGA-BLCA': 'GDC TCGA Bladder Cancer (BLCA)',
    'TCGA-GBM': 'GDC TCGA Glioblastoma (GBM)',
    'TCGA-LIHC': 'GDC TCGA Liver Cancer (LIHC)',
    'TCGA-KIRC': 'GDC TCGA Kidney Clear Cell Carcinoma (KIRC)',
    'TCGA-CESC': 'GDC TCGA Cervical Cancer (CESC)',
    'TCGA-KIRP': 'GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)',
    'TCGA-SARC': 'GDC TCGA Sarcoma (SARC)',
    'TCGA-ESCA': 'GDC TCGA Esophageal Cancer (ESCA)',
    'TCGA-PAAD': 'GDC TCGA Pancreatic Cancer (PAAD)',
    'TCGA-PCPG': 'GDC TCGA Pheochromocytoma & Paraganglioma (PCPG)',
    'TCGA-READ': 'GDC TCGA Rectal Cancer (READ)',
    'TCGA-TGCT': 'GDC TCGA Testicular Cancer (TGCT)',
    'TCGA-LAML': 'GDC TCGA Acute Myeloid Leukemia (LAML)',
    'TCGA-THYM': 'GDC TCGA Thymoma (THYM)',
    'TCGA-ACC': 'GDC TCGA Adrenocortical Cancer (ACC)',
    'TCGA-MESO': 'GDC TCGA Mesothelioma (MESO)',
    'TCGA-UVM': 'GDC TCGA Ocular melanomas (UVM)',
    'TCGA-KICH': 'GDC TCGA Kidney Chromophobe (KICH)',
    'TCGA-UCS': 'GDC TCGA Uterine Carcinosarcoma (UCS)',
    'TCGA-CHOL': 'GDC TCGA Bile Duct Cancer (CHOL)',
    'TCGA-DLBC': 'GDC TCGA Large B-cell Lymphoma (DLBC)',
}

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
    'gistic': 'template.gistic.meta.json',
    'star_counts': 'template.rna.meta.json',
    'methylation27': 'template.methylation.meta.json',
    'methylation450': 'template.methylation.meta.json',
}

# Jinja2 template variables for corresponding "xena_dtype".
METADATA_VARIABLES = {
    'htseq_counts': {'gdc_type': 'HTSeq - Counts'},
    'htseq_fpkm': {'gdc_type': 'HTSeq - FPKM', 'unit': 'fpkm'},
    'htseq_fpkm-uq': {'gdc_type': 'HTSeq - FPKM-UQ', 'unit': 'fpkm-uq'},
    'mirna': {'gdc_type': 'miRNA Expression Quantification'},
    'mirna_isoform': {'gdc_type': 'Isoform Expression Quantification'},
    'cnv': {'gdc_type': 'Copy Number Segment'},
    'masked_cnv': {'gdc_type': 'Masked Copy Number Segment'},
    'muse_snv': {'gdc_type': 'MuSE Variant Aggregation and Masking'},
    'mutect2_snv': {'gdc_type': 'MuTect2 Variant Aggregation and Masking'},
    'somaticsniper_snv': {
        'gdc_type': 'SomaticSniper Variant Aggregation and Masking'
    },
    'varscan2_snv': {'gdc_type': 'VarScan2 Variant Aggregation and Masking'},
    'gistic': {
        'gdc_type': 'GISTIC - focal score by gene',
        'notes': 'The GDC GISTIC copy number dataset is derived from focal '
                 'copy number estimates. Larger chromosomal-level deletions '
                 'may not be not captured in this dataset.',
    },
    'star_counts': {'gdc_type': 'STAR - Counts'},
    'methylation27': {'platform_num': '27'},
    'methylation450': {'platform_num': '450'},
}
valid_dtype = [
    'htseq_counts',
    'htseq_fpkm',
    'htseq_fpkm-uq',
    'mirna',
    'cnv',
    'masked_cnv',
    'muse_snv',
    'mutect2_snv',
    'somaticsniper_snv',
    'varscan2_snv',
    'GDC_phenotype',
    'survival',
    'gistic',
    'star_counts',
    'methylation27',
    'methylation450',
]
