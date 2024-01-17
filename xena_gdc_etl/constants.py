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
    'star_counts': 'template.rna.meta.json',
    'star_tpm' : 'template.rna.meta.json',
    'star_fpkm' : 'template.rna.meta.json',
    'star_fpkm-uq' : 'template.rna.meta.json',
    'mirna': 'template.mirna.meta.json',
    'mirna_isoform': 'template.mirna_isoform.meta.json',
    'segment_cnv_ascat-ngs': 'template.cnv.meta.json', 
    'masked_cnv': 'template.cnv.meta.json',
    'somaticmutation_snv': 'template.snv.meta.json',
    'GDC_phenotype': 'template.phenotype.meta.json',
    'survival': 'template.survival.meta.json',
    'gene-level_ascat-ngs': 'template.ascat-ngs.meta.json',
    'methylation_epic': 'template.methylation.meta.json',
    'methylation27': 'template.methylation.meta.json',
    'methylation450': 'template.methylation.meta.json',
}

# Jinja2 template variables for corresponding "xena_dtype".
METADATA_VARIABLES = {
    'star_counts': {'gdc_type': 'STAR - Counts'},
    'star_tpm' : {'gdc_type' : 'STAR - TPM', 'unit' : 'tpm'}, 
    'star_fpkm' : {'gdc_type' : 'STAR - FPKM', 'unit': 'fpkm'},
    'star_fpkm-uq' : {'gdc_type' : 'STAR - FPKM-UQ', 'unit' : 'fpkm-uq'},
    'mirna': {'gdc_type': 'miRNA Expression Quantification'},
    'mirna_isoform': {'gdc_type': 'Isoform Expression Quantification'},
    'segment_cnv_ascat-ngs': {'gdc_type': 'Copy Number Segment'},
    'masked_cnv': {'gdc_type': 'Masked Copy Number Segment'},
    'somaticmutation_snv': {'gdc_type': 'Ensemble Somatic Variant'},
    'methylation_epic': {'platform_num': 'epic'},
    'methylation27': {'platform_num': '27'},
    'methylation450': {'platform_num': '450'},
}

valid_dtype = [
    'star_counts',
    'star_tpm',
    'star_fpkm',
    'star_fpkm-uq',
    'mirna',
    'mirna_isoform',
    'segment_cnv_ascat-ngs', 
    'masked_cnv',
    'somaticmutation_snv',
    'gene-level_ascat-ngs',
    'methylation_epic',
    'methylation27',
    'methylation450',
    'clinical',
    'survival',
]