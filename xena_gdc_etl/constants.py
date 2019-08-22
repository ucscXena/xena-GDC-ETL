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
    'gistic': {'gdc_type': 'GISTIC - focal score by gene'},
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
    'Xena_phenotype',
    'survival',
    'gistic',
    'star_counts',
    'methylation27',
    'methylation450',
]
CASES_FIELDS_EXPANDS = {
    "CPTAC-3": {
        "fields": [
            "case_id",
            "demographic.cause_of_death",
            "demographic.days_to_birth",
            "demographic.days_to_death",
            "demographic.demographic_id",
            "demographic.ethnicity",
            "demographic.gender",
            "demographic.race",
            "demographic.submitter_id",
            "demographic.vital_status",
            "demographic.year_of_birth",
            "demographic.year_of_death",
            "diagnoses.age_at_diagnosis",
            "diagnoses.ajcc_clinical_m",
            "diagnoses.ajcc_pathologic_m",
            "diagnoses.ajcc_pathologic_n",
            "diagnoses.ajcc_pathologic_stage",
            "diagnoses.ajcc_pathologic_t",
            "diagnoses.ajcc_staging_system_edition",
            "diagnoses.classification_of_tumor",
            "diagnoses.days_to_last_follow_up",
            "diagnoses.days_to_last_known_disease_status",
            "diagnoses.days_to_recurrence",
            "diagnoses.diagnosis_id",
            "diagnoses.last_known_disease_status",
            "diagnoses.lymph_nodes_positive",
            "diagnoses.morphology",
            "diagnoses.primary_diagnosis",
            "diagnoses.prior_malignancy",
            "diagnoses.progression_or_recurrence",
            "diagnoses.residual_disease",
            "diagnoses.site_of_resection_or_biopsy",
            "diagnoses.submitter_id",
            "diagnoses.tissue_or_organ_of_origin",
            "diagnoses.tumor_grade",
            "diagnoses.tumor_largest_dimension_diameter",
            "diagnoses.tumor_stage",
            "diagnoses.year_of_diagnosis",
            "disease_type",
            "exposures.alcohol_history",
            "exposures.alcohol_intensity",
            "exposures.bmi",
            "exposures.cigarettes_per_day",
            "exposures.exposure_id",
            "exposures.height",
            "exposures.pack_years_smoked",
            "exposures.submitter_id",
            "exposures.tobacco_smoking_onset_year",
            "exposures.tobacco_smoking_quit_year",
            "exposures.tobacco_smoking_status",
            "exposures.weight",
            "exposures.years_smoked",
            "primary_site",
            "samples.composition",
            "samples.current_weight",
            "samples.days_to_sample_procurement",
            "samples.freezing_method",
            "samples.initial_weight",
            "samples.is_ffpe",
            "samples.longest_dimension",
            "samples.oct_embedded",
            "samples.preservation_method",
            "samples.sample_id",
            "samples.sample_type",
            "samples.sample_type_id",
            "samples.submitter_id",
            "samples.time_between_excision_and_freezing",
            "samples.tissue_type",
            "samples.tumor_code",
            "samples.tumor_descriptor",
        ],
        "expand": [],
    },
    "GDC-PANCAN": {
        "fields": [
            "project.name",
            "project.project_id",
            "samples.is_ffpe",
            "samples.sample_id",
            "samples.sample_type",
            "samples.sample_type_id",
            "samples.tissue_type",
            "samples.tumor_code",
            "samples.tumor_code_id",
            "samples.tumor_descriptor",
            "tissue_source_site.name",
            "samples.submitter_id",
        ],
        "expand": [
            "demographic",
            "diagnoses",
            "exposures",
            "family_histories",
        ],
    }
}
LIST_FIELDS = {
    "FISH_test_component",
    "FISH_test_component_percentage_value",
    "ablation_performed_indicator",
    "ablation_type",
    "ablation_type_other",
    "additional_pharmaceutical_therapy",
    "additional_radiation_therapy",
    "additional_surgery_locoregional_procedure",
    "additional_surgery_metastatic_procedure",
    "additional_treatment_completion_success_outcome",
    "adjuvant_rad_therapy_prior_admin",
    "agent_total_dose_count",
    "anatomic_neoplasm_subdivision",
    "anatomic_treatment_site",
    "antireflux_treatment_type",
    "assessment_timepoint_category",
    "bcr_ablation_barcode",
    "bcr_ablation_uuid",
    "bcr_drug_barcode",
    "bcr_drug_uuid",
    "bcr_followup_barcode",
    "bcr_followup_uuid",
    "bcr_radiation_barcode",
    "bcr_radiation_uuid",
    "birth_control_pill_history_usage_category",
    "brachytherapy_administered_status",
    "brachytherapy_first_reference_point_administered_total_dose",
    "brachytherapy_method_other_specify_text",
    "brachytherapy_method_type",
    "breast_carcinoma_estrogen_receptor_status",
    "breast_carcinoma_immunohistochemistry_pos_cell_score",
    "breast_carcinoma_progesterone_receptor_status",
    "cancer_diagnosis_cancer_type_icd9_text_name",
    "chemotherapy_negation_radiation_therapy_concurrent_administered_text",
    "chemotherapy_negation_radiation_therapy_concurrent_not_administered_reason",  # noqa: E501
    "chemotherapy_regimen_type",
    "clinical_trail_drug_classification",
    "colorectal_cancer",
    "concurrent_chemotherapy_dose",
    "concurrent_chemotherapy_indicator",
    "course_number",
    "cytogenetic_abnormality",
    "day_of_form_completion",
    "days_to_ablation_performed",
    "days_to_additional_surgery_locoregional_procedure",
    "days_to_additional_surgery_metastatic_procedure",
    "days_to_brachytherapy_begin_occurrence",
    "days_to_brachytherapy_end_occurrence",
    "days_to_chemotherapy_end",
    "days_to_chemotherapy_start",
    "days_to_completion_of_curative_tx",
    "days_to_death",
    "days_to_drug_therapy_end",
    "days_to_drug_therapy_start",
    "days_to_embolization_performed",
    "days_to_first_biochemical_recurrence",
    "days_to_last_followup",
    "days_to_last_known_alive",
    "days_to_new_tumor_event_additional_surgery_procedure",
    "days_to_new_tumor_event_after_initial_treatment",
    "days_to_new_tumor_event_transplant",
    "days_to_performance_status_assessment",
    "days_to_radiation_therapy_end",
    "days_to_radiation_therapy_start",
    "days_to_second_biochemical_recurrence",
    "days_to_stem_cell_transplantation",
    "days_to_third_biochemical_recurrence",
    "death_cause_text",
    "diabetes",
    "diagnostic_ct_result_outcome",
    "diagnostic_mri_result",
    "diagnostic_mri_result_outcome",
    "discontiguous_lesion_count",
    "disease_after_curative_tx",
    "disease_code",
    "dose_frequency_text",
    "drug_name",
    "eastern_cancer_oncology_group",
    "embolization_performed_indicator",
    "embolization_therapy_drug_name",
    "embolization_therapy_drug_name_other",
    "embolization_type",
    "embolizing_agent_utilized",
    "embolizing_agent_utilized_other",
    "er_detection_method_text",
    "er_level_cell_percentage_category",
    "esophageal_tumor_involvement_site",
    "excess_adrenal_hormone_history_type",
    "external_beam_radiation_therapy_administered_paraaortic_region_lymph_node_dose",  # noqa: E501
    "external_beam_radiation_therapy_administered_status",
    "family_cancer_type_txt",
    "family_history_cancer_type",
    "family_history_cancer_type_other",
    "family_member_relationship_type",
    "fdg_or_ct_pet_performed_outcome",
    "first_degree_relative_history_thyroid_gland_carcinoma_diagnosis_relationship_type",  # noqa: E501
    "first_progression_histology_type",
    "first_progression_histology_type_other",
    "first_recurrence_biopsy_confirmed",
    "fluorescence_in_situ_hybridization_diagnostic_procedure_chromosome_17_signal_result_range",  # noqa: E501
    "followup_case_report_form_submission_reason",
    "followup_treatment_success",
    "genetic_abnormality_method",
    "genetic_abnormality_method_other",
    "genetic_abnormality_results",
    "genetic_abnormality_results_other",
    "genetic_abnormality_tested",
    "genetic_abnormality_tested_other",
    "genotyping_results_gene_mutation_not_reported_reason",
    "her2_and_centromere_17_positive_finding_other_measurement_scale_text",
    "her2_erbb_method_calculation_method_text",
    "her2_erbb_pos_finding_cell_percent_category",
    "her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text",  # noqa: E501
    "her2_immunohistochemistry_level_result",
    "her2_neu_and_centromere_17_copy_number_analysis_input_total_number_count",
    "her2_neu_and_centromere_17_copy_number_metastatic_breast_carcinoma_analysis_input_total_number_count",  # noqa: E501
    "her2_neu_breast_carcinoma_copy_analysis_input_total_number",
    "her2_neu_chromosone_17_signal_ratio_value",
    "her2_neu_metastatic_breast_carcinoma_copy_analysis_input_total_number",
    "hist_hepato_carc_fact",
    "histologic_disease_progression_present_indicator",
    "histologic_disease_progression_present_text",
    "histologic_disease_progression_present_type",
    "histological_percentage",
    "histological_type",
    "history_immunological_disease",
    "history_immunosuppresive_rx",
    "horm_ther",
    "human_papillomavirus_type",
    "hypertension",
    "i_131_first_administered_dose",
    "i_131_subsequent_administered_dose",
    "i_131_total_administered_dose",
    "i_131_total_administered_preparation_technique",
    "immunohistochemistry_positive_cell_score",
    "immunophenotype_cytochemistry_percent_positive",
    "immunophenotype_cytochemistry_testing_result",
    "immunophenotypic_analysis_method",
    "immunophenotypic_analysis_results",
    "immunophenotypic_analysis_tested",
    "karnofsky_performance_score",
    "lab_proc_her2_neu_immunohistochemistry_receptor_status",
    "lab_procedure_her2_neu_in_situ_hybrid_outcome_type",
    "lesions_count",
    "loss_expression_of_mismatch_repair_proteins_by_ihc_result",
    "lost_follow_up",
    "lymph_node_involvement_site",
    "lymph_node_location_positive_pathology_name",
    "lymph_node_preoperative_assessment_diagnostic_imaging_type",
    "measure_of_response",
    "metastatic_breast_carcinoma_erbb2_immunohistochemistry_level_result",
    "metastatic_breast_carcinoma_estrogen_receptor_detection_method_text",
    "metastatic_breast_carcinoma_estrogen_receptor_level_cell_percent_category",  # noqa: E501
    "metastatic_breast_carcinoma_estrogen_receptor_status",
    "metastatic_breast_carcinoma_fluorescence_in_situ_hybridization_diagnostic_proc_centromere_17_signal_result_range",  # noqa: E501
    "metastatic_breast_carcinoma_her2_erbb_method_calculation_method_text",
    "metastatic_breast_carcinoma_her2_erbb_pos_finding_cell_percent_category",
    "metastatic_breast_carcinoma_her2_erbb_pos_finding_fluorescence_in_situ_hybridization_calculation_method_text",  # noqa: E501
    "metastatic_breast_carcinoma_her2_neu_chromosone_17_signal_ratio_value",
    "metastatic_breast_carcinoma_immunohistochemistry_er_pos_cell_score",
    "metastatic_breast_carcinoma_immunohistochemistry_er_positive_finding_scale_type",  # noqa: E501
    "metastatic_breast_carcinoma_immunohistochemistry_pr_pos_cell_score",
    "metastatic_breast_carcinoma_immunohistochemistry_progesterone_receptor_positive_finding_scale_type",  # noqa: E501
    "metastatic_breast_carcinoma_lab_proc_her2_neu_immunohistochemistry_receptor_status",  # noqa: E501
    "metastatic_breast_carcinoma_lab_proc_her2_neu_in_situ_hybridization_outcome_type",  # noqa: E501
    "metastatic_breast_carcinoma_pos_finding_her2_erbb2_other_measure_scale_text",  # noqa: E501
    "metastatic_breast_carcinoma_pos_finding_other_scale_measurement_text",
    "metastatic_breast_carcinoma_pos_finding_progesterone_receptor_other_measure_scale_text",  # noqa: E501
    "metastatic_breast_carcinoma_progesterone_receptor_detection_method_text",
    "metastatic_breast_carcinoma_progesterone_receptor_level_cell_percent_category",  # noqa: E501
    "metastatic_breast_carcinoma_progesterone_receptor_status",
    "metastatic_neoplasm_confirmed_diagnosis_method_name",
    "metastatic_neoplasm_initial_diagnosis_anatomic_site",
    "metastatic_site",
    "metastatic_site_at_diagnosis",
    "method_of_curative_tx",
    "mitotane_therapy_adjuvant_setting",
    "molecular_analysis_abnormality_testing_percentage_value",
    "molecular_analysis_abnormality_testing_result",
    "molecular_test_result",
    "month_of_form_completion",
    "new_neoplasm_confirmed_diagnosis_method_name",
    "new_neoplasm_event_occurrence_anatomic_site",
    "new_neoplasm_event_post_initial_therapy_diagnosis_method_text",
    "new_neoplasm_event_post_initial_therapy_diagnosis_method_type",
    "new_neoplasm_event_type",
    "new_neoplasm_histology_other",
    "new_neoplasm_occurrence_anatomic_site_text",
    "new_non_melanoma_event_histologic_type_text",
    "new_primary_melanoma_anatomic_site",
    "new_primary_melanoma_anatomic_site_other_text",
    "new_tumor_cellular_differentiation",
    "new_tumor_event_ablation_embo_tx",
    "new_tumor_event_additional_surgery_procedure",
    "new_tumor_event_after_initial_treatment",
    "new_tumor_event_liver_transplant",
    "new_tumor_metastasis_anatomic_site",
    "new_tumor_metastasis_anatomic_site_other_text",
    "nte_anatomic_site_count",
    "nte_pathologic_tumor_depth",
    "nte_pathologic_tumor_length",
    "nte_pathologic_tumor_width",
    "nte_radiologic_tumor_depth",
    "nte_radiologic_tumor_length",
    "nte_radiologic_tumor_width",
    "nte_tumor_sizes",
    "number_cycles",
    "numfractions",
    "other_chemotherapy_agent_administration_specify",
    "pathologic_tumor_burden",
    "pathologic_tumor_depth",
    "pathologic_tumor_length",
    "pathologic_tumor_width",
    "patient_death_reason",
    "patient_personal_medical_history_thyroid_gland_disorder_name",
    "performance_status_scale_timing",
    "person_neoplasm_cancer_status",
    "pet_scan_results",
    "pgr_detection_method_text",
    "pharm_regimen",
    "pharm_regimen_other",
    "pos_finding_her2_erbb2_other_measurement_scale_text",
    "pos_finding_metastatic_breast_carcinoma_estrogen_receptor_other_measuremenet_scale_text",  # noqa: E501
    "pos_finding_progesterone_receptor_other_measurement_scale_text",
    "positive_finding_estrogen_receptor_other_measurement_scale_text",
    "post_op_ablation_embolization_tx",
    "post_orchi_lymph_node_dissection",
    "post_surgical_procedure_assessment_thyroid_gland_carcinoma_status",
    "postoperative_rx_tx",
    "pregnancies",
    "prescribed_dose",
    "prescribed_dose_units",
    "primary_anatomic_site_count",
    "primary_therapy_outcome_success",
    "prior_tamoxifen_administered_usage_category",
    "progesterone_receptor_level_cell_percent_category",
    "progression_determined_by",
    "progression_determined_by_notes",
    "project_code",
    "radiation_dosage",
    "radiation_therapy",
    "radiation_therapy_administered_dose_text",
    "radiation_therapy_administered_preparation_technique_text",
    "radiation_therapy_not_administered_reason",
    "radiation_therapy_not_administered_specify",
    "radiation_treatment_ongoing",
    "radiation_type",
    "radiation_type_notes",
    "radiologic_tumor_burden",
    "radiologic_tumor_depth",
    "radiologic_tumor_length",
    "radiologic_tumor_width",
    "radiosensitizing_agent_administered_indicator",
    "recurrence_second_surgery_neoplasm_surgical_procedure_name",
    "recurrence_second_surgery_neoplasm_surgical_procedure_name_other",
    "regimen_indication",
    "regimen_indication_notes",
    "regimen_number",
    "relation_testicular_cancer",
    "relative_cancer_type",
    "relative_family_cancer_hx_text",
    "residual_disease_new_tumor_event",
    "residual_disease_post_new_tumor_event_margin_status",
    "residual_tumor",
    "route_of_administration",
    "rt_administered_type",
    "rt_pelvis_administered_total_dose",
    "site_of_additional_surgery_new_tumor_event_mets",
    "smokeless_tobacco_use_age_at_quit",
    "smokeless_tobacco_use_age_at_start",
    "smokeless_tobacco_use_at_diag",
    "smokeless_tobacco_use_per_day",
    "smokeless_tobacco_use_regularly",
    "source_of_patient_death_reason",
    "stem_cell_transplantation",
    "stem_cell_transplantation_type",
    "subsequent_primary_melanoma_during_followup",
    "therapeutic_mitotane_levels_achieved",
    "therapeutic_procedure_new_neoplasm_required_additional_therapy_type",
    "therapy_ongoing",
    "therapy_type",
    "therapy_type_notes",
    "thyroid_gland_carcinoma_involvement_regional_lymph_node_type",
    "thyroid_gland_carcinoma_regional_lymph_node_involvement_anatomic_sites_text",  # noqa: E501
    "total_dose",
    "total_dose_units",
    "treatment_cycles_count",
    "tumor_level",
    "tumor_morphology_percentage",
    "tumor_multifocal",
    "tumor_progression_post_ht",
    "tumor_tissue_site",
    "tumor_tissue_site_other",
    "tx_on_clinical_trial",
    "type_ofprogression_of_disease_ht",
    "units",
    "viral_hepatitis_serology",
    "vital_status",
    "year_of_form_completion",
}
