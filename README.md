# scRNA_Lep

R scripts for the analysis of Single-cell RNA sequencing data of lepromatous leprosy.

The QC step was included in file s1_QC.R;
The normalization and dimension reduction steps were included in file s2_dimension_reduction.R;
The clustering and differential gene expression between cases and controls steps were included in file s3_diff.R.

Files with data required for the analysis: gene_cell_exprs_table_symbol.xls, predicted_label.xls and rRNA.gene.list

There were two study sets, 'Discovery_cohort_PBMC" and "Discovery_cohort_skin". "predicted_label.xls" and "rRNA.gene.list" files were in compressed files input_files_for_Discovery_cohort_PBMC.tar.gz and input_files_for_Discovery_cohort_skin.tar.gz.
The compressed version of "gene_cell_exprs_table_symbol.xls" files exceeds 25MB and cannot be uploaded here. You can ask the author for these files.
 
The  parameter to run the Discovery_cohort_PBMC set were:
 --count gene_cell_exprs_table_symbol.xls --tax human --double_file predicted_label.xls --sample_list P4PBMCJIAN,P5PBMCJIAN,P6PBMCJIAN,P8PBMCJIAN,P9PBMCJIAN,P11PBMCJIAN,P12PBMCJIAN,C1PBMCJIAN,C2PBMCJIAN,ZQPBMCblank,control_FXK_PBMC,control_MMY_PBMC,control_ZXJ_PBMC --outdir Discovery_cohort_PBMC  --resolution 0.2 --pcs_number 30 --group_list case,case,case,case,case,case,case,control,control,control,control,control,control --compare_list case_vs_control --redcell_cutoff 10 --mito_cutoff 10 --nFeature_lowlimit 200 --nFeature_highlimit 10000 --onlypos FALSE --batch TRUE

The  parameter to run the Discovery_cohort_skin set were:
--count gene_cell_exprs_table_symbol.xls --tax human --double_file predicted_label.xls --sample_list XWHskinNormal5,SHPSKINnormal5,WMskinNormal25,P5skin,P8skin,P11skin --outdir Discovery_cohort_skin  --resolution 0.2 --pcs_number 30  --batch TRUE --group_list control,control,control,case,case,case --compare_list case_vs_control --redcell_cutoff 10 --mito_cutoff 10 --nFeature_lowlimit 200 --nFeature_highlimit 10000 --onlypos FALSE --batch TRUE

