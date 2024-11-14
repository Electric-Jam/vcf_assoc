vcf_file="/home/eup009/dnanexus/500k_analysis/res_call/call_output_Case_White_batch1/normalized_merged_variants.vcf.gz"
pheno_file="/home/eup009/dnanexus/500k_analysis/res_call/test_assoc_pheno_cov.txt"
sample_id_col="IID"
pheno_col="L929"
covariate_cols="Age Sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10"
chrom='chr7'
start_bp=74774011
end_bp=74789314
out_dir="/home/eup009/vcf_association/test_files"

# Run Fisher's exact test for association testing
python3 association_test.py \
    --vcf ${vcf_file} \
    --genotype_field GT \
    --chromosome ${chrom} \
    --start_bp ${start_bp} \
    --end_bp ${end_bp} \
    --phenotype ${pheno_file} \
    --sample_id_col ${sample_id_col} \
    --phenotype_col ${pheno_col} \
    --covariate_cols ${covariate_cols} \
    --method fisher \
    --output ${out_dir}/results_fisher_GT.txt

python3 association_test.py \
    --vcf ${vcf_file} \
    --genotype_field PGT \
    --chromosome ${chrom} \
    --start_bp ${start_bp} \
    --end_bp ${end_bp} \
    --phenotype ${pheno_file} \
    --sample_id_col ${sample_id_col} \
    --phenotype_col ${pheno_col} \
    --covariate_cols ${covariate_cols} \
    --method fisher \
    --output ${out_dir}/results_fisher_PGT.txt

python3 association_test.py \
    --vcf ${vcf_file} \
    --genotype_field both \
    --chromosome ${chrom} \
    --start_bp ${start_bp} \
    --end_bp ${end_bp} \
    --phenotype ${pheno_file} \
    --sample_id_col ${sample_id_col} \
    --phenotype_col ${pheno_col} \
    --covariate_cols ${covariate_cols} \
    --method fisher \
    --output ${out_dir}/results_fisher_both.txt



# Run logistic regression for association testing
python3 association_test.py \
    --vcf ${vcf_file} \
    --genotype_field GT \
    --chromosome ${chrom} \
    --start_bp ${start_bp} \
    --end_bp ${end_bp} \
    --phenotype ${pheno_file} \
    --sample_id_col ${sample_id_col} \
    --phenotype_col ${pheno_col} \
    --covariate_cols ${covariate_cols} \
    --method logistic \
    --output ${out_dir}/results_logistic_GT.txt

python3 association_test.py \
    --vcf ${vcf_file} \
    --genotype_field PGT \
    --chromosome ${chrom} \
    --start_bp ${start_bp} \
    --end_bp ${end_bp} \
    --phenotype ${pheno_file} \
    --sample_id_col ${sample_id_col} \
    --phenotype_col ${pheno_col} \
    --covariate_cols ${covariate_cols} \
    --method logistic \
    --output ${out_dir}/results_logistic_PGT.txt

python3 association_test.py \
    --vcf ${vcf_file} \
    --genotype_field both \
    --chromosome ${chrom} \
    --start_bp ${start_bp} \
    --end_bp ${end_bp} \
    --phenotype ${pheno_file} \
    --sample_id_col ${sample_id_col} \
    --phenotype_col ${pheno_col} \
    --covariate_cols ${covariate_cols} \
    --method logistic \
    --output ${out_dir}/results_logistic_both.txt
