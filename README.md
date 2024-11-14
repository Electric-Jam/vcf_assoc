## Usage

### Command-Line Arguments

- `--vcf`: Path to the input VCF file (compressed with gzip, `.vcf.gz`).
- `--genotype_field`: Genotype field to analyze (`GT`, `PGT`, or `both`). Default is `GT`.
- `--chromosome`: (Optional) Chromosome to filter (e.g., `chr1`).
- `--start_bp`: (Optional) Start base pair position for filtering.
- `--end_bp`: (Optional) End base pair position for filtering.
- `--phenotype`: Path to the phenotype file (tab-delimited text file).
- `--sample_id_col`: Column name for sample IDs in the phenotype file.
- `--phenotype_col`: Column name for the binary phenotype in the phenotype file.
- `--covariate_cols`: (Optional) Column names for covariates in the phenotype file (space-separated).
- `--method`: Association test method (`fisher` or `logistic`). Default is `fisher`.
- `--output`: Output file name for the results (e.g., `results.txt`).

### Example Command

```bash
python association_test.py \
    --vcf input.vcf.gz \
    --genotype_field GT \
    --phenotype phenotype.txt \
    --sample_id_col SampleID \
    --phenotype_col CaseControl \
    --covariate_cols Age Sex \
    --method fisher \
    --output results.txt
