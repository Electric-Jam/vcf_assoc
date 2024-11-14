#!/usr/bin/env python3

import argparse
import gzip
import pandas as pd
import numpy as np
from scipy.stats import chisquare, chi2_contingency
import statsmodels.formula.api as smf
import warnings
import logging
import sys
import os
from tqdm import tqdm
import itertools

warnings.filterwarnings('ignore')

def parse_args():
    parser = argparse.ArgumentParser(description='Association test for binary phenotypes using VCF.')
    parser.add_argument('--vcf', required=True, help='Input VCF.gz file')
    parser.add_argument('--genotype_field', choices=['GT', 'PGT', 'both'], default='GT', help='Genotype field to analyze')
    parser.add_argument('--chromosome', help='Chromosome to filter')
    parser.add_argument('--start_bp', type=int, help='Start base pair position')
    parser.add_argument('--end_bp', type=int, help='End base pair position')
    parser.add_argument('--phenotype', required=True, help='Phenotype file (tab-delimited)')
    parser.add_argument('--sample_id_col', required=True, help='Column name for sample ID')
    parser.add_argument('--phenotype_col', required=True, help='Column name for phenotype')
    parser.add_argument('--covariate_cols', nargs='*', help='Column names for covariates')
    parser.add_argument('--method', choices=['fisher', 'logistic'], default='fisher', help='Association test method')
    parser.add_argument('--output', required=True, help='Output file name')
    return parser.parse_args()

def read_phenotype(phenotype_file, sample_id_col, phenotype_col, covariate_cols):
    pheno_df = pd.read_csv(phenotype_file, sep='\t')
    required_cols = [sample_id_col, phenotype_col]
    if covariate_cols:
        required_cols += covariate_cols
    pheno_df = pheno_df[required_cols]
    pheno_df.dropna(subset=[phenotype_col], inplace=True)
    pheno_df[phenotype_col] = pheno_df[phenotype_col].astype(int)
    # Normalize sample IDs
    pheno_df[sample_id_col] = pheno_df[sample_id_col].astype(str).str.strip().str.lower()
    return pheno_df

def perform_hwe(genotypes):
    # Remove missing genotypes
    genotypes = genotypes.dropna()
    if len(genotypes) == 0:
        return np.nan, np.nan

    # Ensure genotypes are integers
    genotypes = genotypes.astype(int)

    # Count genotypes
    genotype_counts = genotypes.value_counts().reindex([0, 1, 2], fill_value=0)
    obs_counts = genotype_counts.values

    # Calculate allele frequencies
    n_individuals = obs_counts.sum()
    n_alleles = 2 * n_individuals
    n_alt_alleles = obs_counts[1] + 2 * obs_counts[2]
    p = n_alt_alleles / n_alleles
    q = 1 - p

    # Expected genotype frequencies under HWE
    expected_freqs = np.array([q ** 2, 2 * p * q, p ** 2])
    expected_counts = expected_freqs * n_individuals

    # Check for zero expected counts
    if np.any(expected_counts < 5):
        # For small expected counts, chi-squared test is not valid
        return np.nan, np.nan

    # Perform chi-squared goodness-of-fit test
    chi2_stat, p_value = chisquare(f_obs=obs_counts, f_exp=expected_counts)

    return p_value, chi2_stat

def serialize_contingency_table(table):
    rows = []
    for genotype in table.index:
        counts = table.loc[genotype]
        counts_str = ','.join([f"{col}:{int(counts[col])}" for col in table.columns])
        row_str = f"{genotype}|{counts_str}"
        rows.append(row_str)
    return ';'.join(rows)

def main():
    args = parse_args()

    # Set up logging to file
    log_filename = os.path.splitext(args.output)[0] + '.log'
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    # Remove any existing handlers
    if logger.hasHandlers():
        logger.handlers.clear()
    # Create file handler
    fh = logging.FileHandler(log_filename, mode='w')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Log the given inputs
    logger.info('Starting association test script with the following parameters:')
    logger.info(f'VCF file: {args.vcf}')
    logger.info(f'Genotype field: {args.genotype_field}')
    logger.info(f'Chromosome: {args.chromosome}')
    logger.info(f'Start position: {args.start_bp}')
    logger.info(f'End position: {args.end_bp}')
    logger.info(f'Phenotype file: {args.phenotype}')
    logger.info(f'Sample ID column: {args.sample_id_col}')
    logger.info(f'Phenotype column: {args.phenotype_col}')
    logger.info(f'Covariate columns: {args.covariate_cols}')
    logger.info(f'Association test method: {args.method}')
    logger.info(f'Output file: {args.output}')
    logger.info('Reading phenotype data...')

    # Read phenotype data
    pheno_df = read_phenotype(args.phenotype, args.sample_id_col, args.phenotype_col, args.covariate_cols)
    logger.info(f'Phenotype data contains {len(pheno_df)} samples.')

    # Initialize counters
    total_variants = 0
    analyzed_variants = 0

    # Open VCF file
    with gzip.open(args.vcf, 'rt') as vcf_file, open(args.output, 'w') as output_file:
        # Write header
        header = ['Chromosome', 'Position', 'Variant_ID', 'Ref', 'Alt', 'MAC_Case', 'MAC_Control',
                  'HWE_p_value', 'Association_p_value', 'Contingency_Table']
        output_file.write('\t'.join(header) + '\n')

        # Read VCF header
        for line in vcf_file:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                vcf_header = line.strip().split('\t')
                samples = vcf_header[9:]
                # Normalize sample IDs in VCF
                samples = [s.strip().lower() for s in samples]
                logger.info('First few sample IDs from VCF:')
                logger.info(samples[:5])
                logger.info('First few sample IDs from phenotype data:')
                logger.info(pheno_df[args.sample_id_col].head().tolist())

                # Build mapping from sample IDs to their indices in the VCF
                sample_id_to_vcf_index = {sample_id: idx for idx, sample_id in enumerate(samples)}
                # Keep only samples present in both phenotype data and VCF
                pheno_df = pheno_df[pheno_df[args.sample_id_col].isin(samples)]
                logger.info(f'Found {len(pheno_df)} overlapping samples between phenotype data and VCF.')
                if len(pheno_df) == 0:
                    logger.error('No overlapping samples found between phenotype data and VCF.')
                    sys.exit(1)
                sample_ids_in_both = pheno_df[args.sample_id_col].tolist()
                # Map sample IDs to their indices in the VCF
                sample_indices = [sample_id_to_vcf_index[sid] for sid in sample_ids_in_both]
                # Reset index of pheno_df for alignment
                pheno_df = pheno_df.reset_index(drop=True)
                break  # Exit header processing
        else:
            logger.error('VCF file does not contain header line.')
            sys.exit(1)

        # Process VCF variants with progress bar
        vcf_file.seek(0)
        # Skip header lines
        for line in vcf_file:
            if line.startswith('#'):
                continue
            else:
                break  # First variant line found
        vcf_file_iter = itertools.chain([line], vcf_file)  # Include the first variant line

        for line in tqdm(vcf_file_iter, desc='Processing variants'):
            total_variants += 1

            fields = line.strip().split('\t')
            chrom, pos, var_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

            # Filter by chromosome and position
            if args.chromosome and chrom != args.chromosome:
                continue
            if args.start_bp and int(pos) < args.start_bp:
                continue
            if args.end_bp and int(pos) > args.end_bp:
                continue

            format_fields = fields[8].split(':')
            genotype_field_indices = [i for i, f in enumerate(format_fields) if f in ['GT', 'PGT']]
            if not genotype_field_indices:
                logger.warning(f'Variant {var_id} at position {pos} skipped: No GT or PGT field.')
                continue

            sample_genotypes = []
            for sid, idx in zip(sample_ids_in_both, sample_indices):
                sample_field = fields[9 + idx].split(':')
                genotype_data = None
                for gt_idx in genotype_field_indices:
                    if args.genotype_field in ['both', format_fields[gt_idx]]:
                        genotype = sample_field[gt_idx]
                        if genotype != '.':
                            alleles = genotype.replace('|', '/').split('/')
                            alleles = [int(a) if a != '.' else np.nan for a in alleles]
                            genotype_data = sum(alleles)
                            break
                if genotype_data is None:
                    genotype_data = np.nan
                sample_genotypes.append(genotype_data)

            genotypes = pd.Series(sample_genotypes)
            combined_df = pheno_df.copy()
            combined_df['Genotype'] = genotypes
            combined_df.dropna(subset=['Genotype'], inplace=True)

            # Skip if no variation
            if combined_df['Genotype'].nunique() == 1:
                logger.warning(f'Variant {var_id} at position {pos} skipped: No variation in genotypes.')
                continue

            # Calculate MAC for cases and controls
            cases = combined_df[combined_df[args.phenotype_col] == 1]
            controls = combined_df[combined_df[args.phenotype_col] == 0]
            mac_case = cases['Genotype'].sum()
            mac_control = controls['Genotype'].sum()

            # Perform HWE test
            hwe_p_value, hwe_chi2 = perform_hwe(combined_df['Genotype'])
            if np.isnan(hwe_p_value):
                logger.info(f'HWE test not performed for variant {var_id} at position {pos} due to insufficient data.')

            # Association test
            contingency_table_str = ''
            if args.method == 'fisher':
                # Build contingency table
                table = pd.crosstab(combined_df['Genotype'], combined_df[args.phenotype_col])
                if table.size == 0:
                    logger.warning(f'Variant {var_id} at position {pos} skipped: Empty contingency table.')
                    continue
                # Serialize contingency table
                contingency_table_str = serialize_contingency_table(table)
                # Use Chi-squared test for independence
                chi2, p_value, dof, expected = chi2_contingency(table)
            else:
                # Logistic regression
                formula = f"{args.phenotype_col} ~ Genotype"
                if args.covariate_cols:
                    covariates = ' + '.join(args.covariate_cols)
                    formula += ' + ' + covariates
                try:
                    model = smf.logit(formula, data=combined_df).fit(disp=0)
                    p_value = model.pvalues['Genotype']
                except Exception as e:
                    logger.warning(f'Variant {var_id} at position {pos} skipped: Logistic regression error: {e}')
                    continue  # Skip variants where logistic regression fails

            # Write results
            output_line = [chrom, pos, var_id, ref, alt, str(mac_case), str(mac_control),
                           f"{hwe_p_value:.4g}" if not np.isnan(hwe_p_value) else 'NA',
                           f"{p_value:.4g}",
                           contingency_table_str]
            output_file.write('\t'.join(output_line) + '\n')
            analyzed_variants += 1

        logger.info(f'Finished processing. Total variants processed: {total_variants}, variants analyzed: {analyzed_variants}.')

if __name__ == '__main__':
    main()
