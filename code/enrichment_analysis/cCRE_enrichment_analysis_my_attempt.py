import sqlite3
from sqlite3 import Error
import pandas as pd
import os

import sys

import datetime
from scipy.stats import hypergeom


# Enrichtment Analysis:
'''
To check if overlap with cCREs is enriched in a celltype (biosample -> file accession).
We need 4 different parameters to fill up the following Contingency table:

                    nt = 15         mt - nt = 485           mt = 500
                    n - nt = 285    m - n - mt + nt = 19215 m - mt = 19500
                    n = 300         m - n = 19700           m = 20000

In this case with example values for a GSEA. 
For our case we need following four values:

m is the total number of distinct cCREs (count ids in database)
mt is the number of cCREs that are annotated in biosample (file accession) we are checking
n is the number of cCREs that are overlapping over all the biosamples
nt is the number of cCREs that are overlapping in the biosample we are making the analysis for

                    +------------------------------------------+                                                                    
                    |         m = all distinct cCREs           |                                                                    
                    |                                          |                                                                    
                    |              +----------------------+    |                                                                    
                    |              |  mt = cCREs in the   |    |                                                                    
                    |              |  biosample           |    |                                                                    
                    |              |                      |    |                                                                    
                    |     +--------|---------+            |    |                                                                    
                    |     |        | nt      |            |    |                                                                    
                    |     |        |         |            |    |                                                                    
                    |     |        +----------------------+    |                                                                    
                    |     |                  |                 |                                                                    
                    |     | n = all          |                 |                                                                    
                    |     | overlapping cCREs|                 |                                                                    
                    |     |                  |                 |                                                                    
                    |     +------------------+                 |                                                                    
                    +------------------------------------------+ 
                    
We can then calculate the p value for observing nt hits given the other parameters of the distribution just by chance
using Fisher's exact test:

In python there are multiple ways:
Using the fisher exact fuction from scipy

                    Contigency Table:
                    nt = 4          mt - nt = 2             mt = 6
                    n - nt = 1      m - n - mt + nt = 3     m - mt = 4
                    n = 5           m - n = 5               m = 10

                    this would correspond to this table 'in code':

                    table = np.array([
                        [4, 2],
                        [1, 3]
                    ])

                    oddsr, p_fisher = fisher_exact(table, alternative='greater')
                    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
                    
Or directly plugging the values into the cumulative density function (or in our case survival function):
                    
                    nt = table[0, 0]
                    mt = table[0, 0] + table[0, 1]
                    n = table[0, 0] + table[1, 0]
                    m = table.sum()
                    
                    p_hyper_geom = hypergeom.sf(nt - 1, m, n, mt)
                    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html

Or it can be calculate by 'manually' computing the cumulative density and substracting it from 1:

                    p_manual = sum([math.comb(n, k)*math.comb(m-n, mt-k)/math.comb(m, mt) for k in range(nt, min(mt, n) + 1)])
                    
'''

if __name__ == '__main__':

    write_console_to_logfile = False

    # Open file to log the outputs:
    log_dir = os.path.join('.', 'ccre_enrichment_logs')
    if not os.path.exists(log_dir):
        os.system(f'mkdir {log_dir}')
    log_file = os.path.join(log_dir, f''
                                     f'{datetime.datetime.now().year}-'
                                     f'{datetime.datetime.now().month}-'
                                     f'{datetime.datetime.now().day}--'
                                     f'{datetime.datetime.now().hour}-'
                                     f'{datetime.datetime.now().minute}-'
                                     f'{datetime.datetime.now().second}.txt')
    if write_console_to_logfile:
        sys.stdout = open(log_file, 'w')

    # Pandas setup for better display of dataframes
    pd.set_option('display.max_rows', 200)
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.width', 1500)
    pd.set_option('max_colwidth', 200)

    # Database with cCRE data
    ccre_db = r'C:\Users\torbe\PycharmProjects\postGWAS_db\cCREs.db'
    cCRE_database_con = None
    try:
        cCRE_database_con = sqlite3.connect(ccre_db)
    except Error as e:
        print(e)


    #  Database with Million Hearts data
    millions_hearts_db = r'C:\Users\torbe\PycharmProjects\postGWAS_db\million_hearts_gwas_and_friends.db'
    million_hearts_con = None
    try:
        million_hearts_con = sqlite3.connect(millions_hearts_db)
    except Error as e:
        print(e)

    # Get the list of 1 Million hearts proxy SNPs and their linked SNPs (within 1000GENOMES:phase_3:ALL)
    snp_query = f'''
        SELECT *
        FROM linked_SNPs_tbl
        
        LEFT JOIN population_tbl
        ON population_tbl.population_id = linked_SNPs_tbl.population_id
        
        WHERE population_tbl.population = "1000GENOMES:phase_3:EUR"
    '''
    print(f'Fetching proxy SNPs and lined variants with r2 >= {(r2_threshhold := 0.6)}.')
    linked_snps = pd.read_sql(sql=snp_query, con= million_hearts_con)
    million_hearts_con.close()

    # Filter for r2 >= 0.6
    linked_snps_r2_bigger_0point6 = linked_snps.loc[linked_snps['r2'] >= r2_threshhold]
    print(linked_snps_r2_bigger_0point6)

    # We get a list of all the file accessions
    file_accession_query = f'''
        SELECT DISTINCT file_accession
        FROM encode_cCRE
    '''
    print('Fetching distinct file accessions:')
    individual_file_accessions = pd.read_sql(sql=file_accession_query, con=cCRE_database_con)
    print(individual_file_accessions)

    # We get a list of all the cCREs
    elements_query = f'''
        SELECT DISTINCT ccre_id
        FROM encode_cCRE
    '''
    print('Fetching distinct cCREs:')
    individual_ccre = pd.read_sql(sql=elements_query, con=cCRE_database_con)
    print(individual_ccre)

    # We get a list of the number of cCREs that are annotated in each file_accession
    elements_per_file_query = f'''
            SELECT file_accession, count(DISTINCT ccre_id) numb_cCREs
            FROM encode_cCRE
            GROUP BY file_accession
        '''
    print('Fetching number distinct cCREs for each biosample:')
    ccre_per_file = pd.read_sql(sql=elements_per_file_query, con=cCRE_database_con)
    print(ccre_per_file)

    # We get a df with all the overlapping cCREs
    print('Fetching overlapping cCREs:')
    df_list = []
    for proxy, proxy_df in linked_snps_r2_bigger_0point6.groupby('proxy_rsid'):
        chr = int(proxy_df['chr'].mean())
        start = proxy_df['pos_hg38'].min()
        end = proxy_df['pos_hg38'].max()

        # Elements have their start between the first and last SNP of the LD strcuture
        # Their end or span it completely
        # And those that reach over all the SNPs
        overlap_query = f'''
            SELECT *
            FROM encode_cCRE 
            WHERE "chr{chr}" = chr
            AND (
                {start} BETWEEN start_hg38_zero_based AND end_hg38_one_based OR
                {end} BETWEEN start_hg38_zero_based AND end_hg38_one_based OR
                (start_hg38_zero_based <= {start} AND {end} <= end_hg38_one_based)
            )
        '''
        ccre_region = pd.read_sql(sql=overlap_query, con=cCRE_database_con)

        df_list.append(ccre_region)

    overlapping_ccres = pd.concat(df_list)

    # Get the number of distinct overlapping cCREs
    distinct_overlapping_ccre = overlapping_ccres['ccre_id'].unique()

    # Get the number of distinct overlapping cCREs in each file asseccion
    number_of_overlaps_dict = {}
    for file_accession in individual_file_accessions['file_accession'].tolist():
        overlaps_for_biosample = overlapping_ccres.loc[
            overlapping_ccres['file_accession'] == file_accession
            ].reset_index(drop=True)
        overlaps_for_biosample.drop_duplicates(subset='ccre_id', inplace=True)
        number_of_overlaps_dict[file_accession] = len(overlaps_for_biosample)

    # Now we have all the values to calculate the enrichment factor and the p-value using the hypergeometric distribution
    summary = {
        'nt': [],
        'mt': [],
        'n': [],
        'm': [],
        'p': [],
        'q': [],
        'enrichment_factor': [],
        'file_accession': []
    }
    for file in individual_file_accessions['file_accession'].tolist():

        # m is the total number of distinct cCREs (count ids in database)
        m = len(individual_ccre)

        # mt is the number of cCREs that are annotated in biosample (file accession) we are checking
        mt = ccre_per_file.loc[ccre_per_file['file_accession'] == file]['numb_cCREs'].iloc[0]

        # n is the number of cCREs that are overlapping over all the biosamples
        n = len(distinct_overlapping_ccre)

        # nt is the number of cCREs that are overlapping in the biosample we are making the analysis for
        nt = number_of_overlaps_dict[file]

        p = hypergeom.sf(nt - 1, m, n, mt)
        enrichment_factor = (nt*m)/(n*mt)
        q = p * len(individual_file_accessions)

        summary['nt'].append(nt)
        summary['mt'].append(mt)
        summary['n'].append(n)
        summary['m'].append(m)
        summary['p'].append(p)
        summary['q'].append(q)
        summary['enrichment_factor'].append(enrichment_factor)
        summary['file_accession'].append(file)

    summary_df = pd.DataFrame(summary)
    summary_df = summary_df.sort_values('q')

    print(summary_df)
    summary_df.to_csv('enrichment_summary.tsv', sep = '\t', index = False)

    if write_console_to_logfile:
        sys.stdout.close()
    quit()