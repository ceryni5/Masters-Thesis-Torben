import pandas as pd
import numpy as np

import sqlite3
import requests

from sqlite3 import Error
from bokeh.models import ColumnDataSource

from datetime import datetime

pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


class Model(object):

    def __init__(self):

        self.database = r'C:\Users\torbe\PycharmProjects\postGWAS_db\million_hearts_gwas_and_friends.db'
        self.gwas_meta = 'gwas_meta_cad'

        ''' General Settings! '''
        self.feedback_text = datetime.now().strftime('%H:%M:%S') + ' Please choose a SNP.'
        self.gene_text = datetime.now().strftime('%H:%M:%S') + ' Please choose a gene.'

        # Important parameters
        self.close_snps = list()

        self.rsid = str()
        self.chrom = int()
        self.snp_pos = int()

        self.window_size = 500_000
        self.w_start = int()
        self.w_end = int()

        # Settings which can be taken from the settings panel
        self.population = '1000GENOMES:phase_3:ALL'


        ''' Column Data Sources '''
        self.gwas_source = ColumnDataSource()
        self.gwas_proxy_source = ColumnDataSource()
        self.gwas_overlay_source = ColumnDataSource()

        self.gene_source = ColumnDataSource()
        self.reg_build_source = ColumnDataSource()
        self.tss_source = ColumnDataSource()

        self.traits_source = ColumnDataSource()

        self.TADs_source = ColumnDataSource()

        self.ABC_source = ColumnDataSource()

        self.Miller_etal_source = ColumnDataSource()

        self.CATlas_source = ColumnDataSource()

        '''Info for selects'''
        db_connection = None
        try:
            db_connection = sqlite3.connect(self.database)

        except Error as e:
            print(f"Can't connect to the database: {e}")
            quit()

        self.all_abc_celltypes = pd.read_sql(
            sql=f'SELECT DISTINCT CellType FROM abc_celltypes_tbl',
            con=db_connection
        )['CellType'].tolist()
        self.all_miller_biosample = pd.read_sql(
            sql=f'SELECT DISTINCT biosample FROM clint_miller_biotypes_tbl',
            con=db_connection
        )[
            'biosample'].tolist()
        self.all_catlas_biosample = pd.read_sql(
            sql=f'SELECT DISTINCT biosample FROM catlas_biotypes_tbl',
            con=db_connection
        )['biosample'].tolist()

        db_connection.close()

    ''' Updates to the widgets '''
    def update_gene(self, gene_input):

        # Establish a connection to the database:
        db_connection = None
        try:
            db_connection = sqlite3.connect(self.database)

        except Error as e:
            print(f"Can't connect to the database: {e}")
            quit()

        # Query for the entered symbol
        gene_query = f'''
            SELECT approved_symbol, hgnc_all_symbols_tbl.hgnc_id
            FROM hgnc_all_symbols_tbl
            LEFT JOIN hgnc_approved_symbols_tbl
            ON hgnc_all_symbols_tbl.hgnc_id = hgnc_approved_symbols_tbl.hgnc_id
            WHERE hgnc_all_symbols_tbl.all_symbols LIKE '{gene_input}';
        '''
        hgnc_info = pd.read_sql(sql=gene_query, con=db_connection)

        # Check if the symbol is an hgnc symbol
        if hgnc_info.empty:
            self.close_snps = []
            self.update_gene_div_text(
                f'❌ Unfortunatly we did not find a hgnc_id for {gene_input}. Please input a valid HGNC gene symbol!')
            db_connection.close()
            return

        hgnc_id = hgnc_info['hgnc_id'].iloc[0]
        hgnc_symbol = hgnc_info['approved_symbol'].iloc[0]

        # Check if the gene is the most likely causal for any of the proxy variants
        proxy_causal_query = f'''
            SELECT rsID
            FROM identified_proxy_SNPs_tbl
            WHERE most_likely_causal_gene_hgnc_id == {hgnc_id}
        '''
        causal_proxies = pd.read_sql(sql=proxy_causal_query, con=db_connection)
        if not causal_proxies.empty:
            self.close_snps = causal_proxies['rsID'].unique()
            self.update_gene_div_text(
                f'✔ {gene_input} was interpreted to be {hgnc_symbol}. It is the most likely causal gene for:')
            db_connection.close()
            return

        # Check if the gene is the nearest gene for any of the proxy variants
        proxy_nearest_query = f'''
                    SELECT rsID
                    FROM identified_proxy_SNPs_tbl
                    WHERE nearest_gene_hgnc_id == {hgnc_id}
                '''
        proxy_nearest = pd.read_sql(sql=proxy_nearest_query, con=db_connection)

        if not proxy_nearest.empty:
            self.close_snps = proxy_nearest['rsID'].unique()
            self.update_gene_div_text(
                f'✔ {gene_input} was interpreted to be {hgnc_symbol}. It is the nearest gene for:')
            db_connection.close()
            return

        # Should this not be the case, we need to get the position of the gene
        gene_list_query = f'''
            SELECT seqid, start_GRCh38, end_GRCh38
            FROM ensembl_genelist_tbl
            WHERE gene_symbol == '{hgnc_symbol}';
        '''
        gene_position = pd.read_sql(sql=gene_list_query, con=db_connection)

        # Check if we got a position back:
        # ToDo: write also e.g. rRNA genes into the database. But this a job I'll probably leave for sb else. - Torben
        if gene_position.empty:
            self.close_snps = []
            self.update_gene_div_text(
                f'❌ {gene_input} was interpreted to be {hgnc_symbol}. Sorry this happend, while you happend to have '
                f'input a valid HGNC ID I cant retrieve its position because the gene product is not protein coding '
                f'and the underlying database at this point only contains protein coding genes. I am aware of this '
                f'and very sorry, this is on my ToDo list...')
            db_connection.close()
            return

        # Check if it might lay on a sex or the MT chrom
        if gene_position['seqid'].iloc[0] in ['X', 'Y', 'MT']:
            self.close_snps = []
            self.update_gene_div_text(
                f'❌ {gene_input} was interpreted to be {hgnc_symbol}. Unfortunatly the gene is on a sex or the MT '
                f'chromosome which are not part of the summary stats.')
            db_connection.close()
            return

        # Get the SNP with the highest -log p in the area (+/-10_000 bp)
        highest_p_snp_query = f'''
            SELECT chr, pos_hg38_one_based, rsid, log10_p_value
            FROM variation
            LEFT JOIN gwas_meta_cad
            ON gwas_meta_cad.variant_id = variation.id
            WHERE chr == {gene_position['seqid'].iloc[0]}
            AND pos_hg38_one_based 
                BETWEEN {gene_position['start_GRCh38'].iloc[0]} AND {gene_position['end_GRCh38'].iloc[0]} 
            ORDER BY log10_p_value ASC
            LIMIT 1
        '''
        highest_p_snp = pd.read_sql(sql=highest_p_snp_query, con=db_connection)

        if not highest_p_snp.empty:
            self.close_snps = highest_p_snp['rsid'].unique()
            self.update_gene_div_text(
                f'✔ {gene_input} was interpreted to be {hgnc_symbol}. The variant with the highest -log(p) in the '
                f'areas is:')
            db_connection.close()
            return

        # How did we even get here.
        self.update_gene_div_text(
            f'❌ {gene_input} was interpreted to be {hgnc_symbol}. WTF how did we even get here. Please excuse my '
            f'laguage but it seems like your gene is not associated or overlapping with any SNP from the summar '
            f'stats. How is this even possible...')
        db_connection.close()
        return

    def update_snp(self, rsid):

        # Establish a connection to the database:
        db_connection = None
        try:
            db_connection = sqlite3.connect(self.database)

        except Error as e:
            print(f"Can't connect to the database: {e}")
            quit()

        # Find the SNP and corresponding information in the database:
        snp_query = f'''
            SELECT chr, pos_hg38_one_based, rsid, log10_p_value
            FROM variation
            LEFT JOIN gwas_meta_cad
            ON gwas_meta_cad.variant_id = variation.id
            WHERE rsid == '{rsid}'
        '''
        snp_df = pd.read_sql(sql=snp_query, con=db_connection)

        if snp_df.empty:
            self.rsid = str()
            self.chrom = int()
            self.snp_pos = int()
            self.w_start = int()
            self.w_end = int()
            self.update_feedback_div_text(
                f'❌ <b>{rsid} does not seem to be part of the summary stats. Please check your input or it might not '
                f'be located on an autosome.</b>')
            return

        else:
            self.rsid = snp_df['rsid'].iloc[0]
            self.chrom = snp_df['chr'].iloc[0]
            self.snp_pos = snp_df['pos_hg38_one_based'].iloc[0]
            self.w_start = int(self.snp_pos - (self.window_size / 2))
            self.w_end = int(self.snp_pos + (self.window_size / 2))

            self.update_feedback_div_text(f'✔ {rsid} is a viable SNP name.')

        db_connection.close()
        return

    def update_feedback_div_text(self, text):

        self.feedback_text = datetime.now().strftime('%H:%M:%S') + f' {text}'

    def update_gene_div_text(self, text):

        self.gene_text = datetime.now().strftime('%H:%M:%S') + f' {text}'


    ''' Updates to data '''
    def update_all_data(self):

        # Establish a connection to the database
        conn_ = None
        try:
            conn_ = sqlite3.connect(self.database)

        except Error as e:
            print(f"Can't connect to the database: {e}")
            quit()

        # Call all the update functions
        self.update_gwas(conn_)
        self.update_genes_regbuild_tss(conn_)
        self.update_gwas_catalog(conn_)
        self.update_TADs(conn_)
        self.update_ABC(conn_)
        self.update_scATAC_seq(conn_)
        self.update_CATlas(conn_)
        conn_.close()

        self.update_feedback_div_text(f'✔ Plot successfully updated for {self.rsid}.')
        self.update_gene_div_text(f'✔ Plot successfully updated for {self.rsid}.')
        return

    def update_population(self, pop):

        self.population = '1000GENOMES:phase_3:' + pop

        # Establish a connection to the database
        conn_ = None
        try:
            conn_ = sqlite3.connect(self.database)

        except Error as e:
            print(f"Can't connect to the database: {e}")
            quit()

        # Call all the update functions
        self.update_gwas(conn_)

        return

    def update_gwas(self, db_connection):

        # Query the database for the GWAS summary stats
        query = f"""
        SELECT v.pos_hg38_one_based, v.rsid, gwas.log10_p_value, v.maf
            FROM {self.gwas_meta} gwas
            INNER JOIN (
                SELECT variation.id, variation.pos_hg38_one_based, variation.rsid, variation.maf
                FROM variation 
                WHERE (variation.chr = '{self.chrom}'
                AND variation.pos_hg38_one_based BETWEEN {self.w_start} AND {self.w_end})) v
            ON gwas.variant_id = v.id
        """
        gwas_df = pd.read_sql(sql=query, con=db_connection)

        # Transfer all the stuff from dataframe into the datasources
        self.gwas_source.data = dict(
            Pos_hg38=gwas_df['pos_hg38_one_based'],
            logP_value=-gwas_df['log10_p_value']
        )

        proxy_df = gwas_df[gwas_df['rsid'] == self.rsid]
        self.gwas_proxy_source.data = dict(
            Pos_hg38=proxy_df['pos_hg38_one_based'],
            logP_value=-proxy_df['log10_p_value'],
        )

        # Query for LD structure and annotation
        ld_query = f"""
            SELECT linked_rsid, chr, pos_hg38, r2, consequence
            FROM linked_SNPs_tbl
            
            LEFT JOIN consequence_tbl
            ON consequence_tbl.consequence_id = linked_SNPs_tbl.consequence_id
            
            LEFT JOIN population_tbl
            ON population_tbl.population_id = linked_SNPs_tbl.population_id
            
            WHERE proxy_rsid = '{self.rsid}'
            AND population = '{self.population}'
        """
        ld_df = pd.read_sql(sql=ld_query, con=db_connection)

        # If the return is empty, fetch the data from ensembl.org
        if ld_df.empty:
            ld_df = self.fetch_lds_from_ensembl()

        # Merge both dataframes and update the datasource for bokeh
        overlay_df = pd.merge(ld_df, gwas_df, how='inner', left_on='linked_rsid',
                              right_on='rsid')  # Inner join to only get SNPs that of the summary stats or linked

        self.gwas_overlay_source.data = dict(
            rsID=overlay_df['rsid'],
            Pos_hg38=overlay_df['pos_hg38_one_based'],
            logP_value=-overlay_df['log10_p_value'],
            r2=overlay_df['r2'],
            most_severe_consequence=overlay_df['consequence'],
            maf=overlay_df['maf']
        )

        return

    def update_genes_regbuild_tss(self, db_connection):

        gene_query = f'''
            SELECT start_GRCh38, end_GRCh38, gene_id, biotype, gene_symbol, strand
            FROM ensembl_genelist_tbl
            LEFT JOIN ensembl_genelist_biotypes_tbl
            ON ensembl_genelist_tbl.biotype_id = ensembl_genelist_biotypes_tbl.biotype_id
            WHERE seqid = {self.chrom} AND (
            
                end_GRCh38 BETWEEN {self.w_start} AND {self.w_end} OR
                start_GRCh38 BETWEEN {self.w_start} AND {self.w_end}
                
            )
        '''
        gene_df = pd.read_sql(sql=gene_query, con=db_connection)

        # Remove 'overhangs' that are not part of the window
        gene_df['start_GRCh38'] = gene_df['start_GRCh38'].apply(lambda x: max(x, self.w_start))
        gene_df['end_GRCh38'] = gene_df['end_GRCh38'].apply(lambda x: min(x, self.w_end))

        # The name is the gene_symbol or the gene_id should no symbol be available
        gene_df['name'] = gene_df.apply(
            lambda row: row['gene_symbol'] if row['gene_symbol'] is not None else row['gene_id'], axis=1)

        # Get the levels where we want to put the genes:
        # Level assignment is done for each biotype:
        biotypes = gene_df['biotype'].unique()
        levels = {'init': -1}
        for idx, biotype in enumerate(biotypes):

            # Get entries of one biotype
            entries = gene_df[gene_df['biotype'] == biotype]

            # Sort entries and get the start and end points of the transcripts
            entries = entries.sort_values(by=['end_GRCh38'], ascending=True)
            starts = entries['start_GRCh38'].values
            ends = entries['end_GRCh38'].values
            names = entries['name'].values
            name_length = [len(x) * 1666 for x in names]

            # Iterate through the genes, check if its start is after the end point of the last gene we put on the level
            distance_correction = 10_000
            levels_whithin_biotype = [0] * len(starts)  # Levels we want to put the genes of this biotype on
            level_ends = [0] * len(starts)  # End point of the last gene we put on this level
            level_ends[0] = -(distance_correction + 1)  # Initilize the first cell,

            for i in range(len(starts)):
                for j in range(len(level_ends)):
                    if starts[i] > (level_ends[j] + distance_correction) or starts[i] > (
                            level_ends[j] + name_length[i]):
                        level_ends[j] = ends[i]
                        levels_whithin_biotype[i] = j
                        break

            # Convert these levels into 'absolut' levels & write these into the levels list outside the loop
            base_level_of_biotype = max(levels.values()) + 1
            for name, level in zip(names, levels_whithin_biotype):
                levels[name] = base_level_of_biotype + level

        # Add level to df
        gene_df['level'] = gene_df['name'].map(levels)

        # Add pos for label to the df
        gene_df['label_pos'] = gene_df['start_GRCh38'] + ((gene_df['end_GRCh38'] - gene_df['start_GRCh38']) / 2)

        # Also fetch data for open target genetics l2g scores!
        l2g_query = f'''
            SELECT approved_symbol, y_proba_full_model
            FROM opentarget_l2g_tbl
            
            LEFT JOIN hgnc_approved_symbols_tbl
            
            ON opentarget_l2g_tbl.hgnc_id = hgnc_approved_symbols_tbl.hgnc_id
            
            WHERE opentarget_l2g_tbl.chrom = {self.chrom}
            AND opentarget_l2g_tbl.pos = {self.snp_pos} 
        '''

        df_l2g_score = pd.read_sql(sql=l2g_query, con=db_connection).groupby(
            by=["approved_symbol"]).median()  # ToDo: why tho?

        if df_l2g_score.empty:  # every thing is gray
            gene_df['color'] = 'gray'
        else:
            gene_df = gene_df.merge(df_l2g_score, how='left', left_on='name', right_on='approved_symbol')
            gene_df['color'] = gene_df['y_proba_full_model'].apply(lambda prob:
                                                                   'red' if prob >= 0.8 else
                                                                   'firebrick' if prob >= 0.6 else
                                                                   'lightcoral' if prob >= 0.3 else
                                                                   'gray'
                                                                   )
        # to CDS
        self.gene_source.data = dict(
            level=gene_df['level'],
            start=gene_df['start_GRCh38'],
            end=gene_df['end_GRCh38'],
            name=gene_df['name'],
            biotype=gene_df['biotype'],
            strand=gene_df['strand'],
            label_pos=gene_df['label_pos'],
            gene_id=gene_df['gene_id'],
            color=gene_df['color']
        )

        # Reg build
        reg_build_query = f'''
            SELECT start_GRCh38, end_GRCh38, ensembl_id, feature_type
            FROM ensembl_reg_build_tbl
            
            LEFT JOIN ensembl_reg_build_features_tbl
            ON ensembl_reg_build_tbl.feature_type_id = ensembl_reg_build_features_tbl.feature_type_id
            
            WHERE seqid = {self.chrom} AND (
                end_GRCh38 BETWEEN {self.w_start} AND {self.w_end} OR
                start_GRCh38 BETWEEN {self.w_start} AND {self.w_end}
            )
        '''
        reg_build_df = pd.read_sql(sql=reg_build_query, con=db_connection)

        # Dicts for mapping color and level of display to the 'description'
        color_dict = {"TSS": 'darkgray',
                      "promoter": '#FF0000',
                      "TF_binding_site": '#CD96CD',
                      "promoter_flanking_region": '#FFA4A4',
                      "open_chromatin_region": 'dimgray',
                      "enhancer": "darkgreen",
                      "CTCF_binding_site": '#40E0D0'}
        counter_dict = {"TSS": 3,
                        "promoter": 5,
                        "TF_binding_site": 4,
                        "promoter_flanking_region": 5,
                        "open_chromatin_region": 0,
                        "enhancer": 2,
                        "CTCF_binding_site": 1}

        # Append df with level and color for plotting
        reg_build_df['level'] = reg_build_df['feature_type'].map(counter_dict)
        reg_build_df['color'] = reg_build_df['feature_type'].map(color_dict)

        self.reg_build_source.data = dict(
            feature_type=reg_build_df['feature_type'],
            level=reg_build_df['level'],
            gene_id=reg_build_df['ensembl_id'],
            end_GRCh38=reg_build_df['end_GRCh38'],
            start_GRCh38=reg_build_df['start_GRCh38'],
            color=reg_build_df['color']
        )

        # TSS
        tss_query = f'''
            SELECT pos_Hg38, transcript, strand, approved_symbol
            FROM tss_tbl
            
            LEFT JOIN hgnc_approved_symbols_tbl
            ON tss_tbl.hgnc_id = hgnc_approved_symbols_tbl.hgnc_id
            
            WHERE chr = 'chr{self.chrom}'
            AND pos_Hg38 BETWEEN {self.w_start} AND {self.w_end};
        '''
        tss_df = pd.read_sql(sql=tss_query, con=db_connection)
        tss_df['level'] = 3
        tss_df['color'] = 'darkgray'

        self.tss_source.data = dict(
            pos_Hg38=tss_df['pos_Hg38'],
            level=tss_df['level'],
            gene_id=tss_df['transcript'],
            strand=tss_df['strand'],
            approved_symbol=tss_df['approved_symbol'],
            color=tss_df['color']
        )

        return

    def update_gwas_catalog(self, db_connection):

        ld_df = pd.DataFrame(self.gwas_overlay_source.data)

        ld_0point6_df = ld_df.loc[ld_df['r2'] >= 0.6][['rsID', 'Pos_hg38']]
        ld_0point6_df.reset_index(inplace=True, drop=True)

        # Fetch associated phenotypes from GWAS catalog table
        # https://stackoverflow.com/questions/30399301/special-character-in-column-name
        query = f"""
                SELECT chr_Hg38, pos_Hg38, PVALUE_MLOG, beta, Odds_Ratio, SNP_ID_CURRENT, [DISEASE/TRAIT], [REPORTED_GENE(S)]
                    FROM gwascatalog_associations_tbl
                    WHERE chr_Hg38 = {self.chrom}
                    AND pos_Hg38 IN ({', '.join(map(str, ld_0point6_df['Pos_hg38'].tolist()))})
                """
        associated_df = pd.read_sql(sql=query, con=db_connection)

        associated_df["marker"] = associated_df['beta'].map(
            lambda x: 'circle' if x == None else 'circle' if x == np.nan else 'triangle' if x > 0 else 'inverted_triangle')
        associated_df["color"] = associated_df['beta'].map(
            lambda x: 'black' if x == None else 'circle' if x == np.nan else '#2166ac' if x > 0 else '#b2182b')
        associated_df["size"] = 5

        self.traits_source.data = dict(
            p_value=associated_df['PVALUE_MLOG'],
            beta=associated_df['beta'],
            Odds_Ratio=associated_df['Odds_Ratio'],
            trait=associated_df['DISEASE/TRAIT'],
            reported_genes=associated_df['REPORTED_GENE(S)'],
            marker=associated_df['marker'],
            color=associated_df['color'],
            size=associated_df['size'],
            chr_Hg38=associated_df['chr_Hg38'],
            pos_Hg38=associated_df['pos_Hg38'],
            rsid=associated_df['SNP_ID_CURRENT']
        )

        return

    def update_TADs(self, db_connection):

        # Fetch TADs
        query = f"""
            SELECT
                tad_tbl.start_hg38,
                tad_tbl.end_hg38,
                tad_sample_tbl.sample

                FROM tad_tbl

                LEFT JOIN tad_sample_tbl
                ON tad_tbl.sample_id = tad_sample_tbl.sample_id

                WHERE (
                    tad_tbl.chr = 'chr{self.chrom}' AND (
                    
                        tad_tbl.end_hg38 BETWEEN {self.w_start} AND {self.w_end} OR
                        tad_tbl.start_hg38 BETWEEN {self.w_start} AND {self.w_end}
                    )
                )
        """
        tad_df = pd.read_sql_query(sql=query, con=db_connection)

        # Remove 'overhangs' that are not part of the window
        tad_df['start_hg38'] = tad_df['start_hg38'].apply(lambda x: max(x, self.w_start))
        tad_df['end_hg38'] = tad_df['end_hg38'].apply(lambda x: min(x, self.w_end))

        self.TADs_source.data = dict(
            start_hg38=tad_df['start_hg38'],
            end_hg38=tad_df['end_hg38'],
            biosample=tad_df['sample'],
        )

        return

    def update_ABC(self, db_connection):

        query = f"""
            SELECT
                abc_tbl.start_hg38,
                abc_tbl.end_hg38,
                abc_tbl.Score,
                abc_celltypes_tbl.CellType,
                abc_targetgenes_tbl.TargetGene

                FROM abc_tbl

                LEFT JOIN abc_celltypes_tbl
                ON abc_celltypes_tbl.CellType_id = abc_tbl.CellType_id

                LEFT JOIN abc_targetgenes_tbl
                ON abc_targetgenes_tbl.TargetGene_id = abc_tbl.TargetGene_id

                WHERE (
                    abc_tbl.chr = 'chr{self.chrom}' AND (
                    
                        abc_tbl.end_hg38 BETWEEN {self.w_start} AND {self.w_end} OR
                        abc_tbl.start_hg38 BETWEEN {self.w_start} AND {self.w_end}
                    
                    )
                    
                )
        """
        abc_df = pd.read_sql_query(sql=query, con=db_connection)

        # Remove 'overhangs' that are not part of the window
        abc_df['start_hg38'] = abc_df['start_hg38'].apply(lambda x: max(x, self.w_start))
        abc_df['end_hg38'] = abc_df['end_hg38'].apply(lambda x: min(x, self.w_end))

        self.ABC_source.data = dict(
            start_hg38=abc_df['start_hg38'],
            end_hg38=abc_df['end_hg38'],
            Score=abc_df['Score'],
            CellType=abc_df['CellType'],
            TargetGene=abc_df['TargetGene'],
        )

        return

    def update_scATAC_seq(self, db_connection):

        # Fetch Miller data
        query = f"""
            SELECT
                clint_miller_tbl.chromStart,
                clint_miller_tbl.chromEnd,
                clint_miller_tbl.score,
                clint_miller_tbl.qValue,
                clint_miller_tbl.peak,
                clint_miller_biotypes_tbl.biosample

                FROM clint_miller_tbl

                LEFT JOIN clint_miller_biotypes_tbl
                ON clint_miller_biotypes_tbl.biosample_id = clint_miller_tbl.biosample_id

                WHERE (
                    clint_miller_tbl.chrom = 'chr{self.chrom}' AND (
                    
                        clint_miller_tbl.chromEnd BETWEEN {self.w_start} AND {self.w_end} OR
                        clint_miller_tbl.chromStart BETWEEN {self.w_start} AND {self.w_end}
                    )
                )
        """
        miller_df = pd.read_sql_query(sql=query, con=db_connection)

        # Remove 'overhangs' that are not part of the window
        miller_df['chromStart'] = miller_df['chromStart'].apply(lambda x: max(x, self.w_start))
        miller_df['chromEnd'] = miller_df['chromEnd'].apply(lambda x: min(x, self.w_end))

        self.Miller_etal_source.data = dict(
            start_hg38=miller_df['chromStart'],
            end_hg38=miller_df['chromEnd'],
            score=miller_df['score'],
            qValue=miller_df['qValue'],
            center=miller_df['chromStart'] + miller_df['peak'],
            biosample=miller_df['biosample'],
        )

        return

    def update_CATlas(self, db_connection):

        # Fetch CATLAS data
        query = f"""
            SELECT
                catlas_tbl.chromStart,
                catlas_tbl.chromEnd,
                catlas_tbl.score,
                catlas_tbl.qValue,
                catlas_tbl.peak,
                catlas_biotypes_tbl.biosample

                FROM catlas_tbl

                LEFT JOIN catlas_biotypes_tbl
                ON catlas_biotypes_tbl.biosample_id = catlas_tbl.biosample_id

                WHERE (
                    catlas_tbl.chrom = 'chr{self.chrom}' AND (
                        catlas_tbl.chromEnd BETWEEN {self.w_start} AND {self.w_end} OR
                        catlas_tbl.chromStart BETWEEN {self.w_start} AND {self.w_end}
                    )
                )
        """
        catlas_df = pd.read_sql_query(sql=query, con=db_connection)

        # Remove 'overhangs' that are not part of the window
        catlas_df['chromStart'] = catlas_df['chromStart'].apply(lambda x: max(x, self.w_start))
        catlas_df['chromEnd'] = catlas_df['chromEnd'].apply(lambda x: min(x, self.w_end))

        self.CATlas_source.data = dict(
            start_hg38=catlas_df['chromStart'],
            end_hg38=catlas_df['chromEnd'],
            score=catlas_df['score'],
            qValue=catlas_df['qValue'],
            center=catlas_df['chromStart'] + catlas_df['peak'],
            biosample=catlas_df['biosample'],
        )

        return

    def fetch_lds_from_ensembl(self):
        """
        Fetch LD structure from ensembl:
        # Computes and returns LD values between the given variant
        # and all other variants in a window centered around the given variant.
        # The window size is set to 500 kb.

        :return: dataFrame with the LDs from our window
        """

        server = "https://rest.ensembl.org"
        # First get info in the SNP:

        var_ext = "/variation/human/%s?" % self.rsid
        r = requests.get(server + var_ext, headers={"Content-Type": "application/json"})

        print('------> rsID: ', self.rsid)

        if not r.ok:
            print(self.rsid + " not found in Ensembl.org try an alias")

        decoded = r.json()

        proxy_chr = decoded["mappings"][0]["seq_region_name"]
        proxy_pos = decoded["mappings"][0]["start"]
        proxy_consequence = decoded['most_severe_consequence']

        ld_ext = "/ld/%s/%s/%s" % ('human', self.rsid, self.population)

        r = requests.get(
            server + ld_ext,
            headers={"Content-Type": "application/json"},
            params={
                'attribs': 1
            }
        )

        if not r.ok:
            print(" not found in Ensembl.org try an alias")
            r.raise_for_status()
        else:
            decoded = r.json()

        proxy_variant, linked_variant, chrom, pos_hg38, r2, consequence, population_list = [self.rsid], \
                                                                                           [self.rsid], \
                                                                                           [proxy_chr], \
                                                                                           [proxy_pos], \
                                                                                           [1.0], \
                                                                                           [proxy_consequence], \
                                                                                           [self.population]

        # Fill these lists:
        for elem in decoded:
            proxy_variant.append(self.rsid)
            linked_variant.append(elem['variation'])
            chrom.append(elem['chr'])
            pos_hg38.append(elem['start'])
            r2.append(float(elem['r2']))
            consequence.append(elem['consequence_type'])
            population_list.append(self.population)

        df = pd.DataFrame({
            'proxy_rsid': proxy_variant,
            'linked_rsid': linked_variant,
            'chr': chrom,
            'pos_hg38': pos_hg38,
            'r2': r2,
            'consequence': consequence,
            'population': population_list,
        })

        return df