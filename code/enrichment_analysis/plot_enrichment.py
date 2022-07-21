import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

import sqlite3
from sqlite3 import Error

plt.rcParams['svg.fonttype'] = 'none'

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

enrichement_df = pd.read_csv(
    'enrichment_summary.tsv',
    sep = '\t'
)

# Translate file_accession to sample
file_accession_query = f'''
    SELECT biosample_term_name, file_accession, biosample_attr_1, biosample_attr_2, biosample_attr_3
    FROM encode_ccRE_meta
'''
file_accession_df = pd.read_sql(file_accession_query, cCRE_database_con)

enrichement_df = enrichement_df.merge(file_accession_df)
sample_color_dict={
    'heart': '#f15025',
    'artery': '#ffc857',
    'lung': '#62c370',
    'liver': '#080357',
    'brain': '#0075a2',
    'immune cell': '#54428e',
    'GI tract': '#820263',
    'other': 'gray'
}
color_sample_dict={y: x for x, y in sample_color_dict.items()}

enrichement_df['biosample_attr_3'] = enrichement_df['biosample_attr_3'].str.replace('gastrointestinal tract', 'GI tract')
enrichement_df['color'] = enrichement_df['biosample_attr_3'].map(sample_color_dict)
enrichement_df['color'] = enrichement_df['color'].fillna('gray')

sig_rows_df = enrichement_df.loc[(enrichement_df['q'] <= 0.05)]

fig, ax = plt.subplots(figsize= (6, 8))


ax.scatter(
    x=enrichement_df["enrichment_factor"],
    y=enrichement_df["q"],
    color=enrichement_df["color"],
    alpha=0.25
)

for color, dff in sig_rows_df.groupby('color'):
    ax.scatter(
        x=dff["enrichment_factor"],
        y=dff["q"],
        c=dff['color'],
        label=f'{color_sample_dict[color]}'
    )

plt.yscale('log')
plt.gca().invert_yaxis()

plt.plot(
    [enrichement_df['enrichment_factor'].min(), enrichement_df['enrichment_factor'].max()],
    [0.05, 0.05],
    ':',
    color='#C0C0C0'
)
plt.plot(
    [2, 2],
    [enrichement_df['q'].min(), enrichement_df['q'].max()],
    ':',
    color='#C0C0C0'
)

sig_rows_df = enrichement_df.loc[(enrichement_df['q'] <= 0.05) & (enrichement_df['enrichment_factor'] >= 3)]

print(enrichement_df.loc[(enrichement_df['q'] <= 0.05)].groupby('biosample_attr_3').count())

print(enrichement_df.loc[(enrichement_df['q'] <= 0.05)])

enrichement_df.loc[(enrichement_df['q'] <= 0.05)].round(1).to_csv('enriched_biosamples.csv', sep=',')

LABELS = []
for idx, row in sig_rows_df.iterrows():
    LABELS.append(ax.text(
        row['enrichment_factor'],
        row['q'],
        row['biosample_term_name'],
        fontsize=8,
        c = 'black'
    ))

adjust_text(
    LABELS,
    expand_points=(2, 2),
    expand_text=(2, 3),
    arrowprops=dict(
        arrowstyle="->",
        lw=1,
        color='#808080'
    ),
    ax=fig.axes[0]
)

# Make everything pretty
# Hide spines
ax.spines["right"].set_color("none")
ax.spines["top"].set_color("none")

# Customize spines color
ax.spines["left"].set_color("#7F7F7F")
ax.spines["left"].set_linewidth(1)
ax.spines["bottom"].set_color("#7F7F7F")
ax.spines["bottom"].set_linewidth(1)

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
ax.set_ylabel("$p_{adj.}$", size=12)

ax.set_xlabel("enrichment factor", size=12)

plt.legend(frameon=False, ncol=2, loc ='lower right', fontsize=10)
plt.show()