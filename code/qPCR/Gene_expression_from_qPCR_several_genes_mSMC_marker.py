import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

plt.rcParams['svg.fonttype'] = 'none'

# Pandas setup for better display of dataframes
pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)


def calc_significance(array1, array2, test="T-test"):
    """ returns sign symbols """
    alpha0 = 0.05
    alpha1 = 0.01
    alpha2 = 0.001
    identical = True
    for elem in array1:
        if elem not in array2:
            identical = False
            break
    if identical:
        return ' ', np.nan

    if test == "T-test":
        if 1.3 > np.std(array1)/np.std(array2) > 0.7:
            equal_std = True
        else:
            equal_std = False

        t_, p_ = stats.ttest_ind(array1, array2, equal_var=True)  # if False = Welch's t test
    else:
        t_, p_ = stats.mannwhitneyu(array1, array2)
        equal_std = 'mannwhitneyu used'

    if alpha1 < p_ <= alpha0:
        return '*', round(p_, 3), equal_std
    elif alpha2 < p_ <= alpha1:
        return '**', round(p_, 3), equal_std
    elif p_ <= alpha2:
        return '***', round(p_, 3), equal_std
    elif p_ <= 0.10:
        return str(round(p_, 3)), round(p_, 3), equal_std
    else:
        return 'n.s.', round(p_, 3), equal_std


def add_sig_bars(df: object, gene_expr: str, list_of_tuples: list):

    conditions = df['sample'].unique()
    conditions_dict = dict(zip(conditions, range(len(conditions))))

    # Calculate significance matrix:
    p_matrix = np.empty(
        dtype='float',
        shape=(len(conditions), len(conditions))
    )
    sig_matrix = np.empty(
        dtype='object',
        shape=(len(conditions), len(conditions))
    )

    for idx1, cond1 in enumerate(conditions):
        for idx2, cond2 in enumerate(conditions):
            array1 = df.loc[df['sample'] == cond1, gene_expr].values
            array2 = df.loc[df['sample'] == cond2, gene_expr].values

            # Remove nan values
            array1 = [x for x in array1 if not np.isnan(x)]
            array2 = [x for x in array2 if not np.isnan(x)]

            result = calc_significance(array1, array2, test=test_to_use)

            p_matrix[idx1][idx2] = result[1]
            sig_matrix[idx1][idx2] = result[0]

    sig_matrix_df = pd.DataFrame(sig_matrix, columns=conditions, index=conditions)
    print(f'Significance matrix for: {gene_expr}')
    print(sig_matrix_df)

    y_of_bar = df[gene_expr].max() # max point in dataframe, to this we'll add the level

    for entry in list_of_tuples:

        cond1 = entry[0]
        cond2 = entry[1]
        level = entry[2]

        plt.text(
            x= (conditions_dict[cond1] + conditions_dict[cond2]) / 2,
            y= y_of_bar + (y_of_bar/20*(level+1)),
            s= f"{sig_matrix_df[cond1][cond2]}",
            ha="center",
            va="center",
            fontdict={"size": 8}
        )
        plt.hlines(
            xmin= min(conditions_dict[cond1], conditions_dict[cond2]) + 0.2,
            xmax= max(conditions_dict[cond1], conditions_dict[cond2]) - 0.2,
            y=y_of_bar + (y_of_bar/20*level),
            linewidth=1,
            color='black'
        )


# Read in data
path = r""
data = "cnn1_mmp9_qpcr_n1&n2&n3&n4&n5.xlsx"
df = pd.read_excel(path + data)

# Define conditions
genes = ['CNN1', 'MMP9']
test_to_use = "T-test"
# test_to_use = '' # default to Mann-Whitney-U Test

# Sort by matrix/plastic:
df[['stimulation', 'base']] = df['sample'].str.split(" ", 1, expand=True,)
df.sort_values(by=['base', 'stimulation'], inplace=True)

df['display_sample'] = df['sample'].str.replace('Plastik', 'Plastic')
df['display_sample'] = df['display_sample'].str.replace(' ', '\n')


conditions = df["display_sample"].unique()

df = df.loc[df['group'] != 'n4'] # n4 was wierd. Expression was overall a lot higher then expected

# Setup the plot
cm = 1/2.54  # centimeters in inches
fig, axes = plt.subplots(len(genes), 1, figsize=(14*cm, 16*cm))
plt.subplots_adjust(wspace=0.7, hspace=0.6, bottom=0.25, left=0.2)
row_counter = 0
my_palette = dict(zip(conditions, ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#8dd3c7', '#ffffb3', '#bebada', '#fb8072']))
my_palette2 = dict(zip(conditions, ['#2b83ba', '#abdda4', '#fdae61', '#d7191c', '#2b83ba', '#abdda4', '#fdae61', '#d7191c']))
# Palettes from: https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=4

for idx, gene_expr in enumerate(df.columns):

    if gene_expr in genes:

        # normalize to the mean of ++ on matrix
        mean = df.loc[(df['stimulation'] == '++') & (df['base'] == 'Matrix'), gene_expr].mean()
        df[gene_expr] = df[gene_expr] / mean

        active_ax = axes[row_counter]
        plt.sca(active_ax)
        sns.barplot(
            ax=active_ax,
            x='display_sample',
            y=gene_expr,
            data=df,
            palette = my_palette2,
            hatch = [None, None, None, None, '//', '//', '//', '//'],
            alpha=0.99,
        )
        sns.swarmplot(
            ax=active_ax,
            x='display_sample',
            y=gene_expr,
            data=df,
            linewidth=0.3,
            size=6,
            palette=my_palette2
        )
        row_counter += 1

        list_of_tuples = [
            # ('sample1', 'sample2', level)
            ('-- Matrix', '+- Matrix', 5),
            ('++ Matrix', '+- Matrix', 5),
            ('-- Matrix', '++ Matrix', 8),

            ('-- Plastik', '+- Plastik', 5),
            ('++ Plastik', '+- Plastik', 5),
            ('-- Plastik', '++ Plastik', 8),

            ('++ Matrix', '++ Plastik', 11),
        ]
        add_sig_bars(df, gene_expr, list_of_tuples)

        plt.sca(active_ax)
        plt.ylabel("Rel. gene expression")
        plt.ylim(0, df[gene_expr].max()*1.7)
        plt.xlabel("")
        plt.title(gene_expr)

        yticks = active_ax.yaxis.get_major_ticks()
        yticks[0].label1.set_visible(False)


# fig.suptitle("Expression of CNN1 \& MMP9 in HAoSMCs", fontsize="x-large")
plt.show()