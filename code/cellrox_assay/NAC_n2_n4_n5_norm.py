import pandas as pd
from matplotlib import pyplot as plt
import os
from PIL import Image
import numpy as np
import glob
import seaborn as sns
from scipy import stats

plt.rcParams['svg.fonttype'] = 'none'

def evaluate_images(green_file, blue_file, exclude_padding = 90, threshhold_green = 90, threshhold_blue = 100):

    # Open images
    green_im = Image.open(green_file)
    blue_im = Image.open(blue_file)

    # Load into arrays
    green_ar = np.asarray(green_im)
    blue_ar = np.asarray(blue_im)

    green_pixel = (green_ar[exclude_padding:-exclude_padding, exclude_padding:-exclude_padding,
                   :] > threshhold_green).sum()
    blue_pixel = (blue_ar[exclude_padding:-exclude_padding, exclude_padding:-exclude_padding,
                  :] > threshhold_blue).sum()

    green_im.close()
    blue_im.close()

    return green_pixel, blue_pixel


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
        return ' '

    if test == "T-test":
        if 1.3 > np.std(array1)/np.std(array2) > 0.7:
            equal_std = True
        else:
            equal_std = False
        print(equal_std)
        t_, p_ = stats.ttest_ind(array1, array2, equal_var=equal_std)  # if False = Welch's t test
    else:
        t_, p_ = stats.mannwhitneyu(array1, array2)
        equal_std = 'mannwhitneyu used'
    print(p_)

    if alpha1 < p_ <= alpha0:
        return '*', equal_std
    elif alpha2 < p_ <= alpha1:
        return '**', equal_std
    elif p_ <= alpha2:
        return '***', equal_std
    elif p_ <= 0.10:
        return str(round(p_, 3)), equal_std
    else:
        return 'n.s.', equal_std


def main(cache_file, test = 'T-test'):

    if not os.path.isfile(cache_file):

        n2_wells = glob.glob(
            r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk/CellRox/' + '2022_03_31_NAC_quench_n2_150min' + r'/B*'
        )
        n4_wells = glob.glob(
            r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk/CellRox/' + '2022_05_04_NAC_quench_n4_150' + r'/B*'
        )
        n5_wells = glob.glob(
            r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk/CellRox/' + '2022_05_12_NAC_quench_n5_180min' + r'/B*'
        )

        wells = n2_wells+n4_wells+n5_wells


        green_images = []
        blue_images = []
        for well in wells:
            green_images.append(glob.glob(well + r'\*CH2.tif')[0])
            blue_images.append(glob.glob(well + r'\*CH1.tif')[0])


        data = dict({
                'condition': [],
                'repeat': [],
                'green_pixel': [],
                'blue_pixel': []
            })
        for i, (green, blue) in enumerate(zip(green_images, blue_images)):

            # Extract green and blue pixel counts
            green_pixel, blue_pixel = evaluate_images(green, blue)
            data['green_pixel'].append(green_pixel)
            data['blue_pixel'].append(blue_pixel)
            data['condition'].append(green.split('\\')[-2].split('_')[2].split('-')[0])
            data['repeat'].append(green.split('\\')[-3].split('_')[5])

        data_df = pd.DataFrame(data)
        repeat_indices = data_df.index.isin([5,7])
        data_df = data_df[~repeat_indices]

        data_df['green/blue'] = data_df['green_pixel']/data_df['blue_pixel']

        data_df.to_csv(
            cache_file,
            sep='\t',
            index=False,
        )

    else:
        data_df = pd.read_csv(
            cache_file,
            sep='\t',
        )


    fig = plt.figure(figsize=(6, 5))

    my_palette = dict(zip(data_df['condition'].unique(), ['#d7191c', '#fdae61', '#ffffbf', '#abdda4', '#2b83ba']))

    '''
    Not normalized!
    '''
    ax1 = fig.add_subplot(1, 2, 1)
    sns.barplot(
        ax=ax1,
        x="condition",
        y='green/blue',
        data=data_df,
        palette='gray',
        errcolor='.2',
    )
    sns.scatterplot(
        ax=ax1,
        x="condition",
        y='green/blue',
        data=data_df,
        palette=my_palette,
        color='black',
        edgecolors='black',
        alpha=1,
        zorder =10
    )

    #calc sig
    posCtrl = data_df.loc[data_df['condition'] == 'posCtrl', 'green/blue'].values
    negCtrl = data_df.loc[data_df['condition'] == 'negCtrl', 'green/blue'].values
    two_mM = data_df.loc[data_df['condition'] == '2mM', 'green/blue'].values
    four_mM = data_df.loc[data_df['condition'] == '4mM', 'green/blue'].values
    eight_mM = data_df.loc[data_df['condition'] == '8mM', 'green/blue'].values

    # significance between Ctrls
    height = data_df['green/blue'].max()
    plt.text(
        x=0.5,
        y=height * 1.22,
        s=f'{calc_significance(posCtrl, negCtrl, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=0.9, y=height * 1.20, linewidth=2, color='black')

    # significance to 2 mM
    plt.text(
        x=1,
        y=height * 1.17,
        s=f'{calc_significance(posCtrl, two_mM, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=1.9, y=height * 1.15, linewidth=2, color='black')

    # significance to 4 mM
    plt.text(
        x=1.5,
        y=height * 1.12,
        s=f'{calc_significance(posCtrl, four_mM, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=2.9, y=height * 1.1, linewidth=2, color='black')

    # significance to 8 mM
    plt.text(
        x=2,
        y=height * 1.07,
        s=f'{calc_significance(posCtrl, eight_mM, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=3.9, y=height * 1.05, linewidth=2, color='black')

    '''
    Normalized
    '''
    mean_signal = data_df.groupby('repeat', as_index=False).sum()

    def normalize_by_mean(row, mean_signals):

        mean = mean_signals[mean_signals['repeat'] == row['repeat']]['green/blue'].iloc[0]
        print(mean)

        row['green/blue_norm'] = row['green/blue'] / mean

        return row

    data_df = data_df.apply(normalize_by_mean, axis = 1, args= (mean_signal, ))

    ax2 = fig.add_subplot(1, 2, 2, sharey=ax1)
    sns.barplot(
        ax=ax2,
        x="condition",
        y='green/blue_norm',
        data=data_df,
        palette='gray',
        errcolor='.2',
    )
    sns.scatterplot(
        ax=ax2,
        x="condition",
        y='green/blue_norm',
        data=data_df,
        palette=my_palette,
        color='black',
        edgecolors='black',
        alpha=1,
        zorder =10
    )

    #calc sig
    posCtrl = data_df.loc[data_df['condition'] == 'posCtrl', 'green/blue_norm'].values
    negCtrl = data_df.loc[data_df['condition'] == 'negCtrl', 'green/blue_norm'].values
    two_mM = data_df.loc[data_df['condition'] == '2mM', 'green/blue_norm'].values
    four_mM = data_df.loc[data_df['condition'] == '4mM', 'green/blue_norm'].values
    eight_mM = data_df.loc[data_df['condition'] == '8mM', 'green/blue_norm'].values

    # significance between Ctrls
    plt.text(
        x=0.5,
        y=height * 1.22,
        s=f'{calc_significance(posCtrl, negCtrl, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=0.9, y=height * 1.20, linewidth=2, color='black')

    # significance to 2 mM
    plt.text(
        x=1,
        y=height * 1.17,
        s=f'{calc_significance(posCtrl, two_mM, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=1.9, y=height * 1.15, linewidth=2, color='black')

    # significance to 4 mM
    plt.text(
        x=1.5,
        y=height * 1.12,
        s=f'{calc_significance(posCtrl, four_mM, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=2.9, y=height * 1.1, linewidth=2, color='black')

    # significance to 8 mM
    plt.text(
        x=2,
        y=height * 1.07,
        s=f'{calc_significance(posCtrl, eight_mM, test=test)[0]}',
        ha="center",
        va="center",
        fontdict={"size": 8}
    )
    plt.hlines(xmin=0.1, xmax=3.9, y=height * 1.05, linewidth=2, color='black')

    plt.ylim(0, height*1.25)
    plt.show()


if __name__ == '__main__':

    cache_file = 'nac_quench_cache_n2_n4_n5.csv'

    # Test to use. 'T-test' uses students t-test, everything else will default to Man Whitney U
    test=""
    main(cache_file, test)