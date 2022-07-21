import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from PIL import Image
import numpy as np
import glob
import os
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
from matplotlib.ticker import FormatStrFormatter

# Pandas setup for better display of dataframes
pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 1500)
pd.set_option('max_colwidth', 200)

plt.rcParams['svg.fonttype'] = 'none'

def generate_colormaps(data_df):

    list_of_colormaps = [
        'Blues',
        'Oranges',
        'Greens',
        'Reds',
        'Purples'
    ]
    colormap_dict = {}
    for idx, time in enumerate(data_df['time'].unique()):

        my_cmap = cmx.ScalarMappable(
            norm=mcolors.Normalize(
                vmin=0,
                vmax=data_df.loc[data_df['time'] == time, 'green/blue'].max()
            ),
            cmap=plt.get_cmap(list_of_colormaps[idx])
        )
        colormap_dict[time] = my_cmap

    return colormap_dict


def map_colors(data_row, colormap_dict):

    active_colormap = colormap_dict[data_row['time']]
    data_row['color'] = active_colormap.to_rgba(data_row['statistic'])

    return data_row


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


def estimate_statistic(data_row, estimator = np.mean, ci = 90, n_boot = 1000):
    # https://github.com/mwaskom/seaborn/blob/v0.8.1/seaborn/categorical.py#L2950

    # calculate a statistic for our grouped data -> in this case the mean
    data_row['statistic'] = estimator(data_row['read_out'])

    # Also calculate a confidence interval, in this case using bootstrapping
    if len(data_row['read_out']) < 2:
        data_row['confint_lower'] = np.nan
        data_row['confint_upper'] = np.nan

    else:
        confint = sns.utils.ci(
            sns.algorithms.bootstrap(
                data_row['read_out'],
                n_boot=n_boot
            ),
            ci
        )

        data_row['confint_lower'] = confint[0]
        data_row['confint_upper'] = confint[1]

    return data_row


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


def main(experiments, test):
    # Parameters for analysis:
    green_threshhold, blue_threshhold = 90, 90
    excluded_padding = 100

    cellROX_path = r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk\CellRox'
    summary_file = rf'titration_analysis_summary.tsv'

    if not os.path.exists(summary_file):

        data = dict({
            'date': [],
            'repeat': [],
            'time': [],
            'comment': [],
            'condition': [],
            'duplicate': [],
            'green_pixel': [],
            'blue_pixel': []
        })

        for repeat in experiments:

            for time in repeat:
                green_files = glob.glob(rf"{cellROX_path}\{time}\*\4x_CH2.tif")
                blue_files = glob.glob(rf"{cellROX_path}\{time}\*\4x_CH1.tif")

                for green_file, blue_file in zip(green_files, blue_files):
                    experiment = green_file.split('\\')[8]
                    well = green_file.split('\\')[9]

                    print(f"Evaluating: {experiment.split('_')[4]}, {experiment.split('_')[5]}, {well.split('_')[2]}")

                    data['date'].append(str.join('-', experiment.split('_')[0:3]))
                    data['repeat'].append(experiment.split('_')[4])
                    data['time'].append(experiment.split('_')[5])
                    data['comment'].append(str.join('_', experiment.split('_')[6:]))
                    data['condition'].append(well.split('_')[2])
                    data['duplicate'].append(well.split('_')[2:])

                    green_pixel, blue_pixel = evaluate_images(green_file, blue_file)
                    data['green_pixel'].append(green_pixel)
                    data['blue_pixel'].append(blue_pixel)

        data_df = pd.DataFrame(data)

        # Rename values
        data_df['condition'] = data_df['condition'].str.replace('+', '', regex=False)
        data_df['condition'] = data_df['condition'].str.replace('-', '0')
        data_df['condition'] = pd.to_numeric(data_df['condition'])

        data_df.to_csv(summary_file, sep="\t", index=None)

    else:
        data_df = pd.read_csv(summary_file, sep = '\t')


    # Kick out 25 and 2.5 ng/mL because we only have one repeat for these values
    data_df = data_df.loc[(data_df['condition'] != 25) & (data_df['condition'] != 2.5)]
    # Kick out data for 800 & 1600 ng/mL because we only have one repeat for these values:
    data_df = data_df.loc[(data_df['condition'] != 800) & (data_df['condition'] != 1600)]
    # Kick out n5:
    # data_df = data_df.loc[(data_df['repeat'] != 'n5') & (data_df['repeat'] != 'n2')]
    data_df = data_df.loc[(data_df['time'] != '300min')]

    data_df['condition'] = data_df['condition'].astype(int)

    # normalize the green via blue
    data_df['green/blue'] = data_df['green_pixel']/data_df['blue_pixel']

    data_df['time'] = data_df['time'].apply(lambda x: int(x[:-3]))
    data_df = data_df.sort_values(by = ['condition', 'time'], ascending=[False, True])

    times = pd.factorize(data_df['time'])
    conditions = pd.factorize(data_df['condition'])
    data_df['time_idx'] = np.array(times[0], dtype=float)
    data_df['cond_idx'] = np.array(conditions[0], dtype=float)

    colormaps_dict = generate_colormaps(data_df)
    color_dict = {}
    for index, elem in enumerate(data_df['time_idx'].unique()):
        color_dict[elem] = list(mcolors.TABLEAU_COLORS)[index]


    # Confidence interval

    mean_signal_df = data_df[['green/blue', 'condition', 'time', 'cond_idx', 'time_idx']]

    mean_signal_df = mean_signal_df.groupby(
        by=['condition', 'time', 'cond_idx', 'time_idx'],
        as_index=False
    ).agg(
        read_out=('green/blue', lambda x: list(x))
    )
    mean_signal_df.reset_index(inplace=True)

    mean_signal_df = mean_signal_df.apply(estimate_statistic, axis=1)

    # Visualize the data
    fig = plt.figure(figsize=(6, 4))
    ax15 = fig.add_subplot(3, 3, (1, 8), projection='3d')
    ax15.view_init(10, 300)

    error_width=0.03
    ax15.bar3d(
        x=mean_signal_df['cond_idx'] - (error_width/2),
        y=mean_signal_df['time_idx'] - (error_width/2),
        z=mean_signal_df['confint_lower'],
        dx=error_width,
        dy=error_width,
        dz=mean_signal_df['confint_upper'] - mean_signal_df['confint_lower'],
        alpha=0.6,
        color='gray'
    )

    ax15.scatter(
        xs=data_df['cond_idx'],
        ys=data_df['time_idx'],
        zs=data_df['green/blue'],
        color='black',
        edgecolors='black',
        alpha=1,
    )

    mean_signal_df = mean_signal_df.apply(map_colors, axis = 1, args=(colormaps_dict, ))

    mean_width = 0.5
    ax15.bar3d(
        x=mean_signal_df['cond_idx'] - (mean_width/2),
        y=mean_signal_df['time_idx'] - (mean_width/2),
        z=0,
        dx=mean_width,
        dy=mean_width,
        dz=mean_signal_df['statistic'],
        alpha=0.6,
        color=mean_signal_df['color'],
        edgecolors='black',
        linewidths=0.3,
    )

    # Label all the stuff
    plt.sca(ax15)
    plt.xticks(ticks = list(set(conditions[0])), labels = conditions[1])
    plt.yticks(ticks = list(set(times[0])), labels = times[1])

    ax15.set_ylabel('Time (min)')
    ax15.set_xlabel('PDGF (ng/mL)')
    ax15.set_zlabel('CellROX Green (RU)')
    ax15.set_zbound(lower = 0)

    axes = []
    for idx, time in enumerate([240, 180, 120]):

        ax = fig.add_subplot(3, 3, 3*(idx + 1))
        data_min_df = data_df.loc[data_df['time'] == time]
        order_min = data_min_df['condition'].unique()
        sns.barplot(
            ax=ax,
            x='condition',
            y='green/blue',
            data=data_min_df,
            color=color_dict[3-idx],
            order = order_min,
            saturation=0.6
            )
        sns.stripplot(
            ax=ax,
            x='condition',
            y='green/blue',
            data=data_min_df,
            linewidth=0.3,
            size=6,
            color=color_dict[3-idx],
            order=order_min
        )
        four_hundred = data_min_df.loc[data_min_df['condition'] == 400, 'green/blue'].values
        two_hundred = data_min_df.loc[data_min_df['condition'] == 200, 'green/blue'].values
        one_hundred = data_min_df.loc[data_min_df['condition'] == 100, 'green/blue'].values
        zero = data_min_df.loc[data_min_df['condition'] == 0, 'green/blue'].values

        plt.sca(ax)

        # significance between 400 and 200:
        height = data_min_df['green/blue'].max()
        plt.text(
            x=0.5,
            y=height * 1.1,
            s=f'{calc_significance(four_hundred, two_hundred, test=test)[0]}',
            ha="center",
            va="center",
            fontdict={"size": 10}
        )
        plt.hlines(xmin=0.1, xmax=0.9, y=height*1.05, linewidth=2, color='black')

        # significance between 200 and 0:
        plt.text(
            x=4,
            y=height * 1.1,
            s=f'{calc_significance(two_hundred, zero, test=test)[0]}',
            ha="center",
            va="center",
            fontdict={"size": 10}
        )
        plt.hlines(xmin=1.1, xmax=6.9, y=height * 1.05, linewidth=2, color='black')

        # significance between 200 and 100:
        plt.text(
            x=1.5,
            y=height * 0.95,
            s=f'{calc_significance(two_hundred, one_hundred, test=test)[0]}',
            ha="center",
            va="center",
            fontdict={"size": 10}
        )
        plt.hlines(xmin=1.1, xmax=1.9, y=height * 0.9, linewidth=2, color='black')

        # significance between 400 and 0:
        plt.text(
            x=3.5,
            y=height * 1.25,
            s=f'{calc_significance(four_hundred, zero, test=test)[0]}',
            ha="center",
            va="center",
            fontdict={"size": 10}
        )
        plt.title(
            label = f'After {time} min:',
            fontdict={"size": 10},
            loc='left'
        )
        plt.hlines(xmin=0.1, xmax=6.9, y=height * 1.2, linewidth=2, color='black')
        plt.ylim(0, height*1.3)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f')) # no decimal points

        # Label all the stuff
        ax.set_ylabel('CellROX Green (RU)')
        ax.set_xlabel('')
        ax.set_xticklabels([])

        axes.append(ax)

    axes[-1].set_xticklabels(data_df['condition'].unique())
    axes[-1].set_xlabel('PDGF (ng/mL)')

    plt.show()

if __name__ == '__main__':
    # List of experiments:
    n1 = [  # n1 (24.11.2021)
        '2021_11_24_titration_n1_60min',
        '2021_11_24_titration_n1_120min',
        '2021_11_24_titration_n1_180min',
        # '2021_11_24_titration_n1_450min',
        # '2021_11_24_titration_n1_nextDay'
    ]

    n2 = [  # n2 (22.12.2021)
        # Careful! The highest concentration measured for this experiment was 200 ng/mL!
        # I messed up my calculations and noticed only after the experiment
        # '2021_12_22_titration_n2_60min_wrong_settings',
        # '2021_12_22_titration_n2_120min_wrong_settings',
        # '2021_12_22_titration_n2_180min_wrong_settings',
        '2021_12_22_titration_n2_180min_corr_settings',
        '2021_12_22_titration_n2_240min_corr_settings',
        '2021_12_22_titration_n2_300min_corr_settings'
    ]

    n3 = [  # n3 (08.03.2022)
        '2022_03_08_titration_n3_60min',
        '2022_03_08_titration_n3_120min',
        '2022_03_08_titration_n3_180min',
        '2022_03_08_titration_n3_240min'
    ]

    n4 = [  # n4 (23.03.2022)
        '2022_03_23_titration_n4_60min',
        '2022_03_23_titration_n4_120min',
        '2022_03_23_titration_n4_180min',
        '2022_03_23_titration_n4_240min'
    ]

    n5 = [
        '2022_04_01_titration_n5_60min',
        '2022_04_01_titration_n5_120min',
        '2022_04_01_titration_n5_180min',
        '2022_04_01_titration_n5_240min'
    ]

    n6 = [
        '2022_05_02_titration_n6_60min',
        '2022_05_02_titration_n6_120min',
        '2022_05_02_titration_n6_180min',
        '2022_05_02_titration_n6_240min'
    ]

    n7 = [
        '2022_05_12_titration_n7_60min',
        '2022_05_12_titration_n7_120min',
        '2022_05_12_titration_n7_180min',
        '2022_05_12_titration_n7_240min'
    ]

    # Test to use. 'T-test' uses students t-test, everything else will default to Man Whitney U
    test = ""
    main([n1, n2, n3, n4, n6, n7], test)