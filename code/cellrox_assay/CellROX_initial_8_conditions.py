# 28.10.21

import os.path
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from PIL import Image
import numpy as np
import glob

plt.rcParams['svg.fonttype'] = 'none'

'''
# List of experiments:

# The first ones:
    # '20211028_1.5h' -> A1 (+++ on matrix) is high, everything basically not visible
    # '20211028_4h' -> A1 (+++ on matrix) is high, everything basically not visible
'''


def get_data(experiment_folder,
             threshhold_g_matrix,
             threshhold_g_plastic,
             threshhold_b_matrix,
             threshhold_b_plastic,
             exclude_padding):
    # Get file names of the green and blue channel
    green_files = glob.glob(
        r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk/CellRox/' + experiment_folder + r'/*/4x_CH2.tif')
    blue_files = glob.glob(
        r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk\CellRox/' + experiment_folder + r'/*/4x_CH1.tif')

    # Prepare lists for sample names and numbers of green/blue pixels
    sample_name_plastic = []
    green_pixels_plastic = []
    blue_pixels_plastic = []
    g_b_plastic = []
    plastic_before = [[0]]
    plastic_after = [[0]]

    sample_name_matrix = []
    green_pixels_matrix = []
    blue_pixels_matrix = []
    g_b_matrix = []
    matrix_before = [[0]]
    matrix_after = [[0]]

    for green_file, blue_file in zip(green_files, blue_files):

        # Open images
        print('Evaluating:', green_file, blue_file)
        green_im = Image.open(green_file)
        blue_im = Image.open(blue_file)

        # Load into arrays
        green_ar = np.asarray(green_im)
        blue_ar = np.asarray(blue_im)

        # Extract the sample name from file name
        if 'atrix' in (sample_name := green_file.split("\\")[-2]):

            sample_name_matrix.append(sample_name.split('_')[2:])

            if '+++' in sample_name or '++200' in sample_name:
                gr_areas = np.where(green_ar[:, :, 1] > threshhold_g_matrix)
                matrix_before = np.copy(green_ar)
                matrix_after = np.copy(green_ar)
                matrix_after[gr_areas] = [200, 0, 0]
                matrix_after[0:exclude_padding, :, :] += 100
                matrix_after[-exclude_padding:-1, :, :] += 100
                matrix_after[exclude_padding:-exclude_padding, 0:exclude_padding, :] += 100
                matrix_after[exclude_padding:-exclude_padding, -exclude_padding:-1, :] += 100

            threshhold_g = threshhold_g_matrix
            threshhold_b = threshhold_b_matrix

            # Count pixel with a higher value in the green/blue channel than the mean
            green_pixel = (green_ar[exclude_padding:-exclude_padding, exclude_padding:-exclude_padding,
                           :] > threshhold_g).sum()
            blue_pixel = (blue_ar[exclude_padding:-exclude_padding, exclude_padding:-exclude_padding,
                          :] > threshhold_b).sum()

            # Append values to lists outside loop
            green_pixels_matrix.append(green_pixel)
            blue_pixels_matrix.append(blue_pixel)
            g_b_matrix.append(green_pixel / blue_pixel)

        elif 'lastic' in (sample_name := green_file.split("\\")[-2]):
            sample_name_plastic.append(sample_name.split('_')[2:])

            if '+++' in sample_name:
                gr_areas = np.where(green_ar[:, :, 1] > threshhold_g_plastic)
                plastic_before = np.copy(green_ar)
                plastic_after = np.copy(green_ar)
                plastic_after[gr_areas] = [200, 0, 0]
                plastic_after[0:exclude_padding, :, :] += 100
                plastic_after[-exclude_padding:-1, :, :] += 100
                plastic_after[exclude_padding:-exclude_padding, 0:exclude_padding, :] += 100
                plastic_after[exclude_padding:-exclude_padding, -exclude_padding:-1, :] += 100

            threshhold_g = threshhold_g_plastic
            threshhold_b = threshhold_b_plastic

            # Count pixel with a higher value in the green/blue channel than the mean
            green_pixel = (green_ar[exclude_padding:-exclude_padding,
                           exclude_padding:-exclude_padding, :] > threshhold_g).sum()
            blue_pixel = (blue_ar[exclude_padding:-exclude_padding,
                          exclude_padding:-exclude_padding, :] > threshhold_b).sum()

            # Append values to lists outside loop
            green_pixels_plastic.append(green_pixel)
            blue_pixels_plastic.append(blue_pixel)
            g_b_plastic.append(green_pixel / blue_pixel)

        green_im.close()
        blue_im.close()

    plastic_df = pd.DataFrame(data={
        'name': sample_name_plastic,
        'green': green_pixels_plastic,
        'blue': blue_pixels_plastic,
        'green/blue': g_b_plastic
    })

    matrix_df = pd.DataFrame(data={
        'name': sample_name_matrix,
        'green': green_pixels_matrix,
        'blue': blue_pixels_matrix,
        'green/blue': g_b_matrix
    })

    return plastic_df, plastic_before, plastic_after, matrix_df, matrix_before, matrix_after

def main():
    experiment_folder = '2021_10_28_1.5h'
    threshhold_g_matrix = 120
    threshhold_g_plastic = 120

    threshhold_b_matrix = 80
    threshhold_b_plastic = 90

    exclude_padding = 100

    plastic_cache_file = f'{experiment_folder}_plastic.tsv'
    matrix_cache_file = f'{experiment_folder}_matrix.tsv'
    if os.path.isfile(plastic_cache_file) and os.path.isfile(matrix_cache_file):
        plastic_df = pd.read_csv(plastic_cache_file, sep='\t')
        matrix_df = pd.read_csv(matrix_cache_file, sep='\t')

    else:
        plastic_df, plastic_before, plastic_after, matrix_df, matrix_before, matrix_after = get_data(experiment_folder,
                                                                                                     threshhold_g_matrix,
                                                                                                     threshhold_g_plastic,
                                                                                                     threshhold_b_matrix,
                                                                                                     threshhold_b_plastic)

        plastic_df.to_csv(
            f'{experiment_folder}_plastic.tsv',
            index=False,
            sep='\t'
        )
        matrix_df.to_csv(
            f'{experiment_folder}_matrix.tsv',
            index=False,
            sep='\t'
        )

    fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(6, 4))

    # Matrix
    bad_indices = matrix_df.index.isin([0,4])
    matrix_df = matrix_df[~bad_indices]

    plt.sca(ax1)
    sns.barplot(
        ax=ax1,
        x='name',
        y='green/blue',
        data=matrix_df,
        palette='gray'
    )
    ax1.set_title('green/blue on Col I matrix')
    ax1.set(ylim=(0, 5))
    ax1.set(xlabel=None)

    print(matrix_df)

    # Plastic
    plt.sca(ax2)
    sns.barplot(
        ax=ax2,
        x='name',
        y='green/blue',
        data=plastic_df,
        color='forestgreen'
    )
    ax2.set_title('green/blue on plastic')
    ax2.set(ylim=(0, 5))
    ax2.set(xlabel='Conditions')
    print(plastic_df)

    plt.show()

if __name__ == '__main__':
    main()

