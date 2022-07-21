from matplotlib import pyplot as plt
from PIL import Image
import numpy as np
import glob

plt.rcParams['svg.fonttype'] = 'none'

n1 = [ # n1 (24.11.2021)
    '2021_11_24_titration_n1_60min',
    '2021_11_24_titration_n1_120min',
    '2021_11_24_titration_n1_180min',
    # '2021_11_24_titration_n1_450min',
    # '2021_11_24_titration_n1_nextDay'
]

n2 = [ # n2 (22.12.2021)
    # Careful! The highest concentration measured for this experiment was 200 ng/uL!
    # I messed up my calculations and noticed only after the experiment
    # '2021_12_22_titration_n2_60min_wrong_settings',
    # '2021_12_22_titration_n2_120min_wrong_settings',
    # '2021_12_22_titration_n2_180min_wrong_settings',
    '2021_12_22_titration_n2_180min_corr_settings',
    '2021_12_22_titration_n2_240min_corr_settings',
    '2021_12_22_titration_n2_300min_corr_settings'
]

n3 = [ # n3 (08.03.2022)
    '2022_03_08_titration_n3_60min',
    '2022_03_08_titration_n3_120min',
    '2022_03_08_titration_n3_180min',
    '2022_03_08_titration_n3_240min'
]

n4 = [ # n4 (23.03.2022)
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

folders_to_visualize = n4
save_file_name = f'titrations_wells_{folders_to_visualize}.svg'

time_points = [entry.split('_')[-1] for entry in folders_to_visualize]

files = glob.glob(
        r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk/CellRox/' + folders_to_visualize[0] + r'/*/4x_CH2.tif')

concentations = [entry.split('\\')[-2].split('_')[-1] for entry in files]

fig1, axes1 = plt.subplots(len(time_points), len(concentations), figsize=(15, 8))

for i, folder in enumerate(folders_to_visualize):

    green_files = glob.glob(
        r'Z:\ICG\Labor\Geraete\1_Keyence\Lisa_Paurat\Torben Falk/CellRox/' + folder + r'/*/4x_CH2.tif'
    )

    for j, file in enumerate(green_files):

        green_im = Image.open(file)
        green_ar = np.asarray(green_im)
        axes1[i, j].imshow(green_ar, interpolation='nearest')
        axes1[i, j].axis('off')

pad=5
for ax, col in zip(axes1[0], concentations):
    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='center', va='baseline')

for ax, row in zip(axes1[:,0], time_points):
    ax.annotate(row, xy=(0, 0.5), xytext=(-pad, 0),
                xycoords='axes fraction', textcoords='offset points',
                size='large', ha='right', va='center', rotation='vertical')

plt.tight_layout()
plt.savefig(save_file_name)
