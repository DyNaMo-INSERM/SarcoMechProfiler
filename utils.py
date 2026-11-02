
import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

color_dict = {'ctrl': 'grey', 'F0262': 'tab:green', 'luminy': 'tab:blue', 'pevk': 'tab:orange',
              'ctrl_larve': 'grey', 'pevk_larva': 'tab:orange', 'sls_short_larva': 'tab:olive',
              'pevk_larva_2026': 'tab:orange', 'ctrl_larve_2026': 'grey', 'sls_short_larva_2026': 'tab:olive'}

E0_min, E0_max = 4, 4.8
beta_min, beta_max = 0.05, 0.14
mask_range = [-0.06, 1.06]


def find_glob_peak(arc_L_arr, high):

    high_pos = np.ones(2)*-1
    for i in range(2):
        tar_high = high[i]
        high_pos[i] = np.abs(arc_L_arr - tar_high).argmin()

    return high_pos


def randomize_sacs(df_all, seed_num=42):
    N_sac = len(df_all['file_sarc'].unique())
    df_random = pd.DataFrame()
    np.random.seed(seed_num)  # fixed seed
    random_binary = np.random.randint(0, 2, size=N_sac)
    for i, (name, df_i) in enumerate(df_all.groupby('file_sarc')):

        df_temp = df_i
        if random_binary[i]:
            df_temp['norm_pk_pk_distance'] = 1 - \
                df_temp['norm_pk_pk_distance']
        df_random = pd.concat([df_random, df_temp], ignore_index=True)
    df_random['norm_pk_pk_distance_rounded'] = df_random['norm_pk_pk_distance'].round(2)
    return df_random

def do_avg_line_profile(df_all_new, coi='E_0', color_i='grey', label='control', ax=None, mask_range=None,bin_percentage=1):
    if mask_range is not None:
        mask_crop = (df_all_new['norm_pk_pk_distance'] > mask_range[0]) & (
            df_all_new['norm_pk_pk_distance'] < mask_range[1])
        df_all_new = df_all_new[mask_crop]
    
    bins = np.arange(mask_range[0], mask_range[1], bin_percentage/100.0)
    df_all_new["x_bin"] = pd.cut(
        df_all_new["norm_pk_pk_distance"], bins=bins)
    
    
    if ax is None:
        plt.figure(figsize=(12, 8))
        ax = plt.gca()

    # Group by bins and compute mean of x and y
    df_mean_std = df_all_new.groupby("x_bin").agg(
        x_mean=("norm_pk_pk_distance_rounded", "mean"),
        y_mean=(coi, "mean"),
        y_std=(coi, "std")

    ).reset_index()

    df_mean_std['sem_YS_myo'] = df_mean_std['y_std'] / \
        np.sqrt(len(df_all_new['date_myo_num'].unique()))

    error_sem = df_mean_std['sem_YS_myo'].to_numpy()
    x = df_mean_std['x_mean']
    y_e = df_mean_std['y_mean'].to_numpy()

    ax.plot(x, y_e, 'k-', color=color_i,
            linewidth=4.5, zorder=2, label=label)

    ax.fill_between(x, y_e-error_sem, y_e+error_sem,
                    alpha=0.3, color=color_i, zorder=2)
    #ax.grid()
    ax.set_xlabel('Normalized sarcomere length')
    ax.set_ylim([E0_min, E0_max] if coi == 'E_0' else [beta_min, beta_max])
    ylabel = r'$\log_{10}[E_0 \, (\mathrm{Pa})]$' if coi == 'E_0' else r'Fluidity index ($\beta$)'

    ax.set_ylabel(ylabel)