
import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

from skimage.feature import match_template
from skimage.transform import rescale, rotate
from scipy import ndimage
from scipy import ndimage as ndi, optimize
import tifffile


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

def rot_ds_to_conf(x, y, ang, up_scale, conf_shape):
    cy, cx = (np.array(conf_shape) - 1) / 2
    # skimage.rotate(img, a) samples output (x,y) from input R(a)·((x,y)-centre)+centre
    t = np.radians(-ang)                 # we called rotate(conf_feat, -ang)
    xr, yr = x - cx, y - cy
    xo = np.cos(t) * xr - np.sin(t) * yr + cx
    yo = np.sin(t) * xr + np.cos(t) * yr + cy
    return xo * up_scale + (up_scale - 1) / 2, yo * up_scale + (up_scale - 1) / 2
def do_overlay(afm, op_rs, shift_coord, mask_bool=True):
    row = int(shift_coord[0])
    col = int(shift_coord[1])

    h, w = afm.shape

    rgb = np.zeros((*op_rs.shape, 3))
    # AFM in red (only inside its window)
    rgb[row:row + h, col:col + w, 0] = afm
    if mask_bool:
        rgb[..., 1] = (op_rs > 0).astype(float)
    else:
        rgb[..., 1] = op_rs
    zeros = np.zeros_like(op_rs)

    r = np.dstack([rgb[..., 0], zeros, zeros]) 
    g = np.dstack([zeros, rgb[..., 1], zeros]) 

    return rgb,r,g

def afm_img_to_confocal_frame(M, afm_channel, d_shape, tif_exp=True, label='none'):
    H, W = afm_channel.shape

    Minv = np.linalg.inv(M)
    Y, X = np.mgrid[:d_shape[0], :d_shape[1]]
    ax_ = Minv[0, 0] * X + Minv[0, 1] * Y + Minv[0, 2]
    ay_ = Minv[1, 0] * X + Minv[1, 1] * Y + Minv[1, 2]
    inside = (ax_ >= -0.5) & (ax_ <= W -
                              0.5) & (ay_ >= -0.5) & (ay_ <= H - 0.5)
    afm_channel_conf = ndi.map_coordinates(
        afm_channel, [ay_, ax_], order=1).astype(np.float32)
    afm_channel_conf[~inside] = 0
    plt.imshow(afm_channel_conf)
    plt.show()
    if tif_exp:
        path = 'output'
        tifffile.imwrite(f"{path}/afm_rot_shear_{label}.tif", afm_channel_conf)
    return afm_channel_conf

def max_rotational_cross_corr(im1, im2, pivot=(0, 0), angle_range=(-5, 5), angle_step=1):
    best_score = 0
    best_angle = 0
    best_pos = [0, 0]
    best = (-np.inf, None)

    # Normalize im1
    # im1 = (im1 - np.mean(im1)) / (np.std(im1) + 1e-8)

    for angle in np.arange(*angle_range, angle_step):
        conf_rot = rotate(im1, -angle, resize=False,
                          preserve_range=True, mode="constant")

        result = match_template(conf_rot.astype(float), im2.astype(float))
        if im2.shape[0] >= conf_rot.shape[0] or im2.shape[1] >= conf_rot.shape[1]:
            continue
        r, c = np.unravel_index(np.argmax(result), result.shape)
        if result[r, c] > best[0]:
            best = (result[r, c], (angle, r, c))
            best_pos[0] = r
            best_pos[1] = c
            best_angle = angle
            best_score = result[r, c]

    return best_angle, best_pos, best_score
