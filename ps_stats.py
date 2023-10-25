#!/usr/bin/env python3

"""
Author : Yubin Yan <yanyb@nwafu.edu.cn>
Purpose: Get statistics of positive selection of retroCNVs
"""

import click
import pandas as pd
import numpy as np
from scipy.stats import ranksums
import glob
from retroCNV_positive_selection import get_retro_site
from retroCNV_positive_selection import tajimaD

@click.command()
@click.option('-v1', '--vcf1', show_default=True, help='A vcf file of retroCNVs')
@click.option('-p', '--pop', show_default=True, help='A text file listing the individuals of a population')


# --------------------------------------------------
def main(vcf1, pop) -> None:
    """ Get statistics of positive selection of retroCNVs """

    out_df = pd.DataFrame(columns=['chr', 'pos', 'N_carriers', 'N_noncarriers',
                                   'true_pi', 'sim_mean_pi', 'sim_sd_pi', 'p_value_pi',
                                   'true_D', 'sim_mean_D', 'sim_sd_D', 'p_value_D',
                                   'true_LD', 'sim_mean_LD', 'sim_sd_LD', 'p_value_LD'])
    for idx, site in enumerate(get_retro_site(vcf1)):
        coor, _, _, N_carriers, N_noncarriers = site[0], site[1], site[2], site[3], site[4]

        dir = pop + '/' + coor

        chr = coor.strip().split('_')[0][3:]
        pos = coor.strip().split('_')[1]

        # LD
        ld = get_df(dir + '/carriers*.geno.ld', type='ld')
        ld = ld.replace('carriers', 'sim', regex=True)

        # Tajima's D
        values = []
        names = []
        for f in glob.glob(dir + '/carriers*.windowed.pi'):  # not sorted!
            win_pi = pd.read_csv(f, sep = '\t')
            values.append(tajimaD(win_pi, 2 * N_carriers))
            names.append(f.split('/')[2].split('.')[0])
        tajima = pd.DataFrame({'ID' : names, 'tajimaD' : values})
        tajima = tajima.replace('carriers', 'sim', regex=True)
        
        df = ld.merge(tajima, on='ID', how='left').sort_values(by='ID').reset_index(drop=True)
        df = df.replace('NA', np.nan)
        
        _, p_tajima = ranksums(df['tajimaD'][0], df['tajimaD'][1:], alternative='less', nan_policy='omit')
        _, p_LD = ranksums(df['LD'][0], df['LD'][1:], alternative='greater', nan_policy='omit')

        true_D, sim_mean_D, sim_sd_D = df['tajimaD'][0], np.mean(df['tajimaD'][1:]), np.std(df['tajimaD'][1:])
        true_LD, sim_mean_LD, sim_sd_LD = df['LD'][0], np.mean(df['LD'][1:]), np.std(df['LD'][1:])
        
        # nucleotide diversity
        if N_noncarriers >= 5:
            carriers_pi = get_df(dir + '/carriers*.windowed.pi', 'pi')
            carriers_pi = carriers_pi.replace('carriers', 'sim', regex=True)

            noncarriers_pi = get_df(dir + '/noncarriers*.windowed.pi', 'pi')
            noncarriers_pi = noncarriers_pi.replace('noncarriers', 'sim', regex=True)

            pi = carriers_pi.merge(noncarriers_pi, on='ID', how='left')
            pi['ratio'] = pi['PI_y'] / pi['PI_x']
            pi = pi.sort_values(by='ID').reset_index(drop=True)

            _, p_pi = ranksums(pi['ratio'][0], pi['ratio'][1:], alternative='greater', nan_policy='omit')
            
            true_pi, sim_mean_pi, sim_sd_pi = pi['ratio'][0], np.mean(pi['ratio'][1:]), np.std(pi['ratio'][1:])
            
        else:
            true_pi, sim_mean_pi, sim_sd_pi, p_pi = 'NA', 'NA', 'NA', 'NA'


        print([chr, pos, N_carriers, N_noncarriers, true_pi, sim_mean_pi, sim_sd_pi, p_pi,
                           true_D, sim_mean_D, sim_sd_D, p_tajima, true_LD, sim_mean_LD, sim_sd_LD, p_LD])

        out_df.loc[idx] = [chr, pos, N_carriers, N_noncarriers, true_pi, sim_mean_pi, sim_sd_pi, p_pi,
                           true_D, sim_mean_D, sim_sd_D, p_tajima, true_LD, sim_mean_LD, sim_sd_LD, p_LD]
    print(out_df)
    out_df.to_csv(pop + '/ps_stats.csv', sep='\t', index=False, encoding='utf-8')


# --------------------------------------------------
def get_df(files, type):
    types = ['pi', 'ld']
    if type not in types:
        raise ValueError("Invalid sim type. Expected one of: %s" % types)

    values = []
    names = []
    for f in glob.glob(files):  # not sorted!
        data = pd.read_csv(f, sep = '\t')
        names.append(f.split('/')[2].split('.')[0])
        if type == 'pi':
            values.append(data['PI'].sum())
        elif type == 'ld':
            values.append(data['R^2'].mean())
    if type == 'pi':
        df = pd.DataFrame({'ID' : names, 'PI' : values})
    elif type == 'ld':
        df = pd.DataFrame({'ID' : names, 'LD' : values})

    return df


# --------------------------------------------------
if __name__ == '__main__':
    main()
