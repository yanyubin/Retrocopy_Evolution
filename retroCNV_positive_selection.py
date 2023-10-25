#!/usr/bin/env python3

"""
Author : Yubin Yan <yanyb@nwafu.edu.cn>
Purpose: detect positive selection of retroCNVs
"""

import click
import re
import os
import shutil
import subprocess
from pybedtools import BedTool
import pandas as pd
import math
import glob

@click.command()
@click.option('-v1', '--vcf1', show_default=True, help='A vcf file of retroCNVs')
@click.option('-v2', '--vcf2', show_default=True, help='A vcf file of population SNPs')
@click.option('-p', '--pop', show_default=True, help='A text file listing the individuals of a population')


# --------------------------------------------------
def main(vcf1, vcf2, pop) -> None:
    """ Detect positive selection of retroCNVs """

    for site in get_retro_site(vcf1):
        coor, carrier, noncarrier, N_carrier, _ = site[0], site[1], site[2], site[3], site[4]

        dir = pop + '/' + coor
        # create a folder for each site
        if os.path.exists(dir):
            choice = input("Directory exists! Type 'y' to overwrite, otherwise terminate: ")
            if choice == 'y':
                shutil.rmtree(dir)
            else:
                break
        os.mkdir(dir)

        # write retroCNV carriers and noncarriers to text files
        with open(dir + '/carriers.list', 'w') as f1:
            f1.write('\n'.join(carrier))
        with open(dir + '/noncarriers.list', 'w') as f2:
            f2.write('\n'.join(noncarrier))

        chr = coor.strip().split('_')[0][3:]
        pos = coor.strip().split('_')[1]

        # observed
        run_vcftools(vcf2, dir, chr, pos, '')

        # simulated
        sim_sites = BedTool().random(l = 1, n = 200, g = './hg19_autosome_size.txt', seed = 2023)
        for idx, sim_site in enumerate(sim_sites):
            run_vcftools(vcf2, dir, sim_site[0][3:], sim_site[2], str(idx + 1))

        tajima = []
        for f in glob.glob(dir + '/carriers*.windowed.pi'):  # not sorted!
            win_pi = pd.read_csv(f, sep = '\t')
            tajima.append(tajimaD(win_pi, 2 * N_carrier))
        df = pd.DataFrame({'tajimaD' : tajima})
        df.to_csv(dir + '/tajimaD_carrier.csv', sep='\t', index=False, encoding='utf-8')


# --------------------------------------------------
def tajimaD(win_pi, N):
    """win_pi: the output of vcftools window pi;
       N: the number of alleles"""

    S = win_pi['N_VARIANTS'].sum()
    pi = win_pi['PI'].sum() * 100000

    if S >= 10 and N >= 5:
        a1 = 0
        for i in range(1, N - 1):
            a1 += 1 / i
        a2 = 0
        for i in range(1, N - 1):
            a2 += 1 / (i ** 2)
        b1 = (N + 1) / (3 * (N - 1))
        b2 = (2 * (N ** 2 + N + 3)) / (9 * N * (N - 1))

        c1 = b1 - 1 / a1
        c2 = b2 - (N + 2) / (a1 * N) + a2 / (a1 ** 2)
        e1 = c1 / a1
        e2 = c2 / (a1 ** 2 + a2)

        D = (pi - S / a1) / math.sqrt(e1 * S + e2 * S * (S - 1))
    else:
        D = 'NA'

    return D


# --------------------------------------------------
def run_vcftools(vcf2, dir, chr, pos, idx):
    # Nucleotide diversity
    subprocess.run(['vcftools', '--vcf', vcf2, '--keep', dir + '/carriers.list', '--chr', chr,
                    '--from-bp', str(int(pos) - 49999), '--to-bp', str(int(pos) + 50000),
                    '--window-pi', '100000', '--out', dir + '/carriers' + idx])
    if os.stat(dir + '/noncarriers.list').st_size > 0:
        subprocess.run(['vcftools', '--vcf', vcf2, '--keep', dir + '/noncarriers.list', '--chr', chr,
                        '--from-bp', str(int(pos) - 49999), '--to-bp', str(int(pos) + 50000),
                        '--window-pi', '100000', '--out', dir + '/noncarriers' + idx])
    # Tajima's D
    """ subprocess.run(['vcftools', '--vcf', vcf2, '--keep', dir + '/carriers.list', '--chr', chr,
                    '--from-bp', str(int(pos) - 49999), '--to-bp', str(int(pos) + 50000),
                    '--TajimaD', '100000', '--out', dir + '/carriers' + idx]) """
    # LD analysis
    subprocess.run(['vcftools', '--vcf', vcf2, '--keep', dir + '/carriers.list', '--chr', chr,
                    '--from-bp', str(int(pos) - 4999), '--to-bp', str(int(pos) + 5000),
                    '--geno-r2', '--ld-window-bp', '10000', '--out', dir + '/carriers' + idx])


# --------------------------------------------------
def get_retro_site(vcf1: str):
    sites = []
    with open(vcf1, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                indv = line.strip().split('\t')
            else:
                record = line.strip().split('\t')
                chrom, coor = record[0], record[1]
                N_carrier = len(re.findall('0/1', ' '.join(record[9:]))) + \
                         len(re.findall('1/1', ' '.join(record[9:]))) * 1
                N_noncarrier = len(re.findall('0/0', ' '.join(record[9:])))
                pop_size = N_carrier + N_noncarrier

                carriers = []
                noncarriers = []
                for idx, x in enumerate(record):
                    if x.startswith(('0/1', '1/1')):
                        carriers.append(indv[idx])
                    elif x.startswith('0/0'):
                        noncarriers.append(indv[idx])

                if N_carrier >= 0.5 * pop_size and N_carrier >= 10:
                    sites.append(((chrom + '_' + coor), carriers, noncarriers, N_carrier, N_noncarrier))
    return sites

# --------------------------------------------------
if __name__ == '__main__':
    main()
