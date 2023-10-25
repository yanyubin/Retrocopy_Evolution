#!/usr/bin/env python3

"""
Author : Yubin Yan <yanyb@nwafu.edu.cn>
Purpose: Get the colocalization frequency of polymorphic retrocopies
"""

import click
import re
import gffutils
from typing import NamedTuple
import pandas as pd


@click.command()
@click.option('-v', '--vcf', help='The retroCNV vcf file')
@click.option('-g', '--gtf', help='The gencode annotation file in gtf/gff3 or db format')


# --------------------------------------------------
def main(vcf: str, gtf: str) -> None:
    """ Get the colocalization frequency of polymorphic retrocopies """

    # create a data frame of retrocopies
    retro_info: list[MyRetro] = []
    with open(vcf, 'r', encoding='utf-8') as retros:
        for retro in retros:
            if retro.startswith('#'):
                continue
            retro_info.append(MyRetro(*get_retro_info(retro)))

    df_retro = pd.DataFrame(retro_info)
    print(df_retro)

    # create a data frame of annotated genes
    gene_info: list[MyGene] = []
    for g in get_annotation(gtf):
        gene_info.append(MyGene(g.seqid, g.start, g.end, *g['gene_name'], g.id[:15]))

    df_gene = pd.DataFrame(gene_info)
    print(df_gene)

    # left join two data frames based on parental gene name
    df = pd.merge(df_retro, df_gene, on='pg_name', how='left')
    print(df)

    # save results to a csv file
    df.to_csv('retroCNV_colocal_indv_final.csv', sep='\t', index=False, encoding='utf-8')


class MyGene(NamedTuple):
    """ Custom tuple to hold gene information """
    pg_chrom: chr
    pg_start: int
    pg_end: int
    pg_name: chr
    pg_id: chr


class MyRetro(NamedTuple):
    """ Custom tuple to hold retroCNV information """
    chrom: chr
    coor: int
    precision: chr
    pg_name: chr
    counts: int


def get_retro_info(record: str) -> MyRetro:
    """ Get information of retroCNVs """

    record = record.strip().split('\t')
    chrom, coor = record[0], record[1]
    pg_name = re.search('PG=(.*);PGTYPE', record[7]).group(1)

    counts = len(re.findall('0/1', ' '.join(record[9:]))) + \
             len(re.findall('1/1', ' '.join(record[9:]))) * 1

    if 'IMPRECISE' in record[7]:
        precision = 'imprecise'
    else:
        precision = 'precise'

    return MyRetro(chrom, coor, precision, pg_name, counts)


def get_annotation(gtf: str):
    """ Get gene information form a gtf file """

    if gtf.endswith(('gtf', 'gtf.gz', 'gff3', 'gff3.gz')):
        db = gffutils.create_db(gtf, dbfn='gencode.v32.annotation.gtf.db', force=False,
                                keep_order=True, merge_strategy='merge',
                                sort_attribute_values=True)
    elif gtf.endswith('db'):
        db = gffutils.FeatureDB(gtf)
    else:
        raise ValueError("This is not a supported annotation file")

    return db.all_features(featuretype='gene')


# --------------------------------------------------
if __name__ == '__main__':
    main()
