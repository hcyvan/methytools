import os.path
import sys
import argparse
import numpy as np
import pandas as pd
import time
import pysam
from pathlib import Path
from .utils import get_param_from_doc


def add_percentile(dmc_file, matrix_bed, target_samples, background_samples=None, dmc_percentile_file=None,
                   minimum_sample_keep_ratio=0.8,
                   verbose=True):
    """
    Add percentile columns the DMCs identified by MOABS:mcomp.

    :param dmc_file: the dmc files generate by mcomp.
    :param matrix_bed: the bed format methylation matrix. This file should bed compressed by bgzip and index in csi format.
    :param target_samples: the names of samples in target group. This names should appear in the header of *matrix_bed*.
    :param background_samples: the names of samples in background group. If not set, samples not in target_samples are considered background_samples.
    :param dmc_percentile_file: Added percentile column to dmc file. If not set, it will print to stdout.
    :param minimum_sample_keep_ratio: The minimum proportion of each group (target group) that must be retained.
                                    A sample will be removed if its value is -1.
    :param verbose: whether to print log to stderr
    :return:
    """
    matrix_bed = Path(matrix_bed)
    matrix_bed_csi = matrix_bed.with_suffix('.gz.csi')
    dmc = pd.read_csv(dmc_file, sep='\t')
    if dmc_percentile_file:
        fo = open(dmc_percentile_file, 'w')
    else:
        fo = sys.stdout
    tabix_file = pysam.TabixFile(str(matrix_bed), index=matrix_bed_csi)
    header = tabix_file.header[0]
    header_items = header.strip().split('\t')
    samples = pd.Series(header_items[3:])
    dmc_new_header = dmc.columns.to_list()
    dmc_new_header += ['percentile']
    # The percentile of low methylated group
    percentiles = [i / 100 for i in reversed(range(0, 101, 1))]
    fo.write('\t'.join(dmc_new_header) + "\n")
    i = 0  # the number of DMCs processed
    total = dmc.shape[0]  # the total number of DMCs
    start = time.time()
    for item in dmc.values:
        fold_change = item[10]
        iterator = tabix_file.fetch(reference=item[0], start=item[1], end=item[2])
        try:
            line = next(iterator)
        except StopIteration:
            continue
        line = line.strip().split('\t')
        full = pd.Series([float(x) for x in line[3:]])
        full[full == -1] = np.nan
        target = full[samples.isin(target_samples)]
        if background_samples is None:
            background = full[~samples.isin(target_samples)]
        else:
            background = full[samples.isin(background_samples)]
        target_keep = target[~target.isna()]
        background_keep = background[~background.isna()]
        if len(target) == 0 or len(background) == 0:
            raise Exception('No samples matched')
        if len(target_keep) / len(target) < minimum_sample_keep_ratio or len(background_keep) / len(
                background) < minimum_sample_keep_ratio:
            continue
        for percentile in percentiles:
            if fold_change < 0:
                left = target_keep.quantile(percentile)
                right = background_keep.quantile(1 - percentile)
            else:
                left = background_keep.quantile(percentile)
                right = target_keep.quantile(1 - percentile)
            if left < right:
                fo.write('\t'.join([str(x) for x in item] + [str(percentile)]) + '\n')
                break
        i += 1
        if i % 10000 == 0:
            passed_seconds = round(time.time() - start)
            processed_r = round(i / total, 3)
            if verbose:
                sys.stderr.write(
                    f'Processed {i}/{total}={processed_r}\tpassed {passed_seconds}s\n')
    fo.close()


def reverse_dmc_file(input_file, output_file=None):
    """
    Modify the orientation of the group comparison
    :param input_file: the input dmc file
    :param output_file: the output dmc file
    :return:
    """
    if output_file:
        fo = open(output_file, 'w')
    else:
        fo = sys.stdout
    with open(input_file, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                fo.write(line)
            else:
                items = line.strip().split('\t')
                difCI = ','.join([str(-float(x)) for x in items[11].split(',')])
                nominalDif = str(-float(items[9]))
                credibleDif = str(-float(items[10]))
                dmcClass = "nc"
                if float(credibleDif) > 0:
                    dmcClass = 'strongHyper'
                elif float(nominalDif) < 0:
                    dmcClass = 'strongHypo'
                items_rev = items[0:3] + items[6:9] + items[3:6] + [nominalDif, credibleDif, difCI, items[12],
                                                                    items[13],
                                                                    dmcClass]
                fo.write('\t'.join(items_rev) + '\n')

    fo.close()


def dmc2dmr(dmc_file, dist=200, minimum=3, dmr_file=None):
    """
    :param dmc_file: the dmc file with header "chrom, start, end, class"
    :param dist: minimum dmc in the dmr
    :param minimum: maximum distance between two dmc
    :param dmr_file:
    :return:
    """
    dmr_chrom = None
    dmr_start = None
    dmr_end = None
    dmr_class = None
    dmr_cpg_num = 0

    if dmr_file:
        out = open(dmr_file, 'w')
    else:
        out = sys.stdout

    for line in open(dmc_file, 'r'):
        dmc = line.strip().split('\t')
        if dmr_cpg_num == 0:
            dmr_chrom = dmc[0]
            dmr_start = dmc[1]
            dmr_end = dmc[2]
            dmr_class = dmc[3]
            dmr_cpg_num = 1
        else:
            if dmr_chrom == dmc[0] and int(dmc[2]) - int(dmr_end) <= dist and dmr_class == dmc[3]:
                dmr_end = dmc[2]
                dmr_cpg_num += 1
            else:
                if dmr_cpg_num >= minimum:
                    out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(dmr_chrom, dmr_start, dmr_end, dmr_class, dmr_cpg_num,
                                                                int(dmr_end) - int(dmr_start)))
                dmr_chrom = dmc[0]
                dmr_start = dmc[1]
                dmr_end = dmc[2]
                dmr_class = dmc[3]
                dmr_cpg_num = 1


def get_args():
    parser = argparse.ArgumentParser(
        description='This tool is employed to further process the results generated by MOABS:mcomp')
    subparsers = parser.add_subparsers(dest='sub', required=True, title='command', description='The available commands',
                                       help='select a sub command to use')
    # ==================================================================================================================
    parser_percentile = subparsers.add_parser('percentile',
                                              help='This command is employed to add percentile columns to the DMC files')
    parser_percentile.add_argument('-i', '--input', required=True,
                                   help=get_param_from_doc('dmc_file', add_percentile))
    parser_percentile.add_argument('-o', '--output', help=get_param_from_doc('dmc_file_filtered', add_percentile))
    parser_percentile.add_argument('-m', '--matrix-bed', required=True,
                                   help=get_param_from_doc('matrix_bed', add_percentile))
    parser_percentile.add_argument('-t', '--target-samples', required=True,
                                   help=get_param_from_doc('target_samples',
                                                           add_percentile) + " These sample names support two format. 1) The samples should split by ',', such as: sample1,sample2,sample3. 2) Store in a file, each line is a sample name")
    parser_percentile.add_argument('-b', '--background-samples',
                                   help=get_param_from_doc('background_samples',
                                                           add_percentile) + " The input format refer to target_sample")
    parser_percentile.add_argument('-k', '--minimum-keep', default=0.8, type=float,
                                   help=get_param_from_doc('minimum_sample_keep_ratio', add_percentile))
    parser_percentile.add_argument('-v', '--verbose', default=True,
                                   help=get_param_from_doc('verbose', add_percentile))
    # ==================================================================================================================
    parser_reverse = subparsers.add_parser('reverse', help='Modify the orientation of the group comparison')
    parser_reverse.add_argument('-i', '--input', required=True, help=get_param_from_doc('input_file', reverse_dmc_file))
    parser_reverse.add_argument('-o', '--output', help=get_param_from_doc('output_file', reverse_dmc_file))
    # ==================================================================================================================
    parser_dmc2dmr = subparsers.add_parser('dmc2dmr', help='Merge DMCs into DMRs')
    parser_dmc2dmr.add_argument('-i', '--input', required=True,
                                help='dmc bed file. "chrom, start, end, class" is necessary')
    parser_dmc2dmr.add_argument('-m', '--minimum', default=3, type=int, help='minimum dmc in the dmr')
    parser_dmc2dmr.add_argument('-d', '--dist', default=200, type=int, help='maximum distance between two dmc')
    parser_dmc2dmr.add_argument('-o', '--output',
                                help='dmr bed file. header: "chrom, start, end, class, dmc_number, length"')

    return parser.parse_args()


def main():
    args = get_args()
    if args.sub == 'dmc2dmr':
        dmc2dmr(args.input, args.dist, args.minimum, args.output)
    elif args.sub == 'percentile':
        if os.path.exists(args.target_samples):
            with open(args.target_samples) as f:
                target_samples = [x.strip() for x in f.readlines()]
        else:
            target_samples = args.target_samples.split(',')
        background_samples = args.background_samples
        if background_samples is not None:
            if os.path.exists(background_samples):
                with open(background_samples) as f:
                    background_samples = [x.strip() for x in f.readlines()]
            else:
                background_samples = background_samples.split(',')

        add_percentile(args.input, args.matrix_bed, target_samples, background_samples, args.output, args.minimum_keep,
                       args.verbose)
    elif args.sub == 'reverse':
        reverse_dmc_file(args.input, args.output)
    else:
        args.print_help()


if __name__ == '__main__':
    main()
