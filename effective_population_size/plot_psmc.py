#!/usr/bin/env python
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import glob
import sys
import argparse


def parse_psmc_file(file_path, target_iteration, last_interval):
    """
    Parses a PSMC output file
    :param file_path: str, path to file
    :param target_iteration: int, target iteration to parse
    :return: (np.array, np.array, float), (time in 2N0, relative population size, theta_0),
            returns (int) (max iteration) if target iteration > max iteration
    """
    # see https://github.com/lh3/psmc for meaning of fields
    t = []
    lambda_k = []
    theta_0 = 0
    with open(file_path, 'r') as infile:
        for line in infile:
            if line.startswith('RD'):
                current_iteration = int(line.strip().split('\t')[1])
            elif line.startswith('TR') and current_iteration == target_iteration:
                theta_0 = float(line.strip().split('\t')[1])
            elif line.startswith('RS') and current_iteration == target_iteration:
                line = line.strip().split('\t')
                t.append(float(line[2]))
                lambda_k.append(float(line[3]))
    infile.close()
    t = np.array(t)
    lambda_k = np.array(lambda_k)
    try:
        # Truncate last atomic interval
        if not last_interval:
            lambda_last = lambda_k[-1]
            t = t[np.where(lambda_k != lambda_last)[0]]
            lambda_k = lambda_k[np.where(lambda_k != lambda_last)[0]]
        return (t, lambda_k, theta_0)
    except IndexError:
        return current_iteration


def scale_results_to_generations(t_k, lambda_k, theta_0, mu, binsize):
    """
    Scale results by N0 to express results in terms of generations
    :param t_k: np.array, TMRCA in 2N0
    :param lambda_k: np.array, relative population size
    :param theta_0: float, scaled mutation rate
    :param mu: float, mutation rate per nucleotide
    :param binsize: int, binsize used to prepare psmcfa
    :return: (np.array, np.array, float), (Time scaled to generations, effective population size, N0)
    """
    # see https://github.com/lh3/psmc
    if not mu:
        # 4mu cancel out
        mu = 1 / 4
    N0 = theta_0 / (4 * mu * binsize)
    # scale to generations
    TMRCA = 2 * N0 * t_k
    # scale relative population size to effective population size
    effective_population_size = N0 * lambda_k
    return TMRCA, effective_population_size, N0


def plot_psmc(T, Ne, N0, mu, minimum_tmrca, maximum_tmrca, target_file=False, scaled_by_generation_time=False,
              generation_time=None, log_yaxis=False, out_file=False, ax=None, color='red', label=None,
              plot_legend=True):
    """
    Plot PSMC
    :param T: np.array, TMRCA
    :param Ne: np.array, Effective population size
    :param N0: int, population size at time 0
    :param mu: float, mutation rate
    :param target_file: boolean, if present PSMC is target file
    :param scaled_by_generation_time: boolean, if time has been scaled by generation time
    :param log_yaxis: boolean, yaxis in log scale
    :param out_file: str, file path to save to
    :param ax: ax object, plot into given ax object
    :param color: str, color of	curve
    :param label: str, species name
    :param plot_legend: boolean, wheter to plot legend or not
    :return: ax, modified ax
    """
    if not ax:
        fig, ax = plt.subplots()
    if mu:
        Ne /= 1000
    if target_file:
        ax.plot(T, Ne, color=color, label=label)
    else:
        ax.plot(T, Ne, color=color, alpha=0.02)

    if mu:
        ax.set_ylabel(r"$N_e \times 10^3$")
    else:
        ax.set_ylabel(r"$4\mu N_e$")

    if mu and not scaled_by_generation_time:
        ax.set_xlabel('T in generations')
    elif mu and scaled_by_generation_time:
        ax.set_xlabel('Years before present (' + r'$\mu=$' + '{:.2e}'.format(mu) +
                      r', $g=$' + '{:.1f})'.format(generation_time))
    else:
        ax.set_xlabel(r'$2\mu T$')
    ax.set_xscale('log')
    if log_yaxis:
        ax.set_yscale('log')
    if not out_file:
        return ax
    if out_file:
        if maximum_tmrca == 0:
            max_x = ax.get_xlim()[1]
        else:
            max_x = maximum_tmrca
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(2.5))
        ax.set_xlim([minimum_tmrca, max_x])
        if plot_legend:
            ax.legend(bbox_to_anchor=(0.5, -.13), loc='upper center', ncol=2)
        fig = plt.gcf()
        fig.savefig(out_file, bbox_inches='tight')
        plt.close()


def parsing_and_plotting_helper(f, iteration, last_interval, mutation_rate, binsize, minimum_tmrca, maximum_tmrca,
                                generation_time, target_file, ax, out_file=None, log_yaxis=False, color='red',
                                label=None, plot_legend=True):
    """
    Helper function for parsing and plotting individual files
    :param f: str, input_file
    :param iteration: int, target iteration to parse
    :param last_interval: boolean, plot last atomic interval
    :param mutation_rate: float, mutation rate per nucleotide
    :param binsize: int, binsize used to prepare psmcfa file
    :param minimum_tmrca: int, minimum TMRCA
    :param maximum_tmrca: int, maximum TMRCA
    :param generation_time: int, generation time in years
    :param target_file: boolean, if file is target file
    :param ax: ax object, plot into given ax object
    :param out_file: str, file path to save to
    :param log_yaxis: boolean, plot y axis in log scale
    :param color: str, color of curve
    :param label: str, species name
    :param plot_legend: boolean, wheter to plot legend or not
    :return: ax, modified ax object
    """
    estimates = parse_psmc_file(f, iteration, last_interval)
    if not isinstance(estimates, int):
        t_k, lambda_k, theta_0 = estimates
    else:
        warnings.warn("PSMC file contains less iteration than the specified target iteration. "
                      "Setting target iteration from {} to {}".format(iteration, estimates))
        t_k, lambda_k, theta_0 = parse_psmc_file(f, estimates, last_interval)
    if t_k.shape[0] > 0:
        T, Ne, N0 = scale_results_to_generations(t_k, lambda_k, theta_0, mutation_rate, binsize)

        scaled_by_generation_time = False
        if generation_time:
            T *= generation_time
            scaled_by_generation_time = True
        # truncate
        Ne = Ne[T >= minimum_tmrca]
        T = T[T >= minimum_tmrca]
        if maximum_tmrca != 0:
            Ne = Ne[T <= maximum_tmrca]
            T = T[T<= maximum_tmrca]
        ax = plot_psmc(T, Ne, N0, mutation_rate, minimum_tmrca, maximum_tmrca,
                       target_file=target_file, scaled_by_generation_time=scaled_by_generation_time,
                       generation_time=generation_time, out_file=out_file, ax=ax, log_yaxis=log_yaxis, color=color,
                       label=label, plot_legend=plot_legend)
    return ax


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_files', nargs='+', required=True,
                        help='Input PSMC files to plot (can be multiple species separated by a space)')
    parser.add_argument('-d', '--input_directories', nargs='+', required=False,
                        help='Plot all .psmc in specified directories as CI')
    parser.add_argument('-l', '--labels', nargs='+', required=True,
                        help='Labels for different PSMC (e.g., species names)')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory where to save plot')
    parser.add_argument('-p', '--prefix', default='psmc', help='Prefix for output file. Default=psmc')
    parser.add_argument('-s', '--binsize', default=100, type=int, help='Bin size used to prepare psmcfa. Default=100')
    parser.add_argument('-m', '--mutation_rate', required=False, type=float, help='Mutation rate per nucleotide '
                                                                                  'per generation')
    parser.add_argument('-x', '--minimum_tmrca', default=0, type=float, help='Minimum TMRCA to plot. [0]')
    parser.add_argument('-X', '--maximum_tmrca', default=0, type=float, help='Maximum TMRCA to plot, 0 means all. [0]')
    parser.add_argument('--last_interval', action='store_true', default=False, help='Show last atomic interval')
    parser.add_argument('-g', '--generation_time', required=False, type=float,
                        help='Generation time to scale into years. If not provided, '
                             'results are expressed in terms of generations')
    parser.add_argument('-n', '--iteration', default=20, type=int, help='Iteration to plot. [20]')
    parser.add_argument('-c', '--colors', default=['red', 'blue', 'orange', 'green'],
                        help="colors to use [['red', 'blue', 'orange', 'green']]")
    parser.add_argument('--log_yaxis', action='store_true', default=False, help='Plot y-axis in log-scale')

    args = parser.parse_args() 
    if args.generation_time and not args.mutation_rate:
        raise AttributeError('If you set a generation time you also need to set a mutation rate')
    ax = None
    if args.input_directories:
        for input_directory in args.input_directories:
            if not input_directory.endswith('/'):
                input_directory += '/'
            input_files = glob.glob(input_directory +'round*.psmc')
            for f in input_files:
                ax = parsing_and_plotting_helper(f, args.iteration, args.last_interval, args.mutation_rate,
                                                 args.binsize, args.minimum_tmrca, args.maximum_tmrca,
                                                 args.generation_time, target_file=False, log_yaxis=args.log_yaxis,
                                                 ax=ax)
    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir += '/'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    out = None
    if len(args.input_files) > len(args.colors):
        raise AttributeError('Please specify more colors')
    for i, (input_file, label, color) in enumerate(zip(args.input_files, args.labels, args.colors)):
        if i == len(args.input_files) - 1:
            out = '{}{}.png'.format(output_dir, args.prefix)
        if len(args.labels) > 1:
            plot_legend = True
        else:
            plot_legend = False
        try:
            ax = parsing_and_plotting_helper(input_file, args.iteration, args.last_interval, args.mutation_rate,
                                             args.binsize, args.minimum_tmrca, args.maximum_tmrca, args.generation_time,
                                             target_file=True, log_yaxis=args.log_yaxis, ax=ax, out_file=out,
                                             color=color, label=label, plot_legend=plot_legend)
        except TypeError:
            # Figure was saved and nothing was returned
            pass


if __name__ == "__main__":
    main(sys.argv[1:])
