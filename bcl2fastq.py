#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function

import matplotlib
matplotlib.use('Agg')

import click
import fileinput
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import shutil
import six
import string
import subprocess as sp
import sys
import tempfile
import time
import warnings
from datetime import datetime
from glob import glob
from matplotlib.cbook import MatplotlibDeprecationWarning
from xml.etree import cElementTree as ET


warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
sns.set_context('paper')
sns.set_style('whitegrid', {'axes.linewidth': 1})
if six.PY2:
    _complement = string.maketrans("ATCG", "TAGC")
else:
    _complement = str.maketrans("ATCG", "TAGC")
complement = lambda seq: seq.translate(_complement)
logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s - %(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


def get_samplesheet(path):
    s = os.path.join(os.path.abspath(path), "SampleSheet.csv")
    logging.info("Using %s", s)
    if not os.path.exists(s):
         raise OSError(2, "No such file", s)
    return s


def get_file_sizes(output_dir):
    total_size = 0
    for f in glob(os.path.join(output_dir, "*.fastq.gz")):
        if os.path.basename(f).startswith("Undetermined_"):
            continue
        total_size += os.path.getsize(f)
    return total_size


def process_samplesheet(samplesheet, new_samplesheet, reverse_complement=False,
    determine=False):
    """Fix hidden characters in sample names and optional reverse complement
    the second index.

    Args:
        samplesheet (str): file path to SampleSheet.csv
        reverse_complement (bool): to reverse complement 'Index2'

    Returns:
        list
    """
    if not determine:
        logging.info("Processing %s", samplesheet)

    samples = []
    start = False
    index2_idx = None

    with open(samplesheet, "rU" if six.PY2 else "r") as ifh, open(new_samplesheet, "w") as ofh:
        for line in ifh:
            toks = line.strip().split(",")
            if not start:
                # table header processing
                if toks[0] == "Sample_ID":
                    start = True
                    if reverse_complement:
                        if "index2" in toks:
                            index2_idx = toks.index("index2")
                        elif "Index2" in toks:
                            index2_idx = toks.index("Index2")
                        else:
                            logging.warn("There is no Index2 to reverse complement")

            elif toks[0]:
                # convert underscores to dashes
                toks[0] = toks[0].replace("_", "-").replace(".", "-")
                toks[1] = toks[0]
                samples.append(toks[0])

                # only adjust on known index
                if reverse_complement and index2_idx:
                    toks[index2_idx] = complement(toks[index2_idx])[::-1]

            # remove blank lines at end of table
            else:
                break
            print(*[t.strip() for t in toks], sep=",", file=ofh)

    run = os.path.basename(os.path.dirname(samplesheet))
    if not determine:
        logging.info('Found %d samples for run %s', len(samples), run)
    return samples


def wait_for_completion(path, no_wait):
    if no_wait:
        return True
    status_xml = os.path.join(path, "RunCompletionStatus.xml")
    sleep_time = 1
    notify = True
    while not os.path.exists(status_xml):
        if notify:
            logging.info("Waiting on run completion.")
            notify = False
        time.sleep(sleep_time)
        if sleep_time < 60:
            sleep_time += 1
    logging.info("Run complete.")
    try:
        doc = ET.parse(status_xml)
        run_status = doc.find("CompletionStatus").text
    except IOError:
        raise IOError(2, "Count not find file", status_xml)
    except AttributeError:
        raise AttributeError(2, "Error parsing file", status_xml)
    return True if run_status == "CompletedAsPlanned" else False


def run_bcl2fastq(runfolder, args, determine=False):
    runlog = os.path.join(runfolder, "bcl2fastq.log")
    cmd = " ".join(map(str, args))
    if determine:
        logging.info("Converting a subset to determine barcodes using: $>%s", cmd)
    else:
        logging.info("Converting .bcl to .fastq using: $>%s", cmd)
    with open(runlog, "w") as fh:
        # bcl2fastq version info...
        sp.check_call("bcl2fastq --version 2>&1 | tail -2 | head -1",
                      stdout=fh,
                      stderr=fh,
                      shell=True)
        sp.check_call(cmd, stdout=fh, stderr=fh, shell=True)
    logging.info(".bcl Conversion successful")


def run_determination_step(runfolder, samplesheet, loading, demultiplexing,
                           processing, writing, barcode_mismatches):
    # set up a temporary working directory
    tmpd = tempfile.mkdtemp(dir=runfolder)

    date = datetime.now().strftime("%Y-%m-%d-%H%M-%S")
    rc_samplesheet = "%s.%s.rc.csv" % (samplesheet, date)
    process_samplesheet(samplesheet, rc_samplesheet, reverse_complement=True,
                        determine=True)

    # run the first set
    cmd_args = ["bcl2fastq", "--tiles", "11105", "--no-lane-splitting",
                "--output-dir", tmpd,
                "--barcode-mismatches", barcode_mismatches,
                "--loading-threads", loading,
                "--demultiplexing-threads", demultiplexing,
                "--processing-threads", processing,
                "--writing-threads", writing,
                "--sample-sheet", rc_samplesheet]
    run_bcl2fastq(tmpd, cmd_args, determine=True)
    rc_file_size = get_file_sizes(tmpd)
    shutil.rmtree(tmpd)

    # try the original
    tmpd = tempfile.mkdtemp(dir=runfolder)
    orig_samplesheet = "%s.%s.orig.csv" % (samplesheet, date)
    samples = process_samplesheet(samplesheet, orig_samplesheet,
                                  reverse_complement=False, determine=True)

    # run the first set
    cmd_args = ["bcl2fastq", "--tiles", "11105", "--no-lane-splitting",
                "--output-dir", tmpd,
                "--barcode-mismatches", barcode_mismatches,
                "--loading-threads", loading,
                "--demultiplexing-threads", demultiplexing,
                "--processing-threads", processing,
                "--writing-threads", writing,
                "--sample-sheet", orig_samplesheet]
    run_bcl2fastq(tmpd, cmd_args, determine=True)
    orig_file_size = get_file_sizes(tmpd)
    shutil.rmtree(tmpd)

    if rc_file_size > orig_file_size:
        logging.info("Using reverse complement of Index2 to demultiplex")
        os.remove(orig_samplesheet)
        return samples, rc_samplesheet
    elif rc_file_size == orig_file_size:
        logging.critical(("The original and reverse complemented barcodes "
                          "yielded the same number of demultiplexed reads."))
        sys.exit(1)
    else:
        logging.info("Using the original barcodes to demultiplex")
        os.remove(rc_samplesheet)
        return samples, orig_samplesheet


def xml_to_df(stats_xml):
    if os.path.exists(stats_xml):
        logging.info("Generating demultiplexing stats file")
        doc = ET.parse(stats_xml)
        root = doc.getroot()
        counts = {}
        for sample in root.iter('Sample'):
            name = sample.get('name')
            if name == "all" or name == "unknown": continue
            counts[name] = {}
            for barcode in sample.iter('Barcode'):
                if barcode.get('name') == "all": continue
                for lane in barcode.iter('Lane'):
                    lane_name = lane.get('number')
                    count = int(lane.findtext('BarcodeCount'))
                    counts[name][lane_name] = count
        return pd.DataFrame(counts)
    else:
        logging.warning('Could not find file %s', stats_xml)
        return None


def barplot_distribution(df, out_file):
    width = max([len(df) / 10, 12])
    f, ax = plt.subplots(figsize=(width, 6))
    df.plot(kind='bar', stacked=True, ax=ax)
    f.savefig(out_file, bbox_inches='tight')
    plt.close()


def Lc(x):
    """Computes the ordinary and generalized Lorenz curve of a list.

    >>> import numpy as np
    >>> t = [1,2,np.nan,7,8]
    >>> p, L, Lg = Lc(t)
    >>> len(p) == len(L) == len(Lg)
    True
    >>> p[1:4]
    array([ 0.25,  0.5 ,  0.75])
    >>> L[1:4] # doctest: +ELLIPSIS
    array([ 0.055...,  0.166...,  0.555...])
    >>> Lg[1:4]
    array([ 0.25,  0.75,  2.5 ])
    >>> t = [1,2,np.nan,7,-8]
    >>> Lc(t) # doctest: +ELLIPSIS
    Traceback (most recent call last):
     ...
    ValueError: x contained negative number
    """
    assert len(x) > 0, "x is empty"
    a = np.array(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.min() < 0:
        raise ValueError("x contained negative number")
    a.sort(kind='mergesort')
    a_len = float(len(a))
    p = np.arange(1, a_len + 1) / a_len
    p = np.append([0], p)
    L = a.cumsum() / a.sum()
    L = np.append([0], L)
    Lg = L * np.mean(a)
    return p, L, Lg


def lorenz_curve(df, out_file):
    p, L, Lg = Lc(df.sum(axis=1).values)
    f, ax = plt.subplots(figsize=(8, 6))
    plt.plot(p, L, axes=ax)
    plt.plot([0,1], axes=ax, color='black', linestyle="--")
    ax.set(title="Distribution of Barcode Mapped Reads")
    f.savefig(out_file, bbox_inches='tight')
    plt.close()


def compile_demultiplex_stats(runfolder, out_dir):
    stats_xml = os.path.join(os.path.abspath(out_dir), "Stats", "DemultiplexingStats.xml")
    df = xml_to_df(stats_xml)
    if df is not None:
        df.sum().to_csv(os.path.join(runfolder, "demultiplexing_stats.csv"))
        try:
            dft = df.transpose().drop("Undetermined", axis=0)
        except ValueError:
            dft = df.transpose()
        barplot_distribution(dft, os.path.join(runfolder, "demultiplexing_distribution.pdf"))
        lorenz_curve(dft, os.path.join(runfolder, "demultiplexing_distribution_curve.pdf"))


@click.command(context_settings=dict(
               help_option_names=['-h', '--help'],
               ignore_unknown_options=True,))
@click.option("--runfolder",
              default=".",
              show_default=True,
              help="path to directory containing run data")
@click.option("--loading",
              default=12,
              type=int,
              show_default=True,
              help="number of threads used for loading BCL data")
@click.option("--demultiplexing",
              default=12,
              type=int,
              show_default=True,
              help="number of threads used for demultiplexing")
@click.option("--processing",
              default=24,
              type=int,
              show_default=True,
              help="number of threads used for processing demultiplexed data")
@click.option("--writing",
              default=12,
              type=int,
              show_default=True,
              help="number of threads used for writing FASTQ data")
@click.option("--barcode-mismatches",
              default=0,
              type=int,
              show_default=True,
              help="number of allowed mismatches per index")
@click.option("--keep-tmp",
              is_flag=True,
              default=False,
              show_default=True,
              help="save Undetermined reads")
@click.option("--reverse-complement",
              is_flag=True,
              default=False,
              show_default=True,
              help="reverse complement index 2 of the sample sheet")
@click.option("--no-wait",
              is_flag=True,
              default=False,
              show_default=True,
              help="process the run without checking its completion status")
@click.option("--overwrite",
              is_flag=True,
              default=False,
              show_default=True,
              help="overwrite existing fastq files in the output directory")
@click.option("--determine",
              is_flag=True,
              default=False,
              show_default=True,
              help="use barcodes in samplesheet as well as the reverse complement of index 2, then demultiplex with best")
@click.option("--no-cleanup",
              is_flag=True,
              default=False,
              show_default=True,
              help="skip all cleaning up -- do not rename fastq output and do not delete undetermined files")
@click.argument('bcl2fastq_args', nargs=-1, type=click.UNPROCESSED)
def bcl2fastq(runfolder, loading, demultiplexing, processing, writing,
              barcode_mismatches, keep_tmp, reverse_complement,
              no_wait, overwrite, determine, no_cleanup, bcl2fastq_args):
    """Runs bcl2fastq2, creating fastqs and concatenating fastqs across lanes.
    Undetermined files are deleted by default.

    Any arguments not matching those outlined below will be sent to the
    `bcl2fastq` call.
    """
    try:
        samplesheet = get_samplesheet(runfolder)
    except OSError:
        logging.critical("Could not find SampleSheet.csv")
        sys.exit(1)

    # will wait on the run to complete without ever checking for samples in
    # the samplesheet
    if not determine:
        # new samplesheet written each time leaving the original unchanged
        date = datetime.now().strftime("%Y-%m-%d-%H%M-%S")
        new_samplesheet = "%s.%s.csv" % (samplesheet, date)
        if os.path.exists(os.path.join(runfolder, "bcl2fastq.log")):
            # this run has already been converted, so don't reverse complement
            # just get the sample names
            if reverse_complement:
                logging.warning("reverse complementing has been skipped as a log file (bcl2fastq.log) was found")
            samples = process_samplesheet(samplesheet, new_samplesheet, False)
        else:
            samples = process_samplesheet(samplesheet, new_samplesheet,
                                          reverse_complement)
        if len(samples) == 0:
            logging.critical(("No samples were found in the SampleSheet. "
                              "Please check its formatting."))
            sys.exit(1)

    # check for 'CompletedAsPlanned' in RunCompletionStatus.xml
    completion_success = wait_for_completion(runfolder, no_wait)
    if not completion_success:
        logging.critical("Run did not complete as planned. Exiting.")
        sys.exit(1)

    # set where we're going to dump the result files
    fastq_dir = os.path.abspath(os.path.join(runfolder, "Data", "Intensities", "BaseCalls"))

    # check original and reverse complement using a few tiles
    if determine:
        samples, new_samplesheet = run_determination_step(runfolder, samplesheet,
                                                          loading, demultiplexing,
                                                          processing, writing,
                                                          barcode_mismatches)

    # run bcl2fastq on the run folder
    cmd_args = ["bcl2fastq", "--sample-sheet", new_samplesheet,
                "--loading-threads", loading,
                "--demultiplexing-threads", demultiplexing,
                "--processing-threads", processing,
                "--writing-threads", writing,
                "--barcode-mismatches", barcode_mismatches,
                "--no-lane-splitting",
                "--runfolder-dir", runfolder] + list(bcl2fastq_args)
    run_bcl2fastq(runfolder, cmd_args)

    # parse DemultiplexingStats.xml into a csv with summary plots
    compile_demultiplex_stats(runfolder, fastq_dir)

    # TODO: deprecate!
    # write file with sample names for downstream parallelization
    with open(os.path.join(fastq_dir, "SAMPLES"), 'w') as ofh:
        print(*samples, sep="\n", file=ofh)

    # cleanup the output directory
    if not no_cleanup:
        logging.info("Cleaning up output directory [%s]" % fastq_dir)
        for f in glob(os.path.join(fastq_dir, "*.fastq*")):
            if not f.endswith(".gz") and not f.endswith(".fastq"):
                continue
            filename = os.path.basename(f)
            if filename.startswith("Undetermined_") and not keep_tmp:
                logging.info("Deleting %s" % filename)
                os.remove(f)
            else:
                try:
                    # AD-332-A10_S1_R1_001.fastq.gz --> AD-332-A10_R1.fastq.gz
                    sample_name, sample_number, read_index, ext = filename.split("_")
                    # munge the file name
                    new_file_name = "%s_%s.%s" % (sample_name, read_index, ext.partition('.')[-1])
                    # prepend the path
                    new_file_name = os.path.join(os.path.dirname(f), new_file_name)
                    if overwrite and os.path.exists(new_file_name):
                        os.remove(new_file_name)
                    os.rename(f, new_file_name)
                except ValueError:
                    logging.warn("Renaming skipped: the output dir contains conflicting FASTQ file for %s" % f)


if __name__ == '__main__':
    bcl2fastq()
