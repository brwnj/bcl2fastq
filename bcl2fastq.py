#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function

import matplotlib
matplotlib.use('Agg')

import click
import fileinput
import os
import pandas as pd
import seaborn as sns
import string
import subprocess as sp
import sys
import time
from datetime import datetime
from glob import glob
from itertools import izip_longest
from xml.etree import cElementTree as ET

sns.set_context('paper')
sns.set_style('whitegrid', {'axes.linewidth': 1})


def log(category, message, *args, **kwargs):
    click.echo('%s: %s' % (click.style(category.ljust(10), fg='cyan'),
        message.replace('{}', click.style('{}', fg='yellow')).format(*args, **kwargs)))


def get_samplesheet(path):
    s = os.path.join(os.path.abspath(path), "SampleSheet.csv")
    log('Info', 'Using {}', s)
    if not os.path.exists(s):
        raise OSError(2, 'No such file', s)
    return s


def process_samplesheet(samplesheet, reverse_complement):
    _complement = string.maketrans("ATCG", "TAGC")
    complement = lambda seq: seq.translate(_complement)
    samples = []
    date = datetime.now().strftime("%Y-%m-%d-%H%M-%S")
    try:
        start = False
        sample_project_idx = 0
        index2_idx = None
        # strip whitespace and rewrite file in place
        for toks in fileinput.input(samplesheet,
                                    mode='rU',
                                    backup='.' + date + '.bak',
                                    inplace=True):
            toks = toks.rstrip("\r\n").split(',')
            if not start:
                # table header processing
                if toks[0] == "Sample_ID":
                    start = True
                    sample_project_idx = toks.index("Sample_Project")
                    if reverse_complement:
                        try:
                            index2_idx = toks.index("index2")
                        except ValueError:
                            index2_idx = toks.index("Index2")
            elif toks[0]:
                # convert underscores to dashes
                toks[0] = toks[0].replace("_", "-").replace(".", "-")
                toks[1] = toks[0]
                # will need to pull sample owner before overwriting
                toks[sample_project_idx] = ""
                samples.append(toks[0])

                # only adjust on known index
                if reverse_complement:
                    toks[index2_idx] = complement(toks[index2_idx])[::-1]
            # remove blank lines at end of table
            else:
                break
            print(",".join([t.strip() for t in toks]))
    finally:
        fileinput.close()
    run = os.path.basename(os.path.dirname(samplesheet))
    log('Info', 'Found {} samples for run {}', len(samples), run)
    return samples


def wait_for_completion(path, no_wait):
    if no_wait:
        return True
    status_xml = os.path.join(path, "RunCompletionStatus.xml")
    sleep_time = 1
    notify = True
    while not os.path.exists(status_xml):
        if notify:
            log('Info', "Waiting on run completion.")
            notify = False
        time.sleep(sleep_time)
        if sleep_time < 60:
            sleep_time += 1
    log('Info', "Run complete.")
    try:
        doc = ET.parse(status_xml)
        run_status = doc.find("CompletionStatus").text
    except IOError:
        raise IOError(2, "Count not find file", status_xml)
    except AttributeError:
        raise AttributeError(2, "Error parsing file", status_xml)
    return True if run_status == "CompletedAsPlanned" else False


def run_bcl2fastq(runfolder, args):
    runlog = os.path.join(runfolder, "bcl2fastq.log")
    cmd = " ".join(map(str, args))
    log("Info", "Converting .bcl to .fastq using: $>{}", cmd)
    with open(runlog, 'w') as fh:
        # bcl2fastq version info...
        sp.check_call("bcl2fastq --version 2>&1 | tail -2 | head -1",
                      stdout=fh,
                      stderr=fh,
                      shell=True)
        sp.check_call(cmd, stdout=fh, stderr=fh, shell=True)
    log("Info", "Conversion successful")
    return True


def compile_demultiplex_stats(runfolder, output_dir):
    stats_xml = os.path.join(os.path.abspath(output_dir), "Stats",
                             "DemultiplexingStats.xml")
    stats_csv = os.path.join(runfolder, "demultiplexing_stats.csv")
    plot_pdf = os.path.join(runfolder, "demultiplexing_distribution.pdf")
    if os.path.exists(stats_xml):
        log("Info", "Generating demultiplexing stats file {}", stats_csv)
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
        counts_df = pd.DataFrame(counts)
        counts_df.sum().to_csv(stats_csv)
        ax = counts_df.T.plot(kind='bar', stacked=True)
        fig = ax.get_figure()
        fig.savefig(plot_pdf, bbox_inches="tight")
    else:
        log('Error', 'Could not find file {}. Demux stats not generated.',
            stats_xml)
    return stats_csv


def build_concat_commands(samples, fastq_dir):
    cmds = []
    for idx, sample in enumerate(samples, start=1):
        for read in ['R1', 'R2']:
            cmd = ['cat']
            result_file = "%s/%s_%s.fastq.gz" % (fastq_dir, sample, read)
            for lane in range(1, 5):
                # build the file paths
                path = "%s/%s_S%d_L00%d_%s_001.fastq.gz" % (fastq_dir, sample,
                                                            idx, lane, read)
                if not os.path.exists(path):
                    sys.exit("Could not find %s. Concatenation failed." % path)
                cmd.append(path)
            cmds.append(" ".join(cmd) + " > " + result_file)
    return cmds


def join_fastqs(samples, fastq_dir, threads=1):
    log("Info", "Joining reads across lanes")
    success = True
    concat_cmds = build_concat_commands(samples, fastq_dir)
    groups = [(sp.Popen(cmd, shell=True) for cmd in concat_cmds)] * threads
    for processes in izip_longest(*groups):
        for p in filter(None, processes):
            p.wait()
            if p.returncode != 0:
                success = False
    return success


def cleanup(patterns):
    log("Info", "Removing intermediate and Undetermined fastq files")
    for p in patterns:
        for f in glob(p):
            os.remove(f)


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
@click.option("--joining",
              default=12,
              type=int,
              show_default=True,
              help="number of threads used for file joining")
@click.option("--keep-tmp",
              is_flag=True,
              default=False,
              show_default=True,
              help="save fastqs across lanes as well as Undetermined")
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
@click.argument('bcl2fastq_args', nargs=-1, type=click.UNPROCESSED)
def bcl2fastq(runfolder, loading, demultiplexing, processing, writing,
              barcode_mismatches, joining, keep_tmp, reverse_complement,
              no_wait, bcl2fastq_args):
    """Runs bcl2fastq2, creating fastqs and concatenating fastqs across lanes.
    Original fastq files and Undetermined files are deleted.
    """
    try:
        samplesheet = get_samplesheet(runfolder)
    except OSError:
        sys.exit("Could not find SampleSheet.csv")
    samples = process_samplesheet(samplesheet, reverse_complement)
    if len(samples) == 0:
        sys.exit(("No samples were found in the SampleSheet. "
                  "Please check its formatting."))
    completion_success = wait_for_completion(runfolder, no_wait)
    if not completion_success:
        sys.exit("Run did not complete as planned. Exiting.")
    fastq_dir = os.path.join(runfolder, "Data", "Intensities", "BaseCalls")
    cmd_args = ["bcl2fastq", "-r", loading, "-d", demultiplexing, "-p",
                processing, "-w", writing, "--barcode-mismatches",
                barcode_mismatches, "-R", runfolder] + list(bcl2fastq_args)
    call_status = run_bcl2fastq(runfolder, cmd_args)
    if not call_status:
        sys.exit("Something went wrong when trying to convert the .bcl files.")
    run_stats = compile_demultiplex_stats(runfolder, fastq_dir)
    # write file with sample names for downstream parallelization
    with open(os.path.join(fastq_dir, "SAMPLES"), 'w') as ofh:
        print(*samples, sep="\n", file=ofh)
    join_status = join_fastqs(samples, fastq_dir, joining)
    if join_status and not keep_tmp:
        cleanup([os.path.join(fastq_dir, "*_S*_L00*_*_001.fastq.gz"),
                 os.path.join(fastq_dir, "Undetermined_*.fastq.gz")])


if __name__ == '__main__':
    bcl2fastq()
