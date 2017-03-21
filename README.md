# NextSeq .bcl Conversion
`bcl_to_fastq` runs bcl2fastq with optional effects to the Sample Sheet and
concatenates reads across lanes into R1 and R2 by sample. By default,
Undetermined and reads across individual lanes are removed on success and
all reads are placed in BaseCalls directory.

Tested on `bcl2fastq2` Conversion Software v2.17.1.14 and Python 2.7 and 3.5.

# Running
```
$ cd /data/nextseq/170111_NS500409_0130_AHHGTMAFXX/SampleSheet.csv
$ bcl_to_fastq --reverse-complement --processing 80
[2017-01-14 19:07:57 - INFO] Using /data/nextseq/170111_NS500409_0130_AHHGTMAFXX/SampleSheet.csv
[2017-01-14 19:07:57 - INFO] Processing /data/nextseq/170111_NS500409_0130_AHHGTMAFXX/SampleSheet.csv
[2017-01-14 19:07:58 - INFO] Found 384 samples for run 170111_NS500409_0130_AHHGTMAFXX
[2017-01-14 19:07:58 - INFO] Run complete.
[2017-01-14 19:07:58 - INFO] Converting .bcl to .fastq using: $>bcl2fastq -r 12 -d 12 -p 80 -w 12 --barcode-mismatches 0 --no-lane-splitting -R .
[2017-01-14 20:54:02 - INFO] .bcl Conversion successful
[2017-01-14 20:54:02 - INFO] Generating demultiplexing stats file
```

# Results
In the run folder, SampleSheet.csv.bak is a backup copy of the original
SampleSheet.csv and is accompanied by:

## bcl2fastq.log

```
$ head bcl2fastq.log
2017-01-14 19:07:58 [f82880] INFO: Create FASTQs for index reads: NO
BCL to FASTQ file converter
bcl2fastq v2.17.1.14
Copyright (c) 2007-2015 Illumina, Inc.

2017-01-14 19:07:58 [16f2880] Command-line invocation: bcl2fastq -r 12 -d 12 -p 80 -w 12 --barcode-mismatches 0 --no-lane-splitting -R .
2017-01-14 19:07:58 [16f2880] INFO: Minimum log level: INFO
2017-01-14 19:07:58 [16f2880] INFO: Sample sheet: './SampleSheet.csv'
2017-01-14 19:07:59 [16f2880] INFO: Runfolder path: '.'
2017-01-14 19:07:59 [16f2880] INFO: Input path: './Data/Intensities/BaseCalls/'
etc...
```

## demultiplexing_stats.csv
```
$ head demultiplexing_stats.csv
AAA003-K10,102570
AAA007-J07,72566
AAA240-I17,146605
AAA240-J05,197833
etc...
```

Fastq files (`<sample>_R?.fastq.gz`) are available in
RunFolder/Data/Intensities/BaseCalls along with a file named SAMPLES which
merely listed the sample IDs that were processed.

# Help
```
$ bcl_to_fastq -h
Usage: bcl_to_fastq [OPTIONS]

  Runs bcl2fastq2, creating fastqs and concatenating fastqs across lanes.
  Original fastq files and Undetermined files are deleted.

Options:
  --runfolder TEXT              path to directory containing run data
                                [default: .]
  --loading INTEGER             number of threads used for loading BCL data
                                [default: 12]
  --demultiplexing INTEGER      number of threads used for demultiplexing
                                [default: 12]
  --processing INTEGER          number of threads used for processing
                                demultiplexed data  [default: 24]
  --writing INTEGER             number of threads used for writing FASTQ data
                                [default: 12]
  --barcode-mismatches INTEGER  number of allowed mismatches per index
                                [default: 0]
  --keep-tmp                    save fastqs across lanes as well as
                                Undetermined  [default: False]
  --reverse-complement          reverse complement index 2 of the sample sheet
                                [default: False]
  --no-wait                     process the run without checking its
                                completion status  [default: False]
  -h, --help                    Show this message and exit.
```

# Requires
+ [click](http://click.pocoo.org/4/)
+ [pandas](http://pandas.pydata.org/)
+ [bcl2fastq2](http://support.illumina.com/downloads/bcl2fastq_conversion_software.html)
+ matplotlib
+ numpy
+ seaborn


# Install
```
git clone git@github.com:brwnj/bcl2fastq.git
cd bcl2fastq
python setup.py install
```
