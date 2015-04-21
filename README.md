# NextSeq .bcl Conversion
`bcl_to_fastq` runs bcl2fastq with optional effects to the Sample Sheet and
concatenates reads across lanes into R1 and R2 by sample. By default,
Undetermined and reads across individual lanes are removed on success and
all reads are placed in BaseCalls directory.

# Running
```
$ cd /mnt/ilmn/150212_NS500409_0021_AH3CHYBGXX
$ bcl_to_fastq --reverse-complement
Info      : Using /mnt/ilmn/150212_NS500409_0021_AH3CHYBGXX/SampleSheet.csv
Info      : Found 48 samples for run 150212_NS500409_0021_AH3CHYBGXX
Info      : Run complete.
Info      : Converting .bcl to .fastq using: $>bcl2fastq -r 12 -d 12 -p 24 -w 12 --barcode-mismatches 0 -R /mnt/ilmn/150212_NS500409_0021_AH3CHYBGXX
Info      : Conversion successful
Info      : Generating demultiplexing stats file /mnt/ilmn/150212_NS500409_0021_AH3CHYBGXX/demultiplexing_stats.csv
Info      : Joining reads across lanes
Info      : Removing intermediate and Undetermined fastq files
```

# Results
In the run folder, SampleSheet.csv.bak is a backup copy of the original
SampleSheet.csv and is accompanied by:

###bcl2fastq.log
```
$ head bcl2fastq.log
bcl2fastq v2.15.0.4
2015-03-09 16:01:59 [7fde98ae2780] Command-line invocation: bcl2fastq -r 12 -d 12 -p 24 -w 12 --barcode-mismatches 0 -R '.'
2015-03-09 16:01:59 [7fde98ae2780] INFO: Minimum log level: INFO
2015-03-09 16:01:59 [7fde98ae2780] INFO: Runfolder path: '.'
2015-03-09 16:01:59 [7fde98ae2780] INFO: Input path: './Data/Intensities/BaseCalls/'
etc...
```

###demultiplexing_stats.csv
```
$ head demultiplexing_stats.csv
AAA003-K10,102570
AAA007-J07,72566
AAA240-I17,146605
AAA240-J05,197833
etc...
```

Joined fastq files (`<sample>_R?.fastq.gz`) are available in
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
  --joining INTEGER             number of threads used for file joining
                                [default: 12]
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
Tested only on Python 2.7.

# Install
```
git clone git@github.com:brwnj/bcl2fastq.git
cd bcl2fastq
python setup.py install
```
