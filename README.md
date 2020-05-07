# T-CNV

- [what's T-CNV](#whats-t-cnv)
- [requirements](#requirements)
- [all options](#all-options)

# what's T-CNV

`T-CNV` is a Python library and command-line software toolkit to infer and visualize copy number 
from high-throughput DNA sequencing data. It is designed for use with hybrid capture, including 
both whole-exome and custom target panels, and short-read sequencing platforms such as Illumina 
and Ion Torrent.

# requirements

`T-CNV` requires pandas, matplotlib, scikit-learn, pysam to be installed. To install, type:

```
pip install pandas matplotlib sklearn pysam
```

# all options
```
optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

information options:
  -t, --tumor    <FILE>     The Tumor Sample BAM file by BQSR or bam file list.
  -n, --normal   <FILE>     The Normal Sample BAM file by BQSR.
  -r, --region   <FILE>     The target region file.
  -f, --fasta    <FILE>     The faidx indexed reference sequence file.
  -s, --sample   <STR>      Output sample name.
  -o, --outdir   <DIR>      Output dir.
  -b, --binSize  <INT>      Stat Depth window. default 40
  -c, --thread   <INT>      How many threads are used to analyze samples. default 10
  -d, --depth    <INT>      The target region average depth threshold. default 30   
  -g, --gcRatio  <FLOAT>    The target region GC Ratio threshold. default 0.8
```
