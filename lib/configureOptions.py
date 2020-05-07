
import os,sys
import argparse

class ConfigureWorkflowOptions():

	def __init__(this):
		document = """
T-CNV is a Python library and command-line software toolkit to infer and visualize copy number 
from high-throughput DNA sequencing data. It is designed for use with hybrid capture, including 
both whole-exome and custom target panels, and short-read sequencing platforms such as Illumina 
and Ion Torrent."""
		this.parser = argparse.ArgumentParser( prog='T-CNV', usage='%(prog)s [-h] [-v,--version]',
				formatter_class=argparse.RawDescriptionHelpFormatter,description = document)
	
	def get_run_options(this):
		document = """
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
"""
		this.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.2')

		OptionGroup = this.parser.add_argument_group(title='information options', description = document.lstrip("\n"))

		OptionGroup.add_argument('-t','--tumor',   type = str, help = argparse.SUPPRESS)

		OptionGroup.add_argument('-n','--normal',  type = str, help = argparse.SUPPRESS)

		OptionGroup.add_argument('-r','--region',  type = str, help = argparse.SUPPRESS)

		OptionGroup.add_argument('-f','--fasta',   type = str, help = argparse.SUPPRESS)

		OptionGroup.add_argument('-s','--sample',  type = str, help = argparse.SUPPRESS)

		OptionGroup.add_argument('-o','--outdir',  type = str, help = argparse.SUPPRESS)

		OptionGroup.add_argument('-c','--thread',  type = int, help = argparse.SUPPRESS, default = 10)

		OptionGroup.add_argument('-b','--binSize', type = int, help = argparse.SUPPRESS, default = 40)

		OptionGroup.add_argument('-d','--depth',   type = int, help = argparse.SUPPRESS, default = 30)

		OptionGroup.add_argument('-g','--gcRatio', type = float, help = argparse.SUPPRESS, default = 0.8)

		workflowOptions = this.parser.parse_args()
		this.checkWorkflowOptions(workflowOptions)

		return workflowOptions

	def checkWorkflowOptions(this,Options):
		## noArgOrError
		noArgOrError = [param for param in [Options.tumor, Options.region, Options.fasta, Options.outdir] if param is None ]
		if len(noArgOrError) > 0:
			this.parser.print_help()
			sys.exit(2)
		
		## outdir
		if not os.path.isdir(Options.outdir):
			try:
				os.makedirs(Options.outdir)
			except Exception as e:
				print("Fatal error: %s. failed to generate %s directory file"%(e,Options.outdir))
				sys.exit(2)
