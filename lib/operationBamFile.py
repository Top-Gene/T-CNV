import numpy as np
import multiprocessing
import os,re,sys,math,pysam
from collections import defaultdict
from baseSharedMethods import createDict

class dealBamFile():

	def filterBamFile(this,bamFile):
		pass

	def getRegionDepth(this, Hash, chorm=None, strat=None, end=None):
		RegionDepth = RegionAverageDepth = 0
		if chorm not in Hash.keys():
			return (RegionDepth,RegionAverageDepth)
		for i in range(int(strat),int(end) + 1):
			if i in Hash[chorm].keys():
				if math.isnan(float(Hash[chorm][i])): continue
				RegionDepth += float(Hash[chorm][i])
		RegionAverageDepth = RegionDepth / (int(end) - int(strat) + 1)
		return (RegionDepth,RegionAverageDepth)

	def getBedCoverage(this, bamFile, bedFile, bp = "20", mpQ = "20"):
		allSiteCover = createDict()
		try:
			result = pysam.depth("-aa", "-q", bp, "-Q", mpQ, "-b", bedFile, bamFile)
			for line in result.split("\n"):
				if re.search(r"^\s*$",line): continue
				cols = line.split("\t")
				allSiteCover[cols[0]][int(cols[1])] = int(cols[2])
		except Exception as e:
			print("error: %s"%e)
			return False
		return allSiteCover

	def getBedCovFile(this, bamFile, bedFile, outFile, bp = "20" , mpQ = "20"):
				try:
					with open(outFile,'w') as fileHandle:
						print(pysam.depth("-aa", "-q", bp, "-Q", mpQ, "-b", bedFile, bamFile),end = "", file=fileHandle)	
				except Exception as e:
					print("error: %s %s"%(bamFile,e))
					return False
				return True

	def getNormalCompare(this, tumorList, bedFile, outdir,thread = 1, bp = "0", mpQ = "0"):
		import warnings
		warnings.filterwarnings(action = 'ignore', category = RuntimeWarning)
		allSiteCover = createDict()
		hashSample = {}
		pool = multiprocessing.Pool(processes=int(thread))
		with open(tumorList,"r") as fileHandle:
			for tumor in fileHandle.readlines():
				tumor = tumor.rstrip("\n")
				sampleName = os.path.basename(tumor).replace(".bam","")
				sampleDepthStatFile = outdir + "/" + sampleName + ".depth.Statistics.txt"
				hashSample[sampleName] = sampleDepthStatFile
				if os.path.isfile(tumor):
					#pass
					pool.apply_async(this.getBedCovFile, (tumor, bedFile, sampleDepthStatFile))
				else:
					print("%s is not exists!!!"%tumor)
					continue

		pool.close()
		pool.join()
		pool.terminate()

		for sampleName in hashSample.keys():
			if not os.path.isfile(hashSample[sampleName]):
				print("Error: the %s depth statistics file not exists!!!"%(hashSample[sampleName]))
				continue
			with open (hashSample[sampleName],'r') as fileHandle:
				for line in fileHandle.readlines():
					line = line.rstrip('\n')
					if re.search(r'^\s*$',line) : continue
					chrom,pos,depth = line.split('\t')
					allSiteCover[sampleName][chrom][int(pos)] = int(depth)
					allSiteCover['Normal'][chrom].setdefault(int(pos),[]).append(int(depth))
		
		if 'Normal' in allSiteCover.keys(): 
			for chorm in allSiteCover['Normal'].keys():
				for pos in allSiteCover['Normal'][chorm].keys():
					if len(allSiteCover['Normal'][chorm][pos]) == 0:
						allSiteCover['Normal'][chorm][pos] = 0
						continue

					std = np.std(allSiteCover['Normal'][chorm][pos],ddof=1)
					mean = np.mean(allSiteCover['Normal'][chorm][pos])
					new_list = []
					for i in allSiteCover['Normal'][chorm][pos]:
						if (i > (mean + 2 * std)) or (i < (mean - 2 * std)):
							continue
						else:
							new_list.append(i)
					median = np.median(new_list)
					allSiteCover['Normal'][chorm][pos] = median
		else:
			print("Fatal error: Contrast Depth Generation Failed!!!")
			sys.exit(2)

		return allSiteCover
				
	def getGC(this,fasta,chorm=None,start=None,end=None):
		try:
			seq = pysam.faidx(fasta,"{}:{}-{}".format(chorm,start,end))
			AllBase = int(end) - int(start) + 1
			G_num = seq.count("G") + seq.count("g")
			C_num = seq.count("C") + seq.count("c")
			return "%.4f"% ((G_num + C_num)/AllBase)
		except:
			return False

	def insert_size(this, bamFile, threads=1):
		result_str = pysam.stats("-@", str(threads), bamFile)
		insert_size = float(re.search('SN\s+insert size average:\s+(\S+)', result_str).group(1))
		insert_std = float(re.search('SN\s+insert size standard deviation:\s+(\S+)', result_str).group(1))
		return (insert_size, insert_std)

	def read_pair(this, bamfile, threads=1, chrom=None, start=None, end=None):
		"""返回read pair

		1. 过滤proper_pair、unmap、duplicate、secondary/supplementary的比对；
		2. 过滤Read1、Read2 不同染色体的比对（可能是translocation）；
		3. 过滤Read1、Read2 比对后方向相同的（可能是Inversion）
			Read1->..........Read2-> 或 <-Read1..........<-Read2；
		4. 保留Read1、Read2 比对后方向不同的比对，
			Read1->..........<-Read2 或 Read2->..........<-Read1，且insert size>3*std的，可能是缺失；
			<-Read1..........Read2-> 或 <-Read2..........Read1->， 可能是扩增；

		"""
		read_dict = defaultdict(lambda: [None, None])
		samfile = pysam.AlignmentFile(bamfile, 'rb', threads=threads)
		for read in samfile.fetch(chrom, start, end):
			# NOTE: bwa mem -M = is_secondary; not -M = is_supplementary
			# NOTE: -F 1294
			if (read.is_proper_pair or read.is_unmapped or read.mate_is_unmapped or read.is_duplicate or read.is_secondary
					or read.is_supplementary or (read.reference_name != read.next_reference_name)
					or read.is_reverse == read.mate_is_reverse):
				continue
			qname = read.query_name
			if qname not in read_dict:
				if read.is_read1:
					read_dict[qname][0] = read
				else:
					read_dict[qname][1] = read
			else:
				if read.is_read1:
					yield read, read_dict[qname][1]
				else:
					yield read_dict[qname][0], read
				del read_dict[qname]

	def yield_discordants(this, bamfile, max_insert_size=0, threads=1, chrom=None, start=None, end=None):
		"""返回可能CNV的比对结果的生成器

		1. 判断类型
		2. 添加CNV开始、结束位置等属性

		Args:
			bamfile: 输入的bam文件, *.bqsr.bam，string
			max_insert_size: mean_insert_size+3*std的值，float
			threads: 线程数，默认为1，int
			chrom: 染色体，为空时遍历整个bam文件，string
			start: 区域开始坐标，为空时遍历整个bam文件，int
			end: 区域结束坐标，为空时遍历整个bam文件，int

		Returns:
			yield class Discordant

		TODO:
			1. 目前只分析read1带有的信息，因此部分比对无法得到read2准确的起始或结束位置；
			2. Raises Error 未完善；

		Raises:
			IOError: An error occurred accessing the bigtable.Table object.
		"""

		class Discordant:
			def __init__(this, read1, read2):
				this.chrom = read1.reference_name
				this.query_name = read1.query_name
				this.read1_start = read1.reference_start
				this.read1_end = read1.reference_end
				this.read2_start = read2.reference_start
				this.read2_end = read2.reference_end
				this.read1_is_reverse = read1.is_reverse
				this.read2_is_reverse = read2.is_reverse
				this.insert_size = abs(read1.template_length)
				this.start = min(this.read1_start, this.read1_end, this.read2_start, this.read2_end)
				this.end = max(this.read1_start, this.read1_end, this.read2_start, this.read2_end)
				if read1.is_reverse is False and read2.is_reverse is True:
					if read1.reference_start < read2.reference_start:
						this.pair_orientation = "F1R2"
					elif read1.reference_start > read2.reference_start:
						this.pair_orientation = "R2F1"
				elif read1.is_reverse is True and read2.is_reverse is False:
					if read1.reference_start < read2.reference_start:
						this.pair_orientation = "R1F2"
					elif read1.reference_start > read2.reference_start:
						this.pair_orientation = "F2R1"

		for read1, read2 in this.read_pair(bamfile, threads, chrom, start, end):
			discordant = Discordant(read1, read2)
			if discordant.pair_orientation in ("F1R2","F2R1") and discordant.insert_size >= max_insert_size:
				discordant.cnv_type = 'Deletion'
				yield discordant
			elif discordant.pair_orientation in ("R2F1","R1F2"):
				discordant.cnv_type = 'Duplication'
				yield discordant
