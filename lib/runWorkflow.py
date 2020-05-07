import random
import pandas as pd
import os,re,sys,math
import multiprocessing
import subprocess
import visualization
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from baseSharedMethods import createDict,computelog2
from operationBamFile import dealBamFile

class DetectByCNVect():

    def __init__(this,workflowOptions):
        this.options = workflowOptions
        this.dealBamObject = dealBamFile()


    def getLog2Ratio(this,regionTumorDepth,regionNormalDepth,TumorAllDepth,NormalAllDepth):
        log2Value = computelog2(regionTumorDepth/regionNormalDepth * (NormalAllDepth / TumorAllDepth))
        return log2Value


    def checkExonCnvByBins(this,tumor,normal,allDepth_t,allDepth_n,statu,chorm=None,start=None,end=None):
        regionSize = int(end) - int(start) + 1
        binsNum = math.ceil(float(regionSize) / this.options.binSize)
        binsInfo = []
        if binsNum < 3 : return statu,binsNum,"NA"

        DelNum = 0
        DupNum = 0
        NonCnvNum = 0

        for num in range(1,binsNum + 1):
            binStart = start + (num -1) * this.options.binSize
            binEnd = binStart + this.options.binSize
            if binEnd >= end:
                binEnd = end

            if (binEnd - binStart) <= this.options.binSize * 0.5:
                continue
            tumorAllDepth,tumorAverageDepth = this.dealBamObject.getRegionDepth(tumor, chorm, binStart, binEnd)
            normalAllDepth,normalAverageDepth = this.dealBamObject.getRegionDepth(normal, chorm, binStart, binEnd)
            if tumorAverageDepth < this.options.depth or normalAverageDepth < this.options.depth: continue
            log2Value = this.getLog2Ratio(tumorAverageDepth,normalAverageDepth,allDepth_t,allDepth_n)
            CnvNum,CnvStatu = this.log2toCnvStatu(log2Value)
            if CnvStatu == "Del":
                DelNum += 1
                binsInfo.append("Del")
            elif CnvStatu == "Dup":
                DupNum += 1
                binsInfo.append("Dup")
            else:
                NonCnvNum += 1
                binsInfo.append("Non-CNV")

        binsInfo = ",".join(binsInfo) if len(binsInfo) > 0 else "NA"

        allBinsNum = DelNum + DupNum + NonCnvNum
        if allBinsNum > 0:
            DelRatio = DelNum / allBinsNum
            DupRatio = DupNum / allBinsNum
            if statu == "Del" and DelRatio > 0.9:
                return "Del",binsNum,binsInfo
            elif statu == "Dup" and DupRatio > 0.9:
                return "Dup",binsNum,binsInfo
            else:
                return "Non-CNV",binsNum,binsInfo
        else:
            return statu,binsNum,binsInfo

    def checkExonCnvByWins(this,tumor,normal,allDepth_t,allDepth_n,statu,chorm=None,start=None,end=None):
        winSize=10
        stepSize=6
        std_threshold = 0.14

        regionSize = int(end) - int(start) + 1
        binsNum = math.floor((regionSize - winSize)/stepSize)+1

        if binsNum < 3: return statu, binsNum

        win_start=int(start)
        win_end = int(start)+winSize
        index=1
        log2result=createDict()

        while win_end <= int(end):
            tumorAllDepth,tumorAverageDepth = this.dealBamObject.getRegionDepth(tumor, chorm, win_start, win_end)
            normalAllDepth,normalAverageDepth = this.dealBamObject.getRegionDepth(normal, chorm, win_start, win_end)
            if tumorAverageDepth < this.options.depth or normalAverageDepth < this.options.depth:
                win_start+=stepSize
                win_end+=stepSize
                continue
            log2Value = this.getLog2Ratio(tumorAverageDepth,normalAverageDepth,allDepth_t,allDepth_n)
            log2result.setdefault('index',[]).append(index)
            log2result.setdefault('log2',[]).append(log2Value)

            index+=1
            win_start+=stepSize
            win_end+=stepSize

        binsNum = len(log2result['log2'])
        if binsNum < 3:
            return statu, binsNum

        tmp_table=pd.DataFrame(log2result)
        std_value = tmp_table['log2'].std()
        if std_value <= std_threshold:
            if all((-0.42 <= x <= 0.32) for x in log2result['log2']):
                return 'Non-CNV',binsNum
            else:
                return statu,binsNum
        else:
            return 'Non-CNV', binsNum

    def gmmClustering(this,df_cnv,Gene,Exon,BinsResult,genelist):

        df_gene = df_cnv[df_cnv["#Gene"].isin(genelist)].loc[:,["#Gene","Exon","log2AdjValue"]].reset_index(drop=True)

        # print(df_gene.head(len(df_gene)))
        data = df_gene.loc[:,["log2AdjValue"]].values

        stand = StandardScaler()
        Features = stand.fit_transform(data)

        gmm = GaussianMixture(n_components=2, covariance_type='tied').fit(Features)
        y_pred = gmm.predict(Features) ## 训练数据

        df_gene["label"] = y_pred

        result = "Non-CNV"

        Gene_label = df_gene[(df_gene["#Gene"] == Gene) & (df_gene["Exon"] == Exon)].reset_index(drop=True)
        Gene_label = Gene_label.ix[0,"label"]

        label0 = df_gene[df_gene["label"] == 0]
        label1 = df_gene[df_gene["label"] == 1]

        label0_log_mean = label0["log2AdjValue"].mean()
        label1_log_mean = label1["log2AdjValue"].mean()

        if abs(0-label0_log_mean) > abs(0-label1_log_mean):
            if Gene_label == 0:
                result = BinsResult
            else:
                result = "Non-CNV"
        else:
            if Gene_label == 0:
                result = "Non-CNV"
            else:
                result = BinsResult
        return result


    def callingCNVByReadDepthMethod(this,tumor,normal,sample,outdir):
        tumorAllDepth = normalAllDepth = 0
        regionDepth = createDict()
        regionGCRatio = createDict()
        with open(this.options.region,"r") as fileHandle:
            for line in fileHandle.readlines():
                if not re.search(r'^chr', line) : continue
                cols = line.rstrip("\n").split("\t")
                region = "{}\t{}\t{}".format(cols[0],cols[1],cols[2])
                try:
                    cols[4] = int(cols[4])
                except ValueError as e:
                    print ("Warning: %s Can't be converted into integers. %s!!!"%(cols[4],e))
                    continue
                AllDepth_t,AverageDepth_t = this.dealBamObject.getRegionDepth(tumor,cols[0],cols[1],cols[2])
                AllDepth_n,AverageDepth_n = this.dealBamObject.getRegionDepth(normal,cols[0],cols[1],cols[2])
                GCRatio = this.dealBamObject.getGC(this.options.fasta,cols[0],cols[1],cols[2])
                tumorAllDepth += AllDepth_t
                normalAllDepth += AllDepth_n
                regionGCRatio[region] = GCRatio
                regionDepth[cols[3]][cols[4]][region]["tumor"] = AverageDepth_t
                regionDepth[cols[3]][cols[4]][region]["normal"] = AverageDepth_n

        log2File = "%s/%s.log2.txt"%(outdir,sample)
        log2adjFile = "%s/%s.log2.adj"%(outdir,sample)

        try:
            outlog2 = open(log2File,'w')
        except IOError as e:
            print("Fatal error: fail to write to %s file %s!!!"%(log2File,e))
            sys.exit(2)
        for Gene in sorted(regionDepth.keys()):
            for Exon in sorted(regionDepth[Gene].keys()):
                for Region in sorted(regionDepth[Gene][Exon].keys()):
                    if regionDepth[Gene][Exon][Region]["tumor"] < this.options.depth or \
                      regionDepth[Gene][Exon][Region]["normal"] < this.options.depth:
                        print("Warning: the region %s average depth below threshold!!!"%(Region))
                        continue
                    if Region not in regionGCRatio.keys():
                        print("Warning: the region %s GC Ratio not available!!!"%Region)
                        continue
                    else:
                        if float(regionGCRatio[Region]) >= this.options.gcRatio:
                            print("Warning: the region %s GC Ratio above threshold!!!"%Region)
                            continue

                    log2Value = this.getLog2Ratio(regionDepth[Gene][Exon][Region]["tumor"],
                      regionDepth[Gene][Exon][Region]["normal"],tumorAllDepth,normalAllDepth)
                    outlog2.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Gene,Exon,Region,regionGCRatio[Region],
                     regionDepth[Gene][Exon][Region]["tumor"],regionDepth[Gene][Exon][Region]["normal"],log2Value))
        outlog2.close()

        ## log2值矫正
        try:
            ##  待完善
            loess = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),'loess_v2.0.R')
            subprocess.check_call("Rscript %s --infile %s --prefix %s"%(loess,log2File,log2adjFile),shell=True)
        except:
            print("Fatal error: fial to adjust log2!!!")
            sys.exit(2)

        df_bed = pd.read_csv(this.options.region,sep="\t",header=None,usecols=[0,3],names=["chr","gene"])
        df_bed.drop_duplicates(keep='first',inplace=True)
        df_bed.reset_index(drop=True,inplace=True)

        df = pd.DataFrame(columns=["#Gene","Exon","Chr","Start","End","gcRatio","TumorAverageDepth","NormalAverageDepth","log2Value","log2AdjValue","CnvNum","CnvStatu","BinsResult","Clustering","Clustering_Info","WinsResult"])

        df_adj = pd.read_csv("%s.result.txt"%log2adjFile,sep="\t",names=["#Gene","Exon","Chr","Start","End","gcRatio","TumorAverageDepth","NormalAverageDepth","log2Value","log2AdjValue"])

        for index in df_adj.index:
            row = df_adj.loc[index,:].to_dict()

            row["CnvNum"],row["CnvStatu"] = this.log2toCnvStatu(row["log2AdjValue"])

            ## 划分 bins 验证CNV
            chorm,start,end = row["Chr"],row["Start"],row["End"]
            row["BinsResult"],binsNum,binsInfo = this.checkExonCnvByBins(tumor,normal,tumorAllDepth,normalAllDepth,
              row["CnvStatu"],chorm,int(start),int(end))

            ## 聚类分析验证
            if row["CnvStatu"] == "Non-CNV":
                row["Clustering"] = "Non-CNV"
                row["Clustering_Info"] = "Exon结果为Non-CNV"
                row['WinsResult'] = "Non-CNV"
                df = df.append(row,ignore_index=True)
                continue

            df_gene = df_adj[df_adj["#Gene"] == row["#Gene"]]
            delNum = len(df_gene[df_gene["log2AdjValue"] <= -0.42])
            dupNum = len(df_gene[df_gene["log2AdjValue"] > 0.32])

            Del_ratio = delNum / len(df_gene)
            Dup_ratio = dupNum / len(df_gene)

            if Del_ratio > 0.9 or Dup_ratio > 0.9:
                row["Clustering"] = row["CnvStatu"]
                row["Clustering_Info"] = "90%以上为CNV"
                row['WinsResult'] = row["CnvStatu"]
                df = df.append(row,ignore_index=True)
                continue

            chorm = df_bed[df_bed["gene"] == row["#Gene"]].reset_index(drop=True,inplace=False)
            chorm = chorm.ix[0,"chr"]

            Gene_list = df_bed[df_bed["chr"] == chorm].reset_index(drop=True,inplace=False)

            ## 染色体基因数目少于3个不进行聚类验证
            if len(Gene_list) < 3:
                row["Clustering"] = row["BinsResult"]
                row["Clustering_Info"] = "基因数目少于3个"
                row['WinsResult'] = row["BinsResult"]
                df = df.append(row,ignore_index=True)
                continue

            # 获取 Gene 的index
            Gene_index = Gene_list[Gene_list["gene"] == row["#Gene"]].index.tolist()[0]

            Gene_1 = []
            Gene_2 = []
            Gene_3 = []
            if Gene_index > 1 and len(Gene_list) >= (Gene_index + 3):
                Gene_1 = Gene_list.iloc[(Gene_index - 2):(Gene_index + 1),:]["gene"]
                Gene_2 = Gene_list.iloc[(Gene_index - 1):(Gene_index + 2),:]["gene"]
                Gene_3 = Gene_list.iloc[Gene_index:(Gene_index + 3),:]["gene"]
            else:
                ## 随机获取两个基因
                allGeneNum = len(Gene_list)
                list1 = this.getRandList(Gene_index,allGeneNum -1)
                Gene_1 = [Gene_list.ix[list1[0],"gene"],Gene_list.ix[Gene_index,"gene"],Gene_list.ix[list1[1],"gene"]]

                list2 = this.getRandList(Gene_index,allGeneNum -1)
                Gene_2 = [Gene_list.ix[list2[0],"gene"],Gene_list.ix[Gene_index,"gene"],Gene_list.ix[list2[1],"gene"]]

                list3 = this.getRandList(Gene_index,allGeneNum -1)
                Gene_3 = [Gene_list.ix[list3[0],"gene"],Gene_list.ix[Gene_index,"gene"],Gene_list.ix[list3[1],"gene"]]

            Gene_1_result = this.gmmClustering(df_adj,row["#Gene"],row["Exon"],row["BinsResult"],list(Gene_1))
            Gene_2_result = this.gmmClustering(df_adj,row["#Gene"],row["Exon"],row["BinsResult"],list(Gene_2))
            Gene_3_result = this.gmmClustering(df_adj,row["#Gene"],row["Exon"],row["BinsResult"],list(Gene_3))

            num = [Gene_1_result,Gene_2_result,Gene_3_result].count(row["BinsResult"])
            if num == 3:
                row["Clustering"] = row["BinsResult"]
                row["Clustering_Info"] = ",".join([Gene_1_result,Gene_2_result,Gene_3_result])
            else:
                row["Clustering"] = "Non-CNV"
                row["Clustering_Info"] = ",".join([Gene_1_result,Gene_2_result,Gene_3_result])

            if row["Clustering"] == "Non-CNV":
                row['WinsResult'] = "Non-CNV"
                df = df.append(row, ignore_index=True)
                continue
            else:
                chorm, start, end = row["Chr"],row["Start"],row["End"]
                row['WinsResult'], binsNum = this.checkExonCnvByWins(tumor, normal, tumorAllDepth, normalAllDepth,row["Clustering"], chorm, int(start), int(end))
                df = df.append(row,ignore_index=True)

        df.to_csv("%s/%s.cnv.result.txt"%(outdir,sample),sep="\t",index=False)
        svgfile = "%s/%s.cnv.svg"%(outdir,sample)
        visualization.Svg("%s/%s.cnv.result.txt"%(outdir,sample), svgfile)

    def getRandList(this,gene_index,allNum):
        List = []
        while True:
            index = random.randint(0,allNum)
            if index != gene_index and index not in List:
                List.append(index)

            if len(List) == 2:
                break
        return List


    def runWorkflow(this):

        if this.options.normal is not None:
            tumor = this.dealBamObject.getBedCoverage(this.options.tumor,this.options.region)
            normal = this.dealBamObject.getBedCoverage(this.options.normal,this.options.region)
            if tumor is not False and normal is not False:
                if this.options.sample is None:
                    this.options.sample = os.path.basename(this.options.tumor).replace(".bam","")
                this.callingCNVByReadDepthMethod(tumor, normal, this.options.sample, this.options.outdir)
            else:
                print("Fatal error: get depth failed!!!")
                sys.exit(2)
        else:
            allSampleCover = createDict()
            allSampleCover = this.dealBamObject.getNormalCompare(this.options.tumor,this.options.region,this.options.outdir,int(this.options.thread))
            normal = allSampleCover['Normal']
            pool = multiprocessing.Pool(processes=int(this.options.thread))
            with open(this.options.tumor,'r') as fileHandle:
                for tumorFullPath in fileHandle.readlines():
                    tumorFullPath = tumorFullPath.rstrip('\n')
                    sampleName = os.path.basename(tumorFullPath).replace(".bam","")
                    if sampleName in allSampleCover.keys():
                        tumor = allSampleCover[sampleName]
                    else:
                        print("error: %s depth Statistics Failed!!!"%(tumorFullPath))
                        continue
                    pool.apply_async(this.callingCNVByReadDepthMethod, (tumor,normal,sampleName,this.options.outdir))
                    this.callingCNVByReadDepthMethod(tumor,normal,sampleName,this.options.outdir)
            pool.close()
            pool.join()
            pool.terminate()


    def log2toCnvStatu(this,logValue):

        if logValue <= -2: # log((0+0.5)/2)
            return 0,"Del"
        elif logValue > -2 and logValue <= -0.42: # log((1+0.5)/2)
            return 1,"Del"
        elif logValue > -0.42 and logValue <= 0.32: # log((2+0.5)/2)
            return 2,"Non-CNV"
        elif logValue > 0.32 and logValue <= 0.81: # log((3+0.5)/2)
            return 3,"Dup"
        elif logValue >= 0.81 and logValue <= 1.17: # log((4+0.5)/2)
            return 4,"Dup"
        elif logValue > 1.17:  # log((5+0.5)/2)
            return 5,"Dup"
        return 2,"Non-CNV"
