
import os,re,sys
import pandas as pd
import matplotlib.pyplot as plt
from baseSharedMethods import createDict

def Svg(infile, svg):

  df = pd.read_csv(infile, sep="\t")
  df = df[~ df["Chr"].isin(["chrX","chrY","chrM"])]
  df["Chr"] = df["Chr"].str.replace("chr","")
  df["Chr"] = df["Chr"].astype(int)
  
  df = df.sort_values(by=["Chr","Start"],axis=0).reset_index(drop=True)
  df["Index"] = df.index

  df_dup = df[df["WinsResult"] == "Dup"]
  df_non_cnv = df[df["WinsResult"] == "Non-CNV"]
  df_del = df[df["WinsResult"] == "Del"]

  plt.figure(figsize=(150,10),dpi=80)
  plt.subplots_adjust(left=0.05, bottom=0.2, right=0.95, top=0.8)

  axs0 = plt.subplot(111,facecolor='#FCFCFC')

  exon_nums = len(df)

  axs0.scatter(df_non_cnv["Index"].tolist(), df_non_cnv["log2AdjValue"].tolist(), marker=".", s=80, color="blue", label="Non-CNV")
  axs0.scatter(df_dup["Index"].tolist(), df_dup["log2AdjValue"].tolist(), marker="*", s=80, color="green", label="Dup")
  axs0.scatter(df_del["Index"].tolist(), df_del["log2AdjValue"].tolist(), marker="p", s=40, color="red", label="Del")

  CNV_y = createDict()
  CNV_text = createDict()
  Hash = createDict()
  Hash_ticks = {}
  for index in df.index:
    row = df.loc[index,:].to_dict()

    df_gene = df[df["#Gene"] == row["#Gene"]]
    df_gene = df_gene.reset_index(drop=True)
    last_index = df_gene.at[len(df_gene)-1,"Index"]
    Chr = row["Chr"]
    Hash_ticks[Chr] = last_index

    if row["#Gene"] not in Hash.keys():
      if row["WinsResult"] == "Non-CNV": continue
      Hash[row["#Gene"]]["Start"] = row["Index"]
      Hash[row["#Gene"]]["End"] = row["Index"]
      Hash[row["#Gene"]]["CNV"] = row["WinsResult"]
      Hash[row["#Gene"]]["annote_x"] = row["Index"] + (Hash[row["#Gene"]]["End"] - Hash[row["#Gene"]]["Start"])/2
    else:
      if row["WinsResult"] == Hash[row["#Gene"]]["CNV"]:
        Hash[row["#Gene"]]["End"] = row["Index"]
        Hash[row["#Gene"]]["annote_x"] = Hash[row["#Gene"]]["Start"] + (Hash[row["#Gene"]]["End"] - Hash[row["#Gene"]]["Start"])/2
      else:
        
        text_x = Hash[row["#Gene"]]["annote_x"]

        if row["#Gene"] in CNV_y.keys() and Hash[row["#Gene"]]["CNV"] in CNV_y[row["#Gene"]].keys():
          if Hash[row["#Gene"]]["CNV"] == "Dup":
            CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]] += 0.1
          else:
            CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]] -= 0.1
          text_y = CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]]  
        else:
          text_y = max([df.at[Hash[row["#Gene"]]["Start"],"log2AdjValue"], df.at[Hash[row["#Gene"]]["End"],"log2AdjValue"]]) 
          text_y = text_y + 0.4 if Hash[row["#Gene"]]["CNV"] == "Dup" else text_y -0.4
          CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]] = text_y

        axs0.annotate(row["#Gene"],xy=(Hash[row["#Gene"]]["Start"],df.at[Hash[row["#Gene"]]["Start"],"log2AdjValue"]),xytext=(text_x,text_y),arrowprops=dict(arrowstyle='->',connectionstyle='arc3',color='red'),fontsize=8)

        Exon_Start = df.at[Hash[row["#Gene"]]["Start"],"Exon"]
        Exon_End = df.at[Hash[row["#Gene"]]["End"],"Exon"]

        if Hash[row["#Gene"]]["Start"] != Hash[row["#Gene"]]["End"]:
          axs0.annotate(row["#Gene"],xy=(Hash[row["#Gene"]]["End"],df.at[Hash[row["#Gene"]]["End"],"log2AdjValue"]),xytext=(text_x,text_y),arrowprops=dict(arrowstyle='->',connectionstyle='arc3',color='red'),fontsize=8)

        if Exon_Start == Exon_End:
          CNV_text[Hash[row["#Gene"]]["CNV"]].setdefault(row["#Gene"],[]).append(Exon_Start)
        else:
          begin = min([Exon_Start,Exon_End])
          end = max([Exon_Start,Exon_End])
          CNV_text[Hash[row["#Gene"]]["CNV"]].setdefault(row["#Gene"],[]).extend(list(range(begin,end+1)))
        
        del Hash[row["#Gene"]]

        if row["WinsResult"] != "Non-CNV":
          Hash[row["#Gene"]]["Start"] = row["Index"]
          Hash[row["#Gene"]]["End"] = row["Index"]
          Hash[row["#Gene"]]["CNV"] = row["WinsResult"]
          Hash[row["#Gene"]]["annote_x"] = row["Index"] + (Hash[row["#Gene"]]["End"] - Hash[row["#Gene"]]["Start"] + 1)/2
        else:
          continue

    if last_index == row["Index"]:
      text_x = Hash[row["#Gene"]]["annote_x"]

      if row["#Gene"] in CNV_y.keys() and Hash[row["#Gene"]]["CNV"] in CNV_y[row["#Gene"]].keys():
        if Hash[row["#Gene"]]["CNV"] == "Dup":
          CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]] += 0.1
        else:
          CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]] -= 0.1
        text_y = CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]]
      else:
        text_y = max([df.at[Hash[row["#Gene"]]["Start"],"log2AdjValue"], df.at[Hash[row["#Gene"]]["End"],"log2AdjValue"]])
        text_y = text_y + 0.4  if Hash[row["#Gene"]]["CNV"] == "Dup" else text_y - 0.4
        CNV_y[row["#Gene"]][Hash[row["#Gene"]]["CNV"]] = text_y

      axs0.annotate(row["#Gene"],xy=(Hash[row["#Gene"]]["Start"],df.at[Hash[row["#Gene"]]["Start"],"log2AdjValue"]),xytext=(text_x + 2,text_y),arrowprops=dict(arrowstyle='->',connectionstyle='arc3',color='red'),fontsize=8)

      Exon_Start = df.at[Hash[row["#Gene"]]["Start"],"Exon"]
      Exon_End = df.at[Hash[row["#Gene"]]["End"],"Exon"]

      if Exon_Start == Exon_End:
        CNV_text[Hash[row["#Gene"]]["CNV"]].setdefault(row["#Gene"],[]).append(Exon_Start)
      else:
        begin = min([Exon_Start,Exon_End])
        end = max([Exon_Start,Exon_End])
        CNV_text[Hash[row["#Gene"]]["CNV"]].setdefault(row["#Gene"],[]).extend(list(range(begin,end+1)))

      if Hash[row["#Gene"]]["Start"] != Hash[row["#Gene"]]["End"]:
        axs0.annotate(row["#Gene"],xy=(Hash[row["#Gene"]]["End"],df.at[Hash[row["#Gene"]]["End"],"log2AdjValue"]),xytext=(text_x + 2,text_y),arrowprops=dict(arrowstyle='->',connectionstyle='arc3',color='red'),fontsize=8)
      del Hash[row["#Gene"]]

  axs0.set_xlim(0,exon_nums)
  axs0.set_ylim(-3,3)

  axs0.plot([0, exon_nums], [0.32, 0.32], linestyle="--",color="green")
  axs0.plot([0, exon_nums], [-0.42, -0.42], linestyle="--",color="red")

  x_coor = []
  x_ticks = []
  for Chr in sorted(Hash_ticks.keys()):
    x_coor.append(Hash_ticks[Chr]) 
    Chr = "chr" + str(Chr)
    x_ticks.append(Chr)

  axs0.set_xticks(x_coor)
  axs0.set_xticklabels(x_ticks,rotation=60)

  axs0.set_ylabel("Exon log2 value")

  axs0.grid(axis='x',linestyle='--',color='#E0E0E0',linewidth=1)

  axs0.legend(ncol=3, loc="upper left")

  DelString = "Candidate Deletion:"
  DupString = "Candidate Duplication:"
  for cnv in CNV_text.keys():
    
    for gene in CNV_text[cnv].keys():
      exon = sorted(CNV_text[cnv][gene])
      exon = [str(i) for i in exon]
      exon = ", ".join(exon)
      if cnv == "Del":
        DelString += "%s Exon [%s] | "%(gene,exon)
      else:
        DupString += "%s Exon [%s] | "%(gene,exon)

  if DelString != "Candidate Deletion:":
    DelString = re.sub(r'\| $','',DelString)
    plt.figtext(0.05,0.04,"%s\n\n"%DelString,color="red",fontsize=16)

  if DupString != "Candidate Duplication:":
    DupString = re.sub(r'\| $','',DupString)
    plt.figtext(0.05,0.08,"%s\n\n"%DupString,color="green",fontsize=16)

  plt.savefig(svg)
