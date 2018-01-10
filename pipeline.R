library(devtools) 
install_git("git://github.com/friendsofstrandseq/MaRyam.git", branch = "master")
library("MaRyam")

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#print(args)
#print(class(args))
args = as.data.frame(strsplit(args, split = "="), stringsAsFactors = F)
#print(args)

binRCfile = args[2,match("binRCfile", as.character(args[1,]))]
BRfile = args[2,match("BRfile", as.character(args[1,]))]
infoFile = args[2,match("infoFile", as.character(args[1,]))]
stateFile = args[2,match("stateFile", as.character(args[1,]))]
outputDir = args[2,match("outputDir", as.character(args[1,]))]
bin.size = as.numeric(args[2,match("bin.size", as.character(args[1,]))])
K = as.numeric(args[2,match("K", as.character(args[1,]))])
maximumCN = as.numeric(args[2,match("maximumCN", as.character(args[1,]))])
haplotypInfo=F
if (any(as.character(args[1,])=="haplotypeInfo")){haplotypInfo = T}

print(paste("binRCfile =", binRCfile))
print(paste("BRfile =", BRfile))
print(paste("infoFile =", infoFile))
print(paste("stateFile =", stateFile))
print(paste("outputDir =", outputDir))
print(paste("bin.size =", bin.size))
print(paste("K =", K))
print(paste("maximumCN =", maximumCN))

binRC = splitChromosomes(changeRCformat(binRCfile, outputDir))
cellTypes = changeCellTypesFormat(stateFile)
NBparams = changeNBparamsFormat(infoFile, K)
p = NBparams[[1]]
r = NBparams[[2]]
segmentsCounts = getSegReadCounts(binRC, BRfile, K, bin.size)

SVcalling_wrapperFunc(bin.size, K, maximumCN, segmentsCounts, r, p, cellTypes, outputDir, hapMode = haplotypInfo)


#dir = "/home/maryam/research/HDhackathon/data/test_small_data/"
#Rscript pipeline.R binRCfile="/home/maryam/research/HDhackathon/data/test_small_data/D2Rfb.100000_fixed.txt.gz" BRfile="/home/maryam/research/HDhackathon/data/test_small_data/few_brs.txt" infoFile="/home/maryam/research/HDhackathon/data/test_small_data/D2Rfb.100000_fixed.info" stateFile="/home/maryam/research/HDhackathon/data/test_small_data/D2Rfb.final.txt" outputDir="/home/maryam/research/HDhackathon/data/test_small_data/" bin.size=100000 K=22 maximumCN=4
# binRCfile="/home/maryam/research/HDhackathon/data/test_small_data/D2Rfb.100000_fixed.txt.gz"
# BRfile="/home/maryam/research/HDhackathon/data/test_small_data/few_brs.txt"
# infoFile="/home/maryam/research/HDhackathon/data/test_small_data/D2Rfb.100000_fixed.info"
# stateFile="/home/maryam/research/HDhackathon/data/test_small_data/D2Rfb.final.txt"
# outputDir="/home/maryam/research/HDhackathon/data/test_small_data/"
# bin.size=100000
# K=22
# maximumCN=4