binRCfile = "/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.100000_fixed.txt.gz"
BRfile = "/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.100000_fixed.many.txt"
infoFile = "/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.100000_fixed.info"
stateFile = "/home/mgharegh/research/HDhackathon/data/skin/D2Rfb.final.txt"
outputDir = "/home/mgharegh/research/HDhackathon/data/skin/"
K = 22
maximumCN = 5
bin.size = 100000
haplotypInfo=F

binRC = splitChromosomes(changeRCformat(binRCfile, outputDir))
cellTypes = changeCellTypesFormat(stateFile)
NBparams = changeNBparamsFormat(infoFile, K)
p = NBparams[[1]]
r = NBparams[[2]]
segmentsCounts = getSegReadCounts(binRC, BRfile, K, bin.size)

SVcalling_wrapperFunc(bin.size, K, maximumCN, segmentsCounts, r, p, cellTypes, outputDir, hapMode = haplotypInfo)
