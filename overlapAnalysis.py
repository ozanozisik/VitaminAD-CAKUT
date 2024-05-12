"""
Created on Mon Jan 25 10:20:47 2021

@author: Ozan
"""

import os
import sys
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import random


def readGmtFileReturnList(gmtPath):
	'''
	Reads GMT file.
	Each line consists of gene set ID, gene set name and genes,
	all tab separated, no header.

	Returns a list of gene set entries, each entry is a list containing
	gene set ID, gene set name, set of genes.
	'''
	geneSetEntries=[]
	try:
		f=open(gmtPath, 'r')
		lines=f.readlines()
		for line in lines:
			tokens=line.strip().split('\t')
			gsID=tokens[0]
			gsName=tokens[1]
			genes=set(tokens[2:])
			geneSetEntries.append([gsID, gsName, genes])
		f.close()
	except IOError:
		print("I/O error while reading gmt file.")
	return geneSetEntries


def readGmtBackground(gmtPath):
	'''
	Reads GMT file.
	Creates background information to be used for enrichment.
	Each line consists of gene set ID, gene set name and genes,
	all tab separated, no header.
	'''
	backgroundGeneSet=set()
	try:
		f=open(gmtPath, 'r')
		lines=f.readlines()
		for line in lines:
			tokens=line.strip().split('\t')
			#gsID=tokens[0]
			#gsName=tokens[1]
			genes=set(tokens[2:])
			backgroundGeneSet.update(genes)
		f.close()
	except IOError:
		print("I/O error while reading gmt file.")
	return backgroundGeneSet



def csvToTex(fileName):
	'''
	Code to convert csv table to tex table.
	Special characters are not handled because in our case there is not any.
	'''

	f=open(fileName, 'r')
	lines=f.readlines()
	f.close()

	fW=open(os.path.splitext(fileName)[0]+'TexVersion.txt', 'w')
	fW.write('\\begin{table}[h!]\n')
	fW.write('\\centering\n')
	fW.write('\\caption{}\n')


	colNum=len(lines[0].strip().split('\t'))
	fW.write('\\begin{tabular}{|')
	for i in range(colNum):
		fW.write('c |')
	fW.write('}\n')
	fW.write('\\hline\n')

	for line in lines:
		#Replaces tab and new line characters with the commands required in tex
		#Code can be added to replace special characters
		#e.g. line=line.replace('>',	'$>$'), line=line.replace('#',	'\#')

		line=line.replace('\t',	' & ')
		line=line.replace('\n',	' \\\\ \hline\n')
		fW.write(line)

	fW.write('\\end{tabular}\n')
	fW.write('\\end{table}\n')

	fW.close()


def createTablesForManuscript(df1, df2, output):
	'''
	In the analyses we used two sources for both vitamin A and vitamin D,
	one from CTD, one from a publication.
	This function creates a table that contains the results for both
	vitamin target lists from two sources.
	'''

	cols=df2.columns
	cols=[c+'2' for c in cols]
	df2.columns=cols
	dfMerged=pd.concat([df1,df2],axis=1)
	dfMerged.drop(['TermSize','QuerySize', 'IntersectionSize', 'DomainSize', 'pValue'], axis=1, inplace=True)
	dfMerged.drop(['TermId2','TermName2', 'TermSize2','QuerySize2', 'IntersectionSize2', 'DomainSize2', 'pValue2'], axis=1, inplace=True)


	toScNot=lambda flt: str('(p-value = {:.2e})'.format(flt))


	dfMerged['pAdj and intersection']=dfMerged['pAdjusted'].map(toScNot) + ' ' + dfMerged['Intersection']
	dfMerged['pAdj2 and intersection2']=dfMerged['pAdjusted2'].map(toScNot) + ' ' +dfMerged['Intersection2']

	dfMerged.drop(['pAdjusted', 'pAdjusted2', 'Intersection', 'Intersection2'], axis=1, inplace=True)

	dfMerged.to_csv(output, sep='\t', float_format='%.2e', index=False)




def analysisHypergeometric(targetGenePaths, geneSetEntries, backgroundSetsDict, outputFolder):
	'''
	For all target lists (e.g. vitamin A target genes from CTD, 
	vitamin A target genes from publication) and for all gene set
	in question, performs hypergeometric test, applies multiple 
	testing correction by Benjamini-Hochberg method.
	'''

	if not os.path.exists(outputFolder):
		os.mkdir(outputFolder)

	results=[]

	for targetGenePath in targetGenePaths:

		df1=pd.read_csv(targetGenePath, header=None)
		targetSet=set(df1[0].values)

		termIDs=[]
		termNames=[]
		termSizes=[]
		querySizes=[]
		intersectionSizes=[]
		domainSizes=[]
		pValues=[]
		intersections=[]

		for gsNo in range(len(geneSetEntries)):
			geneSetEntry=geneSetEntries[gsNo]
			geneSet=geneSetEntry[2]
			backgroundSet=backgroundSetsDict[geneSetEntry[3]]

			#M is the population size
			#n is the number of successes in the population
			#N is the sample size
			#x is the number of drawn “successes”.

			M=len(backgroundSet)
			n=len(geneSet)
			N=len(targetSet.intersection(backgroundSet))#Taking only genes that are also in background
			intersection=list(geneSet.intersection(targetSet))
			x=len(intersection)
			#print(M, n, N, x)

			pval = hypergeom.sf(x-1, M, n, N)
			#print(geneSetEntry[0], geneSetEntry[1], "{:.2e}".format(pval))

			termIDs.append(geneSetEntry[0])
			termNames.append(geneSetEntry[1])
			termSizes.append(n)
			querySizes.append(N)
			intersectionSizes.append(x)
			domainSizes.append(M)
			pValues.append(pval)
			intersection.sort()
			intersections.append(', '.join(intersection))

		reject, pValsAdj, alphacSidak, alphacBonf = multipletests(pValues, alpha=0.05, method='fdr_bh')

		df=pd.DataFrame({'TermId':termIDs,
						 'TermName':termNames,
						 'TermSize':termSizes,
						 'QuerySize':querySizes,
						 'IntersectionSize':intersectionSizes,
						 'DomainSize':domainSizes,
						 'pValue':pValues,
						 'pAdjusted':pValsAdj,
						 'Intersection':intersections
						 })

		df.to_csv(outputFolder+os.path.splitext(os.path.basename(targetGenePath))[0]+'.csv','\t', index=False, float_format="%.2E")

		results.append(df)


	createTablesForManuscript(results[1].copy(), results[0].copy(), outputFolder+'VitA-Merged.csv')
	#csvToTex(outputFolder+'VitA-Merged.csv')

	createTablesForManuscript(results[2].copy(), results[3].copy(), outputFolder+'VitD-Merged.csv')
	#csvToTex(outputFolder+'VitD-Merged.csv')




def analysisRandomizedSampling(targetGenePaths, geneSetEntries, backgroundSetsDict, outputFolder):
	'''
	This is a supplementary analysis.
	Calculates the significance of overlap using ramdomized sampling instead of hypergeometric test.
	'''

	if not os.path.exists(outputFolder):
		os.mkdir(outputFolder)

	results=[]

	bootstrapNum=2000

	for targetGenePath in targetGenePaths:

		df1=pd.read_csv(targetGenePath, header=None)
		targetSet=set(df1[0].values)

		termIDs=[]
		termNames=[]
		termSizes=[]
		querySizes=[]
		intersectionSizes=[]
		domainSizes=[]
		pValues=[]

		for gsNo in range(len(geneSetEntries)):
			geneSetEntry=geneSetEntries[gsNo]
			geneSet=geneSetEntry[2]
			backgroundSet=backgroundSetsDict[geneSetEntry[3]]
			backgroundList=list(backgroundSet)

			intersectionSize=len(geneSet.intersection(targetSet))

			targetSize=len(targetSet.intersection(backgroundSet))#Taking only genes that are also in background, this is consistent with the hypergeometric test.

			betterOrEqual=0

			for i in range(bootstrapNum):

				targetSetBS=set(random.sample(backgroundList, targetSize))
				intersectionBSSize=len(geneSet.intersection(targetSetBS))

				if intersectionBSSize>=intersectionSize:
					betterOrEqual=betterOrEqual+1

			p=betterOrEqual/bootstrapNum


			termIDs.append(geneSetEntry[0])
			termNames.append(geneSetEntry[1])
			termSizes.append(len(geneSet))
			querySizes.append(targetSize)
			intersectionSizes.append(intersectionSize)
			domainSizes.append(len(backgroundSet))
			pValues.append(p)

		reject, pValsAdj, alphacSidak, alphacBonf = multipletests(pValues, alpha=0.05, method='fdr_bh')

		df=pd.DataFrame({'TermId':termIDs,
						 'TermName':termNames,
						 'TermSize':termSizes,
						 'QuerySize':querySizes,
						 'IntersectionSize':intersectionSizes,
						 'DomainSize':domainSizes,
						 'pValue':pValues,
						 'pAdjusted':pValsAdj,
						 })

		df.to_csv(outputFolder+os.path.splitext(os.path.basename(targetGenePath))[0]+'-RS.csv','\t', index=False, float_format="%.2E")

		results.append(df)




def overlapOfCAKUTCausalGenesWithOthers(geneSetEntries, outputFolder):
	'''
	This is a supplementary analysis.
	Finds the overlap of 'CAKUT causal genes' with other gene sets.
	'''
	if not os.path.exists(outputFolder):
		os.mkdir(outputFolder)

	cakutFound=False
	otherGeneSetEntries=[]

	for geneSetEntry in geneSetEntries:
		if geneSetEntry[1]=='CAKUT causal genes':
			cakutGenes=geneSetEntry[2]
			cakutFound=True
		else:
			otherGeneSetEntries.append(geneSetEntry)

	if cakutFound:

		termIDs=[]
		termNames=[]
		termSizes=[]
		querySizes=[]
		intersectionSizes=[]
		intersections=[]

		for geneSetEntry in otherGeneSetEntries:

			intersection=list(cakutGenes.intersection(geneSetEntry[2]))
			intersection.sort()
			intersectionSize=len(intersection)

			termIDs.append(geneSetEntry[0])
			termNames.append(geneSetEntry[1])
			termSizes.append(len(geneSetEntry[2]))
			querySizes.append(len(cakutGenes))
			intersectionSizes.append(intersectionSize)
			intersections.append(', '.join(intersection))

		df=pd.DataFrame({'TermId':termIDs,
						 'TermName':termNames,
						 'TermSize':termSizes,
						 'QuerySize':querySizes,
						 'IntersectionSize':intersectionSizes,
						 'Intersection':intersections
						 })

		df.to_csv(outputFolder+'CAKUTCausalGenesOverlap.csv','\t', index=False)






#GMT file that consists of pathways of interest
gmtFile='Data/PathwaysOfInterest.gmt'

#File that lists the names of the GMT files that contain all the pathways in
#the domain. These GMT files will be used to calculate background gene number
#and to ignore query genes that are not in the domain
backgroundGmtFilesFile='Data/PathwaysOfInterestBackground.txt'

#Files that contain lists of genes to use for enrichment analysis
targetGenePaths=['Data/VitA-Balmer2002-Genes.txt',
				 'Data/VitA-CTD-Genes.txt',
				 'Data/VitD-CTD-Genes.txt',
				 'Data/VitD-Ramagopalan2010.txt']



#Read GMT file that consists of pathways of interest
geneSetEntries=readGmtFileReturnList(gmtFile)

#Read background/domain files list
df1=pd.read_csv(backgroundGmtFilesFile, header=None)
backgroundGmtFilesList=list(df1[0].values)

if(len(geneSetEntries)!=len(backgroundGmtFilesList)):
	print('Number of gene sets and number of background sets do not match')
	sys.exit()


directory=os.path.dirname(backgroundGmtFilesFile)

#Read background/domain GMT files
backgroundSetsDict=dict()
for gsNo in range(len(geneSetEntries)):
	backgroundGmtFile=backgroundGmtFilesList[gsNo]
	geneSetEntries[gsNo].append(backgroundGmtFile)
	backgroundSetsDict[backgroundGmtFile]=readGmtBackground(directory+'/'+backgroundGmtFile)


analysisHypergeometric(targetGenePaths, geneSetEntries, backgroundSetsDict, 'Result/')




####################################################
### Supplementary analyses


## Calculates the significance of overlap using ramdomized sampling instead of hypergeometric test.
#analysisRandomizedSampling(targetGenePaths, geneSetEntries, backgroundSetsDict, 'RandomizedSampling/')

## Finds the overlap of 'CAKUT causal genes' with other gene sets.
#overlapOfCAKUTCausalGenesWithOthers(geneSetEntries, './')


####################################################
