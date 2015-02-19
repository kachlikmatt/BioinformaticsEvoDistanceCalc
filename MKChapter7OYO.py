'''
Name: Matthew Kachlik
Due: 2/19/15
Description: This program Determines Evolutionary distance based off user choice of distance algorithm 
and joining algorithm
'''
#imports used for function calls
import collections
from collections import defaultdict
import sys
import math
'''
	Jukes
	measures evolutionary distance between two nucleotide sequences
'''
def Jukes(seq1,seq2):

	difCtr = float(0)
	lenCtr = float(0)

	#Count differences between sequences ignoring gaps
	for i in range(0, len(seq2)):
		if(seq1[i] != '-' and seq2[i] != '-'):
			lenCtr += 1
			if(seq1[i] != seq2[i]):
				difCtr += 1

	#Calculate and output results

	difCtr = difCtr/lenCtr
	distance = (-3.0/4) * math.log(1-(4.0/3*difCtr))


	return distance
'''
	Kimura
	uses transitions and transversions to add to find 
	the evolutionary distances between species 
'''
def  Kimura(seq1,seq2,aDict):
	s = float(0)
	v = float(0)
	trans = float(0)
	sub = float(0)
	lenCtr = float(0)
	kimura = float(0)
	for i in range(0, len(seq2)):
		if(seq1[i] != '-' and seq2[i] != '-'):
				lenCtr += 1
				#checks the dict whether its a transition or transversion
				if(seq1[i] != seq2[i]):
					if(aDict[seq1[i]][seq2[i]] == 1):
						trans += 1

					elif(aDict[seq1[i]][seq2[i]] == 2):
						sub += 1


	s = float(trans/lenCtr)
	v = float(sub/lenCtr)
	kimura = 1.0/2 * math.log(1/(1-2*s-v)) + 1.0/4*math.log(1/(1-2*v))
	return kimura
'''
	Tamura 
	finds the evolutionary distances between two species 
	using the gc content along with transitions and transversions
	then uses the algorithm to determine evolutionary distances
'''
def Tamura(seq1,seq2,aDict):
	s = float(0)
	v = float(0)
	trans = float(0)
	sub = float(0)
	lenCtr = float(0)
	tamura = float(0)
	gc1 = 0.0
	gc2 = 0.0
	c = float(0)
	for i in range(0, len(seq2)):
		if(seq1[i] != '-' and seq2[i] != '-'):
			lenCtr += 1
			#calculating g c content
			if(seq1[i] == 'G' or seq1[i] == 'C'):
				gc1 += 1
			if(seq2[i] == 'G' or seq2[i] == 'C'):
				gc2 += 1
			#checks the dict whether its a transition or transversion
			if(seq1[i] != seq2[i]):
				if(aDict[seq1[i]][seq2[i]] == 1):
					trans += 1
				elif(aDict[seq1[i]][seq2[i]] == 2):
					sub += 1

	gc1 =  float(gc1/len(seq1))	
	gc2 =  float(gc2/len(seq2))	
	c = gc1 + gc2 - 2 * gc1 * gc2
	s = trans/lenCtr
	v = sub/lenCtr
	tamura = (c*-1.0)*math.log(1-s/c-v) - 1.0/2 *(1-c)*math.log(1-2*v)
	return tamura
'''
shortest
calculates the shortest distance between two clusters and returns the clusters' names
'''
def shortest(clusterDict,shortestI,shortestJ,lowest):
		#finding smallest
	for k,v in clusterDict.items():
		for k2,v2 in v.items():
			if(k != k2):
				if(int(clusterDict[k][k2]) < int(lowest)):
					lowest = int(clusterDict[k][k2])
					shortestI = k
					shortestJ = k2
	return shortestI,shortestJ
'''
	singleLinkage
This method calculates distances between the new cluster
and all other clusters using single linkage 
'''
def singleLinkage(clusterDict,newClusterName,originalDict,clusterNames):
	for cluster in clusterNames:
		smallestD = sys.maxint
		for c1 in cluster:
			for c2 in newClusterName:
				if(originalDict[c1][c2] < smallestD):
					smallestD = originalDict[c1][c2]
		clusterDict[newClusterName][cluster] = smallestD
		clusterDict[cluster][newClusterName] = smallestD	
	return clusterDict
'''
	neighborJoining
	uses the neighbor joing Algorithm to find the optimal clustering
'''
def neighborJoining(clusterDict,originalDict,clusterNames,numSeq):
	
	transition = collections.defaultdict(lambda:defaultdict(float))
	rValues = collections.defaultdict(float)
	distance = int(0)
	newClusterName = ""
	newick = {}
	prevNewCluster = ""
	branchLength = collections.defaultdict(float)
	while(numSeq > 2):

		#finding rvalues for all the arrays by finding the average distance from each cluster
		for c1 in clusterNames:
			distance = 0
			for c2 in clusterNames:
				if(c1 != '' and c2 != '' and c1 != c2):
					distance += float(clusterDict[c1][c2])
			rValues[c1] = distance/float(numSeq-2)	

		#finding transition values
		for c1 in clusterNames:
			for c2 in clusterNames:
				if(c1 != '' and c2 != '' and c1 != c2):
					transition[c1][c2] = float(clusterDict[c1][c2]) - rValues[c1] - rValues[c2]	
		
		#finding the two to merge based on being the closet
		shortestI = clusterNames[0]
		shortestJ =	clusterNames[1]
		lowest = clusterDict[clusterNames[0]][clusterNames[1]]

		
		shortestI,shortestJ = shortest(clusterDict,shortestI,shortestJ,lowest)
		newClusterName = shortestI + shortestJ



		#removing the shortestJ and shortestI from the cluster names so they are not selected again
		clusterNames.remove(shortestI)
		clusterNames.remove(shortestJ)



		

		#adding previous new cluster to the array then calling neighborJoing
		prevClusterNames.append(prevNewCluster)

		#calculating new distance for the merged points and the remaining clusters
		for c1 in clusterNames:
			clusterDict[newClusterName][c1] = (float(clusterDict[shortestI][c1]) + float(clusterDict[shortestJ][c1]) - float(clusterDict[shortestI][shortestJ]))/2.0
			clusterDict[c1][newClusterName] = (float(clusterDict[shortestI][c1]) + float(clusterDict[shortestJ][c1]) - float(clusterDict[shortestI][shortestJ]))/2.0
							
		
			
		#calculating the branch lengths
		branchLength[shortestI] = (clusterDict[shortestI][shortestJ] + rValues[shortestI] - rValues[shortestJ])/2.0
		branchLength[shortestJ] = (clusterDict[shortestI][shortestJ] + rValues[shortestJ] - rValues[shortestI])/2.0
		
		'''
		this if checks if both are in the the newick already
		then adds them together then deletes the smaller
		or will check if each individual part is in there and deleting 
		finally if no of those cases are true it adds the shortestI and shortestJ
		without a deletion
		'''
		if(shortestI in newick and shortestJ in newick):
			newick[newClusterName] = '('+shortestI+':'+str(branchLength[shortestI])+','+shortestJ+":"+str(branchLength[shortestJ])+')'
			del newick[shortestJ]
			del newick[shortestI]
		elif(shortestI in newick):
			newick[newClusterName] = '('+shortestI+':'+str(branchLength[shortestI])+','+shortestJ+":"+str(branchLength[shortestJ])+')'
			del newick[shortestI]
		elif(shortestJ in newick):
			newick[newClusterName] = '('+shortestI+':'+str(branchLength[shortestI])+','+shortestJ+":"+str(branchLength[shortestJ])+')'
			del newick[shortestJ] 
		else:
			newick[newClusterName] = '('+shortestI+':'+str(branchLength[shortestI])+','+shortestJ+":"+str(branchLength[shortestJ])+')'

		#deleting the old clusters
		for k,v in clusterDict.items():
			for k2,v2 in v.items():
				if(k2 == shortestI):
					del clusterDict[k][k2]
				elif(k2 == shortestJ):
					del clusterDict[k][k2]
		for k,v in clusterDict.items():
			if(k == shortestI):
				del clusterDict[k]
			elif(k == shortestJ):
				del clusterDict[k]

		#adding the new cluster
		clusterNames.append(newClusterName)
		prevNewCluster = newClusterName

		print "merging clusters " + shortestI + " " +str(branchLength[shortestI])+","+ shortestJ + " " +str(branchLength[shortestJ])

		numSeq -= 1
		prevClusterNames.remove(shortestI)
		prevClusterNames.remove(shortestJ)
	#checking for a single cluster remaining at the end hardcoded since the newick adds to the end
	if(len(newick) == 1):
		newick[clusterNames[0]] = "(" + clusterNames[0] + ")"
	print "remaining clusters:"
	branchLength[clusterNames[0]] = (clusterDict[clusterNames[0]][clusterNames[1]] + rValues[clusterNames[0]] - rValues[clusterNames[1]])/2.0
	branchLength[clusterNames[1]] = (clusterDict[clusterNames[0]][clusterNames[1]] + rValues[clusterNames[1]] - rValues[clusterNames[0]])/2.0
	print clusterNames[0] +"," + clusterNames[1]
	print "(" +newick[clusterNames[0]]+":"+str(branchLength[clusterNames[0]]) +","+ newick[clusterNames[1]]+":"+str(branchLength[clusterNames[1]])+")"


	'''
	Agglomerative
		uses an agglomerative clustering Algorithm with single linkage to find the optimal
		clustering
	'''
def agglomerativeClustering(clusterDict,originalDict,clusterNames,numSeq):
	#making a dictionary for the newick format to be built
	newick = {}
	prevNewCluster = ""
	while(numSeq > 2):
			
		#finding the first two values to start the lowest value
		shortestI = clusterNames[0]
		shortestJ =	clusterNames[1]
		lowest = clusterDict[clusterNames[0]][clusterNames[1]]

		shortestI,shortestJ = shortest(clusterDict,shortestI,shortestJ,lowest)


		newClusterName = shortestI + shortestJ
		

		#removing the connections from shortestI and shortestJ to start the merge
		for k,v in clusterDict.items():
			for k2,v2 in v.items():
				if(k2 == shortestI):
					del clusterDict[k][k2]
				elif(k2 == shortestJ):
					del clusterDict[k][k2]
			

		#removing the shortestJ and shortestI from the cluster names so they are not selected again
		clusterNames.remove(shortestI)
		clusterNames.remove(shortestJ)


		#making newick format from keys

		'''
		this if checks if both are in the the newick already
		then adds them together then deletes the smaller
		or will check if each individual part is in there and deleting 
		finally if no of those cases are true it adds the shortestI and shortestJ
		without a deletion
		'''
		if(shortestI in newick and shortestJ in newick):
			newick[newClusterName] = '('+newick[shortestI]+','+newick[shortestJ]+')'
			del newick[shortestJ]
			del newick[shortestI]
		elif(shortestI in newick):
			newick[newClusterName] = "("+newick[shortestI]+","+shortestJ+")"
			del newick[shortestI]
		elif(shortestJ in newick):
			newick[newClusterName] = "("+newick[shortestJ]+","+shortestI+")"
			del newick[shortestJ] 
		else:
			newick[newClusterName] = '('+shortestI+','+shortestJ+')'

		#adding previous new cluster to the array then calling neighborJoing
		prevClusterNames.append(prevNewCluster)
		
		
		
		singleLinkage(clusterDict,newClusterName,originalDict,clusterNames)

		#deleting the rest of the shortestI and shortestJ to complete the merge
		for k,v in clusterDict.items():
			if(k == shortestI):
				del clusterDict[k]
			elif(k == shortestJ):
				del clusterDict[k]

		#adding the new cluster and saving it for NJ
		clusterNames.append(newClusterName)
		prevNewCluster = newClusterName
		#printing the merged clusters
		print "merging clusters " + shortestI +","+ shortestJ
		numSeq -= 1
		prevClusterNames.remove(shortestI)
		prevClusterNames.remove(shortestJ)

	print "remaining clusters:"
	print clusterNames[0] +"," + clusterNames[1]
	print newick[clusterNames[0]] + newick[clusterNames[1]]

alignDict = collections.defaultdict(dict)
alignDict = collections.defaultdict(lambda:defaultdict(int))

#transversions
alignDict["G"]["A"] = 1
alignDict["A"]["G"] = 1
alignDict["C"]["T"] = 1
alignDict["T"]["C"] = 1

#transitions 
alignDict["A"]["C"] = 2
alignDict["C"]["A"] = 2
alignDict["G"]["T"] = 2
alignDict["T"]["G"] = 2
alignDict["A"]["T"] = 2
alignDict["T"]["A"] = 2
alignDict["C"]["G"] = 2
alignDict["G"]["C"] = 2
#making arrays to hold the clusterNames previousClusterNames and setting the number of sequences to 0
clusterNames = []
prevClusterNames = []
originalClusterName = []
sequences = []
numSeq = 0
seq = ""
if1 = open(raw_input("Enter the file name in fasta format "),"r")






for line in if1:
	line = line.replace('\n', '')
	line = line.replace('\r','')
	
	if(line[0] == ">"):
		numSeq += 1
		sequences.append(seq) 
		seq = ""
		originalClusterName.append(line)
		clusterNames.append(line)
		prevClusterNames.append(line)
	else:
		seq += line
sequences.append(seq) 
if1.close()
#removing blank sequenc
sequences = filter(None, sequences)
#sequences.pop(0)
print sequences
originalDict = collections.defaultdict(lambda:defaultdict(int))

clusterDict = collections.defaultdict(lambda:defaultdict(int))
'''
making two dictionaries one to be the original to save the value
the other to manipulate 
'''
distanceAlgo = raw_input("Enter which distance algorithm you would like to use\n1)Jukes-Cantor\n2)Tamura\n3)Kimura\n(invalid input will use Kimura)\n")
if(distanceAlgo == "1"):
	for i in range(0,numSeq):
		for j in range(0,numSeq):
			distance = Jukes(sequences[i],sequences[j])
			originalDict[clusterNames[i]][clusterNames[j]] = distance
			clusterDict[clusterNames[i]][clusterNames[j]] = distance
elif(distanceAlgo == "2"):
	for i in range(0,numSeq):
		for j in range(0,numSeq):
			distance = Tamura(sequences[i],sequences[j],alignDict)
			originalDict[clusterNames[i]][clusterNames[j]] = distance
			clusterDict[clusterNames[i]][clusterNames[j]] = distance
else:
	for i in range(0,numSeq):
		for j in range(0,numSeq):
			distance = Kimura(sequences[i],sequences[j],alignDict)
			originalDict[clusterNames[i]][clusterNames[j]] = distance
			clusterDict[clusterNames[i]][clusterNames[j]] = distance
clusteringAlgo = "1"
#prompting user to select which clustering algorithm they want to use
while(int(clusteringAlgo) > 0 and int(clusteringAlgo) < 3):
	clusteringAlgo = raw_input("Enter which clustering algorithm you would like to use\n1)Agglomerative\n2)Neighbor Joining\n3)to quit(this will also quit on an invalid entry)\n")
	if(clusteringAlgo == "1"):
		agglomerativeClustering(clusterDict,originalDict,clusterNames,numSeq)
	elif(clusteringAlgo == "2"):
		neighborJoining(clusterDict,originalDict,clusterNames,numSeq)
	else:
		sys.exit(0)
	#reset for another run
	clusterNames = []
	prevClusterNames = []
	for i in range(0,numSeq):
		clusterNames.append(originalClusterName[i])
		prevClusterNames.append(originalClusterName[i])
		for j in range(0,numSeq):
			clusterDict[originalClusterName[i]][originalClusterName[j]] = originalDict[originalClusterName[i]][originalClusterName[j]]

	