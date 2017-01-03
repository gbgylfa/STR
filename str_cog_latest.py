#lang=python
import os
import csv
from os.path import basename

#script to get COG of /s/ and /sh/ phonemes 
#futher codes instances of /s/ as being 'str' if the following two segments are /tr/
#assumes that phonemes are annotated in all caps (FAAV default)
#does this for all wav files with annotated textgrids of the same name in a root directory
#each wav file its accompanying textgrid must be in the same base directory
#outputs to a csv file 

#######Praat script by Eric Doty####
#####adapted to Python by Gudrun Gylfadottir ######

go("clearinfo")

directory = "/Users/Duna/Box Sync/IHELP/Data"
directoryOut = "/Users/Duna/Documents/Ling/STR/Output"
#create a csv file that will go in the working directory
#if you are running on multiple files and want to create one folder, change "wb" to "a" and it will append to your existing csv file
fileOut = open(directoryOut + "/" + "IHELP.csv", "wb")
#create a writer object
fileWrite = csv.writer(fileOut)
#make header row, if you are adding to an existing csv, comment it out
fileWrite.writerow(["filename","speaker","start","stop","length","word","previous","segment","following","position","cog","maxcog"])


#where your segments are
segTier = 1
#where your words are
wordTier = 2
#how much of the sibilant you want
sibPercentage = 50
# window length (broadband=0.005, narrowband=0.030)
#all sibilants whose 50% length is shorter than this won't be measured
winLen = 0.020
# distance between frames
timeStep = 0.007
# distance between bins
freqStep = 50
# window type
window = "Hamming (raised sine-squared)"
#frequency you filter at?
maxFreq = 15000



def getWavFilePaths(directory):
	'''walks through a directory's subfolders, gets all files with .wav extension, creates list of entire path'''
	wavFiles = []
	for root, dirs, filenames in os.walk(directory, topdown = True):
		#limit to certain folders
		#dirs[:] = [d for d in dirs if "PH8" in d]
		#dirs[:] = [d for d in dirs if "old" not in d]
		filenames = [ f for f in filenames if os.path.splitext(f)[1] in ('.wav', '.WAV') ]
		for filename in filenames:
			wavFiles.append(os.path.join(root, filename))
	return wavFiles

def getName(wavFile):
	'''gets just the filename so we can put it in the csv file later '''
	name = os.path.basename(wavFile)
	baseName = name.rsplit(".",1)[0]
	return baseName

def getSpeaker(name):
	'''gets the speaker's name from the file name by chopping off everything starting at the hypthen'''
	speaker = name.rsplit('-',1)[1]
	return speaker

def getTextGridPath(wavFile):
	'''creates text grid path from given full wav path'''
	path = wavFile.rsplit(".",1)[0]
	textGrid = path + ".TextGrid" 
	return textGrid

def openWavFile(wavFile):
	'''opens the wav file, returns selected Praat object'''
	go("Open long sound file...", wavFile) 
	mySound = selected()
	return mySound

def getWordLabel(wordTier,segTier,segNum):
	'''Gets the word label of a given segment'''
	start = getNum("Get starting point...", segTier, segNum)
	wordInterval = getNum("Get interval at time..." ,wordTier, start+0.000001)
	wordLabel = getString("Get label of interval...", wordTier, wordInterval)
	return wordLabel

def getPoints(segTier,segNum):
	'''takes a tier and a segment number and returns the starting and ending points, length
	and starting and ending points of the middle 50%'''
	sibStart = getNum("Get starting point...", segTier, segNum)
	sibEnd = getNum("Get end point...", segTier, segNum)
	sibLength = sibEnd - sibStart
	measStart = sibStart+(sibLength/4) 
	measEnd = sibEnd-(sibLength/4)
	points = [sibStart,sibEnd,sibLength,measStart,measEnd]
	return points

def getLabels(segNum,segTier,wordTier):
	'''returns a tuple of two lists, one with the segment, one with the word, for the previous segment, the target segment, and the three following segments'''
	segLabels = []
	wordLabels = []
	for num in range(-1,4):
		currSegNum = segNum + num
		try:
			currSegLabel = getString("Get label of interval...", segTier, currSegNum)
			currWordLabel = getWordLabel(wordTier,segTier,currSegNum)
		except:
			currSegLabel = "NA"
			currWordLabel = "NA"
		segLabels.append(currSegLabel)
		wordLabels.append(currWordLabel)
	return segLabels,wordLabels


def IsStr(segLabels,wordLabels):
	'''codes whether the segment is STR'''
	if segLabels[2] == "T":
		if segLabels[3] == "R":
			if wordLabels[3] == wordLabels[1]:
				STR = True
			else: 
				STR = False
		else:
			STR = False
	else:
		STR = False
	return STR


wavFiles = getWavFilePaths(directory)

for wavFile in wavFiles:
	name = getName(wavFile)
	print wavFile
	speaker = getSpeaker(name)
	try:
		mySound = openWavFile(wavFile)
	except:
		print "fail to open wav file"
		continue
	textGrid = getTextGridPath(wavFile)
	try:
		go("Read from file...", textGrid) 
	except:
		print "fail to open text grid"
		continue
	myGrid = selected()
	numSegs = int(getNum("Get number of intervals...",segTier))

	for segNum in xrange(1,numSegs+1):
		select(myGrid)
		segLabel = getString("Get label of interval...", segTier, segNum)
		word = getWordLabel(wordTier,segTier,segNum)
		if segLabel in ("S","SH"):
			if "(" not in word:
				labels = getLabels(segNum,segTier,wordTier)
				wordLabels = labels[1]
				segLabels = labels[0]
				previous = segLabels[0]
				#code for stuff
				if IsStr(segLabels,wordLabels):
					segment = "STR"
					following = segLabels[4]
					if wordLabels[0] != word:
						position = "initial"
					else:
						position = "medial"
				else:
					segment = segLabels[1]
					following = segLabels[2]
					if wordLabels[0] != word:
						position = "initial"
					elif wordLabels[2] != word:
						position = "final"
					else:
						position = "medial"
				###actual measurment part#####		
				points = getPoints(segTier,segNum)
				if (points[2]/2) >= winLen:
					select(mySound)
					try:
						go("Extract part...", points[3], points[4], "rectangular", 1, "no")
					except:
						break
					mySibilant = selected()
					go("Filter (stop Hann band)...", 0, 750, 1) 
					myFilteredSibilant = selected()
					try:
						go("noprogress To Spectrogram...", winLen, maxFreq, timeStep, freqStep, window)
					except: 
						break
					mySpectrogram = selected()
		#get all the frames
					frames = int(getNum("Get number of frames"))
		#get cog for each frame, add total, divide by number of frames to get average
					allCogs = []
					currentTotal = 0
					for frameNum in range(1,frames + 1):
						frameTime = getNum("Get time from frame number...", frameNum)
						go("To Spectrum (slice)...", frameTime)
						mySpectrum=selected()
						cog = getNum("Get centre of gravity...", 2)
						allCogs.append(cog)
						currentTotal=currentTotal+1
						remove(mySpectrum)
						select(mySpectrogram)
					mean = sum(allCogs)/currentTotal
					maxCog = max(allCogs)
					#must be inside if /s/ or /sh/ statement
					remove(mySpectrogram)
					remove(myFilteredSibilant)
					remove(mySibilant)
				else:
					mean = "NA"
					maxCog="NA"
				writeList = [name, speaker, points[0],points[1],points[2],word, previous,segment,following,position,mean,maxCog]
		#write out
				fileWrite.writerow(writeList)
		
		#to check next phoneme in the tier, have to reselect the text grid
		select(myGrid)

#cleanup
	remove(mySound)
	remove(myGrid)


