#!/usr/bin/python
#written by mwalter
# 7 Juli 2016

import sys
sys.path.append('/home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED/')
from functionsrunquery import *
import time

startTime = time.time()


settings={} ## NOTE ALL SETTINGS SAVED AS STRING!!

jobName = sys.argv[1]
workingFolder = "/home/mwalter/project/v0.0.02/projectfiles"
jobFolder = "%s/results/%s"%(workingFolder, jobName)
query1Name = "query1.fa"


settings["RUNNING"]="TRUE"

settings["CONFIGFILENAME"] = "%s/settings.config"%(jobFolder)
settings["JOBFOLDER"] = jobFolder
settings["ADDSSSCRIPTLOCATION"] = "/home/mwalter/hh/scripts/addss.pl"

settings["UNIPROTDATABASELOCATION"]="%s/databases/uniprotHMMs"%(workingFolder)
settings["PDB70DATABASELOCATION"]="%s/databases/pdb70HMMs"%(workingFolder)

settings["UNIPROTDATABASENAME"]= "uniprot20_2016_02"
settings["PDB70DATABASENAME"]= "pdb70"

settings["PICKLELOCATION"]="%s/pickle.pickle"%(workingFolder)


settings["QUERY1NAME"] = query1Name


settings["NUMBERITERATIONSUNIPROT"] = "1"
settings["NUMBERITERATIONSPDB"] = "1"

settings["HHR1OUTPUTLOCATION"] = "%s/%s.hhr" %(jobFolder, query1Name)

settings["S2CLOCATION"]= "%s/databases/s2c" % (workingFolder)


#---
#currently not used:
#settings["INTERFACEPREDICTIONSCRIPT"]= "/home/software/haddock/haddock2.2/tools/contact-chainID"
#usage:  /home/software/haddock/haddock2.2/tools/contact-chainID 4hhb.pdb 6
#----

settings["LOOKUPDATABASELOCATION"]="%s/databases/lookuptable" %(workingFolder)

settings["PDBSEQRESFILE"] = "%s/databases/PDBallFASTA/pdb_seqres.txt"%(workingFolder)


settings["HHRUNIPROTEVALUE"]= "0.0001"
settings["HHRPDBPVALUE"]="40"


#-----------------------------------------
#load settings from query-specific config file:

configFileName = getStringFromSettings("CONFIGFILENAME",settings)
if fileExists(configFileName):
      fillSettings(configFileName, settings) #get all settings from the config textfile
else: throwError("Config file missing!: %s"%(configFileName),settings)

#-----------------------------------------

class NoAlignmentException(Exception):
  pass

print 'Starting process'

#the files that the output will be written to:
humanFile  = "%s/%s.outputHUMAN"%(getStringFromSettings("JOBFOLDER", settings), getStringFromSettings("QUERY1NAME",settings))
computerFile = "%s/%s.outputCOMPUTER"%(getStringFromSettings("JOBFOLDER", settings), getStringFromSettings("QUERY1NAME",settings))
hOutput = open(humanFile, 'w')
cOutput = open(computerFile, 'w')

def hWrite(s):
  hOutput.write(s)

def cWrite(s):
  cOutput.write(s)





#----
#use 1 CPU, make msa with x iterations, add ssinfo, do profile-profile search. results in filename.hhr and log in filename.log

def hhblitsQuery(iFile, settings):
  jobFolder = getStringFromSettings("JOBFOLDER", settings)
  inputFileLocation = "%s/%s"%(jobFolder, iFile)
  UNIPROTEVALUE = getStringFromSettings("HHRUNIPROTEVALUE",settings)
  PDBPVALUE = getStringFromSettings("HHRPDBPVALUE", settings)


  if not fileExists(inputFileLocation):
    throwError("Input query file missing: %s"%(inputFileLocation), settings)
  uniprotLoc = "%s/%s" %(getStringFromSettings("UNIPROTDATABASELOCATION",settings), getStringFromSettings("UNIPROTDATABASENAME",settings))
  pdb70Loc ="%s/%s" %(getStringFromSettings("PDB70DATABASELOCATION",settings), getStringFromSettings("PDB70DATABASENAME",settings))
  numberIterationsUniprot = getStringFromSettings("NUMBERITERATIONSUNIPROT",settings)
  numberIterationsPdb = getStringFromSettings("NUMBERITERATIONSPDB",settings)


  cmd1 = "hhblits -cpu 1 -i %s -d %s -p %s -oa3m %s/%s.a3m -n %s > %s.log 2>&1" %(inputFileLocation, uniprotLoc, UNIPROTEVALUE , jobFolder,iFile, numberIterationsUniprot, inputFileLocation)

  #either show errors and output or write them to the log:
  #cmd2 = "%s %s/%s.a3m >> %s.log 2>&1" %(getStringFromSettings("ADDSSSCRIPTLOCATION",settings),jobFolder, iFile, inputFileLocation)
  cmd2 = "%s %s/%s.a3m" %(getStringFromSettings("ADDSSSCRIPTLOCATION",settings),jobFolder, iFile)

  cmd3 = "hhblits -cpu 1 -i %s/%s.a3m -d %s -E %s -o %s/%s.hhr -n %s >> %s.log 2>&1" %(jobFolder,iFile, pdb70Loc, PDBPVALUE  ,jobFolder,iFile, numberIterationsPdb, inputFileLocation)

  sendToConsole(cmd1)
  sendToConsole(cmd2)
  sendToConsole(cmd3)


hhblitsQuery(getStringFromSettings("QUERY1NAME",settings), settings)



#-----------------------------------------------------------------------
#interpret HHR file:

hits1 = {}

hhrFile1 = getStringFromSettings("HHR1OUTPUTLOCATION",settings)

readHhr(hhrFile1,hits1, settings)



#-----------------------------------------------------------------------

#Great! now we have a dict hits with the hmm names in [0]

#-----------------------------------------------------------------------
#This function can be used to repair corrupted (non-unicode) unpickled dictionaries. Currently seems not needed.

def convert_keys_to_string(dictionary):
  #"""Recursively converts dictionary keys to strings."""
  if not isinstance(dictionary, dict):
    return dictionary
  return dict((str(k), convert_keys_to_string(v))
    for k, v in dictionary.items())





#-----------------------------------------------------------------------
import cPickle

alignmentsDict1 = {}
pickleStartTime = time.time()

def fillAlignmentsDict(hits, alignmentsDict, settings):
  for x in hits:
    nm = hits[x][0]
    lookupTableLoc = getStringFromSettings("LOOKUPDATABASELOCATION",settings)
    #print "\n\n%s\n"%(nm) #bugfixing
    folderLoc = "%s/%s" % (lookupTableLoc, nm)
    if os.path.exists(folderLoc):
      if not nm in alignmentsDict:
        alignmentsDict[nm]=[]
      for fileName in os.listdir(folderLoc):
        if fileName.endswith(".pickle"):
          pklLoc = "%s/%s" % (folderLoc, fileName)
          pkl_file = open(pklLoc, 'rb')
          alignmentsDict[nm].append(cPickle.load(pkl_file))
          pkl_file.close()

      #bugfixing:
      #for i in alignmentsDict[nm]:
        #print i.hsps[0].align_length

fillAlignmentsDict(hits1, alignmentsDict1, settings)

print "So far so good. All pickles loaded in %s seconds" % (str(time.time()-pickleStartTime))



#-----------------------------------------------------------------------

def getComplexNames(ali): #get all complex pdb names from the description in Li's complexfile
  complexNameDict = {}
  cTitle = ali.title.split('>gi|')
  for t in cTitle:
    if t.strip().startswith('gi|'):  #maybe double check why needed...
      t = t[3:]
      #print "t changed to: %s"%(t)
    ts = t.strip()
    if not ts == "":
      cDescr = ts.split()
      cdparts = cDescr[0].split('|')
      complexID = cdparts[2].lower()
      complexChain = cdparts[3].lower()
      if complexID in complexNameDict:
        complexNameDict[complexID].append(complexChain)
      else:
        complexNameDict[complexID]= [complexChain]
  return complexNameDict


#complexNameDict = getComplexNames(alignment.title)-> more intuitive: next line, so send the whole alignment
#complexNameDict = getComplexNames(alignment)

#-------------------------------------

#hits1[x][0] contains the hit name
def getQueryRecord(settings):
  iFile = getStringFromSettings("QUERY1NAME",settings)
  jobFolder = getStringFromSettings("JOBFOLDER", settings)
  inputFileLocation = "%s/%s"%(jobFolder, iFile)
  if not fileExists(inputFileLocation):
    throwError("Input query file missing: %s"%(inputFileLocation), settings)

  handle = open(inputFileLocation, "rU")
  record = SeqIO.read(handle, "fasta")
  handle.close()
  return record


#-------------
PDBALLSEQRES = {}
inputFile = getStringFromSettings("PDBSEQRESFILE",settings)
if not fileExists(inputFile):
  throwError("Input query file missing: %s"%(inputFile), settings)
pdbHandle = open(inputFile, "rU")
for record in SeqIO.parse(pdbHandle, "fasta"):
  PDBALLSEQRES[record.id]=record
pdbHandle.close()


def getSubjectRecord(recordID, settings):
  return PDBALLSEQRES[recordID].seq
#-------------


def calculateIdentity(outputQuery, outputSubject):
  score=0
  total = len(outputQuery)
  if total == 0 :
    print "error! Empty outputQuery %s and outputSubject %s"%(outputQuery, outputSubject)
    return "0" #don't divide by zero
  for x in range(0,total):
    if outputQuery[x]==outputSubject[x]: score+=1
  perc = str((score*10000)/total)
  return "%s.%s"% (perc[:-2],perc[-2:])
  #return str((score*100)/total)

def calculatePositives(outputQuery, outputSubject):
  #right now don't calculate with BLOSUM matrix yet
  return "NaN"



def chopWithoutCounting(seq, howMuch, name):
  seq[name]= seq[name][howMuch:]

def chopWithCounting(seq, counter, howMuch, name):
  counter[name]+=howMuch
  chopWithoutCounting(seq,howMuch, name)

def chop1(seq,counter,name):
  chopWithCounting(seq, counter,1,name)


def returnEntryInBlastFormat(queryRecord, hitsEntry, alignment, accessionNumber, settings):
    hsp = alignment.hsps[0]


    #make 3 objects: one with start and stop nrs, one with sequences, and one with counters
    startStop = {}
    seq = {}
    counter = {}


    #note: can remove second typecast later
    startStop["queryHHRstart"] =   int(hitsEntry[8].split('-')[0].strip())-1
    startStop["queryHHRstop"] =    int(hitsEntry[8].split('-')[1].strip())-1
    startStop["subjectHHRStart"] = int(hitsEntry[9].split('-')[0].strip())-1
    startStop["subjectHHRStop"] =  int(hitsEntry[9].split('-')[1].strip())-1

    #blast enumeration is '1-based'
    startStop["queryBLASTStart"] =   int(hsp.query_start)-1
    startStop["queryBLASTStop"] =    int(hsp.query_end)-1
    startStop["subjectBLASTStart"] = int(hsp.sbjct_start)-1
    startStop["subjectBLASTStop"] =  int(hsp.sbjct_end)-1

    seq["queryOriginal"] =        queryRecord.seq
    seq["queryHHR"] =             hitsEntry[12]
    seq["alignmentHHR"] =         hitsEntry[14]
    seq["subjectHHR"] =           hitsEntry[16]
    seq["subjectHHROriginal"] =   getSubjectRecord(hitsEntry[0], settings)
    seq["queryBLAST"] =           hsp.query
    seq["alignmentBLAST"] =       hsp.match
    seq["subjectBLAST"] =         hsp.sbjct
    subjectBLASTName = "%s_%s"%(alignment.title.split()[0].split('|')[3].lower(), alignment.title.split()[0].split('|')[4])
    seq["subjectBLASToriginal"] = getSubjectRecord(subjectBLASTName, settings)


    counter["queryOriginal"] =   startStop["queryHHRstart"]
    counter["queryHHR"] =        0
    counter["alignmentHHR"] =    0
    counter["subjectHHR"] =      0
    counter["subjectHHROriginal"] = startStop["subjectHHRStart"]
    counter["queryBLAST"] =      0
    counter["alignmentBLAST"] =  0
    counter["subjectBLAST"] =    0
    counter["subjectBLASToriginal"] = startStop["subjectBLASTStart"]

    one = seq["queryOriginal"][0:startStop["queryHHRstart"]]
    two = ""
    three = seq["subjectBLASToriginal"][0:startStop["subjectBLASTStart"]]


    #this is wrong??: depends which case??
    #if (int(startStop["queryHHRstop"]) - int(startStop["queryHHRstart"]) <  int(startStop["queryBLASTStart"]) - int(startStop["subjectHHRStart"])):
      #print "Error! No overlap between alignments!!"


    if int(startStop["queryBLASTStart"]) > int(startStop["subjectHHRStart"]):
      two = seq["subjectHHROriginal"][0:startStop["subjectHHRStart"]]


      difference = int(startStop["queryBLASTStart"]) - int(startStop["subjectHHRStart"])
      if len(seq["queryHHR"].replace("-", ""))<difference:
        print "error! non-overlapping alignment (1)!"
        raise NoAlignmentException
      i=0
      j=0
      k=0
      while i<difference:
        one+= seq["queryHHR"][j]
        two+= seq["subjectHHR"][j]
        if not seq["subjectHHR"][j] == '-':
          i+=1
        if not seq["queryHHR"][j]== '-':
          k+=1
        j+=1

      gapsQ = j-k
      gapsS = j-i

      #for testing:
      print "\ni is now %s and j is now %s" %(i,j)
      print "difference is %s\n" %(difference)
      if not i==j: print "unequal\n"



      counter["queryOriginal"]+= difference + gapsS -gapsQ
      counter["queryHHR"] += difference + gapsS
      counter["alignmentHHR"] += difference + gapsS
      counter["subjectHHR"] += difference + gapsS
      counter["subjectHHROriginal"] += difference

    if int(startStop["subjectHHRStart"]) > int(startStop["queryBLASTStart"]):
      two = seq["subjectHHROriginal"][0:startStop["queryBLASTStart"]]  #doublecheck SUBJECTblaststart?
      difference = int(startStop["subjectHHRStart"]) - int(startStop["queryBLASTStart"])

      if len(seq["queryBLAST"].replace("-", ""))<difference:
        print "error! non-overlapping alignment!(2)" #should be subjectBLAST
        raise NoAlignmentException


      i=0
      j=0
      k=0
      while i<difference:
        two+= seq["queryBLAST"][j]
        three+= seq["subjectBLAST"][j]
        if not seq["queryBLAST"][j] == '-':
          i+=1
        if not seq["subjectBLAST"] == '-':
          k+=1
        j+=1

      gapsQ = j-i
      gapsS = j-k



      #for testing:
      print "\ni is now %s and j is now %s" %(i,j)
      print "difference is %s\n" %(difference)
      if not i==j: print "unequal\n"


      #? why again not?: counter["queryOriginal"]+= difference

      counter["queryBLAST"] += difference +gapsQ
      counter["alignmentBLAST"] += difference + gapsQ
      counter["subjectBLAST"] += difference + gapsQ
      counter["subjectBLASToriginal"] += difference -gapsS +gapsQ

      #for testing:
      print "after adding, query counter is: %s" %(str(counter["queryBLAST"]))
      print "after adding, subject counter is: %s\n" % (str(counter["subjectBLAST"]))



    for i in seq:
      chopWithoutCounting(seq,counter[i],i)

    outputQuery = ""
    outputAlignment = ""
    outputSubject = ""

    #for testing:
    #print "\n---------------------------------------------------------------------------------\n"
    #print "\nNow queryHHR and subjectHHR and queryBLAST and subjectBLAST are:"
    #print seq["queryHHR"]
    #print seq["subjectHHR"]
    #print seq["queryBLAST"]
    #print seq["subjectBLAST"]
    #print "\n\n"
    #print "full HHR alignment: \n"
    #print hitsEntry[12]
    #print hitsEntry[16]
    #print "\n"
    #print "full blast alignment: \n"
    #print hsp.query
    #print hsp.sbjct
    #print "\n"


    #difference = int(startStop["queryBLASTStart"]) - int(startStop["subjectHHRStart"])
    #print "Difference: %s\n" % (difference)
    #print "queryHHRstart: %s"% (startStop["queryHHRstart"]) #fix type capital s later!
    #print "subjectHHRstart: %s"% (startStop["subjectHHRStart"])
    #print "queryBLASTstart: %s" % (startStop["queryBLASTStart"])
    #print "subjectBLASTstart: %s" % (startStop["subjectBLASTStart"])

    #print "\nqueryHHRstop: %s"% (startStop["queryHHRstop"])
    #print "subjectHHRstop: %s"% (startStop["subjectHHRStop"])
    #print "queryBLASTstop: %s" % (startStop["queryBLASTStop"])
    #print "subjectBLASTstop: %s" % (startStop["subjectBLASTStop"])

    #print "\nOriginal sequences (query, intermediate, subject): \n%s\n%s\n%s\n\n" %(queryRecord.seq, getSubjectRecord(hitsEntry[0], settings),getSubjectRecord(subjectBLASTName, settings))
    #print "One: \t%s" %(one)
    #print "Two: \t%s" %(two)
    #print "Three: \t%s" %(three)
    #print "---\n"
    #  end for testing

    while len(seq["queryHHR"])>0 and len(seq["subjectBLAST"])>0 :

      if seq["subjectHHR"][0]=='-' and not seq["queryBLAST"][0]=='-': #if both are '-', just align them.
        outputQuery += seq["queryHHR"][0]
        outputAlignment += '-'
        outputSubject += '-'
        chop1(seq, counter, "queryOriginal")
        chop1(seq, counter, "queryHHR")
        chop1(seq, counter, "alignmentHHR")
        chop1(seq, counter, "subjectHHR")
        chop1(seq, counter, "subjectHHROriginal")

      elif seq["queryBLAST"][0]=='-' and not seq["subjectHHR"][0]=='-':
        outputQuery += '-'
        outputAlignment += '-'
        outputSubject += seq["subjectBLAST"][0]
        chop1(seq, counter, "queryBLAST" )
        chop1(seq, counter, "alignmentBLAST" )
        chop1(seq, counter, "subjectBLAST" )
        chop1(seq, counter, "subjectBLASToriginal" )

      else:
        outputQuery += seq["queryHHR"][0]
        outputAlignment += seq["subjectHHR"][0]
        if not seq["subjectHHR"][0]==seq["queryBLAST"][0]: print "Error! Something went wrong: %s, %s" %(seq["subjectHHR"][0],seq["queryBLAST"][0]) #for testing
        outputSubject += seq["subjectBLAST"][0]
        for a in seq: chop1(seq,counter,a)


    l = max(len(one), len(two), len(three))
    #oops.. we work in python 2.6, so following is 2.7+
    #one = '{:>{}s}'.format(one, l)
    #two = '{:>{}s}'.format(two, l)
    #three = '{:>{}s}'.format(three, l)
    #redone in old format:
    one = '%*s' % ((l), one)
    two = '%*s' % ((l), two)
    three= '%*s' % ((l), three)




    #for testing:
    one+="|||"
    two+="|||"
    three+="|||"


    one += outputQuery
    two += outputAlignment
    three += outputSubject

    #for testing:
    one+="|||"
    two+="|||"
    three+="|||"

    #--

    if int(startStop["subjectHHRStop"]) > int(startStop["queryBLASTStop"]):
      difference = int(startStop["subjectHHRStop"]) - int(startStop["queryBLASTStop"])
      i=0
      j=0
      while i<difference:
        one+= seq["queryHHR"][j]
        two+= seq["subjectHHR"][j]
        if not seq["subjectHHR"][j] == '-':
          i+=1
          chop1(seq, counter, "subjectHHROriginal" )
        if not seq["queryHHR"][j] == '-':
          chop1(seq, counter, "queryOriginal" )
        j+=1


    if int(startStop["queryBLASTStop"]) > int(startStop["subjectHHRStop"]):
      difference = int(startStop["queryBLASTStop"]) - int(startStop["subjectHHRStop"])
      i=0
      j=0
      while i<difference:
        two+= seq["queryBLAST"][j]
        three+= seq["subjectBLAST"][j]
        if not seq["queryBLAST"][j] == '-':
          i+=1
          chop1(seq, counter, "subjectHHROriginal")
        if not seq["subjectBLAST"][j] == '-':
          chop1(seq, counter, "subjectBLASToriginal" )
        j+=1

    one+=seq["queryOriginal"]
    two+=seq["subjectHHROriginal"]
    three+=seq["subjectBLASToriginal"]

    l=max(len(one), len(two), len(three))
    one = '%*s' % ((-l), one)
    two = '%*s' % ((-l), two)
    three = '%*s' % ((-l), three)


    print "\nAlignment complete!!:"
    hWrite("\nAlignment complete!!:\n")
    print "Query:  %s, HHRhit: %s, Subject: %s\n"%(queryRecord.id, hitsEntry[0], alignment.title.split()[0])
    hWrite("Query:  %s, HHRhit: %s, Subject: %s\n"%(queryRecord.id, hitsEntry[0], alignment.title.split()[0]))

    print one
    hWrite(one+"\n")
    print two
    hWrite(two+"\n")
    print three
    hWrite (three+"\n\n")
    print "\n"





    # Fields: subject ids, query length, subject length, query seq, subject seq, q.
    # start, q. end, s. start, s. end, bit score, evalue, alignment length, %
    # identity, % positives, query acc.

    #outputLine = []

    #return outputLine
    #print outputQuery
    #print outputSubject

    #subjectID = "gi|%s|pdb|%s|%s" % (accessionNumber, hitsEntry[0][0:4], hitsEntry[0][-1:]) #NOW SUBJECT ID CONTAINS HHR HIT
    subjectID = alignment.title.split()[0]
    queryLength = len(queryRecord.seq)
    subjectLength = len(hsp.sbjct) #WRONG!! TAKES ONLY THE ALIGNED PART AND INCLUDES GAPS!
    querySeq =  outputQuery
    subjectSeq = outputSubject
    qStart =hitsEntry[8].split('-')[0].strip()
    qEnd = hitsEntry[8].split('-')[1].strip()
    #hhrFoundStart = hitsEntry[9].split('-')[0].strip()
    #hhrFoundEnd= hitsEntry[9].split('-')[1].strip()
    #hhrSentStart = hsp.query_start
    #hhrSentEnd = hsp.query_end
    sStart = hsp.sbjct_start
    sEnd = hsp.sbjct_end
    bitScore = hitsEntry[5]
    #insert one or more HHpred result statistics here
    #choose from:
    #prob(2),E-value(3), P-value(4), SCORE (5)
    eValue = hitsEntry[3]
    #alternative: hsp.expect = blast calculated e-value -> use HHR instead!
    alignmentLength = len(outputQuery) #calculate this from the (possibly short) resulting alignment
    percentIdentity = calculateIdentity(outputQuery, outputSubject)
    percentPositives = calculatePositives(outputQuery, outputSubject)
    queryAcc = queryRecord.id
    #question -> should double gaps (in both query and result) be deleted again? ->slightly problematic -> answer: no.


    return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (subjectID, queryLength , subjectLength, querySeq, subjectSeq, qStart, qEnd, sStart, sEnd, bitScore, eValue, alignmentLength, percentIdentity, percentPositives,queryAcc)




def printOutputInBlastFormat(hits, alignDict, settings):
  outputList = []
  rc = getQueryRecord(settings)
  #add meta header here if wanted
  accessionNumber = 111111 #arbitrary but unique accession ID assigned -> can be removed, not used any more!
  print "test6 %s" %(str(len(hits)))
  for x in hits:
    qb = hits[x][0]
    if qb in alignDict:
      for y in range(0, len(alignDict[qb])):
        align = alignDict[qb][y]
        try:
          outputList.append(returnEntryInBlastFormat(rc, hits[x], align, accessionNumber, settings))
        except NoAlignmentException:
          print "No overlap between alignments Error caught!"
        accessionNumber+=1


  #print each entry in outputList to file:

  #first print Header
  print "\n\nStarting output: \n\n"
  version = "HHPRED&PBLAST"
  iteration = getStringFromSettings("NUMBERITERATIONSUNIPROT", settings) #doublecheck later: should not be NUMBERITERATIONSPDB????
  queryName = "USERQUERY"
  dbPDB70 = getStringFromSettings("PDB70DATABASELOCATION" ,settings)
  nrHits = len(outputList)

  print "# %s"%(version)
  cWrite("# %s\n"%(version))
  print "# Iteration: %s" % (iteration)
  cWrite("# Iteration: %s\n" % (iteration))
  print "# Query: %s" %(rc.id)
  cWrite("# Query: %s\n" %(rc.id))
  print "# Database: %s" %(dbPDB70)
  cWrite("# Database: %s\n" %(dbPDB70))
  print "# Fields: subject ids, query length, subject length, query seq, subject seq, q. start, q. end, s. start, s. end, bit score, evalue, alignment length, % identity, % positives, query acc."
  cWrite("# Fields: subject ids, query length, subject length, query seq, subject seq, q. start, q. end, s. start, s. end, bit score, evalue, alignment length, % identity, % positives, query acc.\n")
  print "# %s hits found" % (nrHits)
  cWrite("# %s hits found\n\n" % (nrHits))
  for i in outputList:
    print i
    cWrite(i+"\n")



printOutputInBlastFormat(hits1, alignmentsDict1, settings)



hOutput.close()
cOutput.close()



print "Total time elapsed: %s"%(str(time.time()-startTime))
print "program done"




