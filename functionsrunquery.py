#written by mwalter
# 8 June 2016


from functionsbasic import *
import sys
sys.path.append('/home/mwalter/modules/pythonmodules/lib64/python')
from Bio import SeqIO
from Bio.Blast import NCBIXML


def getLegend():
  legend = []
  legend.append("PDB ID")               # 0 represents PDB ID
  legend.append("Descr.")          # 1 represents the description
  legend.append("Prob")                 # 2
  legend.append("E-value")              # 3
  legend.append("P-Value")              # 4
  legend.append("Score")                # 5
  legend.append("SS")                   # 6
  legend.append("Cols")                 # 7
  legend.append("Query")            # 8
  legend.append("Templ.")             # 9
  legend.append("Total")                # 10 up to here is from header
  legend.append("11")             # 11 these will still change later /now start of align ->SSPRED
  legend.append("12")             # 12 QUERY
  legend.append("13")             # 13 CONSENSUS
  legend.append("Align")             # 14
  legend.append("15")             # 15
  legend.append("16")             # 16
  legend.append("17")             # 17
  legend.append("18")             # 18
  legend.append("SS")             # 19 CONFIDENCE??
  legend.append("P")             # 20
  #legend.append("E")             # 21
  #legend.append("S")             # 22
  #legend.append("C")             # 23
  #legend.append("24")             # 24
  #legend.append("25")             # 25
  #legend.append("26")             # 26
  #legend.append("27")             # 27
  #legend.append("")             #
  #etc, finish this later.
  return legend

def lineIsNewEntry(str1): #improve this later
  #if ">" in str1:
  if str1[:2]=="No":
    return 1
  return 0


def parseHeader(line, varDict):
  #ignore header for now.
  #if "No Hit" in line:
  if line.startswith(" No Hit"):
    varDict["currentPortion"]+=1


def parseFirstResults(line, hits,varDict):
  #if "No 1" in line(:3):
  if line.strip()=="":
    varDict["currentPortion"]+=1
  #print "moving to part 2\n\n"
  #print "\n\n"
  else:
    #add stuff to list of hits
    array1 = line.split(None,2)
    text1 = array1[2][:23]
    finalArray = [array1[0],array1[1], text1]
    for i in array1[2][24:].split():
      #workaround to fix bug where numbers are so big that there is no whitespace:
      if '(' in i.strip()[1:]:
        composite = i.split('(')
        finalArray.append(composite[0])
        finalArray.append("(%s"%(composite[1]))
      else:
        finalArray.append(i)

    #print "finalArray length: %d"%(len(finalArray))

    hits[int(finalArray[0])] = finalArray[1:]
    for i in range(0,9):
      hits[int(finalArray[0])].append("") #sequence placeholders
    #print "Hits length: %s\t Entry: %s" %(len(hits[int(finalArray[0])]), finalArray[0])
    #print hits[int(finalArray[0])]

def checkEmptyAlignment(currentEntry, hits):
  if hits[currentEntry][12] == "":
    print "Error! Empty HHR entry: %s"%(hits[currentEntry][0])#for testing 26-10-2016 #NOTE: last entry not checked! -> now checked in another portion.



def parseSecondResults(line, hits, varDict):
  currentEntry = varDict["currentEntry"]
  if lineIsNewEntry(line):
    #if currentEntry>0 and  hits[currentEntry][12] == "": print "Error! Empty HHR entry: %s"%(hits[currentEntry][0])#for testing 26-10-2016 #NOTE: last entry not checked!
    if currentEntry>0: checkEmptyAlignment(currentEntry, hits)
    varDict["currentEntry"]+=1
    if(str(varDict["currentEntry"])!=line.split()[1]):
      throwError("Error reading HHM file entry nr: %d" % (varDict["currentEntry"]), settings)  #print line
  else:
    if line[0]=='>':
      #print "%d\t%s" % (currentEntry, line)
      varDict["lineAfterTitle"] = 1
    elif varDict["lineAfterTitle"]:
      varDict["lineAfterTitle"] = 0
      for i in line.split():
        value1 = i.split('=')
        hits[currentEntry].append(value1[1])
    elif not line.strip()=="" :
      if line[0]=='Q':
        pieces = line.split()
        if pieces[1]=="ss_pred":
          hits[currentEntry][11]+=pieces[2]
          varDict["currentLineLength"] = len(pieces[2])
        elif pieces[1]=="Consensus":
          hits[currentEntry][13]+=pieces[3]
        else :
          hits[currentEntry][12]+=pieces[3]
      elif line[0]=='T':
        pieces = line.split()
        if pieces[1]=="ss_pred":
          hits[currentEntry][18]+=pieces[2]
        elif pieces[1]=="Consensus":
          hits[currentEntry][15]+=pieces[3]
        elif pieces[1]=="ss_dssp":
          hits[currentEntry][17]+=pieces[2]
        else :
          hits[currentEntry][16]+=pieces[3]
      #elif line.split(None, 1)[0] == "Confidence":
        #hits[currentEntry][19] += line[10:].strip('\n\t')
      #else :
        #hits[currentEntry][14] += line.strip('\n\t')

      elif line.split(None, 1)[0] == "Confidence":
        hits[currentEntry][19] += line[22:varDict["currentLineLength"]+22].strip('\n\t') #whitespace apparently made of spaces..
      else :
        hits[currentEntry][14] += line[22:varDict["currentLineLength"]+22].strip('\n\t')



  #testing:
  #if (currentEntry > 10):
  #  leg=getLegend()
  #  for i in range(5,10):
  #    for j in range(0,len(hits[i])):
  #      print "%s\t%s"%(leg[j],hits[i][j])
  #    print "\n\n\n"
  #  sys.exit()









def readHhr(fn, hits, settings):
  if not fileExists(fn): throwError("HHR file %s does not exist!"%(fn),settings)

  varDict = {}
  varDict["currentPortion"]= 0
  varDict["currentEntry"] = 0

  varDict["lineAfterTitle"]=0    #boolean
  varDict["currentLineLength"]=0 #length of the last line.

  with open(fn) as hhr:
    for line in hhr:
      if varDict["currentPortion"] == 0: #currently in header
        parseHeader(line, varDict)
      elif varDict["currentPortion"] == 1: #currently in 1st results
        parseFirstResults(line,hits,varDict)
      elif varDict["currentPortion"] == 2: #currently in 1st results
        parseSecondResults(line,hits,varDict)
      elif varDict["currentPortion"] >= 3:
        throwError("Cannot read HHR file %s"%(fn), settings)
  if len(hits)>0: checkEmptyAlignment(varDict["currentEntry"], hits) #also check the last entry! the rest is checked during 'parseSecondResults'
  #print "should close hhr here still or not? i dont think needed"





def read_s2c_mapping(complex_name, chain_name, settings):
  fn = "%s/%s.sc"%(getStringFromSettings("S2CLOCATION", settings) , complex_name)
  mappingList = []
  mappingCounter= 0
  with open(fn) as s2c:
    for line in s2c:
      pts = line.split()
      if pts[0]=="SEQCRD" and pts[1]==chain_name:
        mappingCounter+=1
        if not str(mappingCounter)==pts[5]:
          throwError("Something wrong with S2C Mapping!!!", settings)
        else:
          mappingList.append(pts[6]) ### Be REALLY careful! the first mapping, position 1, is put in place 0 of the list!!!!!
  return mappingList



def readInterfaceFile(fn, chain1, chain2, settings):
  answer = {}
  answer[chain1]={}
  answer[chain2]={}
  with open(fn) as ifacefile:
    for line in ifacefile:
      pts = line.split()
      pts[1]= pts[1].lower()
      pts[4]= pts[4].lower()
      if (pts[1] == chain1 and pts[4]==chain2) or (pts[1]==chain2 and pts[4]==chain1):
        #format: answer[chain][pdsresnr] = [distance1, distance2, etc..]
        if not pts[0] in answer[pts[1]]:
          answer[pts[1]][pts[0]] = [pts[6]]
        else: answer[pts[1]][pts[0]].append(pts[6])
        if not pts[3] in answer[pts[4]]:
          answer[pts[4]][pts[3]] = [pts[6]]
        else: answer[pts[4]][pts[3]].append(pts[6])
  return answer


