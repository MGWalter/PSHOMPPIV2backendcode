#!/usr/bin/python
#written by mwalter
# 7 June 2016


from functionsbasic import *
import fileinput
import sys
sys.path.append('/home/mwalter/modules/pythonmodules/lib64/python')
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pickle




def fixHhsuiteLinks(settings):
  fileToPut = getStringFromSettings("HHLINKSFILE",settings)
  dbLinks={}
  dbLinks["our $execdir"]= getStringFromSettings("PSIPREDBIN",settings)

  dbLinks["our $datadir"]= getStringFromSettings("PSIPREDDATA",settings)
  dbLinks["our $ncbidir"]= getStringFromSettings("NCBIBIN",settings)
# dbLinks["our $pdbdir"]= getStringFromSettings("PDBSTRUCTURES",settings)
  dbLinks["our $dsspdir"]= getStringFromSettings("DSSPDATABASELOCATION",settings)
  dbLinks["our $dssp "]= getStringFromSettings("DSSPBIN",settings) #the space is important, or dsspdir gets overwritten


  for line in fileinput.input(fileToPut, inplace=True):
    for i in dbLinks:
      if line.startswith(i):
        p=line.split(';')
        remainder = p[1]
        line = "%s\t= \"%s\";\t%s"%(i,dbLinks[i],remainder)
    print "%s" % (line),



#-------NOTES---------------------
#databases to grab
#-pdballfasta (done)
#-complexblastdb (done)
#-dssp structures
#-uniprot
#-pdb70 hhms
#-pdb structures ?
#---------------------------------






#---------------------------------------------------------------
#make a blast database out of the complexfastafile
def makeBlastDatabaseFromComplexList(settings):
  #make a blast database out of the complexfastafile
  cfile = getStringFromSettings("COMPLEXFASTAFILE",settings)
  if fileExists(cfile):
    sendToConsole("%s/makeblastdb -in %s -parse_seqids -dbtype prot > %s.info 2>&1"%(getStringFromSettings("BLASTHOME",settings),cfile, cfile))
  else: throwError("COMPLEXFASTAFILE file missing: %s"%(getStringFromSettings("COMPLEXFASTAFILE", settings)), settings)
  if not dirExists(getStringFromSettings("COMPLEXDATABASELOCATION",settings)):
    if getBoolFromSettings("AUTOCREATEMISSINGDATABASELOCATIONS",settings) and dirExists(getStringFromSettings("DEFAULTDATABASELOCATION",settings)):
      sendToConsole("mkdir %s/complexDB" % (getStringFromSettings("DEFAULTDATABASELOCATION",settings)))
      settings["COMPLEXDATABASELOCATION"]= "%s/complexDB" %(getStringFromSettings("DEFAULTDATABASELOCATION",settings))

    else:
      throwError("COMPLEXDATABASELOCATION does not exist: %s"%(getStringFromSettings("COMPLEXDATABASELOCATION", settings)),settings)

  sendToConsole("mv %s.* %s"%(getStringFromSettings("COMPLEXFASTAFILE", settings),getStringFromSettings("COMPLEXDATABASELOCATION",settings)))

#---------------------------------------------------------------
#Download allPDB in fasta format.
def downloadAllPdb(settings):
  pdballfolder = getStringFromSettings("PDBALLFASTASERVERFOLDER", settings)
  pdballfile = getStringFromSettings("PDBALLFASTASERVERFILE", settings)
  if not dirExists(getStringFromSettings("PDBALLFASTADATABASELOCATION", settings)):
    if getBoolFromSettings("AUTOCREATEMISSINGDATABASELOCATIONS",settings) and dirExists(getStringFromSettings("DEFAULTDATABASELOCATION",settings)):
      sendToConsole("mkdir %s/PDBallFASTA" % (getStringFromSettings("DEFAULTDATABASELOCATION",settings)))
      settings["PDBALLFASTADATABASELOCATION"]= "%s/PDBallFASTA"%(getStringFromSettings("DEFAULTDATABASELOCATION",settings))
    else:
      throwError("PDBALLFASTADATABASELOCATION does not exist: %s"%(getStringFromSettings("PDBALLFASTADATABASELOCATION",settings)),settings)

  dbloc = getStringFromSettings("PDBALLFASTADATABASELOCATION",settings)
  if dirExists(dbloc):
    sendToConsole("cd %s;wget %s%s.gz > download.log 2>&1"%(dbloc,pdballfolder, pdballfile))
  if fileExists("%s/%s.gz"%(dbloc,pdballfile)):
    sendToConsole("gunzip %s/%s.gz >> %s/download.log 2>&1"%(dbloc,pdballfile,dbloc))
  else:
    throwError("PDBALL fasta txt file could not be downloaded! Aborting", settings)
#---------------------------------------------------------------
#Download pdb70 HMM files.
def downloadPdb70Hmm(settings):
  pdb70folder = getStringFromSettings("PDB70SERVERFOLDER", settings)
  pdb70file = getStringFromSettings("PDB70SERVERFILE", settings)
  if not dirExists(getStringFromSettings("PDB70DATABASELOCATION", settings)):
    if getBoolFromSettings("AUTOCREATEMISSINGDATABASELOCATIONS",settings) and dirExists(getStringFromSettings("DEFAULTDATABASELOCATION",settings)):
      sendToConsole("mkdir %s/pdb70HMMs" % (getStringFromSettings("DEFAULTDATABASELOCATION",settings)))
      settings["PDB70DATABASELOCATION"]= "%s/pdb70HMMs"%(getStringFromSettings("DEFAULTDATABASELOCATION",settings))
    else:
      throwError("PDB70DATABASELOCATION does not exist: %s"%(getStringFromSettings("PDB70DATABASELOCATION",settings)),settings)

  dbloc = getStringFromSettings("PDB70DATABASELOCATION",settings)
  if dirExists(dbloc):
    sendToConsole("cd %s;wget %s%s > download.log 2>&1"%(dbloc,pdb70folder, pdb70file))
  if fileExists("%s/%s"%(dbloc,pdb70file)):
    sendToConsole("tar zxvf %s/%s -C %s >> %s/download.log 2>&1 "%(dbloc,pdb70file,dbloc,dbloc))
    #clean up the tgz file
    sendToConsole("rm %s/%s"%(dbloc, pdb70file))
    #now move the files one folder up.
    #sendToConsole("mv %s/%s/* %s"%(dbloc, pdb70file[:-4], dbloc)) #apparently pdb70 doesnt need moving, its not in a deeper folder
  else:
    throwError("PDB70 HMM files could not be downloaded! Aborting", settings)

#----------------------------------------------------------------
#Download Uniprot HMM files.
def downloadUniprotHmm(settings):
  uniprotfolder = getStringFromSettings("UNIPROTSERVERFOLDER", settings)
  uniprotfile = getStringFromSettings("UNIPROTSERVERFILE", settings)
  if not dirExists(getStringFromSettings("UNIPROTDATABASELOCATION", settings)):
    if getBoolFromSettings("AUTOCREATEMISSINGDATABASELOCATIONS",settings) and dirExists(getStringFromSettings("DEFAULTDATABASELOCATION",settings)):
      sendToConsole("mkdir %s/uniprotHMMs" % (getStringFromSettings("DEFAULTDATABASELOCATION",settings)))
      settings["UNIPROTDATABASELOCATION"]= "%s/uniprotHMMs"%(getStringFromSettings("DEFAULTDATABASELOCATION",settings))
    else:
      throwError("UNIPROTDATABASELOCATION does not exist: %s"%(getStringFromSettings("UNIPROTDATABASELOCATION",settings)),settings)

  dbloc = getStringFromSettings("UNIPROTDATABASELOCATION",settings)
  if dirExists(dbloc):
    sendToConsole("cd %s;wget %s%s > download.log 2>&1"%(dbloc,uniprotfolder, uniprotfile))
  if fileExists("%s/%s"%(dbloc,uniprotfile)):
    sendToConsole("tar zxvf %s/%s -C %s >> %s/download.log 2>&1"%(dbloc,uniprotfile,dbloc,dbloc))
    #clean up the tgz file
    sendToConsole("rm %s/%s"%(dbloc, uniprotfile))
    #now move the files one folder up.
    sendToConsole("mv %s/%s/* %s"%(dbloc, uniprotfile[:-4], dbloc))
    #!!maybe implement later to remove the empty folder too
  else:
    throwError("UNIPROT HMM files could not be downloaded! Aborting", settings)


#----------------------------------------------------------------
#Download DSSP files
def downloadDsspFiles(settings):
  dsspfolder = getStringFromSettings("DSSPSERVERFOLDER", settings)
#  dsspfile = getStringFromSettings("DSSPSERVERFILE", settings)
  if not dirExists(getStringFromSettings("DSSPDATABASELOCATION", settings)):
    if getBoolFromSettings("AUTOCREATEMISSINGDATABASELOCATIONS",settings) and dirExists(getStringFromSettings("DEFAULTDATABASELOCATION",settings)):
      sendToConsole("mkdir %s/dsspDB" % (getStringFromSettings("DEFAULTDATABASELOCATION",settings)))
      settings["DSSPDATABASELOCATION"]= "%s/dsspDB"%(getStringFromSettings("DEFAULTDATABASELOCATION",settings))
    else:
      throwError("DSSPDATABASELOCATION does not exist: %s"%(getStringFromSettings("DSSPDATABASELOCATION",settings)),settings)

  dbloc = getStringFromSettings("DSSPDATABASELOCATION",settings)
  if dirExists(dbloc):
    sendToConsole("rsync -a rsync://%s %s > %s/download.log 2>&1" %(dsspfolder, dbloc, dbloc))
  #TRY TO BUILD AN ERROR CHECK HERE!
#    sendToConsole("cd %s;wget %s%s.gz>download.log"%(dbloc,dsspfolder, dsspfile))
#  if fileExists("%s/%s.gz"%(dbloc,dsspfile)):
#    sendToConsole("gunzip %s/%s.gz>>download.log"%(dbloc,dsspfile))
#  else:
#    throwError("DSSP files could not be downloaded! Aborting", settings)


#rsync -a rsync://rsync.cmbi.ru.nl/dssp /data/dssp

#----------------------------------------------------------------
#DOWNLOAD PDB STRUCTURE FILES? -> not yet










#----------------------------------------------------------------

#backup any of the existing databases:
def backupAllDatabases(settings):
  fulldb=getStringFromSettings("DEFAULTDATABASELOCATION", settings)
  pdbAllFastaLoc = getStringFromSettings("PDBALLFASTADATABASELOCATION",settings)
  complexLoc = getStringFromSettings("COMPLEXDATABASELOCATION",settings)
  uniprotLoc = getStringFromSettings("UNIPROTDATABASELOCATION",settings)
  pdb70Loc = getStringFromSettings("PDB70DATABASELOCATION",settings)
  dsspLoc = getStringFromSettings("DSSPDATABASELOCATION",settings)
  linkerDB = getStringFromSettings("LINKERDATABASELOCATION",settings)
  #PDB STRUCTURES?
  allSecondary = [pdbAllFastaLoc,complexLoc,uniprotLoc,pdb70Loc,dsspLoc,linkerDB] #also add PDB structures?
  #maybe make stepwise program, such as allSecondary.append(db1), etc, so individual dbs can easily be commented out of backup process.
  #allSecondary= []
  #allSecondary.append(pdbAllFastaLoc)
  #allSecondary.append(complexLoc)
  #allSecondary.append(uniprotLoc)
  #allSecondary.append(pdb70Loc)
  #allSecondary.append(dsspLoc)
  #allSecondary.append(linkerDB)
  #

  #fill in others here.
  if not dirExists(getStringFromSettings("BACKUPLOCATION",settings)):
      sendToConsole("mkdir %s/DBbackups"%(getStringFromSettings("WORKINGFOLDER",settings)))
      settings["BACKUPLOCATION"]="%s/DBbackups"%(getStringFromSettings("WORKINGFOLDER",settings))

  bl = getStringFromSettings("BACKUPLOCATION",settings)

  if dirExists(fulldb): backup(fulldb,bl)
  #also make backups of the linked databases, if they are not in the normal DB folder:
  for s in allSecondary:
    #if dirExists(pdbAllFastaLoc) and not pdbAllFastaLoc.startswith(fulldb): backup(pdbAllFastaLoc,bl)
    if dirExists(s) and not s.startswith(fulldb): backup(s,bl)


def cleanLookupTable(settings):
  dbLocation= getStringFromSettings("LOOKUPDATABASELOCATION", settings)
  b1 = getStringFromSettings("BACKUPLOCATION", settings)
  if dirExists(dbLocation): backup(dbLocation,b1)
  else: sendToConsole("mkdir %s"%(dbLocation))





def makeFileToBlast(settings):
  allFastaDict={}
  allFastaFileName = "%s/%s"%(getStringFromSettings("PDBALLFASTADATABASELOCATION",settings),getStringFromSettings("PDBALLFASTASERVERFILE",settings))
  seqToBlast = {}
  pdb70IndexFileName = "%s/%s"%(getStringFromSettings("PDB70DATABASELOCATION",settings),"pdb70_pdb.index")
  for seqRecord in SeqIO.parse(allFastaFileName, "fasta"):
    allFastaDict[seqRecord.id]= seqRecord.seq


  with open(pdb70IndexFileName) as indexFile:
    #counter =0
    for line in indexFile:
      #counter += 1
      idFromIndex = line[:6]
      if idFromIndex in allFastaDict:
        seqToBlast[idFromIndex]=allFastaDict[idFromIndex]
      #if counter >100: break


  writeLocation = "%s/fileToBlast.fasta" %(getStringFromSettings("PDBALLFASTADATABASELOCATION",settings))
  newFile = open(writeLocation, "w")
  for i in seqToBlast:
      fastaEntry = ">%s\n%s\n"%(i,seqToBlast[i])
      newFile.write(fastaEntry)
  #then print all to file
  newFile.close()










# WORK IN PROGRESS!!

#  uniprotfolder = getStringFromSettings("UNIPROTSERVERFOLDER", settings)
#  uniprotfile = getStringFromSettings("UNIPROTSERVERFILE", settings)
#  if not dirExists(getStringFromSettings("UNIPROTDATABASELOCATION", settings)):
#    if getBoolFromSettings("AUTOCREATEMISSINGDATABASELOCATIONS",settings) and dirExists(getStringFromSettings("DEFAULTDATABASELOCATION",settings)):
#      sendToConsole("mkdir %s/uniprotHMMs" % (getStringFromSettings("DEFAULTDATABASELOCATION",settings)))
#      settings["UNIPROTDATABASELOCATION"]= "%s/uniprotHMMs"%(getStringFromSettings("DEFAULTDATABASELOCATION",settings))
#    else:
#      throwError("UNIPROTDATABASELOCATION does not exist: %s"%(getStringFromSettings("UNIPROTDATABASELOCATION",settings)),settings)

#  dbloc = getStringFromSettings("UNIPROTDATABASELOCATION",settings)
#  if dirExists(dbloc):
#    sendToConsole("cd %s;wget %s%s > download.log 2>&1"%(dbloc,uniprotfolder, uniprotfile))
#    if fileExists("%s/%s"%(dbloc,uniprotfile)):
#      sendToConsole("tar zxvf %s/%s -C %s >> download.log 2>&1"%(dbloc,uniprotfile,dbloc))
#      #clean up the tgz file
#      sendToConsole("rm %s/%s"%(dbloc, uniprotfile))
#      #now move the files one folder up.
#      sendToConsole("mv %s/%s/* %s"%(dbloc, uniprotfile[:-4], dbloc))
#      #!!maybe implement later to remove the empty folder too
#    else:
#      throwError("UNIPROT HMM files could not be downloaded! Aborting", settings)






#blast that file
def blastAgainstComplex(settings):
  blastHome = getStringFromSettings("BLASTHOME",settings)
  #partsHolder = getStringFromSettings("COMPLEXFASTAFILE", settings).split('/')
  #complexFastaFileName = partsHolder[len(partsHolder)-1]
  complexFastaFileName = getStringFromSettings("COMPLEXDATABASENAME", settings)
  dbloc = "%s/%s"%(getStringFromSettings("COMPLEXDATABASELOCATION",settings),complexFastaFileName)
  fileToBlast = "%s/fileToBlast.fasta" %(getStringFromSettings("PDBALLFASTADATABASELOCATION",settings))
  outputLoc = "%s/output.blast" %(getStringFromSettings("PDBALLFASTADATABASELOCATION",settings))
  evalue = getStringFromSettings("BLAST_E-VALUE",settings)
  blastCommand = "%s/blastp -db %s -query %s -out %s -outfmt 5 -evalue %s -num_threads 1"%(blastHome,dbloc, fileToBlast,outputLoc, evalue)
  #include all HSP!!? how? -> already done?

  # use -num_threads 01?
  sendToConsole(blastCommand)
  #"%sblastp -db %s -query %s -evalue %s"%(blasthome,dblocation,   complexfastafile, evalue)




#deprecated soon? -> currently not used
#-> maybe use in next script, or create a new db?
def getComplexNames(ali): #get all complex pdb names from the description in Li's complexfile
  complexNameDict = {}
  cTitle = ali.split('>gi|')
  for t in cTitle:
    if t.strip().startswith('gi|'):  #maybe double check why needed...
      t = t[3:]
      #print "t changed to: %s"%(t)
    ts = t.strip()
    if not ts == "":
      cDescr = ts.split()
      #print cDescr[0]
      #print "\n\n%s\n\n"%(ali)
      cdparts = cDescr[0].split('|')
      complexID = cdparts[2].lower()
      complexChain = cdparts[3].lower()
      if complexID in complexNameDict:
        complexNameDict[complexID].append(complexChain)
      else:
        complexNameDict[complexID]= [complexChain]
  return complexNameDict

#deprecated soon? -> currently not used
def processBlastOutput(resultDict, alignDict, settings):
  blastOutputFile = "%s/output.blast" %(getStringFromSettings("PDBALLFASTADATABASELOCATION",settings))

  result_handle = open(blastOutputFile)
  blast_records = NCBIXML.parse(result_handle)
  for blast_record in blast_records:
    queryName = blast_record.query
    resultDict[queryName]={}
    for alignment in blast_record.alignments:
      alignmentNumber = len(alignDict)
      alignDict[alignmentNumber] = alignment
      complexNameDict = getComplexNames(alignment.title)
      for c in complexNameDict:
        if not c in resultDict[queryName]:
          resultDict[queryName][c]={}
        for chain in complexNameDict[c]:
          if not chain in resultDict[queryName][c]:
            resultDict[queryName][c][chain]= alignmentNumber
            #print "DEBUG: res -\t%s\t%s\t%s\t%s" %(queryName, c, chain, str(alignmentNumber)) #for debugging
          else: print "error! chain already present!"






def checkIfQualifiedHit(alignment, settings):

  #check if there are not more than 1 HSP in the hit.
  if len(alignment.hsps)>1:
    print "Found more than 1 hsp: %s" % (alignment.title) #comment out later
    return False
  elif len(alignment.hsps)<1:
    print "No hsp found!: %s" % (alignment.title)
    return False

  #check if Sequence identity (SID) is high enough
  sidThreshold = int(getStringFromSettings("BLASTSIDTHRESHOLD", settings))
  hsp = alignment.hsps[0] #not a dict but a list!!
  sidHit = int(float(hsp.identities)/float(hsp.align_length)*100)
  if sidHit < sidThreshold:
    print "SID to low: %s for %s" % (sidHit,alignment.title )
    return False
  return True




#def checkIfQualifiedHit(blastEntry, settings):
#  #check if all conditions are met to be a valid hit.
#  #e-value is already checked by blast.
#  if len(blastEntry.alignments)>1:
#    #FIX THIS LATER!! CANNOT IGNORE MULTIPLE HITS ON SAME
#    print "Error: more than one alignment after all!"
#    return False
#  #else: entry=blastEntry.alignments.values()[0]
#  elif len(blastEntry.alignments)>0: entry=blastEntry.alignments[0] #not a dict but a list!!
#  else:
#    print "ERROR, no alignment found."
#    return False
#
#  #check if there are not more than 1 HSP in the hit.
#  if len(entry.hsps)>1:
#    print "Found more than 1 hsp: %s" % (entry.title) #comment out later
#    return False
#  elif len(entry.hsps)<1:
#    print "No hsp found!: %s" % (entry.title)
#    return False
#
#  #check if Sequence identity (SID) is high enough
#  sidThreshold = int(getStringFromSettings("BLASTSIDTHRESHOLD", settings))
#  #hsp = entry.hsps.values()[0]
#  hsp = entry.hsps[0] #not a dict but a list!!
#  sidHit = int(float(hsp.identities)/float(hsp.align_length)*100)
#  if sidHit < sidThreshold:
#    print "SID to low: %s for %s" % (sidHit,entry.title )
#    return False
#  return True


def fillLookupTable(settings):
  lookupTableLocation = getStringFromSettings("LOOKUPDATABASELOCATION", settings)
  blastOutputFile = "%s/output.blast"%(getStringFromSettings("PDBALLFASTADATABASELOCATION",settings))
  result_handle = open(blastOutputFile)
  blast_records = NCBIXML.parse(result_handle)
  for blast_record in blast_records:
    if len(blast_record.alignments)>0:
      queryName = blast_record.query
      recordFolder = "%s/%s"%(lookupTableLocation,queryName)
      sendToConsole("mkdir %s"%(recordFolder))
      #for alignment in blast_record:
      for x in range(0,len(blast_record.alignments)):
        alignment = blast_record.alignments[x]
        if checkIfQualifiedHit(alignment, settings):

          pickleName = "%s/%s.pickle"%(recordFolder,str(x))
          output = open(pickleName, 'wb')
          pickle.dump(alignment, output)
          output.close()

#    if checkIfQualifiedHit(blast_record, settings):#how to get this alignment?
#      queryName = blast_record.query
#      recordFolder = "%s/%s"%(lookupTableLocation,queryName)
#      sendToConsole("mkdir %s"%(recordFolder))
#      pickleName = "%s/pickle.pickle"%(recordFolder)
#      output = open(pickleName, 'wb')
#      pickle.dump(blast_record, output)
#      output.close()






#def pickleDatastructure(item1, item2, settings):
#  #print "packing file ... (to be implemented)"
#  pickleLoc= getStringFromSettings("PICKLELOCATION",settings)
#  output = open(pickleLoc, 'wb')
#  pickle.dump(item1, output)
#  pickle.dump(item2, output)
#  output.close()








