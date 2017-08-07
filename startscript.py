#!/usr/bin/python
#written by mwalter
# 7 June 2016

import os.path
import sys
sys.path.append('/home/mwalter/project/v0.0.02/projectfiles/CODE-PPI-HHPRED')
from functionsstartscript import *

settings={}  ## NOTE ALL SETTINGS SAVED AS STRING!!
workingFolder = "/home/mwalter/project/v0.0.02/projectfiles"


settings["CONFIGFILENAME"]="%s/settings.config"%(workingFolder)
settings["WORKINGFOLDER"]= workingFolder
settings["BLAST_E-VALUE"] = "0.001"
settings["DIAGNOSTICS"] = "FALSE"
settings["FILENAME"]="./testscript1.py"
settings["RUNNING"]="FALSE"
settings["BLASTHOME"]="/home/sbgrid/programs/x86_64-linux/blastplus/2.3.0/bin"
settings["COMPLEXFASTAFILE"]="%s/nr_pdbaa_s2c_new"%(workingFolder)
settings["COMPLEXDATABASELOCATION"]="%s/databases/complexDB"%(workingFolder)
settings["COMPLEXDATABASENAME"]="nr_pdbaa_s2c_new"

settings["AUTOCREATEMISSINGDATABASELOCATIONS"]="TRUE"
# if AUTOCREATEMISSINGDATABASELOCATIONS==True then when other db folders do not exist, the program will put the DB in this folder without throwing error.
]]tings["DEFAULTDATABASELOCATION"]="%s/databases"%(workingFolder)
settings["RECREATEDATABASES"]="FALSE"
settings["PDBALLFASTASERVERFOLDER"]= "ftp://ftp.wwpdb.org/pub/pdb/derived_data/"
settings["PDBALLFASTASERVERFILE"] = "pdb_seqres.txt"
settings["PDBALLFASTADATABASELOCATION"] = "%s/databases/PDBallFASTA"%(workingFolder)
settings["BACKUPPREVIOUSDATABASES"] = "TRUE"
settings["BACKUPLOCATION"]= "%s/DBbackups"%(workingFolder)


settings["HHLINKSFILE"]="/home/mwalter/hh/scripts/HHPaths.pm"

#settings["PSIPREDBIN"]="/home/sbgrid/programs/x86_64-linux/psipred/3.2.1/bin";  # path to PSIPRED V2 binaries
settings["PSIPREDBIN"]="/home/mwalter/software/PSIPRED/bin"  # path to PSIPRED V2 binaries
#settings["PSIPREDDATA"]= "/home/sbgrid/programs/x86_64-linux/psipred/3.2.1/data"; # path to PSIPRED V2 data files
settings["PSIPREDDATA"]= "/home/mwalter/software/PSIPRED/data" # path to PSIPRED V2 data files
settings["NCBIBIN"]= "/home/sbgrid/programs/x86_64-linux/blast/2.2.26/bin"    # path to NCBI binaries (for PSIPRED in addss.pl)
settings["DSSPBIN"]= "/home/sbgrid/programs/x86_64-linux/dssp/2.2.0/bin/mkdssp"  # where is the dssp binary? Used in addss.pl
#settings["PDBSTRUCTURES"]= "/home/mwalter/databases/structuresPDB"            # where are the pdb files? (pdb/divided directory will also work)
settings["DSSPDATABASELOCATION"]= "%s/databases/dsspDB"%(workingFolder)       # where are the dssp files? Used in addss.pl

settings["UNIPROTDATABASELOCATION"]="%s/databases/uniprotHMMs"%(workingFolder)
settings["PDB70DATABASELOCATION"]="%s/databases/pdb70HMMs"%(workingFolder)

settings["PDB70SERVERFOLDER"]="http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/"
settings["PDB70SERVERFILE"]="pdb70_24Mar16.tgz"

settings["UNIPROTSERVERFOLDER"]="http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/"
settings["UNIPROTSERVERFILE"]="uniprot20_2016_02.tgz"

settings["DSSPSERVERFOLDER"]="rsync.cmbi.ru.nl/dssp"


settings["LOOKUPDATABASELOCATION"]="%s/databases/lookuptable" %(workingFolder)

settings["BLASTSIDTHRESHOLD"]= "70" # ust be in int format, for now



configFileName = getStringFromSettings("CONFIGFILENAME",settings)
if fileExists(configFileName):
  fillSettings(configFileName, settings) #get all settings from the config textfile
else: throwError("Config file missing!: %s"%(configFileName),settings)

argList = sys.argv
argList[0]= "FILENAME=%s"%(argList[0])
for argument in argList:
  parseArg(argument, settings) # apply any settings given in command line argument (overriding the defaults)

#to see default settings
if getBoolFromSettings("DIAGNOSTICS", settings):  # if program is called with DIAGNOSTICS=true as an argument, do this:
  printSettings(settings)
  exit()




#backup any of the existing databases:


# if default database directory doesnt exist, make one in workingfolder.
if not dirExists(getStringFromSettings("DEFAULTDATABASELOCATION",settings)):
  sendToConsole("mkdir %s/databases"%(getStringFromSettings("WORKINGFOLDER",settings)))
  settings["DEFAULTDATABASELOCATION"]="%s/databases"%(getStringFromSettings("WORKINGFOLDER",settings))


if getBoolFromSettings("RECREATEDATABASES",settings) and getBoolFromSettings("BACKUPPREVIOUSDATABASES",settings):
  backupAllDatabases(settings)




if getBoolFromSettings("RECREATEDATABASES",settings):

  #---------------------------------------------------------------
  #make a blast database out of the complexfastafile
  makeBlastDatabaseFromComplexList(settings)
  #---------------------------------------------------------------
  #Download allPDB in fasta format.
  downloadAllPdb(settings)
  #---------------------------------------------------------------


  #-dssp structures
  downloadDsspFiles(settings)
  #-uniprot
  downloadUniprotHmm(settings)
  #-pdb70 hhms
  downloadPdb70Hmm(settings)
  #-pdb structures ?


makeFileToBlast(settings) #make a file in PDBALLFASTA to blast against PPIDB complex blast DB.

blastAgainstComplex(settings) #blast HMM master sequences against PPIDB blast DB.

cleanLookupTable(settings)# backup the old lookup table


fillLookupTable(settings)



fixHhsuiteLinks(settings)
#--------------------------------------------------------------
print "program done"

#-----------------------------------------------------


