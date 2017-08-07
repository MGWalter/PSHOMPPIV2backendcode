#!/usr/bin/python
#written by mwalter
#08 June 2016

import subprocess
import os

def fileExists(fileName):
  return os.path.exists(fileName)


def backup(folderName,backupLoc):
  sendToConsole("mv --backup=numbered %s %s/"%(folderName,backupLoc))
  sendToConsole("mkdir %s"%(folderName))

def dirExists(dirName):
  return os.path.isdir(dirName)

def throwSettingsError(argType, settingsDict):
      throwError("%s is not a valid settings syntax! Terminating program."%(argType), settingsDict)
        #later instead of exiting here, throw error to try statement

def parseArg(argument, settingsDict):
  argumentWithoutComments = argument.split("#")
  #argument = argumentWithoutComments[0]
  argData = argumentWithoutComments[0].split('=')
  if argumentWithoutComments[0].strip() != "":
    if len(argData) != 2:
      #throw error
      throwSettingsError(argumentWithoutComments[0],settingsDict)
    else:
      argType = argData[0].strip()
      checkArgExists(argType, settingsDict)
      settingsDict[argType] = argData[1].strip()

def sendToConsole(command):
    subprocess.call(command, shell=True)

def fillSettings(configLoc, settingsDict):
  #read from the config textfile, all default settings.
  with open(configLoc) as settingsFile:
    for line in settingsFile:
      parseArg(line, settingsDict)

def printSettings(settingsDict):
  for s in settingsDict:
    print "%s\t=\t%s"%(s,settingsDict[s])

def getFloatFromSettings(param,settingsDict):  #be very careful with rounding!
  checkArgExists(param, settingsDict)
  try:
    return float(settingsDict[param])
  except ValueError:
    throwError("Setting %s is not a float! Terminating." %(param), settingsDict)

def getBoolFromSettings(param,settingsDict):
  checkArgExists(param, settingsDict)
  if settingsDict[param] in ["1", "TRUE", "true", "True"] :
    return True
  else:
    return False

def getStringFromSettings(param, settingsDict):
  checkArgExists(param, settingsDict)
  return settingsDict[param]

def checkArgExists(argType, settingsDict):
  if not argType in settingsDict:
    throwSettingsError(argType, settingsDict)

def throwError(errMsg, settingsDict):
  #whenever this happens, write to error file and exit program if RUNNING setting is false
  if not "RUNNING" in settingsDict:
    errMsg = "Error: RUNNING variable is missing!!"
    sendToConsole('echo "'+errMsg+'">ERROR.log;')

  elif getBoolFromSettings("RUNNING",settingsDict):
    #create error file at correct location.
    print "For now print error msg while running: \t%s"%(errMsg)
    sendToConsole('echo "'+errMsg+'">ERROR.log;')
  else:
    print "Error:\t%s"%(errMsg)
    sendToConsole('echo "'+errMsg+'">ERROR.log;')
    exit()

