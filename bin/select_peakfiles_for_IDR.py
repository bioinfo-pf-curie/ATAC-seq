#!/usr/bin/env python

##################### Imports ##################################
import argparse
import pandas as pd
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--caller",help = "Precise the peak caller used to produce the files Genrich or Macs2", required = True)
parser.add_argument("--peakfiles", nargs="+",help="Either Genrich or Macs2 peakfiles")
args = parser.parse_args()
#print("")
#print("Caller = : " + args.caller)
#print("Looking for two biggest peakfiles among the " + str(len(args.peakfiles)) + "provided... \n")

filenames = []
numberofpeaks = []

################### Extract number of peaks in each peak file #####################################
if args.caller == "Macs2":
  for f in args.peakfiles:
      try:
        peaks = pd.read_csv(f, sep ="\t",skiprows=1,header=None).loc[:,3].str.replace('.*_peak_', '').str.replace('[a-zA-Z]','').unique()
      except IOError:
        print("Oops! " + f + " is not a valid tabulated file")
        quit()

      filenames.append(f)
      peaks = peaks.astype(np.int)
      peaks = np.append(peaks,0)
      numberofpeaks.append(max(peaks))

elif args.caller == "Genrich":
  for f in args.peakfiles:
      try:
         peaks = pd.read_csv(f, sep ="\t",skiprows=1,header=None).loc[:,3].str.replace('peak_', '').unique()
      except IOError:
        print("Oops! " + f + " is not a valid tabulated file")
        quit()

      filenames.append(f)
      peaks = peaks.astype(np.int)
      peaks = np.append(peaks,0)
      numberofpeaks.append(max(peaks))


########### Select two peakfiles with the biggest number of peaks and print it to shell #############################
selected = np.argsort(numberofpeaks)[-2:]
selected = [filenames[fileindex] for fileindex in selected]
print(selected[0] + ' ' + selected[1])
