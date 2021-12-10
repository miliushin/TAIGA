#!/usr/bin/env python

import sys
import ROOT
import numpy as np

if len(sys.argv) < 2:
	print(" Usage: deltaT.py data_time.txt")
	sys.exit(1)

File = open(sys.argv[1],'r')
outputFile = "out.root"

data = np.genfromtxt(File, delimiter='\t', dtype=np.float)

maxT = 1.
Nt = 100
hist = ROOT.TH1F("time", " ", Nt, 0, maxT)
hist.SetTitle(";#Delta t, s;Nevent, 1/s");

EndTime = data[len(data)-1]
StartTime = data[0]

weight = 1./(EndTime - StartTime)

for i in range(0,len(data)-1):
	delta = data[i+1] - data[i]
	hist.Fill(delta,weight)
	
out_root = ROOT.TFile(outputFile,"RECREATE")
hist.Write("time")
File.close()