#!/usr/bin/env python

import sys
import ROOT


# Read file with energy deposite and PDG-id
if len(sys.argv) < 3:
	print(" Usage: readED.py Energy.txt ParticleID.txt")
	sys.exit(1)

FileED = open(sys.argv[1],'r')
FilePN = open(sys.argv[2],'r')

outputFile = "out.root"

dataED = [map(float, line.split(" ")) for line in FileED]
dataPN = [map(float, line.split(" ")) for line in FilePN]


all_PDG_ID = []
for i in range(0,len(dataPN)):
	for j in range(0,len(dataPN[i])):
		if dataPN[i][j] != 0:
			all_PDG_ID.append(dataPN[i][j])

PDG_ID = sorted(set(all_PDG_ID), key=lambda d: all_PDG_ID.index(d))


total  = ROOT.TH1F( 'total', 'energy deposite, MeV', 100, 0.0, 20 )
can = ROOT.TCanvas( 'can', 'energy deposite', 200, 10, 600, 400 )
leg = ROOT.TLegend(0.6,0.6,0.8,0.8)

histNE = []
Integ  = []
j = 0
for i in range(0,len(PDG_ID)):
	histNE.append(ROOT.TH1F("ED_"+str(i), "energy deposite", 100, 0.0, 20))
	if i!=4:
		j = j + 1
	else:
		j = j + 2
	histNE[i].SetLineColor(j)
	histNE[i].SetLineWidth(4)
	leg.AddEntry(histNE[i],str(PDG_ID[i]),"f")

total.SetTitle("; Total ED, MeV; Nevent")
total.SetFillColor(5)
total.SetLineColor(4)
leg.AddEntry(total,"ALL","f")

for i in range(0,len(dataED)):
	for j in range(0,len(dataED[i])):
		if dataED[i][j] != 0.:
			index = PDG_ID.index(dataPN[i][j])
			histNE[index].Fill(dataED[i][j])
			total.Fill(dataED[i][j])

for i in range(0,len(PDG_ID)):
	Integ.append(histNE[i].Integral())

SortedInteg = sorted(range(len(Integ)), key=lambda k: Integ[k])


total.Draw()
for i in range(0,len(SortedInteg)):
	histNE[int(SortedInteg[i])].Draw("same")

out_root = ROOT.TFile(outputFile,"RECREATE")

leg.Draw()
can.Write("hist")

FileED.close()
FilePN.close()