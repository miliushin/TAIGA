#ifndef MUON_H
#define MUON_H

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>

#include <fstream>

#include "globals.h"

class MUON {
 public:
  MUON(char *name_file);
  ~MUON();

  double GetTime(unsigned char ch[8]);
  int ReadADC(FILE *data_file);
  int GetAmplitude(int ch, double d[Aperture]);

  unsigned int ID;                              // Station indentification
  unsigned int NumBytes;                        // Data packet size
  unsigned int NumSt;
  int data[NUMBER_OF_CHANEL][Aperture];
  unsigned int fin;

 private:
  double TimeEvent;
  double ADC[NUMBER_OF_CHANEL][Aperture];
  TFile *RootFile;

  unsigned long int END_File;
  unsigned long int i_File;

  Int_t NumEvent;  // Number of event

  TGraph *GData;
  TGraph *GDataNew;
  TCanvas *can;
  TTree *TData;
};

#endif