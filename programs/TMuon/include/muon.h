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

  unsigned int ID;                              // Station indentification
  unsigned int NumBytes;                        // Data packet size
  int data[NUMBER_OF_CHANEL][Aperture];         // ADC data
  unsigned int fin;

 private:
  double TimeEvent;
  double ADC[NUMBER_OF_CHANEL][Aperture];
  TFile *RootFile;

  unsigned long int END_File;
  unsigned long int i_File;

  Int_t NumEvent;  // Number of event
  TTree *TData;
};

#endif