#ifndef GRANDE_H
#define RANDE_H

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>

#include <fstream>

#include "globals.h"

class GRANDE {
 public:
  GRANDE(char *name_file);
  ~GRANDE();

  double GetTime(unsigned char ch[8]);
  int    ReadADC(FILE *data_file);

 private:
  unsigned int ID;                                // Station indentification
  unsigned int NumBytes;                          // Data packet size
  unsigned int Delay;                             // The optic line
  unsigned int Claster;
  unsigned int fin;
  double TimeEvent;
  int ADC[NUMBER_OF_CHANNELS][Aperture];          // ADC data
  unsigned long int END_File;
  unsigned long int i_File;

  Int_t NumEvent;                                 // Number of event
  TTree *TData;                                   // Data TTree
  TFile *RootFile;                                // Output file: ddmmyy.root

  long Addr;
  long AddrBase[16]={ 0x80000,0x80800,0x84000,0x84800,
                        0xA0000,0xA0800,0xA4000,0xA4800,
                        0xC0000,0xC0800,0xC4000,0xC4800,
                        0xE0000,0xE0800,0xE4000,0xE4800 };
};

#endif