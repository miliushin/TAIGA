#ifndef DATA_H
#define DATA_H

#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TTree.h>
#include <TSystem.h>
#include <cstdio>



#include "globals.h"

using namespace std;

class DATA {
 public:
  DATA(const char *File);  
  ~DATA();

 private:
  TFile  *RootData;            // File DDMMYY.root with data
  TFile  *RootOut;             // Ouyput file
  TFile  *Root_tmp;
  TTree  *TData[MAX_OF_STAT];  // Max number of stations
  FILE   *FileTMP;

  TGraph *GADC;
  TF1    *FitFunc; 
  TH1F   *HInteg[NUM_OF_CHANNELS]; // Output hist
  TH1F   *HAmpl[NUM_OF_CHANNELS];
  TH1F   *HDeltaT;

  int         NST;       // Number of stations
  int         FirstNST;  // Number of first station

  double      TimeEvent;                       // Tree format data time
  int         ADC[NUM_OF_CHANNELS][Aperture];  // Tree format data ADC


  double      XAperture[Aperture];  // X-axis of ADC
  double      INT_ADC;              // Integral of PMT pulse
  int         AMP_ADC;              // Max amplitude of PMT pulse
  double      TIME_ADC;             // Time of PMT pulse
  double      TimeSt[MAX_OF_EVEN];  // Time Station in sec
  double      DelT;               // Time between events

  double      TimeCounter[NUM_OF_CHANNELS];
  double      IntegCounter[NUM_OF_CHANNELS];
  int         AmplCounter[NUM_OF_CHANNELS];

  int         GetIA_MUON(int data[Aperture]);  // Geting INT_ADC and AMP_ADC
  int         GetIA_GRANDE(int Nch,int data[Aperture]);
  double      FindLevelOneParticle(TH1F *hist); 

};

#endif