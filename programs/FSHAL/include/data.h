#ifndef DATA_H
#define DATA_H

#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TTree.h>

#include "globals.h"

using namespace std;

class DATA {
 public:
  DATA(const char *File, int Type);  // Type = 0 - MUON, 1 - GRANDE
  ~DATA();

 private:
  TFile  *RootData;            // File DDMMYY.root with data
  TFile  *RootOut;             // Ouyput file
  TTree  *TData[MAX_OF_STAT];  // Max number of stations

  TGraph *GADC;
  TH1F   *HInteg[NUM_OF_CHANNELS]; // Output hist
  TH1F   *HAmpl[NUM_OF_CHANNELS];
  TH1F   *HDeltaT;

  int    NST;       // Number of stations
  int    FirstNST;  // Number of first station

  double TimeEvent;                       // Tree format data time
  int    ADC[NUM_OF_CHANNELS][Aperture];  // Tree format data ADC

  double XAperture[Aperture];  // X-axis of ADC
  double INT_ADC;              // Integral of PMT pulse
  int    AMP_ADC;              // Max amplitude of PMT pulse
  double TimeSt[MAX_OF_EVEN];  // Time Station in sec
  double DeltaT;               // Time between events

  void   GetIA(int Nch, int data[Aperture]);  // Geting INT_ADC and AMP_ADC
  void   AddRoot();
};

#endif