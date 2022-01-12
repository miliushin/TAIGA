#include "data.h"

DATA::DATA(const char *File, int Type) {
  if (Type == 0) {
    NST = 3;
    FirstNST = 1;
  } else if (Type == 1) {
    NST = 19;
    FirstNST = 31;
  }

  for (int i = 0; i < Aperture; i++) {
    XAperture[i] = i + 1.;
  }
  INT_ADC = 0.;
  AMP_ADC = 0;

  int len = strlen(File);
  char Date[7];
  Date[6] = 0;
  for (int i = 0; i < 6; i++) {
    Date[i] = File[len - 11 + i];
  }
  char OutName[18];
  OutName[17] = 0;
  sprintf(OutName, "OUTPUT/%s.root", Date);

  char cNch[4];
  for (int i = 0; i < NUM_OF_CHANNELS; i++) {
    sprintf(cNch, "I%d", i);
    HInteg[i] = new TH1F(cNch, " ", 50000, 0, 50000);
    HInteg[i]->SetTitle(";Integral;Nevents");
    sprintf(cNch, "A%d", i);
    HAmpl[i] = new TH1F(cNch, " ", 1000, 0, 1000);
    HAmpl[i]->SetTitle(";Amplitude; Nevents");
  }
  HDeltaT = new TH1F("HDeltaT", " ", 1000, 0, 1);
  HDeltaT->SetTitle(";#Delta t, s; Nevents");

  RootData = new TFile(File, "READ");
  RootOut = new TFile(OutName, "UPDATE");
  printf("Read: %s\n", File);

  //----------Read Data Tree-------------------------------------
  char NstName[4];

  for (int i = FirstNST; i < (NST + FirstNST); i++) {
    sprintf(NstName, "%d", i);
    TData[i - FirstNST] = (TTree *)RootData->Get(NstName);
    TData[i - FirstNST]->SetBranchAddress("TimeEvent", &TimeEvent);
    TData[i - FirstNST]->SetBranchAddress("ADC", ADC);
    // read all entries tree
    Long64_t nentries = TData[i - FirstNST]->GetEntries();
    for (Long64_t j = 0; j < nentries; j++) {
      if (j > MAX_OF_EVEN) {
        printf("ERROR: Events more MAX_OF_EVEN = %d\n", MAX_OF_EVEN);
        exit(1);
      }
      TData[i - FirstNST]->GetEntry(j);
      TimeSt[j] = TimeEvent;
      if (j > 0) {
        DeltaT = TimeSt[j] - TimeSt[j - 1];
        if (DeltaT < 0) DeltaT *= -1.;
        HDeltaT->Fill(DeltaT);
      }

      for (int k = 0; k < NUM_OF_CHANNELS; k++) {
        GetIA(k, ADC[k]);
        HInteg[k]->Fill(INT_ADC);
        HAmpl[k]->Fill(AMP_ADC);
      }
    }
    //----------Write Hist--------------------------------------
    RootOut->cd();
    gDirectory->cd(NstName);
    gDirectory->cd("Amplitude");
    for (int j = 0; j < NUM_OF_CHANNELS; j++) {
      HAmpl[j]->Write();
      HAmpl[j]->Reset();
    }

    gDirectory->cd("../");
    gDirectory->cd("Integral");
    for (int j = 0; j < NUM_OF_CHANNELS; j++) {
      HInteg[j]->Write();
      HInteg[j]->Reset();
    }

    gDirectory->cd("../");
    gDirectory->cd("DeltaT");
    HDeltaT->Write();
    HDeltaT->Reset();
  }
}

DATA::~DATA() {
  for (int i = 0; i < NUM_OF_CHANNELS; i++) {
    delete HAmpl[i];
    delete HInteg[i];
  }
  delete HDeltaT;

  RootOut->Close();
  RootData->Close("R");
}

void DATA::GetIA(int Nch, int data[Aperture]) {
  double New[Aperture];
  int MaxRateNoise;
  int NumberMaxRateNoise = 0;
  int Noise[Aperture];
  int RateNoise[Aperture];
  int MaxAplitude = 0;
  int iMaxAmplitude = 0;

  Noise[0] = data[0];
  RateNoise[0] = 0;
  int iNoise = 0;

  for (int i = 0; i < Aperture; i++) {
    Noise[i] = 0;
    RateNoise[i] = 0;
  }
  bool FlagNoise = false;

  for (int i = 0; i < Aperture; i++) {
    FlagNoise = false;
    for (int j = 0; j <= iNoise; j++) {
      if (data[i] == Noise[j]) {
        RateNoise[j]++;
        FlagNoise = true;
        break;
      }
    }
    if (FlagNoise == false) {
      iNoise++;
      Noise[iNoise] = data[i];
      RateNoise[iNoise] = 1;
    }
  }

  MaxRateNoise = 0;

  for (int i = 0; i <= iNoise; i++) {
    if (RateNoise[i] > MaxRateNoise) {
      MaxRateNoise = RateNoise[i];
      NumberMaxRateNoise = i;
    }
  }

  for (int i = 0; i < Aperture; i++) {
    New[i] = (double)(data[i] - Noise[NumberMaxRateNoise]);

    if (i >= 0 && ((New[i] <= NOISE_RANGE && New[i] >= -NOISE_RANGE))) {
      New[i] = 0;
    }

    if (New[i] > MaxAplitude) {
      MaxAplitude = New[i];
      iMaxAmplitude = i;
    }
  }

  int k = 0;
  bool Flag = false;
  for (int i = 0; i < Aperture; i++) {
    if (New[i] != 0.) {
      for (int j = i; j < Aperture; j++) {
        New[j] = data[j] - Noise[NumberMaxRateNoise];
        if (New[j] <= 0.) {
          New[j] = 0.;
          k = j;
          Flag = true;
          break;
        }
      }
      if (Flag == true) {
        for (int j = i; j >= 0; j--) {
          New[j] = data[j] - Noise[NumberMaxRateNoise];
          if (New[j] <= 0.) {
            New[j] = 0.;
            i = k;
            Flag = false;
            break;
          }
        }
      }
    }
  }

  GADC = new TGraph(Aperture, XAperture, New);

  INT_ADC = GADC->Integral(0, Aperture);
  AMP_ADC = MaxAplitude;

  delete GADC;
}

void DATA::AddRoot() {}