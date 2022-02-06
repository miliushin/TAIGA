#include "TCanvas.h"
#include "TStyle.h"

#include "data.h"

DATA::DATA(const char *File) {

  //gStyle->SetOptFit(0000);
  FileTMP = fopen("point.txt","w");

  for (int i = 0; i < Aperture; i++) {
    XAperture[i] = i + 1.;
  }

  int len = strlen(File);
  char Date[7];
  Date[6] = 0;
  for (int i = 0; i < 6; i++) {
    Date[i] = File[len - 11 + i];
  }
  char OutName[18];
  OutName[17] = 0;
  sprintf(OutName, "OUTPUT/%s.root", Date);

  char OutName_tmp[22];
  OutName_tmp[21] = 0;
  sprintf(OutName_tmp, "OUTPUT/%s_tmp.root", Date);


  char cNch[4];
  for (int i = 0; i < NUM_OF_CHANNELS; i++) {
    sprintf(cNch, "I%d", i);
    HInteg[i] = new TH1F(cNch, " ", 50000, 0, 50000);
    HInteg[i]->SetTitle(";Integral;Nevents");
    sprintf(cNch, "A%d", i);
    HAmpl[i] = new TH1F(cNch, " ", 2048, 0, 2048);
    HAmpl[i]->SetTitle(";Amplitude; Nevents");
  }
  HDeltaT = new TH1F("HDeltaT", " ", 1000, 0, 1);
  HDeltaT->SetTitle(";#Delta t, s; Nevents");

  RootData = new TFile(File, "READ");
  RootOut  = new TFile(OutName, "UPDATE");
  Root_tmp = new TFile(OutName_tmp, "UPDATE");
  printf("Read: %s\n", File);

  //----------Read Data Tree-------------------------------------
  const char *NstName;
  int NumberOfTrees = RootData->GetListOfKeys()->GetSize();

  for (int i = 0; i < NumberOfTrees; i++) {
    NstName = RootData->GetListOfKeys()->At(i)->GetName();

    TTree *Ttmp = (TTree*)Root_tmp->Get(NstName);
    Ttmp->SetBranchAddress("Time", TimeCounter);
    Ttmp->SetBranchAddress("Amplitude", AmplCounter);
    Ttmp->SetBranchAddress("Integral", IntegCounter);


    TData[i] = (TTree *)RootData->Get(NstName);
    printf("%s\n", NstName);

    TData[i]->SetBranchAddress("TimeEvent", &TimeEvent);
    TData[i]->SetBranchAddress("ADC", ADC);
    // read all entries tree

      //int a1 = 0;
      //int a2 = 0;
      //double dif = 0;
      //TH1F *hh = new TH1F("hh","",400,-2,2);


    Long64_t nentries = TData[i]->GetEntries();
    for (Long64_t j = 0; j < nentries; j++) {

      if (j > MAX_OF_EVEN) {
        printf("ERROR: Events more MAX_OF_EVEN = %d\n", MAX_OF_EVEN);
        exit(1);
      }
      TData[i]->GetEntry(j);

      TimeSt[j] = TimeEvent;
      if (j > 0) {
        DeltaT = TimeSt[j] - TimeSt[j - 1];
        if (DeltaT < 0) DeltaT *= -1.;
        HDeltaT->Fill(DeltaT);
      }

      int grande_status = 10;
      bool FlagStat = false;



      for (int k = 0; k < NUM_OF_CHANNELS; k++) {
        AmplCounter[k] = 0;
        IntegCounter[k] = 0.;
        TimeCounter[k] = 0.;
     
        INT_ADC  = 0.;
        AMP_ADC  = 0;
        TIME_ADC = 0.;
        

      
        if(atoi(NstName)<=3) GetIA_MUON(ADC[k]);
        else if(atoi(NstName)>3) {

          if (k==0){
            grande_status = GetIA_GRANDE(k,ADC[k]);
            //a1 = AMP_ADC;
            if (grande_status == 0){
              INT_ADC /= 10.;
              AMP_ADC /= 10;
              FlagStat = true;

            }
          }
          /*
          if (k==1 && FlagStat==true){
            grande_status=GetIA_GRANDE(k,ADC[k]);
            double a2 = double(AMP_ADC);

            if(grande_status==0) {
              dif = (double(a1)/10. - a2)/a2;
            }
          }
          */

          if (k==1 && FlagStat == false){
            grande_status = GetIA_GRANDE(k,ADC[k]);
            if(grande_status==0) FlagStat = true;
          }
          else if(k==1) {FlagStat = false; continue;}
      
          if (k==2){
            grande_status = GetIA_GRANDE(k,ADC[k]);
            if (grande_status == 0){
              INT_ADC /= 10.;
              AMP_ADC /= 10;
              FlagStat = true;
            }
          }
          if (k==3 && FlagStat==false){
            grande_status = GetIA_GRANDE(k,ADC[k]);
            if(grande_status==0) FlagStat = true;
          }
          else if(k==3) {FlagStat = false; continue;}

        }
        //hh ->Fill(dif);

        int jj;

        if(AMP_ADC>0 && atoi(NstName)<=3){
          HInteg[k]->Fill(INT_ADC);
          HAmpl[k]->Fill(AMP_ADC);     
          AmplCounter[k] = AMP_ADC;
          IntegCounter[k] = INT_ADC;
          TimeCounter[k] = TIME_ADC + TimeEvent;
        }
        else if(AMP_ADC>0 && atoi(NstName)>30){
          
          if (k==0) jj = 0;
          if (k==1) jj = 0;
          if (k==2) jj = 1;
          if (k==3) jj = 1;
          if (k<4 && FlagStat==true){
         
            HInteg[jj]->Fill(INT_ADC);
            HAmpl[jj]->Fill(AMP_ADC);     
            AmplCounter[jj] = AMP_ADC;
            IntegCounter[jj] = INT_ADC;
            TimeCounter[jj] = TIME_ADC + TimeEvent;

            if(k==1 || k==3)FlagStat = false;
          }
        }



      }
      Ttmp->Fill();

    }
    /*
    TFile *tfile =new TFile("tt.root","UPDATE");
    hh->Write();
    delete hh;
    tfile->Close();
    */
    Root_tmp ->cd();
    Ttmp->Write("", TObject::kOverwrite);
    delete Ttmp;
    //----------Write Hist--------------------------------------
    RootOut->cd();
    gDirectory->cd(NstName);
    gDirectory->cd("Amplitude");
    for (int j = 0; j < NUM_OF_CHANNELS; j++) {

/*
      int meanA = HAmpl[j]->GetMaximumBin();
      double xmin = double(meanA) - 2.;
      double xmax = double(meanA) + 2.; 
      FitFunc = new TF1("FitFunc", "[2]*TMath::Gaus(x,[0],[1])",xmin,xmax);
      FitFunc->SetParameters(meanA,5,1);
      FitFunc->SetParLimits(0,meanA-1,meanA+1);

      HAmpl[j]->Fit(FitFunc,"Rq");
      fprintf(FileTMP,"%s %d %f\n",NstName,j,FitFunc->GetParameter(0));
*/
      fprintf(FileTMP,"%s %d %f\n",NstName,j,double(HInteg[j]->GetMaximumBin()));

      HAmpl[j]->Write();
      HAmpl[j]->Reset();
      //delete FitFunc;
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

  fclose(FileTMP);
  RootOut->Close();
  RootData->Close("R");
  Root_tmp->Close("R");
}

int DATA::GetIA_MUON(int data[Aperture]) {
  double New[Aperture];
  int MaxRateNoise;
  int NumberMaxRateNoise = 0;
  int Noise[Aperture];
  int RateNoise[Aperture];
  int MaxAmplitude = 0;
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

    if (New[i] > MaxAmplitude) {
      MaxAmplitude = New[i];
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
  


  if (MaxAmplitude<=0 || iMaxAmplitude>990){
    INT_ADC = 0.;
    AMP_ADC = 0;
    TIME_ADC = 0.;
    return 1;
  }

  GADC = new TGraph(Aperture, XAperture, New);
  INT_ADC = GADC->Integral(0, Aperture);
  AMP_ADC = MaxAmplitude; 

  for(int i=iMaxAmplitude;i>=0;i--){
    if(New[i]==0){
      TIME_ADC = double(i)*5.E-9;
      break;
    }
  }  

/*
  double xmin = double(iMaxAmplitude) - 20.;
  double xmax = double(iMaxAmplitude) + 20.; 
  FitFunc = new TF1("FitFunc", "[2]*TMath::Gaus(x,[0],[1])",xmin,xmax);
  FitFunc->SetParameters(iMaxAmplitude,5,1);
  FitFunc->SetParLimits(0,iMaxAmplitude-4,iMaxAmplitude+4);

  char  HistName[2];
  sprintf(HistName,"%d",Nch);

  GADC->Fit("FitFunc","Rq");

  double A1 = FitFunc->GetParameter(2);
  double A2 = FitFunc->GetParameter(0);
  double A3 = FitFunc->GetParameter(1);

  double Func = A2-sqrt(-2.*A3*A3*log(MaxAplitude/3./A1));
*/


/*
  TIME_ADC = FitFunc->GetX(double(MaxAplitude)/3.,xmin,double(iMaxAmplitude));
  if (TIME_ADC!=TIME_ADC){
    TFile *tmp = new TFile("tmp.root", "UPDATE");
     GADC->Write(HistName);
    tmp->Close();
  }
*/
  //delete FitFunc;

  delete GADC;
  return 0;
}
int DATA::GetIA_GRANDE(int Nch, int data[Aperture]) {
  double New[Aperture];
  //double tt[Aperture];
  int iMaxAmplitude = 0;
  int MaxAmplitude = 0.;
  double MeanADC = 0.;
  int EndMean;

  for (int i = 0; i < Aperture; i++) {
    New[i] = 0.;
    //tt[i] = 0.;
  }


  for (int i = 400; i < 600; i++) {
    if (data[i] > MaxAmplitude) {
      MaxAmplitude = data[i];
      iMaxAmplitude = i;
    }
  }

  EndMean = iMaxAmplitude - 20;
  if (EndMean > 0) {
    int k = 0;
    for (int i = 400; i < EndMean; i++) {
      MeanADC += data[i];
      k++;
    }
    MeanADC /= k;
  }

  for (int i = 400; i < 600; i++) {
    New[i] = data[i] - MeanADC;
    //tt[i] = data[i];

    if (New[i] >= -(NOISE_RANGE) && New[i] <= (NOISE_RANGE)) {
      New[i] = 0.;
    }
  }


  MaxAmplitude = MaxAmplitude - MeanADC;
  if (MaxAmplitude >= -(NOISE_RANGE) && MaxAmplitude <= (NOISE_RANGE)){
    MaxAmplitude = 0;
  }

  int k = 0;
  bool Flag = false;
  for (int i = 400; i < 600; i++) {
    if (New[i] != 0) {
      for (int j = i; j < 600; j++) {
        if (New[j] <= 0.) {
          New[j] = 0.;
          k = j;
          Flag = true;
          break;
        }
      }
      if (Flag == true) {
        for (int j = i; j >= 400; j--) {
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


  if (MaxAmplitude<=0){
    INT_ADC = 0.;
    AMP_ADC = 0;
    TIME_ADC = 0.;
    return 1;
  }
  else if (MaxAmplitude>=2030){
    return 2;
  }


  for(int i=iMaxAmplitude;i>=400;i--){
    if(New[i]==0){
      TIME_ADC = double(i)*5.E-9;
      break;
    }
  }  


  GADC = new TGraph(Aperture, XAperture, New);
  INT_ADC = GADC->Integral(0, Aperture);
  AMP_ADC = MaxAmplitude;

/*
  TGraph *GADCt = new TGraph(Aperture, XAperture, tt);
  GADC->SetLineColor(4);
  GADC->SetLineWidth(3);


if(MaxAmplitude<10 && MaxAmplitude>0){
  TCanvas *can = new TCanvas( "can", " ", 200, 10, 600, 400 );
  char  HistName[3];
  sprintf(HistName,"%d",Nch);
  TFile *tmp = new TFile("tmp.root", "UPDATE");
  GADCt->Draw();
   GADC->Draw("same");
   can->Write(HistName);
  tmp->Close();
  delete can;
  delete GADCt;
}
*/
  delete GADC;
  return 0;
}