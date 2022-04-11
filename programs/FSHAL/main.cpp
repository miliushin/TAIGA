#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>

#include "data.h"
#include "coin3D.h"
#include "axis.h"

using namespace std;

string DAYS[MAX_OF_DAYS];

string GetDay(const char *name) {
  int len = strlen(name);
  char Date[7];

  for (int i = 0; i < 30; i++) {
    if (name[len - i] == 'r' && name[len - i - 1] == '.') {
      for (int k = 0; k < 20; k++) {
        if (name[len - i - 2 - k] == '/') {
          for (int j = 0; j < 6; j++) {
            Date[j] = name[len - i - 1 - k + j];
          }
          Date[6] = 0;
          break;
        }
      }
      break;
    }
  }
  return Date;
}

void MakeRoot(const char *name) {

  char cNumSt[3];
  cNumSt[2] = 0;

  char OutName[18];
  OutName[17] = 0;

  char OutName_tmp[22];
  OutName_tmp[21] = 0;
  double Time[12];
  int    Ampl[12];
  double Integ[12];

  sprintf(OutName, "OUTPUT/%s.root", name);           // Result file
  sprintf(OutName_tmp, "OUTPUT/%s_tmp.root", name);   // temporary file

  TFile *OutRoot = new TFile(OutName, "RECREATE");
  TFile *OutRoot_tmp = new TFile(OutName_tmp, "RECREATE");

  int j = 0;
  for (int i = 1; i < MAX_OF_STAT + 1; i++) {
    if (i<4) j = i;
    else     j = i + 27;
    OutRoot->cd();
    sprintf(cNumSt, "%d", j);
    gDirectory->mkdir(cNumSt);
    gDirectory->cd(cNumSt);
    gDirectory->mkdir("Amplitude");
    gDirectory->mkdir("Integral");
    gDirectory->mkdir("DeltaT");
  }
  //OutRoot->cd();
  //gDirectory->mkdir("Angles");
  OutRoot->Close();
     // OutRoot_tmp->cd();

     // gDirectory->mkdir("pulse");

  j = 0;
  for (int i = 1; i < MAX_OF_STAT + 1; i++){
    if (i<4) j = i;
    else     j = i + 27;
    OutRoot_tmp->cd();
    sprintf(cNumSt, "%d", j);
    TTree *Ttmp = new TTree(cNumSt, "TAIGA-GRANDE/MUON data");
    Ttmp->Branch("Time", Time, "Time[12]/D");
    Ttmp->Branch("Amplitude", Ampl, "Ampl[12]/I");
    Ttmp->Branch("Integral", Integ, "Integ[12]/D");
    Ttmp->Write();
    delete Ttmp;

/*
    sprintf(cNumSt, "%d", j);
    gDirectory->mkdir(cNumSt);
    gDirectory->cd(cNumSt);

    gDirectory->mkdir("pulse");
    gDirectory->cd("pulse");
*/

  }
  OutRoot_tmp->Close();


  delete OutRoot;
  delete OutRoot_tmp;
}

int main(int argc, char **argv) {
  string file;
  string day;
  ifstream FileList;
  int iDay = 0;
  int Type;
  bool FlagDay;

  char OutName[18];
  OutName[17] = 0;
  char OutName_tmp[22];
  OutName_tmp[21] = 0;

  // FileList opening
  if (argc > 1) {
    FileList.open(argv[1]);

    printf("FileList: %s opened\n", argv[1]);
  } else {
    printf("For RUN:  ./FSHAL FileList\n");

    exit(1);
  }


  system("mkdir OUTPUT");

  while (getline(FileList, file)) {
    FlagDay = false;
    day = GetDay(file.c_str());

    if (iDay == 0) {
      DAYS[iDay] = day;
      MakeRoot(day.c_str());
      iDay++;
    } else {
      for (int i = 0; i < iDay; i++) {
        if (DAYS[i] == day) {
          FlagDay = true;
          break;
        }
      }
      if (FlagDay == false) {
        DAYS[iDay] = day;
        MakeRoot(day.c_str());
        iDay++;
      }
    }
    sprintf(OutName, "OUTPUT/%s.root", day.c_str()); 
    sprintf(OutName_tmp, "OUTPUT/%s_tmp.root", day.c_str()); 

    DATA *date = new DATA(file.c_str());
    delete date;

    Coin3D *coin = new Coin3D(OutName_tmp);
    delete coin;

    Axis *ax = new Axis(OutName,OutName_tmp);
    delete ax;

  }
  return 0;
}