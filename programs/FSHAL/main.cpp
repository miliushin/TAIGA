#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>

#include "data.h"

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

void MakeRoot(const char *name, int Type) {
  int NumSt;
  int FirstSt;
  if (Type == 0) {
    NumSt = 3;
    FirstSt = 1;
  } else if (Type == 1) {
    NumSt = 19;
    FirstSt = 31;
  }
  char cNumSt[3];
  char channel[3];

  char OutName[18];
  OutName[17] = 0;
  sprintf(OutName, "OUTPUT/%s.root", name);
  TFile *OutRoot = new TFile(OutName, "RECREATE");

  for (int i = 0; i < NumSt; i++) {
    OutRoot->cd();
    sprintf(cNumSt, "%d", i + FirstSt);
    gDirectory->mkdir(cNumSt);
    gDirectory->cd(cNumSt);

    gDirectory->mkdir("Amplitude");
    gDirectory->cd("Amplitude");

    gDirectory->cd("../");
    gDirectory->mkdir("Integral");
    gDirectory->cd("Integral");

    gDirectory->cd("../");
    gDirectory->mkdir("DeltaT");
  }

  OutRoot->Close();
  delete OutRoot;
}

int main(int argc, char **argv) {
  FILE *file_data;
  string file;
  string day;
  ifstream FileList;
  int iDay = 0;
  int Type;
  bool FlagDay;

  // FileList opening
  if (argc > 2) {
    FileList.open(argv[1]);

    printf("FileList: %s opened\n", argv[1]);
  } else {
    printf("For RUN:  ./FSHAL FileList Type\n");
    printf("Where Type = 0 - MUON, 1 - GRANDE\n");

    exit(1);
  }
  Type = atoi(argv[2]);
  if (Type == 0 || Type == 1) {
  } else {
    printf("ERROR: Type=%d\n", Type);
    exit(1);
  }

  system("mkdir OUTPUT");

  while (getline(FileList, file)) {
    FlagDay = false;
    day = GetDay(file.c_str());

    if (iDay == 0) {
      DAYS[iDay] = day;
      MakeRoot(day.c_str(), Type);
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
        MakeRoot(day.c_str(), Type);
        iDay++;
      }
    }

    DATA *date = new DATA(file.c_str(), Type);
    delete date;
  }
  return 0;
}