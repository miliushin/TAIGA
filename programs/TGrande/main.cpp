#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>

#include "grande.h"

using namespace std;

string DAYS[MAX_OF_DAYS];

string GetDay(const char *name) {
  int len = strlen(name);
  char Date[7];

  for (int i = 0; i < 30; i++) {
    if (name[len - i] == 'S') {
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
  char OutName[12];
  char NstName[4];
  OutName[11] = 0;
  NstName[3] = 0;
  sprintf(OutName, "%s.root", name);
  TFile *OutRoot = new TFile(OutName, "RECREATE");
  TTree *TData[NUMBER_OF_STATIONS];

  Int_t NubmerOfEvent;
  Double_t TimeEvent;
  Int_t ADC[NUMBER_OF_CHANNELS][Aperture];
  Double_t METEO[3];
  Int_t Ncounts = Aperture;

  for (int i = 31; i < 31 + NUMBER_OF_STATIONS; i++) {
    sprintf(NstName, "%d", i);
    TData[i - 31] = new TTree(NstName, "TAIGA-GRANDE data");
    TData[i - 31]->Branch("NumberOfEvent", &NubmerOfEvent, "NumberOfEvent/I");
    TData[i - 31]->Branch("TimeEvent", &TimeEvent, "TimeEvent/D");
    TData[i - 31]->Branch("ADC", ADC, "ADC[12][1024]/I");
    TData[i - 31]->Write();
  }

  for (int i = 0; i < NUMBER_OF_STATIONS; i++) {
    delete TData[i];
  }
  OutRoot->Close();

  delete OutRoot;
}

int main(int argc, char **argv) {
  FILE *file_data;
  string file;
  ifstream FileList;
  int iDay = 0;
  bool FlagDay;

  // FileList opening
  if (argc > 1) {
    FileList.open(argv[1]);

    printf("FileList: %s opened\n", argv[1]);
  } else {
    printf("FileList not found\n");
    exit(1);
  }

  while (getline(FileList, file)) {
    FlagDay = false;
    file_data = fopen(file.c_str(), "r");
    printf("File: %s", file.c_str());

    if (file_data != NULL) {
      printf(" opened succefully\n");
    } else {
      printf(" open ERROR\n");
      continue;
    }

    // Read file
    char name_file[255];
    strcpy(name_file, file.c_str());
    if (iDay == 0) {
      DAYS[iDay] = GetDay(name_file);
      MakeRoot(DAYS[iDay].c_str());
      iDay++;
    } else {
      for (int i = 0; i < iDay; i++) {
        if (DAYS[i] == GetDay(name_file)) {
          FlagDay = true;
          break;
        }
      }
      if (FlagDay == false) {
        DAYS[iDay] = GetDay(name_file);
        MakeRoot(DAYS[iDay].c_str());
        iDay++;
      }
    }

    GRANDE *grande = new GRANDE(name_file);
    if (FlagDay == false) grande->ReadMeteo(file);
    grande->ReadADC(file_data);
    delete grande;
    fclose(file_data);
  }
  return 0;
}