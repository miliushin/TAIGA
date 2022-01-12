#include "muon.h"

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TROOT.h>

#include <cstdlib>

using namespace std;

MUON::MUON(char* name_file) {
  int len = strlen(name_file);
  char Nst[4];
  char Date[7];

  Nst[0] = name_file[len - 7];
  Nst[1] = name_file[len - 6];
  Nst[2] = name_file[len - 5];
  Nst[3] = 0;
  char NstName[4];
  sprintf(NstName, "%d", atoi(Nst));

  for (int i = 0; i < 20; i++) {
    if (name_file[len - i] == 'M') {
      for (int j = 0; j < 6; j++) {
        Date[j] = name_file[len - i - 1 - 6 + j];
      }
      Date[6] = 0;

      break;
    }
  }

  char OutName[12];
  OutName[11] = 0;
  sprintf(OutName, "%s.root", Date);

  int NCounts = Aperture;
  RootFile = new TFile(OutName, "UPDATE");
  TData = (TTree*)RootFile->Get(NstName);
  TData->SetBranchAddress("NumberOfEvent", &NumEvent);
  TData->SetBranchAddress("TimeEvent", &TimeEvent);
  TData->SetBranchAddress("ADC", ADC);
  // TData->Branch("METEO",ADC,"ADC/D");
}

MUON::~MUON() {
  TData->Write("", TObject::kOverwrite);
  delete TData;
  RootFile->Close();
  delete RootFile;
}

double MUON::GetTime(unsigned char ctt[8]) {
  int data_time[4];
  short unsigned int h, m, s, mls, mks, dns;

  for (int i = 0; i < 4; i++) {
    data_time[i / 2] = ctt[i * 2 + 1] * 256 + ctt[i * 2];
  }

  data_time[0] = ctt[1] * 256 + ctt[0];
  data_time[1] = ctt[3] * 256 + ctt[2];
  data_time[2] = ctt[5] * 256 + ctt[4];
  data_time[3] = ctt[6] * 256 + ctt[6];

  dns = (data_time[0] & 0x7f) * 10;
  mks = (data_time[0] & 0xff80) >> 7;
  mks |= (data_time[1] & 1) << 9;
  mls = (data_time[1] & 0x7fe) >> 1;
  s = (data_time[1] & 0xf800) >> 11;
  s |= (data_time[2] & 1) << 5;
  m = (data_time[2] & 0x7e) >> 1;
  h = (data_time[2] & 0xf80) >> 7;

  double time = (double)dns / 1000000000.0 + (double)mks / 1000000.0 +
                (double)mls / 1000.0 + (double)s + (double)m * 60.0 +
                (double)h * 3600.0;

  return (time);
}

int MUON::ReadADC(FILE* data_file) {
  unsigned char buf[SIZE_OF_HEADER];
  unsigned char data_read[0x8000];
  int dat = 0;

  fseek(data_file, 0, SEEK_END);
  END_File = ftell(data_file);
  i_File = 0;
  fseek(data_file, 0, SEEK_SET);
  int Nevent;
  while (i_File < END_File) {
    // Read header of file
    if (!(fread(buf, SIZE_OF_HEADER, 1, data_file))) {
      printf("ERROR: read file ");
      return -1;
    }
    Nevent = 0;

    ID = buf[1] * 256 + buf[0];
    NumBytes = buf[3] * 256 + buf[2];

    if (NumBytes > 0x8000) {
      printf("ERROR:  Too much bytes = 0x%X\n", NumBytes);
    }

    NumEvent =
        buf[7] * 256 * 256 * 256 + buf[6] * 256 * 256 + buf[5] * 256 + buf[4];

    TimeEvent = GetTime(buf + 12);

    fread(data_read, NumBytes, 1, data_file);

    for (int ich = 0; ich < NUMBER_OF_CHANNELS; ich++) {
      for (int ibin = 0; ibin < Aperture; ibin++) {
        data[ich][ibin] = 0;
      }
    }

    int dib = 0;
    int dch = 0;
    int dch_add = 0;
    int fbin = 0;

    for (int ibin = 0; ibin < NumBytes - 4; ibin += 2) {
      dat = data_read[ibin + 1] * 256 + data_read[ibin];
      int ddd = data_read[ibin + 12] * 256 + data_read[ibin + 11];
      dch_add = (dat >> 12);

      if (dch_add != dch) {
        dib = 0;
        dch = dch_add;
      }

      data[dch][dib] = dat & 0xFFF;
      data[dch][dib] -= 2048;

      fbin = ibin + 2;
      dib++;
    }

    fin = data_read[fbin + 3] * 256 * 256 * 256 +
          data_read[fbin + 2] * 256 * 256 + data_read[fbin + 1] * 256 +
          data_read[fbin + 0];

    if (fin != 0xFFFFFFFF) {
      printf("ERROR::  uncorrect finich of package\n");
      return -1;
    }

    for (int i = 0; i < NUMBER_OF_CHANNELS; i++) {
      for (int j = 0; j < Aperture; j++) {
        ADC[i][j] = data[i][j];
      }
    }

    TData->Fill();

    i_File = ftell(data_file);
  }
  return 1;
}
