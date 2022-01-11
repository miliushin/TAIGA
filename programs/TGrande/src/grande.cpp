#include "grande.h"

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TROOT.h>

#include <cstdlib>

using namespace std;

GRANDE::GRANDE(char* name_file) {
  int len = strlen(name_file);
  char Nst[4];
  char Date[7];

  Nst[0] = name_file[len - 7];
  Nst[1] = name_file[len - 6];
  Nst[2] = name_file[len - 5];
  Nst[3] = 0;
  char NstName[4];
  sprintf(NstName, "%d", atoi(Nst));


  for (int i = 0; i < 30; i++) {
    if (name_file[len - i] == 'S') {
      for (int k = 0; k < 20; k++) {
        if (name_file[len - i - 2 - k] == '/'){
          for (int j = 0; j < 6; j++) {
            Date[j] = name_file[len - i - 1 - k + j];
           }
           Date[6] = 0;
          break;
        }
      }
      break;
    }
  }

  char OutName[12];
  OutName[11] = 0;
  sprintf(OutName, "%s.root", Date);

  RootFile = new TFile(OutName, "UPDATE");
  TData = (TTree*)RootFile->Get(NstName);
  TData->SetBranchAddress("NumberOfEvent", &NumEvent);
  TData->SetBranchAddress("TimeEvent", &TimeEvent);
  TData->SetBranchAddress("ADC", ADC);
  // TData->Branch("METEO",ADC,"ADC/D");
}

GRANDE::~GRANDE() {
  TData->Write("", TObject::kOverwrite);
  delete TData;
  RootFile->Close();
  delete RootFile;
}

double GRANDE::GetTime(unsigned char ctt[8]) {
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


int GRANDE::ReadADC(FILE* data_file) {
  unsigned char buf[SIZE_OF_HEADER];
  unsigned char data_read[20000];
  int dat = 0;

  fseek(data_file, 0, SEEK_END);
  END_File = ftell(data_file);
  i_File = 0;
  fseek(data_file, 0, SEEK_SET);
  int Nevent;
  while (i_File < END_File) {
    //--------------------------------------- Read header of file
    if (!(fread(buf, SIZE_OF_HEADER, 1, data_file))) {
      printf("ERROR: read file ");
      return -1;
    }
    Nevent = 0;

    ID = buf[1] * 256 + buf[0];

    switch (ID) {
      case 3032:
        Claster = buf[20];
        break;
      case 3033:
        Claster = buf[20];
        break;
      case 3034:
        Claster = buf[20];
        break;
      case 3035:
        Claster = buf[20];
        break;
      default:
        printf("WROND ID: %i\n", ID);
        return -1;
    }

    NumBytes = buf[3] * 256 + buf[2];
    NumEvent =
        buf[7] * 256 * 256 * 256 + buf[6] * 256 * 256 + buf[5] * 256 + buf[4];

    TimeEvent = GetTime(buf + 12);
    Delay = (buf[23] * 256 + buf[22]) * 5;
    TimeEvent += (double)Delay / 1000000000.0;


    for (int ich = 0; ich < NUMBER_OF_CHANNELS; ich++) {
      for (int ibin = 0; ibin < Aperture; ibin++) {
        ADC[ich][ibin] = 0;
      }
    }

    //-------------------------------------- END Read header of file

    int NumPoints = 0;
    int index;
    int j = 0;
    int ch = 0;
    unsigned char data[20000];
    //---------------------Read Data---------------------------------
    fread(data_read, NumBytes, 1, data_file);
    for (int i = 0; i < 8; i++) {
      Addr = (unsigned long)data_read[j + 3] * 256L * 256L * 256L +
             (unsigned long)data_read[j + 2] * 256L * 256L +
             (unsigned long)data_read[j + 1] * 256L + (unsigned long)data_read[j + 0];

      Addr &= 0x7fffefff;
      NumPoints = (unsigned int)data_read[j + 5] * 256 + data_read[j + 4];



      ch = (unsigned int)data_read[j + 6] & 0xf;
      if (ch==12) break;                                 // only 12 channels of GRANDE now
      int ipst = (Addr - AddrBase[ch])/2;
      int ipfn = ipst + NumPoints;
      if(ipst<0 || ipfn>1024){
        printf("ERROR:   Addr=0x%lx,  NumPoints=%i   ipst=%i ipfn=%i \n",Addr,NumPoints,ipst,ipfn);
        ipst = 400;
      } 

      for (int k = 0; k < NumPoints * 2; k += 2) {
        index = j + k + 8;
        ADC[ch][ipst + k / 2] =
            (data_read[index + 1] * 256 + data_read[index]) - 2048;
      }

      j += (NumPoints * 2 + 8);

    }
    if(ch == 11) TData->Fill();
    i_File = ftell(data_file);
  }

  return 1;
}
