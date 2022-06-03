//
//  Analysis.hpp
//  muon_rs
//
//  Created by Mikhail Iliushin on 18.05.2022.
//

#ifndef Analysis_hpp
#define Analysis_hpp

#include <stdio.h>

#include <fstream>

#include "globals.h"

class Analysis {
 public:
  Analysis(int opt);  // opt = Mode
  ~Analysis();

  int ResaveDataMUON(FILE *FileMu, FILE *fileMuNew, int Nstation);
  int ReadTimeGRANDE(FILE *File_GRANDE_tim);
  int ReadTimeHiSCORE(FILE *File_HiSCORE_tim);
  void ReadTimeMUON(FILE *FileMu, int Nstation);
  void ResetNeventMUON();
  void GetRSPATH(char *path);
  int GetMode();
  // int ReadTimeIACT(FILE *File_IACT_tim);

 private:
  int Mode;  // Mode: 0 - GRANDE, 1 - HiSCORE
  double TimeDelayMUON[NUMBER_OF_MUON_ST];
  double **TimeMUON;
  double **TimeGRANDE;
  unsigned int NeventGRANDE[NUMBER_OF_GRANDE_ST];
  unsigned int NeventMUON[NUMBER_OF_MUON_ST];
  unsigned int
      NeventMUON_NEW[NUMBER_OF_MUON_ST];  // Number of events after coin
  char *rs_path;

  double GetTimeMUON(unsigned char ch[8]);
  int FindCoin_MUON_GRANDE(
      double Time);  // coincidences between the MUON and the GRANDE st.
  int FindCoin_MUON_MUON(
      double Time,
      int Nstation);  // coincidences between a MUON and another MUON st.
};

#endif /* Analysis_hpp */
