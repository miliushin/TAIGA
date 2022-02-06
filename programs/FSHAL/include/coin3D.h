#ifndef COIN3D_H
#define COIN3D_H


#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TTree.h>
#include <TSystem.h>

#include "globals.h"

class Coin3D {

public:
  Coin3D(const char *File1,const char *File2);
  ~Coin3D();

private:
  TFile  *OutRoot;
  TFile  *Root_tmp;
  FILE   *FileTMP;
  FILE   *FileXYZ;
  TH1F   *HTheta;
  TH1F   *HPhi;

  double TimeCounter[NUM_OF_CHANNELS];
  double IntegCounter[NUM_OF_CHANNELS];
  //int    AmplCounter[NUM_OF_CHANNELS];

  int    Nevent[MAX_OF_STAT];
  //int    Nst3D[MAX_OF_3DSTAT][MAX_OF_3DSTAT][MAX_OF_3DSTAT];
  
  double **TimeTop;
  double **NpartTop;

  double ***TimeCoin;
  double **NpartCoin;
  int    NCoin[MAX_OF_3DSTAT];
  int    nSt1, nSt2, nSt3;

  double Coef[MAX_OF_STAT][NUM_OF_CHANNELS];

  double sX[MAX_OF_STAT];
  double sY[MAX_OF_STAT];
  double sZ[MAX_OF_STAT];

  double Phi;
  double Theta;
  double ERR_Phi;
  double ERR_Theta;

  void   FindCoin(int iComb,int St1, int St2, int St3);
  double FindThetaPhi(int iComb,int iCoin);
  int    GetNumStat(int iComb);


};

#endif
