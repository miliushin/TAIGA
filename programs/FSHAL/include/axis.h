#ifndef AXIS_H
#define AXIS_H


#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TSystem.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "globals.h"


class Axis {

public:
  Axis(const char *File1, const char *File2);
  ~Axis();

private:
  TFile  *ResultRoot;
  TFile  *TmpRoot;


  TTree  *ResultTree;
  TTree  *TmpTree;

  TF2    *LDF1;
  TF2    *LDF2;
  TF2    *LDF3;
  TF1    *LDF4;
  TF1    *FDensityE_Teor;

  FILE   *FileXYZ;

  double sX[MAX_OF_STAT];
  double sY[MAX_OF_STAT];
  double sZ[MAX_OF_STAT];

  double time;
  double Phi;
  double Theta;
  double ERR_Phi;
  double ERR_Theta;
  double Nall_Top;
  double rho_200;
  double S;
  double NparticalTop[MAX_OF_STAT];
  int    NDet[MAX_OF_STAT];
  double x0, y0;
  double err_x0, err_y0;
  double *Ne;
  double *err_Ne;
};
#endif