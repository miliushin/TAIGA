#include "coin3D.h"

Coin3D::Coin3D(const char *File1,const char *File2){

  // Read the single-paticle level
  int i1, i2;
  double i3;
  OutRoot = new TFile(File1,"UPDATE");

  FileTMP = fopen("point.txt","r");
  FileXYZ = fopen("xyz.conf","r");
  HTheta = new TH1F("theta","",90,0.,90.);
  HPhi   = new TH1F("phi","",360,0.,360.);
 
  TimeTop  = new double*[MAX_OF_STAT];
  NpartTop = new double*[MAX_OF_STAT]; 

  TimeCoin  = new double**[3];
  NpartCoin = new double*[MAX_OF_3DSTAT];

  for(int i=0;i<3;i++){
    TimeCoin[i] = new double*[MAX_OF_3DSTAT];
    for(int j=0;j<MAX_OF_3DSTAT;j++){
      TimeCoin[i][j] = new double[MAX_OF_EVEN_COIN];
      for(int k=0;k<MAX_OF_EVEN_COIN;k++){
        TimeCoin[i][j][k] = 0.;
      }
    }
  }

  for(int i=0;i<MAX_OF_STAT;i++){
    Nevent[i] = 0;
    TimeTop[i] = new double[MAX_OF_EVEN];
    NpartTop[i] = new double[MAX_OF_EVEN];
    for(int j=0;j<MAX_OF_EVEN;j++){
      TimeTop[i][j] = 0.;
      NpartTop[i][j] = 0.;
    }
  }

  for(int i=0;i<MAX_OF_3DSTAT;i++){
    NCoin[i] = 0;
    NpartCoin[i] = new double[MAX_OF_EVEN_COIN];
    for(int j=0;j<MAX_OF_EVEN_COIN;j++){
      NpartCoin[i][j] = 0.;
    }
  }

  while (!feof(FileTMP)) {
    fscanf(FileTMP, "%d %d %lf\n",&i1,&i2,&i3);
    if (i1>3) i1 -= 28;
    else i1 -= 1; 
    Coef[i1][i2] = i3;
  }

  // Read station coordinates
  double xx, yy, zz;
  while (!feof(FileXYZ)) {
    fscanf(FileXYZ, "%d %lf %lf %lf\n",&i1,&xx,&yy,&zz);
    if (i1>3) i1 -= 28;
    else i1 -= 1; 
    sX[i1] = xx;
    sY[i1] = yy;
    sZ[i1] = zz;
  }

  Root_tmp = new TFile(File2,"READ");

  //----------Read tmp Tree-------------------------------------
  const char *NstName;
  int Nst;
  int NumberOfTrees = Root_tmp->GetListOfKeys()->GetSize();

  for (int i = 0; i < NumberOfTrees; i++) {
    NstName = Root_tmp->GetListOfKeys()->At(i)->GetName();
    Nst = atoi(NstName);
    if(Nst>3) Nst -= 28;
    else Nst -= 1; 
    TTree *Ttmp = (TTree*)Root_tmp->Get(NstName);
    Ttmp->SetBranchAddress("Time", TimeCounter);
    //Ttmp->SetBranchAddress("Amplitude", AmplCounter);
    Ttmp->SetBranchAddress("Integral", IntegCounter);
    Nevent[Nst] = Ttmp->GetEntries();
    printf("St:%d events:%d\n ",Nst,Nevent[Nst]);
    for (Long64_t j = 0; j < Nevent[Nst]; j++) {
      Ttmp->GetEntry(j);
      int MaxTopCh;
      // if TAIGA-GRANDE
      if(Nst>2)
        MaxTopCh = 2;
      // if TAIGA-MUON
      else
        MaxTopCh = 8;

      double MinTime = 0.;
      int iTime =0;
      for(int l=0;l<MaxTopCh;l++){
        if(TimeCounter[l]>0){
          iTime++;
          TimeTop[Nst][j] += TimeCounter[l];
          /*
          if(MinTime==0.){
            MinTime = TimeCounter[l];
          }
          else if (TimeCounter[l]<MinTime){
            MinTime = TimeCounter[l];
          }
          */
          NpartTop[Nst][j] += IntegCounter[l]/Coef[Nst][l];


        }
      }
      //TimeTop[Nst][j] = MinTime;
      TimeTop[Nst][j] /=double(iTime);

    }
  }

  // Find all 3D combinations of stations
  int nComb = 0;
  for(int i=0;i<MAX_OF_STAT-2;i++){
    for(int j=0;j<MAX_OF_STAT-2-i;j++){
      for(int k=0;k<MAX_OF_STAT-2-i-j;k++){
        int ii = i;
        int jj = j+i+1;
        int kk = k+j+i+2; 
        if(Nevent[ii]>0 && Nevent[jj]>0 && Nevent[kk]>0){
          FindCoin(nComb,ii,jj,kk);
          //printf("Ncoin: %d\n",NCoin[nComb]);

        }
        nComb++;
      }
    }
  }

  bool FlagCoin = false;
  for(int i=0;i<nComb;i++){
    int max_i = 0;
    int max_j = 0;
    for(int j=0;j<NCoin[i];j++){
      double MaxNpart = 0.;
      for(int ii=0;ii<nComb;ii++){
        for(int jj=0;jj<NCoin[ii];jj++){
          if(i==ii) continue;
          double delta = TimeCoin[0][i][j] - TimeCoin[0][ii][jj];
          if (delta<=DELTAT && delta>=-1.*DELTAT){
            FlagCoin = true;
            if(NpartCoin[i][j]>=NpartCoin[ii][jj]){
              MaxNpart = NpartCoin[i][j];
              max_i = i;
              max_j = j;
            }
            else{
              MaxNpart = NpartCoin[ii][jj];
              max_i = ii;
              max_j = jj;
            }
          }
        }
      }
      if(FlagCoin==false){
        MaxNpart = NpartCoin[i][j];
        max_i = i;
        max_j = j;
      }
      else{
        FlagCoin=false;
      }
      if(MaxNpart>0.){
        if(GetNumStat(max_i)==1){
          if(FindThetaPhi(max_i,max_j)==1){
            //printf("%f %f\n",Theta,Phi);
            HTheta->Fill(Theta);
            HPhi ->Fill(Phi);
          }
        }
        else{
          printf("ERROR: GetNumStat()\n");
        }
      }
    }
  }

}

Coin3D::~Coin3D(){
  OutRoot->cd();
  gDirectory->cd("Angles");
  HTheta->Write("theta");
  HPhi->Write("phi");
  OutRoot->Close();


  for(int i=0;i<MAX_OF_STAT;i++){
    delete [] TimeTop[i];
    delete [] NpartTop[i];
  }
  delete [] TimeTop;
  delete [] NpartTop;


  for(int i=0;i<MAX_OF_3DSTAT;i++){
    delete [] NpartCoin[i];
  }
  delete [] NpartCoin;

  for(int i=0;i<3;i++){
    for(int j=0;j<MAX_OF_3DSTAT;j++){
      delete [] TimeCoin[i][j];
    }
    delete [] TimeCoin[i];
  }
  delete [] TimeCoin;

  Root_tmp->Close();
  fclose(FileTMP);
  system("rm point.txt");

}


void Coin3D::FindCoin(int iComb,int St1, int St2, int St3) {

  double delta;

  for(int i=0;i<Nevent[St1];i++){
    if  (TimeTop[St1][i]==0.) continue;
    for(int j=0;j<Nevent[St2];j++){
      if  (TimeTop[St2][j]==0.) continue;
      delta = TimeTop[St1][i] - TimeTop[St2][j];
      if(delta >= -DELTAT && delta <= DELTAT){
        for(int k=0;k<Nevent[St3];k++){
          if  (TimeTop[St3][k]==0.) continue;
          delta = TimeTop[St1][i] - TimeTop[St3][k];
          if(delta >= -DELTAT && delta <= DELTAT){
            TimeCoin[0][iComb][NCoin[iComb]] = TimeTop[St1][i];
            TimeCoin[1][iComb][NCoin[iComb]] = TimeTop[St2][j];
            TimeCoin[2][iComb][NCoin[iComb]] = TimeTop[St3][k];
            NpartCoin[iComb][NCoin[iComb]] = NpartTop[St1][i] + NpartTop[St2][j] + NpartTop[St3][k];
            NCoin[iComb]++;
            if(NCoin[iComb]>MAX_OF_EVEN_COIN){
              printf("ERROR: Number of coin more MAX_OF_EVEN_COIN\n");
              break;
            }
           // break;
          }
        }
        //break;
      }
    }
  }

}

double Coin3D::FindThetaPhi(int iComb, int iCoin) {
  double A;
  double B;
  double C;
  double DIS;
  double d;

  double alpha, betta;
  double gamma1, gamma2, gamma;
  double Theta1, Theta2;

  double ERR_B, ERR_C, ERR_gamma, ERR_alpha, ERR_betta;

  double X13, X12, Y13, Y12, Z12, Z13;
  double t12, t13;

  double X1 = sX[nSt1];
  double X2 = sX[nSt2];
  double X3 = sX[nSt3];

  double Y1 = sY[nSt1];
  double Y2 = sY[nSt2];
  double Y3 = sY[nSt3];

  double Z1 = sZ[nSt1];
  double Z2 = sZ[nSt2];
  double Z3 = sZ[nSt3];

  X13 = X1 - X3;
  Y13 = Y1 - Y3;
  Z13 = Z1 - Z3;
  X12 = X1 - X2;
  Y12 = Y1 - Y2;
  Z12 = Z1 - Z2;

  t12 = TimeCoin[0][iComb][iCoin] - TimeCoin[1][iComb][iCoin];
  t13 = TimeCoin[0][iComb][iCoin] - TimeCoin[2][iComb][iCoin];


  d = X13 * Y12 - X12 * Y13;

  A = pow(X13 * Z12 - X12 * Z13, 2.) + pow(Y12 * Z13 - Y13 * Z12, 2.) + d * d;
  B = -2. * VC * (X13 * t12 - X12 * t13) * (X13 * Z12 - X12 * Z13) -
      2. * VC * (Y12 * t13 - Y13 * t12) * (Y12 * Z13 - Y13 * Z12);
  C = pow(VC * (X13 * t12 - X12 * t13), 2.) +
      pow(VC * (Y12 * t13 - Y13 * t12), 2.) - d * d;
  DIS = B * B - 4. * A * C;
  if (DIS < 0) return 0;

  gamma1 = (-B + sqrt(DIS)) / 2. / A;
  gamma2 = (-B - sqrt(DIS)) / 2. / A;

  Theta1 = acos(gamma1) / Pi * 180.;
  Theta2 = acos(gamma2) / Pi * 180.;

  gamma = (-B + sqrt(DIS)) / 2. / A;

  alpha = (1. / d) *
          ((VC * t13 - gamma * Z13) * Y12 - (VC * t12 - gamma * Z12) * Y13);
  betta = (1. / d) *
          ((VC * t12 - gamma * Z12) * X13 - (VC * t13 - gamma * Z13) * X12);

  Phi = atan2(betta, alpha) / Pi * 180.;
  Theta = acos(gamma) / Pi * 180.;
  if (Phi < 0.) Phi = Phi + 360.; 
  // Find errors

  ERR_B = ERR_T * 2. * VC *
          sqrt(pow(X13 * Z12 - X12 * Z13, 2.) * (X13 * X13 + X12 * X12) +
               pow(Y12 * Z13 - Y13 * Z12, 2.) * (Y12 * Y12 + Y13 * Y13));
  ERR_C = 2. * VC * VC * ERR_T *
          sqrt(pow(X13 * t12 - X12 * t13, 2.) * (X13 * X13 + X12 * X12) +
               pow(Y12 * t13 - Y13 * t12, 2.) * (Y12 * Y12 + Y13 * Y13));
  ERR_gamma =
      (1. / 2. / sqrt(DIS) / A) *
      sqrt((2. * B * B - 4. * A * C) * ERR_B * ERR_B + pow(2. * A * ERR_C, 2.));
  ERR_Theta = ERR_gamma / sqrt(1. - gamma * gamma) / Pi * 180.;

  ERR_alpha = sqrt(pow(ERR_T * VC, 2.) * (Y12 * Y12 + Y13 * Y13) +
                   pow(ERR_gamma * (Y12 * Z13 - Y13 * Z12), 2.)) /
              d;
  ERR_betta = sqrt(pow(ERR_T * VC, 2.) * (X13 * X13 + X12 * X12) +
                   pow(ERR_gamma * (X13 * Z12 - X12 * Z13), 2.)) /
              d;
  ERR_Phi = sqrt(pow(ERR_betta / alpha, 2.) +
                 pow(betta * ERR_alpha / alpha / alpha, 2.)) /
            (1. + pow(betta / alpha, 2.)) / Pi * 180.;

  // printf("B: %f %f \n",B,ERR_B);
  // printf("C: %f %f \n",C,ERR_C);
  // printf("gamma: %f %f \n",gamma,ERR_gamma);
  return 1;
}


int Coin3D::GetNumStat(int iCoin){

  int nComb = 0;
  for(int i=0;i<MAX_OF_STAT-2;i++){
    for(int j=0;j<MAX_OF_STAT-2-i;j++){
      for(int k=0;k<MAX_OF_STAT-2-i-j;k++){
        int ii = i;
        int jj = j+i+1;
        int kk = k+j+i+2;
        if (iCoin == nComb){
          nSt1 = ii;
          nSt2 = jj;
          nSt3 = kk;
          return 1;
        }
        nComb++;
      }
    }
  }

  return 0;
}