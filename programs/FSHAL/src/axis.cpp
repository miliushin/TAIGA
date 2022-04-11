
#include "axis.h"

#include <TGraph.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TGraph2D.h>
#include <TH2F.h>
Double_t fLDF(Double_t *x, Double_t *par) {
  Double_t CS = TMath::Gamma(4.5-par[0])/(2.*3.1415926535*par[1]*par[1]*TMath::Gamma(par[0])*TMath::Gamma(4.5-2.*par[0]));
  return CS*par[2]*pow(sqrt(pow(x[0]-par[3],2.)+pow(x[1]-par[4],2.))/par[1],par[0]-2.)*
   pow(1.+sqrt(pow(x[0]-par[3],2.)+pow(x[1]-par[4],2.))/par[1],par[0]-4.5);
} 

Double_t fLDF1(Double_t *x, Double_t *par) {
  Double_t CS = TMath::Gamma(4.5-par[0])/(2.*Pi*par[1]*par[1]*TMath::Gamma(par[0])*TMath::Gamma(4.5-2.*par[0]));
  return log10(CS*par[2]*pow(x[0]/par[1],par[0]-2.)*
    pow(1.+x[0]/par[1],par[0]-4.5));
} 


Axis::Axis(const char *File1, const char *File2){

  TH1F *hist = new TH1F("hist", "", 500, 0.0, 500);

  int ii;

  //LDF = new TF2("LDF","tgamma(4.5-[0])/tgamma([0])/tgamma(4.5-2.*[0])/2./3.14159265358979323846/[1]/[1] \
*[2]*pow(sqrt(pow(x-[3],2.)+pow(y-[4],2.))  \
    /[1],[0]-2.)*pow(1.+sqrt(pow(x-[3],2.)+pow(y-[4],2.))/[1],[0]-4.5)",-1000.,1000.,-1000.,1000.);

  gStyle->SetOptStat(111111);
  
  LDF1 = new TF2("LDF1",fLDF,-1000.,1000.,-1000.,1000.,5);
  LDF2 = new TF2("LDF2",fLDF,-1000.,1000.,-1000.,1000.,5);
  LDF3 = new TF2("LDF3",fLDF,-1000.,1000.,-1000.,1000.,5);
  LDF4 = new TF2("LDF4",fLDF,-1000.,1000.,-1000.,1000.,5);

  TH2F *Hxy = new TH2F("hxy","",100,-1000,1000,100,-1000,1000);

//LDF3 = new TF2("LDF3","TMath::Gamma(4.5-[0])/(2.*3.1415926535*[1]*[1]*TMath::Gamma([0])*TMath::Gamma(4.5-2.*[0]))* \
  pow(sqrt(pow(x-[3],2.)+pow(y-[4],2.))/[1],[0]-2.)* \
  pow(1.+sqrt(pow(x-[3],2.)+pow(y-[4],2.))/[1],[0]-4.5)",-1000.,1000.,-1000.,1000.);

/*
  LDF1 = new TF1("LDF1",fLDF1,-1000.,1000.,3);
  LDF2 = new TF1("LDF2",fLDF1,-1000.,1000.,3);
  LDF3 = new TF1("LDF3",fLDF1,-1000.,1000.,3);
*/
  // Read station coordinates
  FileXYZ = fopen("xyz.conf","r");
  double xx, yy, zz;
  while (!feof(FileXYZ)) {
    fscanf(FileXYZ, "%d %lf %lf %lf\n",&ii,&xx,&yy,&zz);
    if (ii>3) ii -= 28;
    else ii -= 1; 
    sX[ii] = xx;
    sY[ii] = yy;
    sZ[ii] = zz;
  }
  fclose(FileXYZ);



  // Read tmp Tree
  TmpRoot = new TFile(File2,"UPDATE");
  TmpTree = (TTree*)TmpRoot->Get("result");
  TmpTree->SetBranchAddress("time",&time);
  TmpTree->SetBranchAddress("theta", &Theta);
  TmpTree->SetBranchAddress("phi", &Phi);
  TmpTree->SetBranchAddress("err_theta", &ERR_Theta);
  TmpTree->SetBranchAddress("err_phi", &ERR_Phi);
  TmpTree->SetBranchAddress("NparticalTop", NparticalTop);
  TmpTree->SetBranchAddress("NDet", NDet);

  // Make result tree for data
  ResultRoot = new TFile(File1,"UPDATE");
  ResultTree = new TTree("result","Result Tree");
  ResultTree->Branch("time",&time, "time/D");
  ResultTree->Branch("theta", &Theta, "Theta/D");
  ResultTree->Branch("err_theta", &ERR_Theta, "ERR_Theta/D");
  ResultTree->Branch("phi", &Phi, "Phi/D");
  ResultTree->Branch("err_phi", &ERR_Phi, "ERR_Phi/D");
  ResultTree->Branch("x0", &x0, "x0/D");
  ResultTree->Branch("err_x0", &err_x0, "err_x0/D");
  ResultTree->Branch("y0", &y0, "y0/D");
  ResultTree->Branch("err_y0", &err_y0, "err_y0/D");
  ResultTree->Branch("NparticalTop", &Nall_Top, "Nall_Top/D");
  ResultTree->Branch("rho_200", &rho_200, "rho_200/D");
  ResultTree->Branch("S",&S,"S/D");

//  ResultTree->Branch("NDet", NDet, "NDet[22]/I");

  double SDensity;
  double SDet;

  Long64_t nentries = TmpTree->GetEntries();
  Ne = new double[nentries];
  err_Ne = new double[nentries];


  TFile *FF = new TFile("errX0.root","recreate");

  double ni;
  double ri, err_ri;
  double fi, err_fi;
  double Sfi, err_Sfi;
  double Rm = 80.;
  int    MaxDet;
  int    ERR_CHI = 0;
  for (Long64_t i = 0; i < nentries; i++) {
    S = 1.1;
    MaxDet = 0;
  	ni = 0.;
    Sfi = 0.;
    err_Sfi = 0.;
    ri = 0.;
    fi = 0.;
    err_fi = 0.;
    err_ri = 0.;
	  SDensity = 0.;
	  x0 = 0.;
  	y0 = 0.;
    err_x0 = 0.;
    err_y0 = 0.;


    double err_x1 = 0.;
    double err_x2 = 0.;
    double err_y1 = 0.;
    double err_y2 = 0.;



    TmpTree->GetEntry(i);

  	for(int j=0;j<MAX_OF_STAT;j++){
  	  if (NDet[j]==0) continue;
  	  if (j>2) SDet = SQUARE_GRANDE*6.;
  	  else               SDet = SQUARE_MUON;

      x0 += NparticalTop[j]*sX[j]/(SDet*NDet[j]);
      y0 += NparticalTop[j]*sY[j]/(SDet*NDet[j]);
      SDensity += NparticalTop[j]/(SDet*NDet[j]);  

      ni +=  NparticalTop[j];
      MaxDet++;
    }

    x0 /=SDensity;
    y0 /=SDensity;



    int iNe = 0;
    for(int j=0;j<MAX_OF_STAT;j++){
      if (NDet[j]==0) continue;
      if (j>2) SDet = SQUARE_GRANDE*6.;
      else               SDet = SQUARE_MUON;
  
      err_x1 += pow(sX[j]/SDet/NDet[j],2.)*NparticalTop[j];
      err_y1 += pow(sY[j]/SDet/NDet[j],2.)*NparticalTop[j];

      err_x2 += NparticalTop[j]/pow(SDet*NDet[j],2.);
      err_y2 += NparticalTop[j]/pow(SDet*NDet[j],2.);
    }

    err_x1 = sqrt(err_x1);
    err_y1 = sqrt(err_y1);
    err_x2 = sqrt(err_x2);
    err_y2 = sqrt(err_y2);

    err_x0 = sqrt(pow(err_x1/SDensity,2.)+pow(err_x2*x0/SDensity,2.));
    err_y0 = sqrt(pow(err_y1/SDensity,2.)+pow(err_y2*y0/SDensity,2.));


    double* rNe = new double[MaxDet];
    double* err_rNe = new double[MaxDet];

    double* xNe = new double[MaxDet];
    double* yNe = new double[MaxDet];

    double* DensityE = new double[MaxDet];
    double* ERRDensityE = new double[MaxDet];

    // 1st level
    iNe = 0;
  	for(int j=0;j<MAX_OF_STAT;j++){

  	  if (NDet[j]==0) continue;
  	  if (j>2) SDet = SQUARE_GRANDE*6.;
  	  else               SDet = SQUARE_MUON;

      ri = sqrt(pow(sX[j]-x0,2.)+pow(sY[j]-y0,2.));
      err_ri = 1./ri*sqrt(pow((sX[j]-x0)*err_x0,2.) + pow((sY[j]-y0)*err_y0,2.));
      fi = tgamma(4.5-S)/tgamma(S)/tgamma(4.5-2.*S)/(2.*Pi*Rm*Rm)*pow(ri/Rm,S-2.)*pow(1.+ri/Rm,S-4.5);
      err_fi = tgamma(4.5-S)/tgamma(S)/tgamma(4.5-2.*S)/(2.*Pi*Rm*Rm)*err_ri/Rm*
               sqrt(pow(S-2.,2.)*pow(ri/Rm,2.*S-6.)*pow(1.+ri/Rm,2.*S-9.) +
               pow(S-4.5,2.)*pow(ri/Rm,2.*S-4.)*pow(1.+ri/Rm,2.*S-11.));

      Sfi += fi*SDet*NDet[j]*cos(Theta*Pi/180.);
      err_Sfi += pow(err_fi*SDet*NDet[j]*cos(Theta*Pi/180.),2.) + pow(fi*SDet*NDet[j]*sin(Theta*Pi/180.)*
        ERR_Theta*Pi/180.,2.);
      xNe[iNe] = sX[j];
      yNe[iNe] = sY[j];
      rNe[iNe] = ri;
      err_rNe[iNe] = err_ri;
      DensityE[iNe] = NparticalTop[j]/(SDet*NDet[j]);
      ERRDensityE[iNe] = sqrt(NparticalTop[j]*pow(1./(SDet*NDet[j]),2.));

      iNe++;

  	}
    if (Sfi!=0) {
      Ne[i] = ni/Sfi;
      err_Ne[i] = sqrt(ni/pow(Sfi,2.) + err_Sfi*pow(ni/Sfi/Sfi,2.));
    }
    else {
      Ne[i] = 0.;
      err_Ne[i] =0.;
    }
    //printf("Ne: %f Err: %f\n",Ne[i],err_Ne[i]);
 
    // 1st-fit
    LDF1->FixParameter(0,S);
    LDF1->FixParameter(1,Rm);
    LDF1->FixParameter(2,Ne[i]);
    LDF1->SetParameter(3,x0);    //LDF1->SetParLimits(3,x0-err_x0,x0+err_x0);
    LDF1->SetParameter(4,y0);    //LDF1->SetParLimits(4,y0-err_y0,y0+err_y0);

    TGraph2DErrors *GNe1 = new TGraph2DErrors(MaxDet,xNe,yNe,DensityE,0,0,ERRDensityE);

    hist->Fill(err_x0);


    if(i==1)GNe1->Fit("LDF1");
    else GNe1->Fit("LDF1","q");
    x0 = LDF1->GetParameter(3); err_x0 = LDF1->GetParError(3);
    y0 = LDF1->GetParameter(4); err_y0 = LDF1->GetParError(4);

    // 2 level
    iNe = 0;
    Sfi = 0.;
    err_Sfi = 0.;
    for(int j=0;j<MAX_OF_STAT;j++){

      if (NDet[j]==0) continue;
      if (j>2) SDet = SQUARE_GRANDE*6.;
      else               SDet = SQUARE_MUON;

      ri = sqrt(pow(sX[j]-x0,2.)+pow(sY[j]-y0,2.));
      err_ri = 1./ri*sqrt(pow((sX[j]-x0)*err_x0,2.) + pow((sY[j]-y0)*err_y0,2.));


      fi = tgamma(4.5-S)/tgamma(S)/tgamma(4.5-2.*S)/(2.*Pi*Rm*Rm)*pow(ri/Rm,S-2.)*pow(1.+ri/Rm,S-4.5);

      err_fi = tgamma(4.5-S)/tgamma(S)/tgamma(4.5-2.*S)/(2.*Pi*Rm*Rm)*err_ri/Rm*
               sqrt(pow(S-2.,2.)*pow(ri/Rm,2.*S-6.)*pow(1.+ri/Rm,2.*S-9.) +
               pow(S-4.5,2.)*pow(ri/Rm,2.*S-4.)*pow(1.+ri/Rm,2.*S-11.));


      Sfi += fi*SDet*NDet[j]*cos(Theta*Pi/180.);
      err_Sfi += pow(err_fi*SDet*NDet[j]*cos(Theta*Pi/180.),2.) + pow(fi*SDet*NDet[j]*sin(Theta*Pi/180.)*
        ERR_Theta*Pi/180.,2.);

      iNe++;
    }  
    if (Sfi!=0){
      Ne[i] = ni/Sfi;
      err_Ne[i] = sqrt(ni/pow(Sfi,2.) + err_Sfi*pow(ni/Sfi/Sfi,2.));

    }
    else Ne[i] = 0.;

    //printf("Ne2: %f %f\n",Ne[i],err_Ne[i]);


    // 2d-fit
    LDF3->SetParameter(0,S);      LDF3->SetParLimits(0,0.2,3.);
    LDF3->FixParameter(1,Rm);
    LDF3->FixParameter(2,Ne[i]);
    LDF3->SetParameter(3,x0);     //LDF3->SetParLimits(3,x0-err_x0,x0+err_x0);
    LDF3->SetParameter(4,y0);     //LDF3->SetParLimits(4,y0-err_y0,y0+err_y0);
    
    if(i==1)GNe1->Fit("LDF3");
    else GNe1->Fit("LDF3","q");
   
    S  = LDF3->GetParameter(0);
    x0 = LDF3->GetParameter(3);  err_x0 = LDF3->GetParError(3);
    y0 = LDF3->GetParameter(4);  err_y0 = LDF3->GetParError(4);


    iNe = 0;
    Sfi = 0.;
    err_Sfi = 0.;
    for(int j=0;j<MAX_OF_STAT;j++){

      if (NDet[j]==0) continue;
      if (j>2) SDet = SQUARE_GRANDE*6.;
      else               SDet = SQUARE_MUON;

      ri = sqrt(pow(sX[j]-x0,2.)+pow(sY[j]-y0,2.));
      err_ri = 1./ri*sqrt(pow((sX[j]-x0)*err_x0,2.) + pow((sY[j]-y0)*err_y0,2.));


      fi = tgamma(4.5-S)/tgamma(S)/tgamma(4.5-2.*S)/(2.*Pi*Rm*Rm)*pow(ri/Rm,S-2.)*pow(1.+ri/Rm,S-4.5);

      err_fi = tgamma(4.5-S)/tgamma(S)/tgamma(4.5-2.*S)/(2.*Pi*Rm*Rm)*err_ri/Rm*
               sqrt(pow(S-2.,2.)*pow(ri/Rm,2.*S-6.)*pow(1.+ri/Rm,2.*S-9.) +
               pow(S-4.5,2.)*pow(ri/Rm,2.*S-4.)*pow(1.+ri/Rm,2.*S-11.));


      Sfi += fi*SDet*NDet[j]*cos(Theta*Pi/180.);
      err_Sfi += pow(err_fi*SDet*NDet[j]*cos(Theta*Pi/180.),2.) + pow(fi*SDet*NDet[j]*sin(Theta*Pi/180.)*
        ERR_Theta*Pi/180.,2.);

      iNe++;
    }  
    if (Sfi!=0){
      Ne[i] = ni/Sfi;
      err_Ne[i] = sqrt(ni/pow(Sfi,2.) + err_Sfi*pow(ni/Sfi/Sfi,2.));

    }
    else Ne[i] = 0.;


    // 2d-fit
    LDF4->SetParameter(0,S);      LDF3->SetParLimits(0,0.2,3.);
    LDF4->FixParameter(1,Rm);
    LDF4->SetParameter(2,Ne[i]);
    LDF4->SetParameter(3,x0);     //LDF3->SetParLimits(3,x0-err_x0,x0+err_x0);
    LDF4->SetParameter(4,y0);     //LDF3->SetParLimits(4,y0-err_y0,y0+err_y0);
 
    GNe1->Fit("LDF4","q");

    S     = LDF4->GetParameter(0);
    Ne[i] = LDF4->GetParameter(2);
    x0    = LDF4->GetParameter(3);  //err_x0 = LDF4->GetParError(3);
    y0    = LDF4->GetParameter(4);  //err_y0 = LDF4->GetParError(4);

    // 3 level
    iNe = 0;
    Sfi = 0.;

    double *rr = new double[MaxDet];
    double *dd = new double[MaxDet];
    for(int j=0;j<MAX_OF_STAT;j++){

      if (NDet[j]==0) continue;
      if (j>2) SDet = SQUARE_GRANDE*6.;
      else               SDet = SQUARE_MUON;

      ri = sqrt(pow(sX[j]-x0,2.)+pow(sY[j]-y0,2.));
      err_ri = 1./ri*sqrt(pow((sX[j]-x0)*err_x0,2.) + pow((sY[j]-y0)*err_y0,2.));

      fi = tgamma(4.5-S)/tgamma(S)/tgamma(4.5-2.*S)/(2.*Pi*Rm*Rm)*pow(ri/Rm,S-2.)*pow(1.+ri/Rm,S-4.5);

      Sfi += fi*SDet*NDet[j]*cos(Theta*Pi/180.);

      rr[iNe] = ri;
      dd[iNe] = NparticalTop[j]/(SDet*NDet[j]);

      rNe[iNe] = log10(ri);

      err_rNe[iNe] = 1./ri/log(10.)*err_ri;
      DensityE[iNe] = log10(NparticalTop[j]/(SDet*NDet[j]));
      ERRDensityE[iNe] = sqrt(pow(1./sqrt(NparticalTop[j])/log(10.),2.));
      iNe++;
    }  



    //if (Sfi!=0) Ne[i] = ni/Sfi;
    //else Ne[i] = 0.;

    TCanvas *can = new TCanvas( "can", " ", 200, 10, 600, 400 );   
    TGraphErrors *GNe = new TGraphErrors(MaxDet,rNe,DensityE,err_rNe,ERRDensityE);
    

    double chi =0.;
    for(int j=0;j<MaxDet;j++){
      chi += pow((dd[j]-(TMath::Gamma(4.5-S)/(2.*Pi*Rm*Rm*
        TMath::Gamma(S)*TMath::Gamma(4.5-2.*S))*Ne[i]*pow(rr[j]/Rm,S-2.)*pow(1.+rr[j]/Rm,S-4.5))),2.)/
      (TMath::Gamma(4.5-S)/(2.*Pi*Rm*Rm*
        TMath::Gamma(S)*TMath::Gamma(4.5-2.*S))*Ne[i]*pow(rr[j]/Rm,S-2.)*pow(1.+rr[j]/Rm,S-4.5)); 
    }
    chi /=MaxDet;
    if (chi>1.) ERR_CHI++;
    else{
      Nall_Top = Ne[i];
      rho_200 = TMath::Gamma(4.5-S)/(2.*Pi*Rm*Rm*
        TMath::Gamma(S)*TMath::Gamma(4.5-2.*S))*Ne[i]*pow(200./Rm,S-2.)*pow(1.+200./Rm,S-4.5);
  
     // printf("%.9f %f %f %f %f\n",time,Theta,Phi,x0,y0);

      ResultTree->Fill();
      Hxy->Fill(x0,y0);
    }

   
    
    int nGR = 100;
    double YGR[nGR];
    double XGR[nGR];

    for (int iGR=1;iGR<nGR;iGR++){
      XGR[iGR] = iGR*10.;
      YGR[iGR] = log10(TMath::Gamma(4.5-S)/(2.*Pi*Rm*Rm*
        TMath::Gamma(S)*TMath::Gamma(4.5-2.*S))*Ne[i]*pow(XGR[iGR]/Rm,S-2.)*pow(1.+XGR[iGR]/Rm,S-4.5));
      XGR[iGR] = log10(XGR[iGR]);
    }

    TGraph *DensityGR = new TGraph(nGR,XGR,YGR);


    //printf("%f %f %f\n",S,Rm,Ne[i]);
    GNe->Draw("AP*");
    DensityGR->Draw("C");
    can->Write("Ne");
    
    
    delete DensityGR;
    delete GNe1;
    delete GNe;
    delete can;
    delete [] rNe;
    delete [] err_rNe;
    delete [] xNe;
    delete [] yNe;
    delete [] DensityE;
    delete [] ERRDensityE;
    delete [] rr;
    delete [] dd;
  }
  printf("%d %d\n",nentries, ERR_CHI);
  hist->Write();
  Hxy->Write("xy");
  FF->Close();


}

Axis::~Axis(){
  delete [] Ne;
  delete [] err_Ne;
  TmpRoot->Close();

  ResultRoot->cd();
  ResultTree->Write();
  ResultRoot->Close();

}