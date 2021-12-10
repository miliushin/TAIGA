#include <fstream>
#include <cstdio>
#include <iostream>
#include <cstring>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH1.h>

#include <TFile.h>

using namespace std;


int main(int argc, char ** argv) {

	double delT = 0.001;
	ifstream FileMuon;
	ifstream FileHiscore;

	string timeH;
	string timeM;

	int N = 10000000;
	double *TimeMuon = new double[N];
	double *TimeHiscore = new double[N];
	double *delta_time = new double[N];

	double StartTime;
	double EndTime;
	double weight;
	bool FirstEvent = false;


	TH1F *hist = new TH1F("time", " ", 100000, -delT, delT);
	hist->SetTitle(";#Delta t, s;Nevent, 1/s");

	// FileList opening
	if(argc>2){
		FileMuon.open(argv[1]);
		printf("File: %s opened\n",argv[1]);

		FileHiscore.open(argv[2]);
		printf("File: %s opened\n",argv[2]);
	}
	else {
		printf("File not found\n");
		exit(1);
	}


	double delta = 0;
	int i = 0;
	while (getline(FileMuon, timeM)){
		TimeMuon[i] = atof(timeM.c_str())-11.;
		i++;
	}
	int iMuon = i;

	i= 0;
	while (getline(FileHiscore, timeH)){
		TimeHiscore[i] = atof(timeH.c_str());
		i++;
	}
	int iHiscore = i;

	int iDelta = 0;
	for (i=0;i<iMuon;i++){
		if(i==1000) break;
		printf("%d %d\n",i,iMuon);
		for(int j=0;j<iHiscore;j++){
			delta = TimeMuon[i] - TimeHiscore[j];
			if(delta>=-delT && delta<=delT){
				
				delta_time[iDelta++] = delta;
				
				if(FirstEvent==false){
					StartTime = TimeMuon[i];
					FirstEvent = true;
				}
				EndTime = TimeMuon[i];
			}
		}
	}
	weight = EndTime - StartTime;
	weight = 1./weight;
	for(i=0;i<iDelta;i++){
		hist->Fill(delta_time[i],weight);
	}

	TFile *RootFile = new TFile("out.root","recreate");

	hist->Write("time");
	RootFile->Close();

	delete []TimeMuon;
	delete []TimeHiscore;
	delete []delta_time;

	return 0;
}