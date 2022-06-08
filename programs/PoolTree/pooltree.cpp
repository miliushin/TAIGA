#include <fstream>
#include <cstdio>
#include <iostream>
#include <cstring>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;

int main(int argc, char ** argv) {

  if (argc !=3) {
    printf("Usage: ./Pool muon.root grande.root");
    exit(1);
  }

  TFile  *File1 = new TFile(argv[1],"READ");
  TFile  *File2 = new TFile(argv[2],"UPDATE");

  char Nst[2];
  Nst[1] = 0;

  for(int i=1;i<4;i++){
  	sprintf(Nst,"%d",i);
    TTree *tree = (TTree*)File1->Get(Nst);
    TTree *treeclone = tree->CloneTree();
    File2->cd();
    treeclone->Write("",TObject::kOverwrite);
    delete tree;
    delete treeclone;
  }

  File1->Close();
  File2->Close();

  return 0;
}