//
//  main.cpp
//  muon_rs
//
//  Created by Mikhail Iliushin on 12.05.2022.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <dirent.h>

#include "Analysis.hpp"

using namespace std;



void ReadGRANDE(string data, Analysis* MuonResave){
    char DATA_DIR[200];
    char ST_DIR[200];
    char FTIME[200];
    string GRANDE_PATH;
    DIR *main_dir;      // main dir of GRANDE
    DIR *data_dir;      // data dir for day
    DIR *st_dir;        // station dir
    struct dirent *main_ent;
    struct dirent *data_ent;
    struct dirent *st_ent;
    bool FlagDataDir;
    bool FlagTim;
    FILE *File_GRANDE_tim;
    ifstream GRANDE_LIST("grande.dir");  // path to GRANDE data directory
    
    // ----------------  read grande
    if (!GRANDE_LIST){
        printf("ERROR: grande.dir is not open\n");
        exit(1);
    }
    else{
        while(getline(GRANDE_LIST,GRANDE_PATH)){
            main_dir = opendir(GRANDE_PATH.c_str());
            // Read main GRANDE directory
            while ((main_ent = readdir (main_dir)) != NULL) {
                if(main_ent->d_type==DT_DIR){
                    FlagDataDir = false;
                    memset(DATA_DIR,'\0',200);
                    sprintf(DATA_DIR, "%s",main_ent->d_name);
                    for(int i=0;i<6;i++){
                        if(DATA_DIR[i]!=data[i]){
                            FlagDataDir = true;
                            break;
                        }
                    }
                    if(FlagDataDir==false){
                        sprintf(DATA_DIR, "%s%s",GRANDE_PATH.c_str(),main_ent->d_name);
                        data_dir = opendir(DATA_DIR);
                        // Read GRANDE data dir for day
                        while ((data_ent = readdir(data_dir)) != NULL){
                            if(data_ent->d_type==DT_DIR){
                                memset(ST_DIR,'\0',200);
                                sprintf(ST_DIR,"%s",data_ent->d_name);
                                if(ST_DIR[0]=='S' && ST_DIR[1]=='t'){
                                    sprintf(ST_DIR,"%s/%s",DATA_DIR,data_ent->d_name);
                                    // Read GRANDE station dir
                                    st_dir  = opendir(ST_DIR);
                                    FlagTim = false;
                                    while((st_ent = readdir(st_dir)) != NULL){
                                        if (st_ent->d_type==DT_REG) {
                                            memset(FTIME,'\0',200);
                                            sprintf(FTIME,"%s",st_ent->d_name);
                                            for(int i=0;i<20;i++){
                                                if(FTIME[i]=='.' && FTIME[i+1]=='t' && FTIME[i+2]=='i'){
                                                    sprintf(FTIME,"%s/%s",ST_DIR,st_ent->d_name);
                                                    File_GRANDE_tim = fopen(FTIME,"r");
                                                    MuonResave->ReadTimeGRANDE(File_GRANDE_tim);
                                                    fclose(File_GRANDE_tim);
                                                    FlagTim = true;
                                                    break;
                                                }
                                            }
                                            if (FlagTim==true)
                                                break;
                                        }
                                    }
                                    closedir(st_dir);
                                }
                            }
                        }
                        closedir(data_dir);
                    }
                }
            }
            closedir(main_dir);
        }
    }
    GRANDE_LIST.close();
    // ----------------  end  read grande
}

void ReadMUON(string data, Analysis* MuonResave){
    char Nst[3];  // char number of muon station
    char DATA_DIR[200];
    char ST_DIR[200];
    char RS_DIR[200];      // new dir for reasave data
    char mk_dir[200];
    char fname[200];        // file data name
    char RSfname[200];
    DIR *main_dir;          // main dir of MOUN
    DIR *data_dir;          // data dir for day
    DIR *st_dir;            // station dir
    struct dirent *main_ent;
    struct dirent *data_ent;
    struct dirent *st_ent;

    bool FlagDataDir;
    string MUON_PATH;
    ifstream MUON_DIR("muon.dir");      // path to muon data directory
    FILE *MuFile;
    FILE *RSMuFile;

    
    
    
    if (!MUON_DIR){
        printf("ERROR: muon.dir is not open\n");
        exit(1);
    }
    else{
        while (getline(MUON_DIR, MUON_PATH)) {
            
            // Open data directory
            main_dir = opendir(MUON_PATH.c_str());
            
            if(!main_dir){
                printf("ERROR: DIR %s is not open\n",MUON_PATH.c_str());
                exit(1);
            }
            
            // Read main MUON directory
            while ((main_ent = readdir (main_dir)) != NULL) {
                if(main_ent->d_type==DT_DIR){
                    FlagDataDir = false;
                    memset(DATA_DIR,'\0',200);
                    sprintf(DATA_DIR, "%s",main_ent->d_name);
                    for(int i=0;i<6;i++){
                        if(DATA_DIR[i]!=data[i]){
                            FlagDataDir = true;
                            break;
                        }
                    }
                    if (FlagDataDir==false){
                        sprintf(DATA_DIR, "%s%s",MUON_PATH.c_str(),main_ent->d_name);
                        memset(RS_DIR,'\0',200);
                        if(MuonResave->GetMode()==0)
                            sprintf(RS_DIR, "%s.rsg/",DATA_DIR);
                        else if(MuonResave->GetMode()==1)
                            sprintf(RS_DIR, "%s.rsh/",DATA_DIR);
                        MuonResave->GetRSPATH(RS_DIR);
                        MuonResave->ResetNeventMUON();
                        memset(mk_dir,'\0',200);
                        sprintf(mk_dir, "mkdir %s", RS_DIR);
                        system(mk_dir);   // make dir .rsh(.rsg)
                        data_dir = opendir(DATA_DIR);
                        // Read MUON data dir for day
                        while((data_ent = readdir(data_dir)) != NULL){
                            if(strcmp(data_ent->d_name, ".") == 0 ||
                               strcmp(data_ent->d_name, "..") == 0)
                                continue;
                            if (data_ent->d_type==DT_DIR){
                                memset(ST_DIR,'\0',200);
                                sprintf(ST_DIR, "%s",data_ent->d_name);
                                if(ST_DIR[0]=='M' && ST_DIR[1]=='u'){
                                    Nst[0] = ST_DIR[2];
                                    Nst[1] = ST_DIR[3];
                                    Nst[2] = 0;
                                    
                                    sprintf(ST_DIR, "%s/%s",DATA_DIR,data_ent->d_name);
                                    memset(mk_dir,'\0',200);
                                    sprintf(mk_dir,"mkdir %s%s",RS_DIR,data_ent->d_name);
                                    system(mk_dir);
                                    // Read muon station dir
                                    st_dir = opendir(ST_DIR);
                                    while ( (st_ent = readdir(st_dir)) != NULL) {
                                        if (st_ent->d_type==DT_REG){
                                            memset(fname, '\0', 200);
                                            sprintf(fname,"%s/%s",ST_DIR,st_ent->d_name);
                                            if ((MuFile = fopen(fname,"rb")) == NULL){
                                                printf("OPEN ERROR:\n %s\n",fname);
                                            }
                                            else{
                                                printf("READ:\n %s\n",fname);
                                                MuonResave->ReadTimeMUON(MuFile, atoi(Nst));
                                                fclose(MuFile);
                                            }
                                        }
                                    }
                                    closedir(st_dir);
                                }
                                else{
                                    memset(mk_dir,'\0',200);
                                    sprintf(mk_dir, "cp -r %s/%s %s/%s",DATA_DIR,data_ent->d_name,RS_DIR,data_ent->d_name);
                                    system(mk_dir);
                                }
                            }
                            else{
                                memset(mk_dir,'\0',200);
                                sprintf(mk_dir, "cp -r %s/%s %s/%s",DATA_DIR,data_ent->d_name,RS_DIR,data_ent->d_name);
                                system(mk_dir);
                            }
                        }
                        rewinddir(data_dir);
                        // Read MUON data dir for day
                        while((data_ent = readdir(data_dir)) != NULL){
                            if (data_ent->d_type==DT_DIR){
                                memset(ST_DIR,'\0',200);
                                sprintf(ST_DIR, "%s",data_ent->d_name);
                                if(ST_DIR[0]=='M' && ST_DIR[1]=='u'){
                                    Nst[0] = ST_DIR[2];
                                    Nst[1] = ST_DIR[3];
                                    Nst[2] = 0;
                                    
                                    sprintf(ST_DIR, "%s/%s",DATA_DIR,data_ent->d_name);
                                    // Read muon station dir
                                    st_dir = opendir(ST_DIR);
                                    while ( (st_ent = readdir(st_dir)) != NULL) {
                                        if (st_ent->d_type==DT_REG){
                                            memset(fname, '\0', 200);
                                            sprintf(fname,"%s/%s",ST_DIR,st_ent->d_name);
                                            memset(RSfname, '\0', 200);
                                            sprintf(RSfname,"%s%s/%s",RS_DIR,data_ent->d_name,st_ent->d_name);
                                            if ((MuFile = fopen(fname,"rb")) == NULL){
                                                printf("OPEN ERROR:\n %s\n",fname);
                                            }
                                            else{
                                                
                                                if( (RSMuFile = fopen(RSfname,"wb")) == NULL){
                                                    printf("WRITE ERROR: \n %s\n", RSfname);
                                                }
                                                else{
                                                    printf("RESAVE:\n %s\n",RSfname);
                                                    MuonResave->ResaveDataMUON(MuFile, RSMuFile, atoi(Nst));
                                                }
                                            }
                                        }
                                    }
                                    closedir(st_dir);
                                }
                            }

                        }
                        
                        closedir(data_dir);
                    }
            
            
                }
            }
            closedir(main_dir);
        }
    }
            
}


int main(int argc, const char * argv[]) {
    
    string data;
    ifstream DATA_LIST("data.list");    // directory for analysis
    
    Analysis *MuonResave = new Analysis(0);

    if(!DATA_LIST){
        printf("ERROR: data.list is not open\n");
        exit(1);
    }
    else{
        while(getline(DATA_LIST,data)){
            if(MuonResave->GetMode()==0)
                ReadGRANDE(data, MuonResave);
            ReadMUON(data, MuonResave);   // read and resave
        }
    }

    delete MuonResave;
    
    return 0;
}
