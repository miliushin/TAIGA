//
//  Analysis.cpp
//  muon_rs
//
//  Created by Mikhail Iliushin on 18.05.2022.
//

#include "Analysis.hpp"

#include <stdio.h>

// convert char to int
double chd(const char *ch, int ich) {
  char buf[1];
  strncpy(buf, ch + ich, 1);
  return atof(buf);
}

Analysis::Analysis(int opt) {
  Mode = opt;
  ResetNeventMUON();
  if (Mode == 0) {
    TimeGRANDE = new double *[NUMBER_OF_GRANDE_ST];
    for (int i = 0; i < NUMBER_OF_GRANDE_ST; i++) {
      NeventGRANDE[i] = 0;
      TimeGRANDE[i] = new double[MAX_OF_EVENTS];
    }
  }
  TimeMUON = new double *[NUMBER_OF_MUON_ST];
  for (int i = 0; i < NUMBER_OF_MUON_ST; i++) {
    TimeMUON[i] = new double[MAX_OF_EVENTS];
  }

  // Reading delay time of muon stations
  int Nst;
  double Time;
  FILE *FileDelay = fopen("delay", "r");
  if (FileDelay == NULL) {
    printf("ERROR: file delay is not open\n");
  } else {
    while (!feof(FileDelay)) {
      fscanf(FileDelay, "%d %lf\n", &Nst, &Time);
      TimeDelayMUON[Nst - 1] = Time;
    }
  }
  fclose(FileDelay);
}

Analysis::~Analysis() {
  // ----------- write status file -----------------
  char status_file_name[200];
  memset(status_file_name, '\0', 200);
  sprintf(status_file_name, "%s/status.txt", rs_path);
  FILE *FStatus;
  if ((FStatus = fopen(status_file_name, "w")) != NULL) {
    fprintf(FStatus, "         Number of events at MUON stations:\n\n");
    fprintf(FStatus, "     St    Before              After\n");
    for (int i = 0; i < NUMBER_OF_MUON_ST; i++) {
      double PERevent = NeventMUON[i] / NeventMUON_NEW[i] * 100.;
      fprintf(FStatus, "     %d    %d(100%%)           %d(%f%%)\n", i,
              NeventMUON[i], NeventMUON_NEW[i], PERevent);
    }
    fclose(FStatus);
  } else {
    printf("ERROR: can't write status file:\n %s\n", status_file_name);
  }
  // ----------- end write status file--------------
  if (Mode == 0) {
    for (int i = 0; i < NUMBER_OF_GRANDE_ST; i++) {
      delete[] TimeGRANDE[i];
    }
    delete[] TimeGRANDE;
  }

  for (int i = 0; i < NUMBER_OF_MUON_ST; i++) {
    delete[] TimeMUON[i];
  }
  delete[] TimeMUON;
  printf("-----------------------------------------\n");
  printf("\n      Reaserve completed   \n");
  printf("-----------------------------------------\n");
}

void Analysis::ResetNeventMUON() {
  for (int i = 0; i < NUMBER_OF_MUON_ST; i++) {
    NeventMUON[i] = 0;
    NeventMUON_NEW[i] = 0;
  }
}

double Analysis::GetTimeMUON(unsigned char ch[8]) {
  int data_time[4];
  short unsigned int h, m, s, mls, mks, dns;

  for (int i = 0; i < 4; i++) {
    data_time[i / 2] = ch[i * 2 + 1] * 256 + ch[i * 2];
  }

  data_time[0] = ch[1] * 256 + ch[0];
  data_time[1] = ch[3] * 256 + ch[2];
  data_time[2] = ch[5] * 256 + ch[4];
  data_time[3] = ch[6] * 256 + ch[6];

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

int Analysis::ResaveDataMUON(FILE *FileMu, FILE *FileMuNew, int Nstation) {
  unsigned char buf[SIZE_OF_HEADER];
  unsigned long int END_File;
  unsigned long int i_File;
  double TimeEvent;
  unsigned int NumBytes;  // Data packet size
  unsigned char data_read[0x8000];

  fseek(FileMu, 0, SEEK_END);
  END_File = ftell(FileMu);
  i_File = 0;
  fseek(FileMu, 0, SEEK_SET);

  while (i_File < END_File) {
    // Read header of file
    if (!(fread(buf, SIZE_OF_HEADER, 1, FileMu))) {
      printf("ERROR: read file of muon data\n");
      return -1;
    }
    TimeEvent = GetTimeMUON(buf + 12) +
                TimeDelayMUON[Nstation - 1];  // - 10. before 10.12.2022

    NumBytes = buf[3] * 256 + buf[2];

    if (NumBytes > 0x8000) {
      printf("ERROR:  Too much bytes = 0x%X\n", NumBytes);
    }

    if (!fread(data_read, NumBytes, 1, FileMu)) {
      printf("ERROR: read data block\n");
    }
    //  if GRANDE coin
    if (Mode == 0) {
      if (FindCoin_MUON_GRANDE(TimeEvent) == 1) {
        fwrite(buf, SIZE_OF_HEADER, 1, FileMuNew);
        fwrite(data_read, NumBytes, 1, FileMuNew);
        NeventMUON_NEW[Nstation - 1]++;
        continue;
      } else if (FindCoin_MUON_MUON(TimeEvent, Nstation) == 1) {
        fwrite(buf, SIZE_OF_HEADER, 1, FileMuNew);
        fwrite(data_read, NumBytes, 1, FileMuNew);
        NeventMUON_NEW[Nstation - 1]++;
        continue;
      }
    }
    i_File = ftell(FileMu);
  }
  return 0;
}

int Analysis::ReadTimeGRANDE(FILE *File_GRANDE_tim) {
  int Nst, i1, i2, i3;
  char strTime[21];
  char strDelay[20];
  char strDelay_tmp[20];

  while (!feof(File_GRANDE_tim)) {
    memset(strDelay, '\0', 20);
    fscanf(File_GRANDE_tim, "%d%d%d%d%s%s", &Nst, &i1, &i2, &i3, strTime,
           strDelay);

    if (Nst < 31 || Nst > 49) continue;

    double h = 10. * chd(strTime, 0) + chd(strTime, 1);
    double m = 10. * chd(strTime, 3) + chd(strTime, 4);
    double s = 10. * chd(strTime, 6) + chd(strTime, 7);
    double mls =
        100. * chd(strTime, 9) + 10. * chd(strTime, 10) + chd(strTime, 11);
    double mks =
        100. * chd(strTime, 13) + 10. * chd(strTime, 14) + chd(strTime, 15);
    double ns =
        100. * chd(strTime, 17) + 10. * chd(strTime, 18) + chd(strTime, 19);

    for (int i = 0; i < 20; i++) {
      if (strDelay[i] == '(') {
        strDelay_tmp[i] = '\0';
        break;
      }
      strDelay_tmp[i] = strDelay[i];
    }

    TimeGRANDE[Nst - 31][NeventGRANDE[Nst - 31]] =
        3600. * h + 60. * m + s + mls / 1000. + mks / 1000000. +
        ns / 1000000000. + atof(strDelay_tmp);
    NeventGRANDE[Nst - 31]++;
  }

  return 0;
}

int Analysis::FindCoin_MUON_GRANDE(double Time) {
  double dif;
  int Ncoin = 0;
  for (int i = 0; i < NUMBER_OF_GRANDE_ST; i++) {
    for (int j = 0; j < NeventGRANDE[i]; j++) {
      dif = Time - TimeGRANDE[i][j];
      if (dif < DELT && dif > (-1. * DELT)) {
        Ncoin++;
        if (Ncoin == NCOIN_GRANDE) return 1;
      }
    }
  }
  return 0;
}

int Analysis::FindCoin_MUON_MUON(double Time, int Nstation) {
  double dif;
  int Ncoin = 0;
  for (int i = 0; i < NUMBER_OF_MUON_ST; i++) {
    if (i == (Nstation - 1)) continue;
    for (int j = 0; j < NeventMUON[i]; j++) {
      dif = Time - TimeMUON[i][j];
      if (dif < DELT && dif > (-1. * DELT)) {
        Ncoin++;
        if (Ncoin >= NCOIN_GRANDE) return 1;
      }
    }
  }
  return 0;
}

void Analysis::ReadTimeMUON(FILE *FileMu, int Nstation) {
  unsigned char buf[SIZE_OF_HEADER];
  unsigned long int END_File;
  unsigned long int i_File;
  unsigned int NumBytes;
  unsigned char data_read[0x8000];

  fseek(FileMu, 0, SEEK_END);
  END_File = ftell(FileMu);
  i_File = 0;
  fseek(FileMu, 0, SEEK_SET);

  while (i_File < END_File) {
    // Read header of file
    if (!(fread(buf, SIZE_OF_HEADER, 1, FileMu))) {
      printf("ERROR: read file of muon data\n");
      break;
    }

    TimeMUON[Nstation - 1][NeventMUON[Nstation - 1]] =
        GetTimeMUON(buf + 12) +
        TimeDelayMUON[Nstation - 1];  // - 10. before 10.12.2022

    NeventMUON[Nstation - 1]++;
    if (NeventMUON[Nstation - 1] >= MAX_OF_EVENTS) {
      printf("ERROR:Events more MAX_OF_EVENTS\n");
      exit(1);
    }
    NumBytes = buf[3] * 256 + buf[2];

    if (NumBytes > 0x8000) {
      printf("ERROR:  Too much bytes = 0x%X\n", NumBytes);
    }

    if (!fread(data_read, NumBytes, 1, FileMu)) {
      printf("ERROR: read data block\n");
    }

    i_File = ftell(FileMu);
  }
  printf("Nevent: %d\n", NeventMUON[Nstation - 1]);
}

void Analysis::GetRSPATH(char *path) {
  rs_path = new char[strlen(path) + 1];
  strcpy(rs_path, path);
}

int Analysis::GetMode() { return Mode; }
