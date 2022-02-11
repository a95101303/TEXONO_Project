#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void Manu_T_read()
{
  Char_t file_path[1000];
  Char_t file_path1[1000];
  Char_t file_path2[1000];
  Char_t file_path3[1000];
  Int_t number_of_file=1;
  Bool_t end_of_set=true;
  FILE *raw1,*raw2,*raw3,*pfile;
  Double_t T[4];
  Double_t RET; // RUN End Time
  
  FILE *file1 =fopen("220105_1221_PhysOFF_Det70x70_ThrdM522_Det50x50_Temp188_ThrdM565_FADC_RAW_Data.txt","w");
  sprintf(file_path,"/SynRaid5/p104_rawdata_twoDet_part4/220105/220105_1221_PhysOFF_Det70x70_ThrdM522_Det50x50_Temp188_ThrdM565_FADC_RAW_Data");
  
  
  for(number_of_file=1; number_of_file<1207; number_of_file++)
     {
       sprintf(file_path1,"%s_%d.bin",file_path,number_of_file);
       sprintf(file_path3,"%s_%d.bin",file_path,number_of_file+1);
       
       raw1=fopen(file_path1,"rb");
       pfile=fopen(file_path3,"rb");
       
       fseek(raw1,-5*sizeof(Double_t),SEEK_END);
       fread(&T,sizeof(Double_t),4,raw1);
       fread(&RET,sizeof(Double_t),1,raw1);
       
       printf("%d      %f      %f      %f      %f       %f\n",number_of_file,T[0],T[1],T[2],T[3],RET);
       fprintf(file1,"%d      %f      %f      %f      %f     %f\n", number_of_file,T[0],T[1],T[2],T[3],RET);
     }
}
