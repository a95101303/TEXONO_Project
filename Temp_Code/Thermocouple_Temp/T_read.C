#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//int main(int argc,char *argv[])
//
//void T_read(Int_t argc,Char_t *argv[])
void T_read()
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

  //if (argc>1)
  // strcpy(file_path, argv[1]); 
  //if (argc>2) { init_run=atoi(argv[2]); }
  //if (argc>3) { final_run=atoi(argv[3]); }

  // sprintf(file_path,"/SynRaid3/p104_rawdata_twoDet_part2/210424/210424_1445_PhysOn_Det70x70_ThrdM520_Det50x50_ThrdM570_FADC_RAW_Data");

   sprintf(file_path,"/SynRaid1/p104_rawdata/200911/200911_1753_PhysOn_Det70x70NaI3ThrdP685_Det50x50NaI5Thrd0_FADC_RAW_Data");

   
   
  //
  for(number_of_file=1; number_of_file<20; number_of_file++)
  //while(end_of_set)
  {
    sprintf(file_path1,"%s_%d.bin",file_path,number_of_file);
    sprintf(file_path3,"%s_%d.bin",file_path,number_of_file+1);

    raw1=fopen(file_path1,"rb");
    pfile=fopen(file_path3,"rb");

    fseek(raw1,-5*sizeof(Double_t),SEEK_END);
    fread(&T,sizeof(Double_t),4,raw1);
    fread(&RET,sizeof(Double_t),1,raw1);
    printf("%d %f  %f  %f  %f  %f\n",number_of_file,T[0],T[1],T[2],T[3],RET);
   	


    // number_of_file++;

    /*
    if(pfile!=NULL)
    { 
	    
      fclose(raw3);
                  
    }
    else
    {
       end_of_set=false;
   fclose(raw3);
      
    }

    fclose(raw1);
    */
    
  }





}
