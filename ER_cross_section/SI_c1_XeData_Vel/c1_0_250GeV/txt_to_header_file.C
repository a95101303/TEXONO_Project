#include "reference_file.h"
void txt_to_header_file()
{ 
  const int numfile = 11;
  char name[numfile][200];
  FILE *filename = fopen("filename_header.txt","r");
 float E1[3000], dcs1[3000];
  
  for(int kk = 0; kk<numfile; kk++)
    {
      fscanf(filename,"%s \n",name[kk]);
      printf("%s \n",name[kk]);
   
  
      char txtfile[100];
      sprintf(txtfile,"%s.txt", name[kk]);
  
      FILE *fin = fopen(txtfile,"r");
      char next;
      int rows = 0;
  
      while((next = fgetc(fin)) != EOF)
	{
	  if(next == 'E') 
	    rows++;
	}
      rewind(fin);
      
      const int data_pt = rows;
      printf("Number  %d  \n",data_pt);
     
  
      for(int ll =0; ll<data_pt;ll++)
	{
	  fscanf(fin, "%e   %e \n", &E1[ll], &dcs1[ll]);
	  printf("%e   %e \n", E1[ll], dcs1[ll]);
	}
  
  
      char outfile[100];
      sprintf(outfile,"%s.h", name[kk]);
  
  
      FILE *fout = fopen(outfile,"w");
  
      fprintf(fout,"double %s[data_bin][2] = {\n",name[kk]);
      for(int ll =0; ll<156;ll++)
	{
	  if(ll<rows)
	    {
	      fprintf(fout,"{ %f, %e},\n",E1[ll], dcs1[ll]);
	    }
	  else
	    {
	      fprintf(fout,"{ %f, %e},\n",reference[ll][0], 0.0);  

	    }

	}
    
      fprintf(fout,"};");



    }


}
