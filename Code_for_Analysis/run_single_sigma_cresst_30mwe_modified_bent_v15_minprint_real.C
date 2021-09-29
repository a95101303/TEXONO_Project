#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_cresst_30mwe_modified_bent_v15_minprint_real(double mx=2.0)//For KS
{//First Stage
    int Event_Number=2500;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE.C");
    string Mass_Point[16]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05"};
    
    //==================2GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,31,1e-37)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,33,1e-35)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,35,1e-33)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,36,1e-32)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,37,1e-31)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,45,1e-30)",Event_Number));
        
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,46,3e-30)",Event_Number));
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,47,5e-30)",Event_Number));
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,48,7e-30)",Event_Number));
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,49,9e-30)",Event_Number));
         
        
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,55,1e-29)",Event_Number));
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,56,2e-29)",Event_Number));
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,57,3e-29)",Event_Number));
      //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,0,%i,58,4e-29)",Event_Number));
         

    //==================0.16GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,30,1e-38)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,31,1e-37)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,33,1e-35)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,34,1e-34)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,35,1e-33)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,36,1e-32)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,37,1e-31)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,38,1e-30)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,39,1e-29)",Event_Number));
    
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,40,3e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,41,5e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,42,7e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,43,9e-29)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,44,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,45,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,46,3e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,47,4e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,48,5e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,49,6e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_CRESST_30MWE(1,9,%i,50,7e-28)",Event_Number));
}
//============Earth===========


