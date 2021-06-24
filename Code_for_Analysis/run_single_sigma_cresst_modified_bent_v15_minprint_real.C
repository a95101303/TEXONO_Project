#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_cresst_modified_bent_v15_minprint_real(double mx=2.0)
{
    int Event_Number=1000;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST.C");
    
    string Mass_Point[5]={"20","10","2","0P2","0P05"};//First Round
    //Prove the length of the air
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,2,%i,1,3e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,2,%i,2,7e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,2,%i,3,8e-29)",Event_Number));

     
}
//============Earth===========


