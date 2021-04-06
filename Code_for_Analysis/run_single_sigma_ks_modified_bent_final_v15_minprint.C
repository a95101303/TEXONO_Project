#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_ks_modified_bent_final_v15_minprint(double mx=2.0)
{
    int Event_Number=5000;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Final.C");
    
    string Mass_Point[3]={"20","2","0P2"};
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Final(0,%i,1,1e-30)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Final(1,%i,1,1.1e-28)",Event_Number));//2
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Final(2,%i,1,1e-27)",Event_Number));//0.2

}
