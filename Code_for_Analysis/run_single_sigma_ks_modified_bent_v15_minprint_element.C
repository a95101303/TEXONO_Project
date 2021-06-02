#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_ks_modified_bent_v15_minprint_element(double mx=2.0)
{
    int Event_Number=1000;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material.C");
    
    string Mass_Point[4]={"15","2","0P2","0P05"};//First Round

    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material(1,0,%i,1,8e-31,0)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material(1,0,%i,1,8e-31,1)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material(1,0,%i,1,8e-31,2)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material(1,0,%i,1,8e-31,3)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material(1,0,%i,1,8e-31,4)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material(1,0,%i,1,8e-31,5)",Event_Number));

}


