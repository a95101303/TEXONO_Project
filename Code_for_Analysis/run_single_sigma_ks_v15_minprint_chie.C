#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_ks_v15_minprint_chie(double mx=2.0)
{
    int Event_Number=1;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_chie.C");
    
    //double WIMP_Mass_Array[15]={3,4,5,10,20};

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(0,%i,100,2e-28)",Event_Number));//V
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(1,%i,100,1.5e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(2,%i,100,5e-28)",Event_Number));
    // gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(3,%i,100,3e-22)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(4,%i,100,5e-27)",Event_Number));

}
