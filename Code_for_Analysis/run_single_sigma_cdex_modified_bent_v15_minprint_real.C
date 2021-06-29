#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_cdex_modified_bent_v15_minprint_real(double mx=2.0)//For KS
{
    //int Event_Number=1;
    int Event_Number=1000;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_CDEX_Run_MODIFIED_BentR.C");
    
    string Mass_Point[16]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05"};

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_CDEX_Run_MODIFIED_BentR(1,0,%i,2,2e-31)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_CDEX_Run_MODIFIED_BentR(1,0,%i,3,3e-31)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_CDEX_Run_MODIFIED_BentR(1,0,%i,4,4e-31)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_CDEX_Run_MODIFIED_BentR(0,0,%i,0,5e-31)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_CDEX_Run_MODIFIED_BentR(1,0,%i,0,5e-31)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_CDEX_Run_MODIFIED_BentR(0,0,%i,0,5e-36)",Event_Number));

}


