#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_ks_modified_bent_v15_minprint_30mwe_real(double mx=2.0)//For KS
{
    //int Event_Number=1;
    int Event_Number=1;
    //int Bent_or_not, int Index_Mass, int Simulated_Event_Number, int Index_Sigma,
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_BentR_30MWE.C");
    
    string Mass_Point[16]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05"};
    //string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P11"};
    
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_BentR_30MWE(1,5,%i,1)",Event_Number));
}


