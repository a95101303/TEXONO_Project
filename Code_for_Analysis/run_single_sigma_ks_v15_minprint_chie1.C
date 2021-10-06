#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_ks_v15_minprint_chie1(double mx=2.0)
{
    int Event_Number=1;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_chie1.C");
    
    //double WIMP_Mass_Array[15]={3,4,5,10,20};

    //m_A'=0
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie1(1,%i,100,2e-24)",Event_Number));
    //m_A'=1e6
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(0,%i,100,2e-36)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(1,%i,100,1e-33)",Event_Number));

}
