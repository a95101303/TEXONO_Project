#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_ks_modified_bent_v15_minprint(double mx=2.0)
{
    int Event_Number=1000;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent.C");
    
    string Mass_Point[4]={"20","10","2","0P2"};
      
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,22,1e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,23,2e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,24,3e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,25,4e-29)",Event_Number));
    
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,22,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,23,3e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,24,4e-28)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,26,7e-29)",Event_Number));
    /*
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,27,7e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,28,8e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,29,9e-28)",Event_Number));

    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,3,%i,25,2e-27)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,3,%i,26,3e-27)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,3,%i,27,4e-27)",Event_Number));
     */
    
    /*
     //20GeV
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,23,2e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,22,3e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,23,4e-29)",Event_Number));

    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,21,2e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,22,3e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,23,4e-29)",Event_Number));
     //10GeV
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,23,4e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,24,5e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,25,6e-29)",Event_Number));
     //2GeV
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,21,1e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,22,2e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,22,2e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,23,3e-28)",Event_Number));

    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,21,7e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,3,%i,21,7e-28)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,22,8e-28)",Event_Number));
     */
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,3,%i,24,1e-27)",Event_Number));

    
     
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,20,8e-30)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,21,9e-30)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,22,1e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,23,2e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,24,3e-29)",Event_Number));

        
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,21,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,22,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,23,3e-28)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,20,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,21,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,22,3e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,23,4e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,24,5e-28)",Event_Number));
     
}
//============Earth===========
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,1,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,1,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,2,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,2,1e-32)",Event_Number));

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,2,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,2,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,2,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,2,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,3,%i,2,1e-30)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,2,1e-30)",Event_Number));

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,0,%i,1,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,1,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,2,%i,1,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,1,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,3,%i,1,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,1,1e-32)",Event_Number));

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,20,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,20,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,21,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,21,1e-32)",Event_Number));
//============Air===========

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,20,5e-29)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,21,6e-29)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,22,7e-29)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,23,8e-29)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,1,%i,24,9e-29)",Event_Number));
 

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(2,%i,1,1e-34)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(2,%i,2,1e-32)",Event_Number));

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(3,%i,1,1e-32)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(3,%i,2,1e-30)",Event_Number));

/*
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,%i,1,8e-30)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(0,%i,2,9e-30)",Event_Number));
 */
    
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,20,8e-30)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,21,9e-30)",Event_Number));
/*
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,22,1e-29)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,23,2e-29)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,0,%i,24,3e-29)",Event_Number));

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,20,3e-29)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,21,4e-29)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,22,5e-29)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,23,6e-29)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,1,%i,24,7e-29)",Event_Number));

 gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,19,1e-28)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,20,2e-28)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,21,3e-28)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,22,4e-28)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,23,5e-28)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,24,6e-28)",Event_Number));
    
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,20,1e-28)",Event_Number));
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,21,2e-28)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,22,3e-28)",Event_Number));
gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,23,4e-28)",Event_Number));
 */
//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,19,1e-28)",Event_Number));

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,2,%i,22,4e-28)",Event_Number));

//gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent(1,3,%i,24,5e-28)",Event_Number));


