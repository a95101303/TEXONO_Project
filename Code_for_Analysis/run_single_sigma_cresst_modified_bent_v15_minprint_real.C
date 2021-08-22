#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_cresst_modified_bent_v15_minprint_real(double mx=2.0)//For KS
{//First Stage
    int Event_Number=5000;
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST.C");
    string Mass_Point[16]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05"};
    //Prove the length of the air
    //==================0.06GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,14,%i,0,6e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,14,%i,1,8e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,14,%i,2,9e-27)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,2,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,3,3e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,4,4e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,5,5e-28)",Event_Number));
    
    //==================0.6GeV==================

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,5,%i,6,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,5,%i,6,1e-28)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,5,%i,7,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,5,%i,9,4e-28)",Event_Number));

    //==================0.06GeV==================
    
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,1,8e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,2,9e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,3,1e-26)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,14,%i,4,8e-29)",Event_Number));
     

    //==================0.07GeV==================
    /*
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,13,%i,1,5e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,13,%i,2,6e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,13,%i,3,7e-29)",Event_Number));
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,13,%i,4,8e-29)",Event_Number));
     */
    //==================0.08GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,12,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,12,%i,4,2e-28)",Event_Number));
   // gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,12,%i,5,3e-28)",Event_Number));
    //==================0.09GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,11,%i,2,9e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,11,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,11,%i,4,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,11,%i,6,4e-28)",Event_Number));

    //==================0.1GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,10,%i,4,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,10,%i,5,3e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,10,%i,6,4e-28)",Event_Number));

    //==================0.2GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,9,%i,4,2e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,9,%i,5,3e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,9,%i,6,4e-27)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,9,%i,7,5e-27)",Event_Number));
    
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,9,%i,0,5e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,9,%i,1,6e-28)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,9,%i,3,8e-28)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,9,%i,5,1e-27)",Event_Number));

    //==================2GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,0,5e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,1,7e-29)",Event_Number));

    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,1,3e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,2,4e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,3,5e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,4,6e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,5,7e-29)",Event_Number));
    
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,0,%i,6,8e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,0,%i,7,9e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,0,%i,8,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(0,0,%i,9,2e-28)",Event_Number));

    //For Brem
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,0,%i,88,6e-29)",Event_Number));

    //==================1GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,1,%i,1,8e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,1,%i,2,9e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,1,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,1,%i,4,2e-28)",Event_Number));
    //==================0.9GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,2,%i,2,9e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,2,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,2,%i,4,2e-28)",Event_Number));
    //==================0.8GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,3,%i,2,9e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,3,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,3,%i,4,2e-28)",Event_Number));
    //==================0.7GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,4,%i,2,9e-29)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,4,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,4,%i,4,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,4,%i,5,3e-28)",Event_Number));
    //==================0.6GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,5,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,5,%i,4,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,5,%i,5,3e-28)",Event_Number));
    //==================0.5GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,6,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,6,%i,4,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,6,%i,5,3e-28)",Event_Number));
    //==================0.4GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,7,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,7,%i,4,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,7,%i,5,3e-28)",Event_Number));
    //==================0.3GeV==================
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,8,%i,3,1e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,8,%i,4,2e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,8,%i,5,3e-28)",Event_Number));
    //gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Bent_From_CRESST(1,8,%i,6,4e-28)",Event_Number));


}
//============Earth===========


