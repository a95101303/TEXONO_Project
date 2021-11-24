#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <sys/stat.h>
#include "TROOT.h"

void run_single_sigma_ks_v15_minprint_chie_Low(double mx=2.0)
{
    int LLL=0;
    int File=1;//Xe_c1[0],Xe_d1[1],Ge_c1[2],Ge_d1[3]
    //vector<double> WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};//Xenon
    vector<double> WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07};

    for(int kkk=LLL; kkk<LLL+1; kkk++)
    {
    gROOT->ProcessLine(".L Hist_SetLimit_Plot_v2_Standard_ER.C");
    gROOT->ProcessLine(Form("Hist_SetLimit_Plot_v2_Standard_ER(%i,%i)",kkk,File));
        cout << "WIMP_mx_Array: " << WIMP_mx_Array[kkk] << endl;
    }

}
