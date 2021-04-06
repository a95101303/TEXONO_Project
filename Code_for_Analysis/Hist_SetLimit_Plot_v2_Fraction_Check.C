#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

#include "VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"
#include "velocity_distribution_2000_Ave.h"
#include "dsigma_dT2.h"

//200eV threshold ==> 1.009keV
//1.009keV        ==> 201.183(km/s) (10GeV)
//200eV threshold ==> 1.009keV
//1.009keV        ==> 277.443(km/s) (2.34GeV)
int Number_of_Digits(double a)
{
    int count=0;
    while(a/10>1)
    {
        a/=10;
        count++;
    }
    return(count);
}

void Hist_SetLimit_Plot_v2_Fraction_Check()
{
                
    double Vecolity[2000];double Possiblity[2000];
    //
    for(int j=0;j<2000;j++)
    {
        float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2]);
        Vecolity[j] = v;
    }
    //
    double sum;double iteration=0;
    for(int j=0;j<2000;j++)
    {
        if(Vecolity[1999-j]>=277.443)
        {
            sum = sum + velo_dist_Ave[1999-j][3]*1e5;
            if(Number_of_Digits(iteration)!=Number_of_Digits(sum))
            {
                cout << "Vecolity[1999-j]: " << Vecolity[1999-j] << endl;
                cout << "sum: " << sum << endl;
            }
            iteration = sum;
        }
    }
    cout << "sum: " << sum << endl;



}

