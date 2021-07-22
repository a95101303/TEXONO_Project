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
#include "TGraph2D.h"

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/velocity_distribution_2000_Ave.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/dsigma_dT2.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/dsigma_dT2_Bent_TEXONO.h"

//200eV threshold ==> 1.009keV
//1.009keV        ==> 201.183(km/s)
void Halo_Model()
{
    TCanvas *c3 = new TCanvas("c3");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.1);
    
    double Vecolity[2000];double Possiblity[2000];
    double Aft_Vecolity[2000];

    double sum; for(int j=0;j<2000;j++){sum = sum + velo_dist_Ave[j][3];}
    for(int j=0;j<2000;j++)
    {
        float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2]);
        Vecolity[j] = v;
        Possiblity[j] = (1/(sum))*velo_dist_Ave[j][3];
    }
    TGraph *Flux = new TGraph(2000,Vecolity,Possiblity);
    TH1F   *Flux_HIST = new TH1F("Flux_HIST","Flux_HIST",2000,0,791);
    Flux_HIST->SetLineColor(1);
    Flux_HIST->SetLineWidth(5);

    Flux->SetTitle("");
    Flux->GetXaxis()->SetTitle("V_{#chi}");
    Flux->GetYaxis()->SetTitle("A.U.");
    Flux->GetYaxis()->SetTitleFont(20);
    Flux->GetYaxis()->SetLabelSize(0.03);

    Flux->Draw("AL");
    
    c3->Print("Halo_Model.pdf");

    return 0;
}

