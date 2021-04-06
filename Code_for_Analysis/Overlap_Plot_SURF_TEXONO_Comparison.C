#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"

void Overlap_Plot_SURF_TEXONO_Comparison()
{
    TCanvas *c4 = new TCanvas("c4");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    char fname[100];string Sfname;
    char fname1[100];string Sfname1;

    sprintf(fname,"/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/10GeV/NU_Line.root ");
    sprintf(fname1,"/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/2_TEXONO_Flux/10GeV/NU_Line.root");

    
    TFile *fin=new TFile(fname);
    TFile *fin1=new TFile(fname1);

    TGraph *SURF=(TGraph*)fin ->Get("Threshold_Plot");
    TGraph *KS  =(TGraph*)fin1->Get("Threshold_Plot");
    TF1 *Linear_Line = new TF1("Linear_Line","x",1e-42,1e-27);

    SURF->SetLineColor(2);
    KS->SetLineColor(3);
    Linear_Line->SetLineColor(4);

    TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->AddEntry("","M_{#chi}=10GeV (Migdal Effect)","");
    leg->AddEntry(SURF,"CDEX","lP");
    leg->AddEntry(KS,  "KS","lP");
    
    KS->Draw("AL");
    SURF->Draw("Lsame");
    Linear_Line->Draw("Lsame");
    leg->Draw();

    gPad->SetLogx();
    gPad->SetLogy();
    c4->Print("CDEX_KS_Comparison.pdf");
}

