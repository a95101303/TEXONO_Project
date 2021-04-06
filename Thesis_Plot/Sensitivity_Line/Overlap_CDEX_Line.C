#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/dsigma_dT2.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/velocity_distribution_2000_Ave.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/cpkkd_calculation_New.h"

void Overlap_CDEX_Line()
{
    TCanvas *c4 = new TCanvas("c4");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    char fname[100];string Sfname;
    char fname1[100];string Sfname1;

    sprintf(fname,"/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/10GeV/NU_Line.root ");
    sprintf(fname1,"/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/2_TEXONO_Flux/0P8GeV/MD_Line.root");

    
    TFile *fin=new TFile(fname);
    TFile *fin1=new TFile(fname1);

    TGraph *SURF=(TGraph*)fin ->Get("Threshold_Plot");
    TGraph *KS  =(TGraph*)fin1->Get("Threshold_Plot");
    
    TF1 *Linear_Line = new TF1("Linear_Line","x",1e-42,1e-27);

    SURF->SetLineColor(2);
    SURF->SetLineWidth(2);
    KS->SetLineColor(3);
    KS->SetLineWidth(2);
    Linear_Line->SetLineColor(4);
    Linear_Line->SetLineWidth(2);

    TLegend *leg = new TLegend(0.1,0.6,0.4,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->AddEntry("","M_{#chi}=0.8GeV","");

    KS->GetXaxis()->SetTitle("real #sigma_{SI} (cm^{2})");
    KS->GetYaxis()->SetTitle("sensitivities #sigma_{SI} (cm^{2})");
    KS->GetYaxis()->SetTitleSize(0.04);
    KS->GetXaxis()->SetTitleSize(0.04);

    TLine *MD_Line_STFlux;
    MD_Line_STFlux= new TLine(0,2.722e-35,1,2.722e-35);
    MD_Line_STFlux->SetLineColor(4);
    MD_Line_STFlux->SetLineStyle(8);
    MD_Line_STFlux->SetLineWidth(3);

    TLine *Vacuum_Constrained_Line;
    Vacuum_Constrained_Line= new TLine(2.722e-35,1e-42,2.722e-35,2.722e-35);
    Vacuum_Constrained_Line->SetLineColor(4);
    Vacuum_Constrained_Line->SetLineStyle(8);
    Vacuum_Constrained_Line->SetLineWidth(2);
    
    TLine *Earth_Constrained_Line_Left;
    Earth_Constrained_Line_Left= new TLine(3.3e-35,1e-42,3.3e-35,2.722e-35);
    Earth_Constrained_Line_Left->SetLineColor(2);
    Earth_Constrained_Line_Left->SetLineStyle(8);
    Earth_Constrained_Line_Left->SetLineWidth(2);

    TLine *Earth_Constrained_Line_Right;
    Earth_Constrained_Line_Right= new TLine(5e-28,1e-42,5e-28,2.722e-35);
    Earth_Constrained_Line_Right->SetLineColor(2);
    Earth_Constrained_Line_Right->SetLineStyle(8);
    Earth_Constrained_Line_Right->SetLineWidth(2);

    leg->AddEntry(Linear_Line,  "Vacuum-constrained line(X=Y)","lP");
    leg->AddEntry(KS,  "Earth-constrained line","lP");

    leg->AddEntry(MD_Line_STFlux,  "Vacuum-constrained #sigma_{SI}","l");
    leg->AddEntry(Earth_Constrained_Line_Right,  "Earth-constrained #sigma_{SI}(Final result)","l");

    KS->Draw("AL");
    Vacuum_Constrained_Line->Draw("Lsame");
    Earth_Constrained_Line_Left->Draw("Lsame");
    Earth_Constrained_Line_Right->Draw("Lsame");
    Linear_Line->Draw("Lsame");
    MD_Line_STFlux->Draw("Lsame");
    leg->Draw();

    gPad->SetLogx();
    gPad->SetLogy();
    c4->Print("KS_Line_overlap.pdf");
}

