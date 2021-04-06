#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/dsigma_dT2.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/velocity_distribution_2000_Ave.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/cpkkd_calculation_New.h"

void  Overlap_TEXONO_stage()
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
    SURF->SetLineWidth(2);
    KS->SetLineColor(4);
    KS->SetLineWidth(2);
    Linear_Line->SetLineColor(4);
    Linear_Line->SetLineStyle(8);
    Linear_Line->SetLineWidth(2);

    TLegend *leg = new TLegend(0.1,0.7,0.4,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->AddEntry("","M_{#chi}=10GeV","");
    //leg->AddEntry(SURF,"Earth-constrained line for CDEX","lP");
    leg->AddEntry(KS,  "Earth-constrained line","lP");
    leg->AddEntry(Linear_Line,"Vacuum-constrained line(X=Y)","lP");

    KS->GetXaxis()->SetTitle("real #sigma_{SI} (cm^{2})");
    KS->GetYaxis()->SetTitle("sensitivities #sigma_{SI} (cm^{2})");
    KS->GetYaxis()->SetTitleSize(0.04);
    KS->GetXaxis()->SetTitleSize(0.04);

    TLine *Earth_Constrained_Line_Left_1;
    Earth_Constrained_Line_Left_1= new TLine(1e-42,1e-42,1e-42,1e-40);
    Earth_Constrained_Line_Left_1->SetLineColor(2);
    Earth_Constrained_Line_Left_1->SetLineStyle(8);
    Earth_Constrained_Line_Left_1->SetLineWidth(4);

    TLine *Earth_Constrained_Line_Right_1;
    Earth_Constrained_Line_Right_1= new TLine(1e-36,1e-42,1e-36,1e-40);
    Earth_Constrained_Line_Right_1->SetLineColor(2);
    Earth_Constrained_Line_Right_1->SetLineStyle(8);
    Earth_Constrained_Line_Right_1->SetLineWidth(4);

    TLine *Earth_Constrained_Line_Right_2;
    Earth_Constrained_Line_Right_2= new TLine(5e-31,1e-42,5e-31,1e-40);
    Earth_Constrained_Line_Right_2->SetLineColor(2);
    Earth_Constrained_Line_Right_2->SetLineStyle(8);
    Earth_Constrained_Line_Right_2->SetLineWidth(4);

    TLine *Earth_Constrained_Line_Right_3;
    Earth_Constrained_Line_Right_3= new TLine(1e-28,1e-42,1e-28,1e-40);
    Earth_Constrained_Line_Right_3->SetLineColor(2);
    Earth_Constrained_Line_Right_3->SetLineStyle(8);
    Earth_Constrained_Line_Right_3->SetLineWidth(4);

    TLatex *tex = new TLatex(5e-41,3e-42,"Stage-1");
    tex->SetTextColor(3);
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->SetLineWidth(2);

    TLatex *tex1 = new TLatex(5e-35,3e-42,"Stage-2");
    tex1->SetTextColor(2);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.05);
    tex1->SetLineWidth(2);

    TLatex *tex2 = new TLatex(6e-31,3e-42,"Stage-3");
    tex2->SetTextColor(1);
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.05);
    tex2->SetLineWidth(2);

    KS->Draw("AL");
    Earth_Constrained_Line_Left_1->Draw("Lsame");
    Earth_Constrained_Line_Right_1->Draw("Lsame");
    Earth_Constrained_Line_Right_2->Draw("Lsame");
    Earth_Constrained_Line_Right_3->Draw("Lsame");
    tex->Draw("same");
    tex1->Draw("same");
    tex2->Draw("same");

    Linear_Line->Draw("Lsame");
    leg->Draw();

    gPad->SetLogx();
    gPad->SetLogy();
    c4->Print("KS_Line.pdf");
}

