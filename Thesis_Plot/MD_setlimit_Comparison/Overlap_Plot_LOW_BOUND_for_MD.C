#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLine.h"

void Overlap_Plot_LOW_BOUND_for_MD()
{
            //====================Input=====================//
            TGraphErrors *Data;TGraph *MD_1;TGraph *MD_2;
    
    TFile *finMD = TFile::Open("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/LOW_BOUND/ROOT/Low_Mass/MD_1_TEXONO.root");
            MD_1=(TGraph*)finMD->Get("RS");
            MD_1->SetLineColor(2);
            MD_1->SetMarkerColor(2);
            MD_1->SetLineStyle(2);
            MD_1->SetLineWidth(4);

    TFile *finBR = TFile::Open("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/LOW_BOUND/ROOT/Low_Mass/MD_14_TEXONO.root");
            MD_2=(TGraph*)finBR->Get("RS");
            MD_2->SetLineColor(4);
            MD_2->SetMarkerColor(4);
            MD_2->SetLineStyle(2);
            MD_2->SetLineWidth(4);

            Data=(TGraphErrors*)finBR->Get("TEXONOData");
            Data->SetLineColor(1);
            Data->SetMarkerColor(1);
            Data->SetLineStyle(1);
            Data->SetLineWidth(1);

            Data->GetXaxis()->SetRangeUser(0,2.4);
            Data->GetYaxis()->SetRangeUser(1,80);

            Data->GetXaxis()->SetTitle("E_{det}[keVee]");
            Data->GetXaxis()->CenterTitle();
            Data->GetYaxis()->CenterTitle();

            TCanvas *c3 = new TCanvas("c3");
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);
            gStyle->SetTitleSize(0.01,"XY");
            gStyle->SetTitleFont(50,"XY");
            gStyle->SetLegendFont(62);

            TLegend *leg = new TLegend(0.2,0.6,0.4,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.06);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            leg->AddEntry(MD_1,"M_{#chi}=1GeV/c^{2}, #sigma^{SI}_{#chi N}:1.511 #times 10^{-35} cm^{2}","l");
            leg->AddEntry(MD_2,"M_{#chi}=0.05GeV/c^{2}, #sigma^{SI}_{#chi N}:1.59 #times 10^{-31} cm^{2}","l");

            Data->Draw("AP");
            MD_1->Draw("LPsame");
            MD_2->Draw("LPsame");
                
            leg->Draw();
            c3->Print("TEXONO_DATA_Comparison_MD.pdf");
}
    

