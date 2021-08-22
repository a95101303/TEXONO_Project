#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLine.h"

void Overlap_Plot_Flux()
{
        TFile *ROOT_FILE  = TFile::Open("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent_MAT/0P6GeV/6_STS_Bent.root");

        TFile *ROOT_FILE1 = TFile::Open("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent_MAT/0P6GeV/6_STS_Bent_Comparison.root");

            //====================Input=====================//
        TH1F *Bent_Flux;TH1F *Non_Bent_Flux;
        Bent_Flux=(TH1F*)ROOT_FILE->Get("Flux_HIST_Aft_Collision_EARTH");
        Bent_Flux->SetTitle("Velocity distirbutioin(SHIELDING)");
        Bent_Flux->SetLineColor(3);
        Bent_Flux->SetMarkerColor(3);
        Bent_Flux->SetLineStyle(1);
        Bent_Flux->SetLineWidth(3);
        Bent_Flux->GetXaxis()->SetRangeUser(0,800);
        Bent_Flux->GetYaxis()->SetRangeUser(0,1);
        Bent_Flux->Scale(1./Bent_Flux->Integral());
        
        Non_Bent_Flux=(TH1F*)ROOT_FILE1->Get("Flux_HIST_Aft_Collision_EARTH");
        Non_Bent_Flux->SetLineColor(2);
        Non_Bent_Flux->SetMarkerColor(2);
        Non_Bent_Flux->SetLineStyle(1);
        Non_Bent_Flux->SetLineWidth(2);
        //Non_Bent_Flux->GetXaxis()->SetLimits(0,5);
        //Non_Bent_Flux->GetYaxis()->SetRangeUser(0,1);
        Non_Bent_Flux->Scale(1./Non_Bent_Flux->Integral());

        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

        TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(Bent_Flux,"MAT(Bending)","l");
        leg->AddEntry(Non_Bent_Flux,"CAT(Straight)","l");

        Bent_Flux->Draw("HIST");
        Non_Bent_Flux->Draw("HISTsame");
        leg->Draw();
        c3->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent_MAT/0P6GeV/Flux.pdf");
}
    

