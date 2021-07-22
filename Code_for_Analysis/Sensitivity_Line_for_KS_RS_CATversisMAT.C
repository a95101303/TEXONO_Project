#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
#include "TLegend.h"

void Sensitivity_Line_for_KS_RS_CATversisMAT()
{
        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

        double X[1]={1};double Y[1]={2};
        TGraph *Tuned_Frame = new TGraph(1,X,Y);
        Tuned_Frame->GetXaxis()->SetRangeUser(0,1);
        Tuned_Frame->GetYaxis()->SetRangeUser(1e7,1e10);

        cout << Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/2GeV/Recoil_Spectrum/NU_%i_STS_Bent.root",FILE) << endl;
        TGraph *The_Vacuum_Case;
        TFile *The_Vacuum_Case_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/2GeV/Recoil_Spectrum/NU_%i_STS_Bent.root",FILE));
            The_Vacuum_Case=(TGraph*)The_Vacuum_Case_FILE->Get("ER_Spectrum_Bef");
            The_Vacuum_Case->SetLineColor(2);
            The_Vacuum_Case->SetLineStyle(1);
            The_Vacuum_Case->GetXaxis()->SetLimits(0,1);
            The_Vacuum_Case->GetYaxis()->SetLimits(1e7,1e10);

        cout << "YES4" << endl;

        TGraph *The_Straight_Case;
        TFile *The_Straight_Case_FILE = TFile::Open(" /Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent_MAT/2GeV/6_STS_Bent_Comparison.root");
            The_Straight_Case=(TGraph*)The_Straight_Case_FILE->Get("ER_Spectrum_Aft");
            The_Straight_Case->SetLineColor(3);
            The_Straight_Case->SetLineStyle(1);
        
        cout << "YES4" << endl;
 
        TGraph *The_Bending_Case;
        TFile *The_Bending_Case_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/2GeV/Recoil_Spectrum/NU_%i_STS_Bent.root",FILE));
            The_Bending_Case=(TGraph*)The_Bending_Case_FILE->Get("ER_Spectrum_Aft");
            The_Bending_Case->SetLineColor(4);
            The_Bending_Case->SetLineStyle(1);

    cout << "YES4" << endl;


            TLegend *leg = new TLegend(0.1,0.1,0.4,0.3);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);

            The_Vacuum_Case->Draw("ALP");
            The_Straight_Case->Draw("LPsame");
            The_Bending_Case->Draw("LPsame");
            
            leg->Draw();
            leg->AddEntry(The_Vacuum_Case,"The Vacuum Case","l");
            leg->AddEntry(The_Straight_Case,"The Straight Line Case","l");
            leg->AddEntry(The_Bending_Case,"The Bending Line Case","l");

        cout << "YES4" << endl;

            c3->SetLogy();
            c3->SetLogx();
            cout << "YES4" << endl;

        c3->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/Recoil_Spectrum/All_%i_2GeV_STS.pdf",FILE));
    
}
    

