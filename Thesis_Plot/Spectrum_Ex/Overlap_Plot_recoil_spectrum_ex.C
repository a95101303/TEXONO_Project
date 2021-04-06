#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"

void Overlap_Plot_recoil_spectrum_ex()//CDEX:Threshold=160eV, TEXONO:Threshold=200eV
{
        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

        //==============Mass==============//
        TFile *ROOT_FILE1 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Thesis_Plot/Spectrum_Ex/NU.root"));
        TGraph *NU = (TGraph*)ROOT_FILE1->Get("Graph");
    
        NU->GetXaxis()->SetTitle("Energy[keVee]");
        NU->GetYaxis()->SetTitle("Count[Evts/kg/day]");
        NU->GetXaxis()->SetLimits(1e-2,1e1);
        NU->GetYaxis()->SetRangeUser(1e-7,1e+7);
        //NU->SetName("NU");
        NU->SetLineColor(2);
        NU->SetMarkerColor(2);
        NU->SetLineWidth(5);

        TFile *ROOT_FILE2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Thesis_Plot/Spectrum_Ex/MD.root"));
        TGraph *MD = (TGraph*)ROOT_FILE2->Get("Graph");
        MD->SetLineColor(3);
        MD->SetMarkerColor(3);
        TFile *ROOT_FILE3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Thesis_Plot/Spectrum_Ex/BR.root"));
        TGraph *BR = (TGraph*)ROOT_FILE3->Get("Graph");
        BR->SetLineColor(4);
        BR->SetMarkerColor(4);

        //==============Input STandard Flux as well as the "attenuated" flux==============
        //=======================Recoil Spectrum set for three processes==============================

        //Energy recoil Spectrum
        TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        //leg->AddEntry("",Form("#sigma_{SI}:%.2f #times 10^{-%.f}, M_{#chi}=%.2fGeV",CNFNV(0,sigma_si),CNFNV(1,sigma_si),mx),"");
        leg->AddEntry(NU,"NR","l");
        leg->AddEntry(MD,"MD","l");
        leg->AddEntry(BR,"BR","l");

        NU->Draw("ALP");
        MD->Draw("LPsame");
        BR->Draw("LPsame");
        leg->Draw();
        c3->SetLogy();
        c3->SetLogx();
    
        c3->Print("Recoil_Spectrum_ex.pdf");

}
    

