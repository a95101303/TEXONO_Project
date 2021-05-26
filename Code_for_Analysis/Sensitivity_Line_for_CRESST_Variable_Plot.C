#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"

void Sensitivity_Line_for_CRESST_Variable_Plot()
{
    string Exp_Flux[2]={"3_CRESST_Flux","4_CRESST_Flux"};
    string Mass_Point[3]={"20","2","0P2"};

    for(int Exp=1; Exp<2; Exp++)
    {
        for(int Mass_INT=1; Mass_INT<2; Mass_INT++)
        {
            TH1F *Bent_Line; TH1F *Straight_Line;
            
            TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/Hist_49_STS_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                Bent_Line=(TH1F*)fin2->Get("Ratio_of_Earth_Origin_After");
                Bent_Line->SetLineColor(2);
                Bent_Line->SetLineStyle(1);
            Bent_Line->GetYaxis()->SetRangeUser(0,1);
            Bent_Line->GetXaxis()->SetRangeUser(0,100);
                Bent_Line->Scale(1/(Bent_Line->GetEntries()));

            TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/Hist_49_STS_Bent_Comparison.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                Straight_Line=(TH1F*)fin3->Get("Ratio_of_Earth_Origin_After");
                Straight_Line->SetLineColor(3);
                Straight_Line->SetLineStyle(1);
                Straight_Line->Scale(1/(Straight_Line->GetEntries()));

            TCanvas *c3 = new TCanvas("c3");
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);
            cout << "YES4" << endl;

            TLegend *leg = new TLegend(0.1,0.1,0.4,0.3);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);

        //Bent_Line->GetXaxis()->SetRangeUser(0,0.1);
        //Bent_Line->GetYaxis()->SetRangeUser(0,200);

            Bent_Line->Draw("HIST");
            //Straight_Line->Draw("HISTsame");
            
            leg->Draw();
            //leg->AddEntry(Bent_Line ,"Bent","L");
            //leg->AddEntry(Straight_Line  ,"Straight","L");
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/Ratio_of_Earth_Origin_After.pdf",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
        }
    }
    
}
    

