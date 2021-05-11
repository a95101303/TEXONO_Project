#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"

void Sensitivity_Line_for_CRESST_Overlap()
{
    string Exp_Flux[1]={"3_CRESST_Flux","4_CRESST_Flux"};
    string Mass_Point[4]={"20","10","2","0P2"};

    for(int Exp=0; Exp<1; Exp++)
    {
        for(int Mass_INT=1; Mass_INT<2; Mass_INT++)
        {
            TH1F *Bent_Line; TH1F *Straight_Line;
            
            TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/20_STS_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                Bent_Line=(TH1F*)fin2->Get("Collision_Time_Hist_Earth");
                Bent_Line->SetLineColor(2);
                Bent_Line->SetLineStyle(1);

            TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/20_STS_Bent_Comparison.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                Straight_Line=(TH1F*)fin3->Get("Collision_Time_Hist_Earth");
                Straight_Line->SetLineColor(3);
                Straight_Line->SetLineStyle(1);
        
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
            Straight_Line->Draw("HISTsame");
            
            leg->Draw();
            leg->AddEntry(Bent_Line ,"Bent","L");
            leg->AddEntry(Straight_Line  ,"Straight","L");
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV_Collision_Time.png",Mass_Point[Mass_INT].c_str()));
        }
    }
    
}
    

