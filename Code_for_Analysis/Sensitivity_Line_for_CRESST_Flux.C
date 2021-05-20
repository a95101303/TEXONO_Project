#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLine.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/dsigma_dT2.h"

void Sensitivity_Line_for_CRESST_Flux()
{
    string Exp_Flux[1]={"4_CRESST_Flux"};
    string Mass_Point[1]={"2"};
    
        for(int Mass_INT=0; Mass_INT<1; Mass_INT++)
        {
            for(int FILE=49; FILE<50; FILE++)
            {
                cout << "Mass: " << Mass_Point[Mass_INT].c_str() << endl;
                cout << "FILE: " << FILE << endl;
        TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/%sGeV/49_STS_Bent_Comparison.root",Mass_Point[Mass_INT].c_str(),FILE));
                cout << "Mass: " << Mass_Point[Mass_INT].c_str() << endl;
                cout << "FILE: " << FILE << endl;

  TFile *ROOT_FILE_Straight = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/%sGeV/%i_STS_Bent_Comparison.root",Mass_Point[Mass_INT].c_str(),FILE));
                cout << "Mass: " << Mass_Point[Mass_INT].c_str() << endl;
                cout << "FILE: " << FILE << endl;

    TFile *ROOT_FILE_Bent = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/%sGeV/%i_STS_Bent.root",Mass_Point[Mass_INT].c_str(),FILE));
                cout << "Mass: " << Mass_Point[Mass_INT].c_str() << endl;
                cout << "FILE: " << FILE << endl;

                //====================Input=====================//
            TH1F *Bef_Flux;TH1F *Aft_Flux_S;TH1F *Aft_Flux_B;
            Bef_Flux=(TH1F*)ROOT_FILE->Get("Flux_HIST_Random");
            Bef_Flux->SetLineColor(3);
            Bef_Flux->SetMarkerColor(3);
            Bef_Flux->SetLineStyle(1);
            Bef_Flux->SetLineWidth(3);
            Bef_Flux->SetTitle("");
            Bef_Flux->GetXaxis()->SetLimits(0,800);
            Bef_Flux->GetYaxis()->SetLimits(0,2);
            Bef_Flux->Scale(1./2500.);
            cout << "Mass: " << Mass_Point[Mass_INT].c_str() << endl;
            cout << "FILE: " << FILE << endl;

            Aft_Flux_S=(TH1F*)ROOT_FILE_Straight->Get("Flux_HIST_Aft_Collision_EARTH");
            Aft_Flux_S->SetLineColor(2);
            Aft_Flux_S->SetMarkerColor(2);
            Aft_Flux_S->SetLineStyle(1);
            Aft_Flux_S->SetLineWidth(2);
            Aft_Flux_S->SetTitle("");
            Aft_Flux_S->GetXaxis()->SetLimits(0,800);
            Aft_Flux_S->GetYaxis()->SetLimits(0,2);
            Aft_Flux_S->Scale(1./2500.);
                cout << "Mass: " << Mass_Point[Mass_INT].c_str() << endl;
                cout << "FILE: " << FILE << endl;

            Aft_Flux_B=(TH1F*)ROOT_FILE_Bent->Get("Flux_HIST_Aft_Collision_EARTH");
            Aft_Flux_B->SetLineColor(4);
            Aft_Flux_B->SetMarkerColor(2);
            Aft_Flux_B->SetLineStyle(1);
            Aft_Flux_B->SetLineWidth(2);
            Aft_Flux_B->SetTitle("");
            Aft_Flux_B->GetXaxis()->SetLimits(0,800);
            Aft_Flux_B->GetYaxis()->SetLimits(0,2);
            Aft_Flux_B->Scale(1./2500.);
                cout << "Mass: " << Mass_Point[Mass_INT].c_str() << endl;
                cout << "FILE: " << FILE << endl;

            TCanvas *c3 = new TCanvas("c3");
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);


            TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            //leg->AddEntry("",Form("#sigma_{SI}:%.2f #times 10^{-%.f}, M_{#chi}=%.2fGeV",CNFNV(0,sigma_si),CNFNV(1,sigma_si),mx),"");
            leg->AddEntry(Bef_Flux,"Vacuum Case","l");
            leg->AddEntry(Aft_Flux_S,"Straight Line","l");
            leg->AddEntry(Aft_Flux_B,"Bending Line","l");
                cout << "EMass: " << Mass_Point[Mass_INT].c_str() << endl;
                cout << "FILE: " << FILE << endl;

            Bef_Flux->Draw("HIST");
            Aft_Flux_S->Draw("HISTsame");
            Aft_Flux_B->Draw("HISTsame");
            leg->Draw();
                
            c3->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/%sGeV/Velocity_Distribution/%i_Comparison.pdf",Mass_Point[Mass_INT].c_str(),FILE));
            }
        }
    
}
    

