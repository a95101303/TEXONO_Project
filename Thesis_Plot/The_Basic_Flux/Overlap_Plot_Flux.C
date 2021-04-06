#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLine.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/dsigma_dT2.h"

void Overlap_Plot_Flux()
{
    string Exp_Name[2]={"CDEX","TEXONO"};
    string Exp_Flux[2]={"1_CDEX_Flux","2_TEXONO_Flux"};
    string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    
    for(int Exp=1; Exp<2; Exp++)
    {
        for(int Mass_INT=10; Mass_INT<11; Mass_INT++)
        {
            for(int FILE=1; FILE<40; FILE++)
            {
                TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/2_TEXONO_Flux/%sGeV/%i.root",Mass_Point[Mass_INT].c_str(),FILE));
                TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
                Double_t mx,sigma_si;
                T1_TREE->SetBranchAddress("mx",&mx);T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
                T1_TREE->GetEntry(0);

                //====================Input=====================//
            TH1F *Bef_Flux;TH1F *Aft_Flux;
            Bef_Flux=(TH1F*)ROOT_FILE->Get("Flux_HIST_Random");
            Bef_Flux->SetLineColor(3);
            Bef_Flux->SetMarkerColor(3);
            Bef_Flux->SetLineStyle(1);
            Bef_Flux->SetLineWidth(3);
            Bef_Flux->SetTitle("");
            Bef_Flux->GetXaxis()->SetLimits(0,800);
            Bef_Flux->GetYaxis()->SetLimits(0,2);
            Bef_Flux->Scale(1./2500.);
            
            Aft_Flux=(TH1F*)ROOT_FILE->Get("Flux_HIST_Aft_Collision_EARTH");
            Aft_Flux->SetLineColor(2);
            Aft_Flux->SetMarkerColor(2);
            Aft_Flux->SetLineStyle(1);
            Aft_Flux->SetLineWidth(2);
            Aft_Flux->SetTitle("");
            Aft_Flux->GetXaxis()->SetLimits(0,800);
                Aft_Flux->GetYaxis()->SetLimits(0,2);
            Aft_Flux->Scale(1./2500.);

            TCanvas *c3 = new TCanvas("c3");
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);

            TLine *MD_Line_STFlux= new TLine(134,0,134,1);
            MD_Line_STFlux->SetLineStyle(9);
            MD_Line_STFlux->SetLineColor(45);
            MD_Line_STFlux->SetLineWidth(5);

            TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            //leg->AddEntry("",Form("#sigma_{SI}:%.2f #times 10^{-%.f}, M_{#chi}=%.2fGeV",CNFNV(0,sigma_si),CNFNV(1,sigma_si),mx),"");
            leg->AddEntry(Bef_Flux,"Vacuum Case","l");
            leg->AddEntry(Aft_Flux,"Earth Effect Case","l");
            leg->AddEntry(MD_Line_STFlux,"Threshold","l");

            Bef_Flux->Draw("HIST");
            Aft_Flux->Draw("HISTsame");
            MD_Line_STFlux->Draw("Lsame");
            leg->Draw();
                
            c3->Print(Form("%i_Used.pdf",FILE));
            }
        }
    }
}
    

