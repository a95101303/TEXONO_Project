#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
//Consider the first bin
#include "Data_FILE/CDEX_real_Migdal_Lower.h"
#include "Data_FILE/CDEX_real_Brem_Lower.h"
#include "Data_FILE/TEXONO_real_Migdal_Lower.h"
//Consider the third Bin, MD>=0.4GeV, BR>=0.2GeV
#include "Data_FILE/TEXONO_real_Brem_Lower_1st.h"
#include "Data_FILE/TEXONO_real_Brem_Lower_3rd.h"
//The boundary cases(0.1~0.2GeV)
#include "Data_FILE/TEXONO_real_Brem_Lower_1st_Bound.h"
#include "Data_FILE/TEXONO_real_Brem_Lower_3rd_Bound.h"

void Overlap_Plot_NU_for_CRESST_Bent_Hist()
{
    string Exp_Name[2]={"CDEX","TEXONO"};
    string Exp_Flux[1]={"3_CRESST_Flux"};
    string Mass_Point[3]={"20","2","0P2"};
    //string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P1"};
    //string Mass_Point[12]={"20","19","17","15","13","11","9","7","5","4","3","2P35"};
    //string Mass_Point[10]={"0P1"};
    for(int Exp=0; Exp<1; Exp++)
    {
        for(int Mass_INT=0; Mass_INT<3; Mass_INT++)
        {
            TH1F *Bent_Line; TH1F *Straight_Line;
            
TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/20_STS_Bent_Prove.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                Bent_Line=(TH1F*)fin2->Get("Collision_Time_Hist_Earth");
                Bent_Line->SetLineColor(2);
                Bent_Line->SetLineStyle(1);

TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/20_STS_Bent_Comparison_Prove.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
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
    

