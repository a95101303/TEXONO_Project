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
//The boundary of NU, MD and BR
#include "Data_FILE/TEXONO_real_NU_HM.h"
#include "Data_FILE/TEXONO_real_MD_HM.h"
#include "Data_FILE/TEXONO_real_BR_HM.h"

void Overlap_Plot_NU_MD_BR_Comparison_1()
{
    string Exp_Name[2]={"CDEX","TEXONO"};
    string Exp_Flux[2]={"1_CDEX_Flux","2_TEXONO_Flux"};
    //string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    //string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P1"};
    string Mass_Point[12]={"20","19","17","15","13","11","9","7","5","4","3","2P35"};
    //string Mass_Point[10]={"0P1"};
    for(int Exp=1; Exp<2; Exp++)
    {
        for(int Mass_INT=0; Mass_INT<12; Mass_INT++)
        {
            TGraph *NU; TGraph *MD;TGraph *BR;
    TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            NU=(TGraph*)fin2->Get("Threshold_Plot");
            NU->SetLineColor(2);
            NU->SetLineWidth(2);
    TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/MD_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            MD=(TGraph*)fin3->Get("Threshold_Plot");
            MD->SetLineColor(3);
            MD->SetLineWidth(2);
    TFile *fin4 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/BR_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            BR=(TGraph*)fin4->Get("Threshold_Plot");
            BR->SetLineColor(4);
            BR->SetLineWidth(2);

             

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

            NU->Draw("ALP");
            MD->Draw("LPsame");
            BR->Draw("LPsame");

            TLine *NU_Line;TLine *MD_Line;TLine *BR_Line;
            
            NU_Line = new TLine(0,TEXONO_real_NU_HM[Mass_INT][1],1,TEXONO_real_NU_HM[Mass_INT][1]);
            NU_Line->SetLineColor(2);
            NU_Line->SetLineStyle(2);
            NU_Line->Draw("Lsame");

            MD_Line = new TLine(0,TEXONO_real_MD_HM[Mass_INT][1],1,TEXONO_real_MD_HM[Mass_INT][1]);
            MD_Line->SetLineColor(3);
            MD_Line->SetLineStyle(2);
            MD_Line->SetLineWidth(3);
            MD_Line->Draw("Lsame");

            BR_Line = new TLine(0,TEXONO_real_BR_HM[Mass_INT][1],1,TEXONO_real_BR_HM[Mass_INT][1]);
            BR_Line->SetLineColor(4);
            BR_Line->SetLineStyle(2);
            BR_Line->Draw("Lsame");

            leg->AddEntry("",Form("%sGeV",Mass_Point[Mass_INT].c_str()),"");
            leg->AddEntry(NU_Line,"NU_Lower_Boundary","l");
            leg->AddEntry(MD_Line,"MD_Lower_Boundary","l");
            leg->AddEntry(BR_Line,"BR_Lower_Boundary","l");

            leg->Draw();

            TF1 *Linear_Line = new TF1("Linear_Line","x",1e-40,1e-28);
            Linear_Line->SetLineColor(6);
            Linear_Line->Draw("Lsame");
            
            cout << "YES7" << endl;

            cout << "YES8" << endl;
            c3->SetLogy();
            cout << "YES9" << endl;
            c3->SetLogx();
            cout << "YES10" << endl;
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Comparison_for_three_processes/%s/All_%sGeV_Try.pdf",Exp_Name[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
        }
    }
}
    

