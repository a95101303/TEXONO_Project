#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
//Consider the first bin
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Data_FILE/CDEX_real_Migdal_Lower.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Data_FILE/CDEX_real_Brem_Lower.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Data_FILE/TEXONO_real_Migdal_Lower.h"
//Consider the third Bin, MD>=0.4GeV, BR>=0.2GeV
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Data_FILE/TEXONO_real_Brem_Lower_1st.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Data_FILE/TEXONO_real_Brem_Lower_3rd.h"
//The boundary cases(0.1~0.2GeV)
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Data_FILE/TEXONO_real_Brem_Lower_1st_Bound.h"
#include "/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Data_FILE/TEXONO_real_Brem_Lower_3rd_Bound.h"

void Overlap_Plot_NU_MD_BR_Comparison()
{
    string Exp_Name[2]={"CDEX","TEXONO"};
    string Exp_Flux[2]={"1_CDEX_Flux","2_TEXONO_Flux"};
    //string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    string Mass_Point[2]={"20","0P2"};
    //string Mass_Point[10]={"0P1"};
    for(int Exp=1; Exp<2; Exp++)
    {
        for(int Mass_INT=0; Mass_INT<1; Mass_INT++)
        {
            TGraph *NU; TGraph *MD;TGraph *BR;TGraph *BR_3;
            
    TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                BR=(TGraph*)fin2->Get("Threshold_Plot");
                BR->SetLineColor(2);
                BR->SetLineWidth(2);

    TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_Line_STS.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                BR_3=(TGraph*)fin3->Get("Threshold_Plot");
                BR_3->SetLineColor(4);
                BR_3->SetLineWidth(2);
            
            TCanvas *c3 = new TCanvas("c3");
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);

            TLegend *leg = new TLegend(0.1,0.7,0.4,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);

            BR->Draw("ALP");
            BR_3->Draw("LPsame");

            TF1 *Linear_Line = new TF1("Linear_Line","x",1e-40,1e-28);
            Linear_Line->SetLineColor(3);
            Linear_Line->Draw("Lsame");

            BR->GetXaxis()->SetTitle("real #sigma_{SI} (cm^{2})");
            BR->GetYaxis()->SetTitle("sensitivities #sigma_{SI} (cm^{2})");

            leg->AddEntry("",Form("M_{#chi}=0.2GeV",Mass_Point[Mass_INT].c_str()),"");
            leg->AddEntry(Linear_Line,"Vacuum-constrained line","l");
            leg->AddEntry(BR,"Earth-constrained line with method1","l");
            leg->AddEntry(BR_3,"Earth-constrained line with method2","l");

            leg->Draw();

            c3->SetLogy();
            c3->SetLogx();
            
            c3->Print(Form("All_%sGeV_STS_TEXONO_NU.pdf",Mass_Point[Mass_INT].c_str()));
        }
    }
}
    

