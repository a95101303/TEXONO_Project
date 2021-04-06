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

void Overlap_Plot_NU_for_CRESST_Bent()
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
            TGraph *NU_STS; TGraph *NU_STS_Bent;TGraph *NU_STS_Earth;TGraph *NU_STS_Earth_Bent;
            
        TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_STS.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS=(TGraph*)fin2->Get("Threshold_Plot");
                NU_STS->SetLineColor(2);
                NU_STS->SetLineStyle(1);

        TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_STS_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS_Bent=(TGraph*)fin3->Get("Threshold_Plot");
                NU_STS_Bent->SetLineColor(2);
                NU_STS_Bent->SetLineStyle(5);
        TFile *fin4 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_STS_Earth.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS_Earth=(TGraph*)fin4->Get("Threshold_Plot");
                NU_STS_Earth->SetLineColor(4);
                NU_STS_Earth->SetLineStyle(1);

TFile *fin5 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_STS_Earth_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS_Earth_Bent=(TGraph*)fin5->Get("Threshold_Plot");
                NU_STS_Earth_Bent->SetLineColor(4);
                NU_STS_Earth_Bent->SetLineStyle(5);

             

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

            /*
            MD->Draw("ALP");
             */
            NU_STS->Draw("ALP");
            NU_STS_Bent->Draw("LPsame");
            NU_STS_Earth->Draw("LPsame");
            NU_STS_Earth_Bent->Draw("LPsame");
            
            leg->Draw();
            leg->AddEntry("",Form("%sGeV",Mass_Point[Mass_INT].c_str()),"l");
            leg->AddEntry("","Solid Line: Without bending","l");
            leg->AddEntry("","Dash  Line: With    bending","l");

            TF1 *Linear_Line = new TF1("Linear_Line","x",1e-40,1e-28);
            Linear_Line->SetLineColor(8);
            Linear_Line->Draw("Lsame");
            cout << "YES7" << endl;

            cout << "YES8" << endl;
            c3->SetLogy();
            cout << "YES9" << endl;
            c3->SetLogx();
            cout << "YES10" << endl;
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/All_%sGeV_STS.png",Mass_Point[Mass_INT].c_str()));
        }
    }
}
    

