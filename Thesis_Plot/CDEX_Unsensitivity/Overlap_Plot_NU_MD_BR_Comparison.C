#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
//Consider the first bin
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Data_FILE/MD_CDEX/CDEX_real_Migdal_Lower.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Data_FILE/CDEX_real_Brem_Lower.h"

void Overlap_Plot_NU_MD_BR_Comparison()
{
    string Exp_Name[2]={"CDEX","TEXONO"};
    string Exp_Flux[2]={"1_CDEX_Flux","2_TEXONO_Flux"};
    string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    //string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P1"};
    //string Mass_Point[12]={"20","19","17","15","13","11","9","7","5","4","3","2P35"};
    //string Mass_Point[10]={"0P1"};
    for(int Exp=0; Exp<1; Exp++)
    {
        for(int Mass_INT=11; Mass_INT<12; Mass_INT++)
        {
            TGraph *NU; TGraph *MD;TGraph *BR;TGraph *BR_3;
            
    TFile *fin1 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/MD_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            MD=(TGraph*)fin1->Get("Threshold_Plot_Bin_Possion");
            MD->GetXaxis()->SetRangeUser(1e-42,1e-27);
            MD->GetYaxis()->SetRangeUser(1e-42,1e-27);
            MD->GetXaxis()->SetTitle("real #sigma_{SI} (cm^{2})");
            MD->GetYaxis()->SetTitle("sensitivities #sigma_{SI} (cm^{2})");

            MD->SetLineColor(1);
     TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/BR_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
             BR=(TGraph*)fin2->Get("Threshold_Plot_Bin_Possion");
             BR->SetLineColor(2);

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

            c3->SetLogy();
            c3->SetLogx();

            TF1 *Linear_Line = new TF1("Linear_Line","x",1e-40,1e-28);
            Linear_Line->SetLineColor(4);

            MD->Draw("ALP");
            BR->Draw("LPsame");
            Linear_Line->Draw("Lsame");


            TLine *MD_Line_STFlux;TLine *BR_Line_STFlux;
            TLine *BR_Line_STFlux_1st;TLine *BR_Line_STFlux_3rd;
            TLine *BR_Line_STFlux_1st_Bound;TLine *BR_Line_STFlux_3rd_Bound;
            
            
            if(Exp==0)
            {
                MD_Line_STFlux= new TLine(0,CDEX_real_Migdal_Lower[Mass_INT][1],1e-27,CDEX_real_Migdal_Lower[Mass_INT][1]);
                MD_Line_STFlux->SetLineColor(1);
                MD_Line_STFlux->SetLineStyle(5);
                MD_Line_STFlux->Draw("Lsame");

                BR_Line_STFlux = new TLine(0,CDEX_real_Brem_Lower[Mass_INT][1],1e-27,CDEX_real_Brem_Lower[Mass_INT][1]);
                BR_Line_STFlux->SetLineColor(2);
                BR_Line_STFlux->SetLineStyle(5);
                BR_Line_STFlux->Draw("Lsame");
            }
            leg->AddEntry("",Form("M_{#chi}=0.09GeV",Mass_Point[Mass_INT].c_str()),"");
            leg->AddEntry(MD,"Earth-constrained MD","l");
            leg->AddEntry(BR,"Earth-constrained Brem","l");
            leg->AddEntry(MD_Line_STFlux,"Vaccum-constrained MD","l");
            leg->AddEntry(BR_Line_STFlux,"Vaccum-constrained Brem","l");
            leg->Draw();
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Thesis_Plot/CDEX_Unsensitivity/All_%sGeV_Try.pdf",Mass_Point[Mass_INT].c_str()));

            /*
            leg->AddEntry(BR,"MD","l");
            leg->AddEntry(BR_3,"BR","l");

            leg->Draw();

            TF1 *Linear_Line = new TF1("Linear_Line","x",1e-40,1e-28);
            Linear_Line->SetLineColor(4);
            Linear_Line->Draw("Lsame");
            cout << "YES7" << endl;

            cout << "YES8" << endl;
            c3->SetLogy();
            cout << "YES9" << endl;
            c3->SetLogx();
            cout << "YES10" << endl;
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Thesis_Plot/CDEX_Unsensitivity/All_%sGeV_Try.pdf",Mass_Point[Mass_INT].c_str()));
             */
        }
    }
}
    

