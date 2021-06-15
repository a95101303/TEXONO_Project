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

void Overlap_Plot_NU_MD_BR_Comparison()
{
    string Exp_Name[2]={"CDEX","TEXONO"};
    string Exp_Flux[2]={"1_CDEX_Flux","2_TEXONO_Flux"};
    string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    //string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P1"};
    //string Mass_Point[12]={"20","19","17","15","13","11","9","7","5","4","3","2P35"};
    //string Mass_Point[10]={"0P1"};
    for(int Exp=; Exp<2; Exp++)
    {
        for(int Mass_INT=9; Mass_INT<10; Mass_INT++)
        {
            TGraph *NU; TGraph *MD;TGraph *BR;TGraph *BR_3;
            /*
    TFile *fin = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/NU_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            NU=(TH1F*)fin->Get("Threshold_Plot");
            NU->SetLineColor(2);
             */
            /*
    TFile *fin1 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/MD_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            MD=(TGraph*)fin1->Get("Threshold_Plot");
            MD->SetLineColor(1);
             */
            /*
    TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/BR_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            BR=(TGraph*)fin2->Get("Threshold_Plot");
            BR->SetLineColor(2);
    TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/BR_Line_Bin3.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
            BR_3=(TGraph*)fin3->Get("Threshold_Plot");
            BR_3->SetLineColor(3);
            */
            
        TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/BR_Line.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                BR=(TGraph*)fin2->Get("Threshold_Plot");
                BR->SetLineColor(2);
        TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/BR_Line_STS.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                BR_3=(TGraph*)fin3->Get("Threshold_Plot");
                BR_3->SetLineColor(3);

             

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
            BR->Draw("ALP");
            BR_3->Draw("LPsame");


            /*
            TLine *Fast_NU = new TLine(2e-32,0,2e-32,1);
            Fast_NU->SetLineColor(42);
            Fast_NU->Draw("Lsame");
            TLine *Fast_MD = new TLine(3.4e-32,0,3.4e-32,1);
            Fast_MD->SetLineColor(46);
            Fast_MD->Draw("Lsame");
             */
            TLine *MD_Line_STFlux;TLine *BR_Line_STFlux;
            TLine *BR_Line_STFlux_1st;TLine *BR_Line_STFlux_3rd;
            TLine *BR_Line_STFlux_1st_Bound;TLine *BR_Line_STFlux_3rd_Bound;
            
            if(Exp==0)
            {
                MD_Line_STFlux= new TLine(0,CDEX_real_Migdal_Down[Mass_INT][1],1,CDEX_real_Migdal_Down[Mass_INT][1]);
                MD_Line_STFlux->SetLineColor(45);
                //MD_Line_STFlux->Draw("Lsame");

                BR_Line_STFlux = new TLine(0,CDEX_real_Brem_Lower[Mass_INT][1],1,CDEX_real_Brem_Lower[Mass_INT][1]);
                BR_Line_STFlux->SetLineColor(41);
                //BR_Line_STFlux->Draw("Lsame");
            }
            if(Exp==1)
            {
                /*
                MD_Line_STFlux= new TLine(0,TEXONO_real_Migdal_Lower[Mass_INT][1],1,TEXONO_real_Migdal_Lower[Mass_INT][1]);
                MD_Line_STFlux->SetLineColor(1);
                MD_Line_STFlux->SetLineStyle(2);
                MD_Line_STFlux->Draw("Lsame");
                 */
                /*
                BR_Line_STFlux = new TLine(0,TEXONO_real_Brem_Lower[Mass_INT][1],1,TEXONO_real_Brem_Lower[Mass_INT][1]);
                BR_Line_STFlux->SetLineColor(3);
                BR_Line_STFlux->Draw("Lsame");
                */
                /*
                BR_Line_STFlux_1st= new TLine(0,TEXONO_real_Brem_Lower_1st[Mass_INT][1],1,TEXONO_real_Brem_Lower_1st[Mass_INT][1]);
                BR_Line_STFlux_1st->SetLineColor(2);
                BR_Line_STFlux_1st->SetLineStyle(2);
                BR_Line_STFlux_1st->Draw("Lsame");

                BR_Line_STFlux_3rd = new TLine(0,TEXONO_real_Brem_Lower_3rd[Mass_INT][1],1,TEXONO_real_Brem_Lower_3rd[Mass_INT][1]);
                BR_Line_STFlux_3rd->SetLineColor(3);
                BR_Line_STFlux_3rd->SetLineStyle(2);
                BR_Line_STFlux_3rd->Draw("Lsame");
                 */
                BR_Line_STFlux_1st_Bound= new TLine(0,TEXONO_real_Brem_Lower_1st_Bound[Mass_INT][1],1,TEXONO_real_Brem_Lower_1st_Bound[Mass_INT][1]);
                BR_Line_STFlux_1st_Bound->SetLineColor(2);
                BR_Line_STFlux_1st_Bound->SetLineStyle(2);
                //BR_Line_STFlux_1st_Bound->Draw("Lsame");

                BR_Line_STFlux_3rd_Bound = new TLine(0,TEXONO_real_Brem_Lower_3rd_Bound[Mass_INT][1],1,TEXONO_real_Brem_Lower_3rd_Bound[Mass_INT][1]);
                BR_Line_STFlux_3rd_Bound->SetLineColor(3);
                BR_Line_STFlux_3rd_Bound->SetLineStyle(2);
                //BR_Line_STFlux_3rd_Bound->Draw("Lsame");


            }

            //leg->AddEntry(MD,"MD","l");
            //leg->AddEntry(MD_Line_STFlux,"MD before the attenuation(1st)","l");

            /*
            leg->AddEntry(BR,"BR(1st)","l");
            leg->AddEntry(BR_3,"BR(3rd)","l");
            
            leg->AddEntry(BR_Line_STFlux_1st,"BR before the attenuation(1st)","l");
            leg->AddEntry(BR_Line_STFlux_3rd,"BR before the attenuation(3rd)","l");
             
            leg->AddEntry(BR_Line_STFlux_1st_Bound,"BR before the attenuation(1st)","l");
            leg->AddEntry(BR_Line_STFlux_3rd_Bound,"BR before the attenuation(3rd)","l");
             */
            leg->AddEntry(BR,"BR_M1","l");
            leg->AddEntry(BR_3,"BR_M2","l");

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
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/Comparison_for_three_processes/%s/All_%sGeV_Try.pdf",Exp_Name[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
        }
    }
}
    

