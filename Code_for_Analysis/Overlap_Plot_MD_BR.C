#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLine.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "Data_FILE/CDEX_real_Migdal_Lower.h"
#include "Data_FILE/CDEX_real_Brem_Lower.h"
#include "Data_FILE/TEXONO_real_Migdal_Lower.h"
#include "Data_FILE/TEXONO_real_Brem_Lower.h"
#include "dsigma_dT2.h"

void Overlap_Plot_MD_BR()
{
    string Exp_Name[2]={"CDEX","TEXONO"};
    string Exp_Flux[2]={"1_CDEX_Flux","2_TEXONO_Flux"};
    string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    for(int Exp=1; Exp<2; Exp++)
    {
        for(int Mass_INT=9; Mass_INT<10; Mass_INT++)
        {
            for(int FILE=1; FILE<40; FILE++)
            {
                TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/2_TEXONO_Flux/%sGeV/%i.root",Mass_Point[Mass_INT].c_str(),FILE));
                TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
                Double_t mx,sigma_si;
                T1_TREE->SetBranchAddress("mx",&mx);T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
                T1_TREE->GetEntry(0);

    string path = Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/Recoil_Spectrum/MD_%i.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str(),FILE);
                ifstream fin(path);
    string path1 = Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/Recoil_Spectrum/BR_%i.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str(),FILE);
                ifstream fin1(path1);

            if(fin.is_open() and fin1.is_open()){//Open

                //====================Input=====================//
            TGraphErrors *Data;TGraph *MD_B;TGraph *MD_A;TGraph *BR_B;TGraph *BR_A;
    TFile *finMD = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/Recoil_Spectrum/MD_%i.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str(),FILE));
            MD_B=(TGraph*)finMD->Get("ER_Spectrum_Bef");
            MD_B->SetLineColor(2);
            MD_B->SetMarkerColor(2);
            MD_B->SetLineStyle(1);
            MD_B->SetLineWidth(4);
                MD_B->GetXaxis()->SetRangeUser(0,1e2);
            //MD_B->GetXaxis()->SetLimits(0,2.4);
            MD_B->GetYaxis()->SetRangeUser(1e-8,1e+15);

            MD_A=(TGraph*)finMD->Get("ER_Spectrum_Aft");
            MD_A->SetLineColor(4);
            MD_A->SetMarkerColor(4);
            MD_A->SetLineStyle(1);
            MD_A->SetLineWidth(4);
                
    TFile *finBR = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/Recoil_Spectrum/BR_%i.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str(),FILE));
            BR_B=(TGraph*)finBR->Get("ER_Spectrum_Bef");
            BR_B->SetLineColor(2);
            BR_B->SetMarkerColor(2);
            BR_B->SetLineStyle(2);
            BR_B->SetLineWidth(4);

            BR_A=(TGraph*)finBR->Get("ER_Spectrum_Aft");
            BR_A->SetLineColor(4);
            BR_A->SetMarkerColor(4);
            BR_A->SetLineStyle(2);
            BR_A->SetLineWidth(4);

            double *RE_DATA_1    =Hist_SetLimit_Plot_v2_Extract_Peak(0);
            double *RE_Rate_1    =Hist_SetLimit_Plot_v2_Extract_Peak(1);
            double *RE_DATA_Err_1=Hist_SetLimit_Plot_v2_Extract_Peak(2);
            double *RE_Rate_Err_1=Hist_SetLimit_Plot_v2_Extract_Peak(3);

                
                for(int kkk=0; kkk<257; kkk++)
                {
                    RE_DATA_Err_1[kkk]=0;
                }
                 
            TGraphErrors *TEXONOData = new TGraphErrors(257,RE_DATA_1,RE_Rate_1,RE_DATA_Err_1,RE_Rate_Err_1);
            TEXONOData->SetName("TEXONOData");

            TCanvas *c3 = new TCanvas("c3");
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);

            TLegend *leg = new TLegend(0.4,0.6,0.6,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            //leg->AddEntry("",Form("#sigma_{SI}:%.2f #times 10^{-%.f}, M_{#chi}=%.2fGeV",CNFNV(0,sigma_si),CNFNV(1,sigma_si),mx),"");
            //leg->AddEntry(MD_B,"Migdal Effect before the attenuation","l");
            //leg->AddEntry(MD_A,"Migdal Effect after the attenuation","l");
            //leg->AddEntry(BR_B,"Brem before the attenuation","l");
            //leg->AddEntry(BR_A,"Brem after the attenuation","l");

                c3->SetLogy();
                c3->SetLogx();
                MD_B->SetTitle("");
                MD_B->GetXaxis()->SetTitle("Energy (keVee)");
                MD_B->GetYaxis()->SetTitle("Count (kg^{-1} keV^{-1} day^{-1})");
            MD_B->Draw("ALP");
            MD_A->Draw("LPsame");
            BR_B->Draw("LPsame");
            BR_A->Draw("LPsame");
            TEXONOData->Draw("LPsame");

            leg->Draw();
            c3->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/%s/%sGeV/Recoil_Spectrum/MD_BR_%i.png",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str(),FILE));
            }
            }
        }
    }
}
    

