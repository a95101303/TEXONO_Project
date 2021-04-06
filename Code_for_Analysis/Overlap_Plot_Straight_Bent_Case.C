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

void Overlap_Plot_Straight_Bent_Case()
{
        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

        TLegend *leg = new TLegend(0.1,0.1,0.4,0.3);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);

        TH2F *Straight_Bent_comparison;
    
        TFile *fin2 = TFile::Open("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/2_TEXONO_Flux/20GeV/27_STS_Bent.root");
        Straight_Bent_comparison=(TH2F*)fin2->Get("Extended_Length_AIR_Bent_vs_Straight_Line");
        Straight_Bent_comparison->SetLineColor(2);

        Straight_Bent_comparison->GetXaxis()->SetRange(0,1000);
        Straight_Bent_comparison->GetYaxis()->SetRange(0,1000);

        Straight_Bent_comparison->Draw("COLZ");

        TF1 *Linear_Line = new TF1("Linear_Line","x",0,1000);
        Linear_Line->SetLineColor(4);
    
        Linear_Line->Draw("Lsame");

        c3->Print("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/2_TEXONO_Flux/20GeV/Straight_Bent_2D.pdf");
}
    

