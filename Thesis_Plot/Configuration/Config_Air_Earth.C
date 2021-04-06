#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <math.h>


#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

void Config_Air_Earth()
{
    double atm_table[20][6] = {
        0, 15.00, 9.807, 10.13, 1.225, 1.789,
        1000, 8.50, 9.804, 8.988, 1.112, 1.758,
        2000, 2.00, 9.801, 7.950, 1.007, 1.726,
        3000, -4.49, 9.797, 7.012, 0.9093, 1.694,
        4000, -10.98, 9.794, 6.166, 0.8194, 1.661,
        5000, -17.47, 9.791, 5.405, 0.7364, 1.628,
        6000, -23.96, 9.788, 4.722, 0.6601, 1.595,
        7000, -30.45, 9.785, 4.111, 0.5900, 1.561,
        8000, -36.94, 9.782, 3.565, 0.5258, 1.527,
        9000, -43.42, 9.779, 3.080, 0.4671, 1.493,
        10000, -49.90, 9.776, 2.650, 0.4135, 1.458,
        15000, -56.50, 9.761, 1.211, 0.1948, 1.422,
        20000, -56.50, 9.745, 0.5529, 0.08891, 1.422,
        25000, -51.60, 9.730, 0.2549, 0.04008, 1.448,
        30000, -46.64, 9.715, 0.1197, 0.01841, 1.475,
        40000, -22.80, 9.684, 0.0287, 0.003996, 1.601,
        50000, -2.5, 9.654, 0.007978, 0.001027, 1.70,
        60000, -26.13, 9.624, 0.002196, 0.0003097, 1.584,
        70000, -53.57, 9.594, 0.00052, 0.00008283, 1.438,
        80000, -74.51, 9.564, 0.00011, 0.00001846, 1.321
    };

    double earth_table[14][6] = {
        0     , 13.0885, 0.0,  0.0   , 0.0, 0.0,
        1221.5, 12.760, 0.0, -8.8381, 0.0, 12.893569,
        3480.0, 9.9, -1.2638, -3.6426, -5.5281, 10.900696,
        3630.0, 5.50, -6.4761, 5.5283, -3.0807, 5.528405,
        5600.0, 4.44, -6.4761, 5.5283, -3.0807, 4.912992,
        5701.0, 4.38, -6.4761, 5.5283, -3.0807, 4.411899,
        5771.0, 3.97, -1.4836, 0.0, 0.0, 3.983938,
        5971.0, 3.72, -8.0298, 0.0, 0.0, 3.848353,
        6151.0, 3.43, -3.8045, 0.0, 0.0, 3.488987,
        6291.0, 3.37, 0.6924, 0.0, 0.0, 3.367155,
        6346.6, 3.38, 0.6924, 0.0, 0.0, 3.377736,
        6356.0, 2.9, 0.0, 0.0, 0.0, 2.900000,
        6368.0, 2.6, 0.0, 0.0, 0.0, 2.600000,
        6371.0, 1.02, 0.0, 0.0, 0.0, 1.020000
    };

    TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.035,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);
         
    double Deep_Earth[14];double Density_Earth[14];
    double Deep_ATM[20];  double Density_ATM[20];

    for(int kkk=0; kkk<20; kkk++)
    {
        Deep_ATM[kkk]    = atm_table[kkk][0];
        Density_ATM[kkk] = atm_table[kkk][4];
    }
    for(int kkk=0; kkk<14; kkk++)
    {
        Deep_Earth[kkk]    = earth_table[kkk][0];
        Density_Earth[kkk] = earth_table[kkk][1];
    }


    TGraph *ATM_Model = new TGraph(20,Deep_ATM,Density_ATM);
    ATM_Model->SetLineColor(2);
    ATM_Model->SetTitle("");
    ATM_Model->GetXaxis()->SetTitle("R(km)");
    ATM_Model->GetYaxis()->SetTitle("kg/m^{3}");
    ATM_Model->GetXaxis()->SetRangeUser(0,6372);
    ATM_Model->GetYaxis()->SetRangeUser(0,15);
    ATM_Model->SetMarkerColor(2);
    ATM_Model->SetMarkerStyle(8);
    //ATM_Model->Draw("ALP");

    
    TGraph *Earth_Model = new TGraph(14,Deep_Earth,Density_Earth);
    Earth_Model->SetLineColor(2);
    Earth_Model->SetTitle("");
    Earth_Model->GetXaxis()->SetTitle("R(km)");
    Earth_Model->GetYaxis()->SetTitle("kg/m^{3}");
    Earth_Model->GetXaxis()->SetRangeUser(0,80000);
    Earth_Model->GetYaxis()->SetRangeUser(0,15);
    Earth_Model->SetMarkerColor(2);
    Earth_Model->SetMarkerStyle(8);
    Earth_Model->Draw("ALP");
    
    TLegend *leg= new TLegend(0.6,0.1,0.9,0.4);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.06);
    leg->SetBorderSize(0);
    leg->SetTextFont(20);

    leg->Draw();
    c1->SetLogy();
    c1->Print("Earth_Model.pdf");
}
