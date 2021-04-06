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

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"

vector<double> Rebin_two(vector<double> ValueX,vector<double> ValueY, vector<double> Sigma)
{
    vector<double> Re_V;
    double Event_Number[2] = { 1/(Sigma[0]*Sigma[0]),1/(Sigma[1]*Sigma[1]) };
    double Total_Event     = Event_Number[0] + Event_Number[1] ;
    double Weighted_ValueX = (Event_Number[0]*ValueX[0]+Event_Number[1]*ValueX[1])/(Total_Event);
    double Weighted_ValueY = (Event_Number[0]*ValueY[0]+Event_Number[1]*ValueY[1])/(Total_Event);
    double Weighted_Error =   std::sqrt(1/(Total_Event));
    Re_V.push_back(Weighted_ValueX);
    Re_V.push_back(Weighted_ValueY);
    Re_V.push_back(Weighted_Error);
    return Re_V;
}


void Hist_SetLimit_Plot_Line()
{

const double kms1_to_cmday1 = 100000.0*86400.0; // km/s to cm/day
const double kms1_to_c = 1000.0/2.99792458e8; // km/s to natural unit
const double MeV1_to_keV1 = 1.0/1000.0; // per MeV to per keV
double sum;
    
for(int j=0;j<2000;j++)
{
    sum = sum + velo_dist_Ave[j][3];
}
TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.05,"XY");
gStyle->SetTitleFont(62,"XY");
gStyle->SetLegendFont(62);
const int Number_mx_Candidate = 20;
double A = AGe;
const int reso_T = 1000;
double T[Number_mx_Candidate][reso_T];
double T_QF[Number_mx_Candidate][reso_T];
double recoil[Number_mx_Candidate][reso_T];
double v;
double WIMP_mx_TEXONO[Number_mx_Candidate];
double WIMP_mx_CDEX[Number_mx_Candidate];
double WIMP_max_T[Number_mx_Candidate];
    // Neutrino_Mass_series
double Initial_point[5] = {1e-1,1e0,1e1,1e2,1e3};
    
    WIMP_mx_TEXONO[0]=2.34;//Threshold
    WIMP_mx_CDEX[0]=2.07;//Threshold

    for(int j=1 ; j<20 ; j++)
    {
        WIMP_mx_TEXONO[j] = j+2;
        WIMP_mx_CDEX[j] = j+2;
    }
    
/*for(int kkk=0 ; kkk<5 ; kkk++)
{
    for(int jjj=1 ; jjj<10 ; jjj++)
    {
        WIMP_mx[kkk*10+jjj-kkk-1] = (Initial_point[kkk]*jjj);
        cout << "kkk*10+jjj-kkk-1; " << kkk*10+jjj-kkk-1 << endl;
    }
}*/

//double mx = 50.0; // GeV
//========For Data
const int Data_element = 257;
double qf0, qf1, qf2, qf3, qf4, alpha=1.0;
qf0 = alpha*0.19816;
qf1 = alpha*0.05052;
qf2 = alpha*0.00378;
qf3 = alpha*0.00192;
qf4 = alpha*0.0016;
double Scaling_Factor_TEXONO_1[Number_mx_Candidate];
double Scaling_Factor_TEXONO[Number_mx_Candidate];
double Scaling_Factor_CDEX[Number_mx_Candidate];

//=========
// Data_implemented_no_rebin
/*
for(int jjj=0 ; jjj<257 ; jjj++)
{
    RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
    RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
    RE_DATA_Err[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][3]*1.64458/0.994;
    RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2];
}*/
// Data_implemented_rebi n_2

double RE_DATA[Data_element]; double RE_Rate[Data_element];
double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
    /*
for(int jjj=0 ; jjj<257 ; jjj++)
{
    if(jjj==0){
        RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0]-0.025;
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994;}
    else if(jjj==256){
        RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0]+0.025;
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994;}
    else{
        RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err[jjj]= 0.025;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994;}
}*/
    
    double *RE_DATA_1    =Hist_SetLimit_Plot_v2_Extract_Peak(0);
    double *RE_Rate_1    =Hist_SetLimit_Plot_v2_Extract_Peak(1);
    double *RE_DATA_Err_1=Hist_SetLimit_Plot_v2_Extract_Peak(2);
    double *RE_Rate_Err_1=Hist_SetLimit_Plot_v2_Extract_Peak(3);
        

    for(int jjj=0 ; jjj<257 ; jjj++)
    {
        if(jjj==0){
            RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
            RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
            RE_DATA_Err[jjj]= 0;
            RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994;}
        else if(jjj==256){
            RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
            RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
            RE_DATA_Err[jjj]= 0;
            RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994;}
        else{
            RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
            RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
            RE_DATA_Err[jjj]= 0.025;
            RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994;}
    }
    double RE_DATA_CDEX[4]= {0.185,0.235,0.285,0.325};
    double RE_Rate_CDEX[4]= {8.2,7.6,7.4,7.9};
    double RE_DATA_Err_CDEX[4]= {0,0,0,0};
    double RE_Rate_Err_CDEX[4]= {3,3.5,1.5,2};

    for(int jjj=0; jjj<Number_mx_Candidate; jjj++)
    {
        //double *Factor_TEXONO=0;
        cout << "WIMP_mx_TEXONO[jjj]: " << WIMP_mx_TEXONO[jjj] << endl;
        cout << "WIMP_mx_CDEX[jjj]: " << WIMP_mx_CDEX[jjj] << endl;
        double *Factor_TEXONO=0;double *Factor_TEXONO_1=0;double *Factor_CDEX=0;
        
        Factor_TEXONO = cpkkd_calculation_Scaling_Factor(257,RE_DATA,RE_Rate,RE_DATA_Err,RE_Rate_Err,WIMP_mx_TEXONO[jjj]);
        Scaling_Factor_TEXONO[jjj]=Factor_TEXONO[0]*1e-40;
        //cout << "Factor_TEXONO[0]: " << Factor_TEXONO[0]*1e-40 << endl;
        //cout << "BIN_DET_TEXONO: " << Factor_TEXONO[1] << endl;

        Factor_TEXONO_1 = cpkkd_calculation_Scaling_Factor(257,RE_DATA_1,RE_Rate_1,RE_DATA_Err_1,RE_Rate_Err_1,WIMP_mx_TEXONO[jjj]);
        Scaling_Factor_TEXONO_1[jjj]=Factor_TEXONO_1[0]*1e-40;
        cout << "Factor_TEXONO_1[0]: " << Factor_TEXONO_1[0]*1e-40 << endl;
        cout << "BIN_DET_TEXONO_1: " << Factor_TEXONO_1[1] << endl;
        
        Factor_CDEX = cpkkd_calculation_Scaling_Factor(4,RE_DATA_CDEX,RE_Rate_CDEX,RE_DATA_Err_CDEX,RE_Rate_Err_CDEX,WIMP_mx_CDEX[jjj]);
        Scaling_Factor_CDEX[jjj]  =Factor_CDEX[0]*1e-40;
        //cout << "Factor_CDEX[0]: " << Factor_CDEX[0]*1e-40 << endl;
        //cout << "BIN_DET_CDEX: " << Factor_CDEX[1] << endl;
    }

     
    
    
        //===============================================================

    cout << "3: " << endl;
char fout_name[100];
sprintf(fout_name,"Set_Limits_Plot_for_ALL.root");
TFile *fout=new TFile(fout_name,"recreate");
    cout << "4: " << endl;

TGraph *Set_Limits_plot = new TGraph(Number_mx_Candidate,WIMP_mx_TEXONO,Scaling_Factor_TEXONO);
Set_Limits_plot->SetName("TEXONO_Without_Subtracting_the_KL_shells");
Set_Limits_plot->SetLineColor(2);
Set_Limits_plot->SetTitle("Set_Limits_plot");
Set_Limits_plot->GetXaxis()->SetTitle("m_{#chi} (GeV/c^{2})");
Set_Limits_plot->GetYaxis()->SetTitle("#sigma_{SI} (cm^{2})");
    Set_Limits_plot->GetXaxis()->SetRangeUser(1.9,22);
Set_Limits_plot->GetYaxis()->SetRangeUser(1e-42,1e-30);
Set_Limits_plot->SetMarkerColor(2);
Set_Limits_plot->SetMarkerStyle(8);

    TGraph *Set_Limits_plot_2 = new TGraph(Number_mx_Candidate,WIMP_mx_TEXONO,Scaling_Factor_TEXONO_1);
    Set_Limits_plot_2->SetName("TEXONO_With_Subtracting_the_KL_shells");
    Set_Limits_plot_2->SetLineColor(4);
    Set_Limits_plot_2->SetMarkerColor(4);
    Set_Limits_plot_2->SetMarkerStyle(8);

    TGraph *Set_Limits_plot_1 = new TGraph(Number_mx_Candidate,WIMP_mx_CDEX,Scaling_Factor_CDEX);
    Set_Limits_plot_1->SetName("CDEX-1a");
    Set_Limits_plot_1->SetLineColor(3);
    Set_Limits_plot_1->SetMarkerColor(3);
    Set_Limits_plot_1->SetMarkerStyle(8);


    cout << "5: " << endl;

TLegend *leg= new TLegend(0.1,0.5,0.6,0.9);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.06);
leg->SetBorderSize(0);
leg->SetTextFont(20);
    
leg->AddEntry("","Standard Flux of WIMP","");
leg->AddEntry("","TEXONO(2020)","");
leg->AddEntry(Set_Limits_plot,"Without Subtracting the K/L shells","lp");
leg->AddEntry(Set_Limits_plot_2,"With Subtracting the K/L shells","lp");
leg->AddEntry(Set_Limits_plot_1,"CDEX-1a","lp");

    
    cout << "6: " << endl;

Set_Limits_plot->Draw("aplE");
    Set_Limits_plot_2->Draw("plEsame");
    Set_Limits_plot_1->Draw("plEsame");

leg->Draw();
    cout << "7: " << endl;

Set_Limits_plot->Write();
    Set_Limits_plot_1->Write();
    Set_Limits_plot_2->Write();

    c1->SetLogy();
    c1->SetLogx();
c1->Print("Set_Limits_Plot_for_ALL.png");
    cout << "8: " << endl;

    
    //
//68% (1sigma) => 90% (2sigma)
}
