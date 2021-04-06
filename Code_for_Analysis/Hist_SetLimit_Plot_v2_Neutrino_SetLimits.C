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

#include "/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"
#include "/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/dsigma_dT2.h"
#include "/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/cpkkd_chi_N_v3.h"
#include "/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/velocity_distribution_2000_Ave.h"
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


void Hist_SetLimit_Plot_v2_Neutrino_SetLimits()
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
const int Number_mx_Candidate = 45;
double A = AGe;
const int reso_T = 1000;
double T[Number_mx_Candidate][reso_T];
double T_QF[Number_mx_Candidate][reso_T];
double recoil[Number_mx_Candidate][reso_T];
double v;
double WIMP_mx[Number_mx_Candidate];
double WIMP_max_T[Number_mx_Candidate];
double Limitation[Number_mx_Candidate];
    
// Neutrino_Mass_series
double Initial_point[5] = {1e-1,1e0,1e1,1e2,1e3};
for(int kkk=0 ; kkk<5 ; kkk++)
{
    for(int jjj=1 ; jjj<10 ; jjj++)
    {
        if(kkk==0 and (jjj==1 or jjj==2)){WIMP_mx[kkk*10+jjj-kkk-1] = 0.3;}
        else{WIMP_mx[kkk*10+jjj-kkk-1] = (Initial_point[kkk]*jjj);}
        
        //cout << "kkk*10+jjj-kkk-1; " << kkk*10+jjj-kkk-1 << endl;
    }
}

// Find the maximum on the dark matter energy recoil
    for(int jjj=0; jjj<Number_mx_Candidate ; jjj++)
{
    WIMP_max_T[jjj] = 1000.0*max_recoil_A(WIMP_mx[jjj], 779.135*1000.0/2.99792458e8, A)+1.8; //keV
}
    
for(int jjj=0 ; jjj<100 ; jjj++)
{
        //cout << "2+0.01*jjj: " << 2+0.01*jjj << endl;
        //cout << "TQF(WIMP_max_T[jjj]): " << TQF(1000.0*max_recoil_A(2+0.01*jjj, 779.135*1000.0/2.99792458e8, A)) << endl;
}
//double mx = 50.0; // GeV
//========For Data
const int Data_element = 257;
const int Neutrino_Data_element = 5000;
//Dark matter spectrum
double RE_DATA[Data_element]    ; double RE_Rate[Data_element];
double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
//Neutrino Spectrum
double RE_Nu_DATA[Neutrino_Data_element]    ; double RE_Nu_Rate[Neutrino_Data_element];
double RE_Nu_DATA_Err[Neutrino_Data_element]; double RE_Nu_Rate_Err[Neutrino_Data_element];
    double RE_Nu_Bin_Number[Neutrino_Data_element];
    
double qf0, qf1, qf2, qf3, qf4, alpha=1.0;
qf0 = alpha*0.19816;
qf1 = alpha*0.05052;
qf2 = alpha*0.00378;
qf3 = alpha*0.00192;
qf4 = alpha*0.0016;

//===============================
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
//===============================
TGraph *Neutrino = new TGraph("Try_1.txt", "%lg  %lg");
Neutrino->SetName("Neutrino");
Neutrino->SetLineWidth(5);
Neutrino->SetLineColor(8);
//===============================
Double_t xp=0; Double_t yp=0;
//cout << "TGraph: " <<     Neutrino->GetPoint(1,xp,yp);
//cout << "xp: " << xp << "yp: " << yp << endl;
//===============================// Set the data for Neutrino
for(int jjjj=0 ; jjjj<Neutrino_Data_element ; jjjj++)
{
    double xp=0; double yp=0;
    Neutrino->GetPoint(jjjj,xp,yp);
    //cout << "==============================================" << endl;
    //cout << "Real_xp: " << xp << endl;
    //cout << "Real_yp: " << yp << endl;
    //RE_Nu_DATA[jjjj] = xp;
    //RE_Nu_Rate[jjjj] = yp;
    //cout << "==============================================" << endl;
    
    for(int jjj=0 ; jjj<2000 ; jjj++)
    {
        if( (xp>=0.0000+jjj*0.05) and (xp<0.0500+jjj*0.05) )
        {
            cout << "xp: " << 0.025+jjj*0.05 << endl;
            cout << "yp: " << yp << endl;
            RE_Nu_DATA[jjj] = 0.025+jjj*0.05;
            RE_Nu_Rate[jjj] = RE_Nu_Rate[jjj] + yp;
            cout << "RE_Nu_Rate[jjj]: " << RE_Nu_Rate[jjj] << endl;
            RE_Nu_Bin_Number[jjj] = RE_Nu_Bin_Number[jjj]+1;
        }
    }
}
    for(int jjj=0 ; jjj<2000 ; jjj++)
    {
        if(RE_Nu_Bin_Number[jjj]!=0){RE_Nu_Rate[jjj] = RE_Nu_Rate[jjj]/RE_Nu_Bin_Number[jjj];}
        else{RE_Nu_Rate[jjj]=0;}
    }
    
//===============================//
    
for(int jjj=0 ; jjj<257 ; jjj++)
{
    if(jjj==0){
        RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*(1.64458/0.994);}
    else if(jjj==256){
        RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*(1.64458/0.994);}
    else{
        RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*(1.64458/0.994);}
}
    
//===============================
int Combine[6] = {12,13,14,23,24,34};
double Constant[6];
double ConstantX[6];

    //Min_Point_found
    //===============================================================
        //Theoretical line of the rate of WIMP Method 2
        double recoilX[Number_mx_Candidate][reso_T];
        for(int jjj=0 ; jjj<Number_mx_Candidate ; jjj++)
        {
            for(int i=0;i<reso_T;i++)
            {
                T[jjj][i] = ((double)i+0.5)*((WIMP_max_T[jjj])/(double)reso_T); // keV
                T_QF[jjj][i] = TQF(T[jjj][i]);
                recoilX[jjj][i] = 0.0;
                for(int j=0;j<2000;j++)
                {
                    float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_c;
                    float v_cm_day = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_cmday1;
                    
                    if((max_recoil_A_keV(WIMP_mx[jjj], v, A))>T[jjj][i])
                    {
                        recoilX[jjj][i] = recoilX[jjj][i] + rate_scale_QF(T[jjj][i])*fdsigma_dT_keV(WIMP_mx[jjj], 1e-40, v, A, T[jjj][i])
                        *N_atom_Ge_1kg*(rohx/WIMP_mx[jjj])*v_cm_day*(1/(sum))*velo_dist_Ave[j][3];
                    }
                }
            }
    
        }
    //===============================================================
        //Theoretical line of the rate of WIMP Method 1
    for(int jjj=0 ; jjj<Number_mx_Candidate ; jjj++)
    {
        for(int i=0;i<reso_T;i++)
        {
            recoil[jjj][i] = 5*cpkkd_chi_N(WIMP_mx[jjj], T[jjj][i], 243.8125, 152.5, A);
        }
        
    }
    //===============================================================
    int dm_spec_resolution=1000;
    double cpkkd_expect_e_reso0[Number_mx_Candidate][dm_spec_resolution];
    double cpkkd_expect_e_reso1[Number_mx_Candidate][dm_spec_resolution];
    double cpkkd_expect_e_reso_C1A0[dm_spec_resolution];
    double cpkkd_expect_e_reso_C1A1[dm_spec_resolution];

    double diff0[dm_spec_resolution];
    double diff[dm_spec_resolution];
    double diff_C1A[dm_spec_resolution];

    double sig_E; double dEx;
    double a0_dE = 33.4992/1000.0;
    double a1_dE = 13.2145/1000.0;
    
    double Factor[Number_mx_Candidate][Neutrino_Data_element];
    double Factor1[Number_mx_Candidate][1000];
    double Factor2[Number_mx_Candidate][1000];

    float  Boundary=0;
    //===============================================================
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        double pEx = T_QF[kkk][0];
        //Find out the Full energy spectrum with the resolution of the detector based on the theory
        for(int jjj=0 ; jjj<dm_spec_resolution ; jjj++)
        {
            for(int j=0 ; j<dm_spec_resolution ; j++)
            {
                //Find the theoretical line of the rate of WIMP
                if (recoilX[kkk][j]>0)
                {
                    sig_E = a0_dE + a1_dE*sqrt(T_QF[kkk][j]);
                    
                    if(j==0) { dEx = T_QF[kkk][0]; } else { dEx = T_QF[kkk][j] - pEx; }
                    
                    if( (recoil[kkk][j]>0)&&(T_QF[kkk][j]>=0.7e-3) )
                    {
                        Factor1[kkk][jjj] = Factor1[kkk][jjj] + recoil[kkk][j]*exp(-pow((T_QF[kkk][j]-T_QF[kkk][jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                    }
                    pEx = T_QF[kkk][j];
                }
            }
        }
        //Find out the Full energy spectrum with the resolution of the detector based on the detector
        for(int jjj=0 ; jjj<Neutrino_Data_element ; jjj++)
        {
            for(int j=0 ; j<dm_spec_resolution ; j++)
            {
                //Find the theoretical line of the rate of WIMP
                if (recoilX[kkk][j]>0)
                {
                    sig_E = a0_dE + a1_dE*sqrt(T_QF[kkk][j]);
                    
                    if(j==0) { dEx = T_QF[kkk][0]; } else { dEx = T_QF[kkk][j] - pEx; }
                    
                    if( (recoil[kkk][j]>0)&&(T_QF[kkk][j]>=0.7e-3) )
                    {
                        Factor[kkk][jjj] = Factor[kkk][jjj] + recoil[kkk][j]*exp(-pow((T_QF[kkk][j]-RE_Nu_DATA[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                    }
                    pEx = T_QF[kkk][j];
                }
            }
            //cout << "Factor[3][jjj]: " << Factor[3][jjj] << endl;
        }
        
    }
    //Find out the scaling factor
    //===============================================================
    double Initial_Factor[Number_mx_Candidate];
    double Ratio[Number_mx_Candidate];
    
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        Initial_Factor[kkk] = Factor[kkk][0]/(RE_Nu_Rate[0]);
    }
        
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        for(int iii=0 ; iii<Neutrino_Data_element ; iii++)
        {
                double Ratio1 = (Factor[kkk][iii]/(RE_Nu_Rate[iii]));
                if(iii>0 && abs(Ratio1)>Initial_Factor[kkk] && abs(1/Ratio1)!=0)
                {
                    Initial_Factor[kkk] = Ratio1;
                }
        }
    }
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        //cout << "1e-40*(1/Initial_Factor[kkk]): " << 1e-40*(1/Initial_Factor[kkk]) << endl;
        Limitation[kkk] = 1e-40*(1/Initial_Factor[kkk]);
    }
        //===============================================================
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        //cout << "1/Ratio[kkk]: " << 1/Ratio[kkk] << endl;
        for(int jjj=0 ; jjj<1000 ; jjj++)
        {
            //cout << "==============================================" << endl;
            //cout << "T_QF[jjj][i] : " << T_QF[kkk][jjj]  << endl;
            //cout << "Factor1[kkk][jjj]: " << Factor1[kkk][jjj] << endl;
            Factor2[kkk][jjj] = Factor1[kkk][jjj]*(1/Initial_Factor[kkk]);
            //cout << "==============================================" << endl;
        }
    }

for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
{
    cout << "WIMP_mx: " << WIMP_mx[kkk] << endl;
    cout << "Limitation[kkk]: " << Limitation[kkk] << endl;
}

//===============================================================

char fout_name[100];
sprintf(fout_name,"Neutrino_floor.root");
TFile *fout=new TFile(fout_name,"recreate");

TGraph *Set_Limits_plot = new TGraph(45,WIMP_mx,Limitation);
Set_Limits_plot->SetLineColor(2);
Set_Limits_plot->SetMarkerStyle(21);
Set_Limits_plot->SetMarkerSize(1.2);
Set_Limits_plot->SetTitle("Set_Limits_plot");
Set_Limits_plot->GetXaxis()->SetTitle("WIMP mass[GeV]");
Set_Limits_plot->GetYaxis()->SetTitle("");
Set_Limits_plot->GetXaxis()->SetRangeUser(0,1e4);
Set_Limits_plot->GetYaxis()->SetRangeUser(1e-49,1e-42);

TLegend *leg= new TLegend(0.5,0.5,0.9,0.9);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.06);
leg->SetBorderSize(0);
leg->SetTextFont(22);
leg->AddEntry(Set_Limits_plot,"Neutrino Floor","lp");


Set_Limits_plot->Draw("apl");
leg->Draw();

Set_Limits_plot->Write();
c1->SetLogy();
c1->SetLogx();
c1->Print("Neutrino_floor.png");
//68% (1sigma) => 90% (2sigma)
}
