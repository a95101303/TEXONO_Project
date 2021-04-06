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
#include "/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/Hist_SetLimit_Plot_v2_Extract_Peak_Debug_1.h"
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


void Hist_SetLimit_Plot_v2_Standard()
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
const int Number_mx_Candidate = 5;
double A = AGe;
const int reso_T = 1000;
double T[Number_mx_Candidate][reso_T];
double T_QF[Number_mx_Candidate][reso_T];
double recoil[Number_mx_Candidate][reso_T];
double v;
double WIMP_mx[Number_mx_Candidate]={0.1785,1,2,3,4};
double WIMP_max_T[Number_mx_Candidate];
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
const int Neutrino_Data_element = 2000;
double RE_DATA[Data_element]; double RE_Rate[Data_element]; double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
double RE_DATA_KL_extracted[Data_element]; double RE_Rate_KL_extracted[Data_element]; double RE_DATA_Err_KL_extracted[Data_element]; double RE_Rate_Err_KL_extracted[Data_element];
double qf0, qf1, qf2, qf3, qf4, alpha=1.0;
qf0 = alpha*0.19816;
qf1 = alpha*0.05052;
qf2 = alpha*0.00378;
qf3 = alpha*0.00192;
qf4 = alpha*0.0016;

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
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994;}
}


    double Data_RE_CDEX[4]= {0.185,0.235,0.285,0.325};
    double Data_RATE_CDEX[4]= {8.2,7.6,7.4,7.9};
    double Data_RE_Err[4]= {0,0,0,0};
    double Data_RE_Rate_Err[4]= {3,3.5,1.5,2};
/*
for(int jjj=0 ; jjj<4 ; jjj++)
{
 if(jjj<4){
        RE_DATA[jjj]= Data_RE_CDEX[jjj];
        RE_Rate[jjj]= Data_RATE_CDEX[jjj];
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = Data_RE_Rate_Err[jjj]*1.64458/0.994;}
    else{
        RE_DATA[jjj]= 0;
        RE_Rate[jjj]= 0;
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = 0;}
}*/

    
    
    
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
            recoil[jjj][i] = rate_scale_QF(T[jjj][i])*cpkkd_chi_N(WIMP_mx[jjj], T[jjj][i], 243.8125, 152.5, A);
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

    double Factor[Number_mx_Candidate][257];
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
                    
                    if( (recoilX[kkk][j]>0)&&(T_QF[kkk][j]>=0.7e-3) )
                    {
                        Factor1[kkk][jjj] = Factor1[kkk][jjj] + recoilX[kkk][j]*exp(-pow((T_QF[kkk][j]-T_QF[kkk][jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                    }
                    pEx = T_QF[kkk][j];
                }
            }
        }
        //Find out the Full energy spectrum with the resolution of the detector based on the detector
        for(int jjj=0 ; jjj<Data_element ; jjj++)
        {
            for(int j=0 ; j<dm_spec_resolution ; j++)
            {
                //Find the theoretical line of the rate of WIMP
                if (recoilX[kkk][j]>0)
                {
                    sig_E = a0_dE + a1_dE*sqrt(T_QF[kkk][j]);
                    
                    if(j==0) { dEx = T_QF[kkk][0]; } else { dEx = T_QF[kkk][j] - pEx; }
                    
                    if( (recoilX[kkk][j]>0)&&(T_QF[kkk][j]>=0.7e-3) )
                    {
                        Factor[kkk][jjj] = Factor[kkk][jjj] + recoilX[kkk][j]*exp(-pow((T_QF[kkk][j]-RE_DATA[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
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
        Initial_Factor[kkk] = Factor[kkk][0]/(RE_Rate[0]+RE_Rate_Err[0]);
    }
        
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        for(int iii=0 ; iii<Data_element ; iii++)
        {
                double Ratio1 = (Factor[kkk][iii]/(RE_Rate[iii]+RE_Rate_Err[iii]));
                if(iii>0 && abs(Ratio1)>Initial_Factor[kkk] && abs(1/Ratio1)!=0)
                {
                    Initial_Factor[kkk] = Ratio1;
                }
        }
    }
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        cout << "1/Ratio[kkk]: " << 1/Initial_Factor[kkk] << endl;

    }
        //===============================================================
    for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
    {
        //cout << "1/Ratio[kkk]: " << 1/Ratio[kkk] << endl;
        for(int jjj=0 ; jjj<1000 ; jjj++)
        {
            //Factor2[kkk][jjj] = Factor1[kkk][jjj]*(1/Initial_Factor[kkk]);
            Factor2[kkk][jjj] = Factor1[kkk][jjj];
            if(kkk==3)
            {
                //cout << "T[jjj][i]: " << T[kkk][jjj] << endl;
                //cout << "T_QF[kkk][jjj]: " << T_QF[kkk][jjj] << endl;
                //cout << "Factor2[kkk][jjj]: " << Factor2[kkk][jjj] << endl;
                //cout << "Factor1[kkk][jjj]: " << Factor1[kkk][jjj] << endl;
            }
        }
    }
    
        //===============================================================
TGraphErrors *cdexdata = new TGraphErrors(257,RE_DATA,RE_Rate,RE_DATA_Err,RE_Rate_Err);
cdexdata->SetName("cdexdata");
cdexdata->SetLineColor(2);
cdexdata->GetXaxis()->SetRangeUser(0,0.6);
cdexdata->GetYaxis()->SetRangeUser(0,200);
cdexdata->SetTitle("The rate of WIMP");
cdexdata->GetXaxis()->SetTitle("Recoil Energy[keV]");
cdexdata->GetYaxis()->SetTitle("Rate[Count/day*kg*keV]");
cdexdata->SetMarkerStyle(21);
    
TGraphErrors *CDEX_DATA = new TGraphErrors(4,Data_RE_CDEX,Data_RATE_CDEX,Data_RE_Err,Data_RE_Rate_Err);
CDEX_DATA->SetName("CDEX_DATA");
CDEX_DATA->SetLineColor(3);
CDEX_DATA->SetMarkerStyle(21);

TGraph *gcpkkdX1 = new TGraph(reso_T,T_QF[0],Factor2[0]);
gcpkkdX1->SetName("gcpkkdX1");
gcpkkdX1->SetLineColor(4);
TGraph *gcpkkdX2 = new TGraph(reso_T,T_QF[1],Factor2[1]);
gcpkkdX2->SetName("gcpkkdX2");
gcpkkdX2->SetLineColor(5);
TGraph *gcpkkdX3 = new TGraph(reso_T,T_QF[2],Factor2[2]);
gcpkkdX3->SetName("gcpkkdX3");
gcpkkdX3->SetLineColor(6);
TGraph *gcpkkdX4 = new TGraph(reso_T,T_QF[3],Factor2[3]);
gcpkkdX4->SetName("gcpkkdX4");
gcpkkdX4->SetLineColor(7);
TGraph *gcpkkdX5 = new TGraph(reso_T,T_QF[4],Factor2[4]);
gcpkkdX5->SetName("gcpkkdX5");
gcpkkdX5->SetLineColor(8);

    
TLine *line1 = new TLine(0,Data_RATE_CDEX[0]+Data_RE_Rate_Err[0],Data_RE_CDEX[0],Data_RATE_CDEX[0]+Data_RE_Rate_Err[0]);
TLine *line2 = new TLine(0,RE_Rate[0]+RE_Rate_Err[0],RE_DATA[0],RE_Rate[0]+RE_Rate_Err[0]);
    
TLegend *leg= new TLegend(0.5,0.5,0.9,0.9);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.03);
leg->SetBorderSize(0);
leg->SetTextFont(22);
    
double SigmaSI[Number_mx_Candidate];
double SigmaSI_Scale_Number[Number_mx_Candidate];
int    SigmaSI_Scale[Number_mx_Candidate];

for(int kkk=0 ; kkk<Number_mx_Candidate ; kkk++)
{
    SigmaSI[kkk] = (1/Initial_Factor[kkk])*(1e-40);
    cout << "SigmaSI[kkk]: " << SigmaSI[kkk] << endl;
}

for(int jjj=0 ; jjj<5 ; jjj++)
{
    for(int kkk=0 ; kkk<45 ; kkk++)
    {
        if(log10(SigmaSI[jjj]*pow(10,kkk))>0 and log10(SigmaSI[jjj]*pow(10,kkk))<1)
        {
            SigmaSI_Scale[jjj] = kkk;
            SigmaSI_Scale_Number[jjj] = SigmaSI[jjj]*pow(10,kkk);
            cout << "SigmaSI_Scale[jjj]: " << SigmaSI_Scale[jjj] << "SigmaSI_Scale_Number[jjj]: " << SigmaSI_Scale_Number[jjj] << endl;
        }
    }
}
cout << "SigmaSI_Scale_Number[0]: " << SigmaSI_Scale_Number[0] << endl;
cout << "SigmaSI_Scale[0]: " << SigmaSI_Scale[0] << endl;

cout << "(1/Initial_Factor[0])*1e-40): " << (1/Initial_Factor[0])*1e-40 << endl;
cout << "First_point: " << RE_Rate[0]+RE_Rate_Err[0] << endl;
leg->AddEntry(cdexdata,"TEXONO","lp");
leg->AddEntry(CDEX_DATA,"CDEX","lp");
//leg->AddEntry(gcpkkdX1,Form("m_{#chi}=0.1785GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[0],SigmaSI_Scale[0]),"lp");
//leg->AddEntry(gcpkkdX2,Form("m_{#chi}=1GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[1],SigmaSI_Scale[1]),"lp");
//leg->AddEntry(gcpkkdX3,Form("m_{#chi}=2GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[2],SigmaSI_Scale[2]),"lp");
leg->AddEntry(gcpkkdX4,Form("m_{#chi}=3GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[3],SigmaSI_Scale[3]),"lp");
//leg->AddEntry(gcpkkdX5,Form("m_{#chi}=4GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[4],SigmaSI_Scale[4]),"lp");

cdexdata->Draw("aplE");
CDEX_DATA->Draw("plsame");
//gcpkkdX1->Draw("pl3same");
//gcpkkdX2->Draw("pl3same");
//gcpkkdX3->Draw("pl3same");
gcpkkdX4->Draw("pl3same");
//gcpkkdX5->Draw("pl3same");
line1->Draw("lsame");
line2->Draw("lsame");
leg->Draw();
 
  
    //
c1->Print("TEXONO_CDEX.png");
//68% (1sigma) => 90% (2sigma)
}
