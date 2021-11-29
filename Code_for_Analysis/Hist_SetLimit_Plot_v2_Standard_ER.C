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

#include "VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"
#include "dsigma_dT2.h"
#include "cpkkd_chi_N_v3.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"
#include "VrV_le_Xenon1T.h"

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


void Hist_SetLimit_Plot_v2_Standard_ER(int Index, int File)
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
const int reso_T = 2510;
double T[Number_mx_Candidate][reso_T];
double T_QF[Number_mx_Candidate][reso_T];
double recoil[Number_mx_Candidate][reso_T];
double v;
double WIMP_mx[Number_mx_Candidate]={11,13,15,17,19};
double WIMP_max_T[Number_mx_Candidate];
    

    
    vector<double> Cross_section;
    vector<double> WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07};

    //=========
    // Data_implemented_no_rebin
    /*
     const int Data_element = 257;
     double RE_DATA[Data_element]; double RE_Rate[Data_element]; double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
     double RE_DATA_KL_extracted[Data_element]; double RE_Rate_KL_extracted[Data_element]; double RE_DATA_Err_KL_extracted[Data_element]; double RE_Rate_Err_KL_extracted[Data_element];

    for(int jjj=0 ; jjj<257 ; jjj++)
    {
        RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][3]*1.64458/0.994;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2];
    }*/
    // Data_implemented_rebi n_2
    //============================================================================
    //For TEXONO
    /*
     const int Data_element = 257;
     double RE_DATA[Data_element]; double RE_Rate[Data_element]; double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
     double RE_DATA_KL_extracted[Data_element]; double RE_Rate_KL_extracted[Data_element]; double RE_DATA_Err_KL_extracted[Data_element]; double RE_Rate_Err_KL_extracted[Data_element];

    for(int jjj=0 ; jjj<257 ; jjj++)
    {
        if(jjj==0){
            RE_DATA[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0]-0.025;
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
     */
    //============================================================================
    //For CDEX
    
         const int Data_element = 4;
         double RE_DATA[Data_element]; double RE_Rate[Data_element]; double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
         double RE_DATA_KL_extracted[Data_element]; double RE_Rate_KL_extracted[Data_element]; double RE_DATA_Err_KL_extracted[Data_element]; double RE_Rate_Err_KL_extracted[Data_element];

        double Data_RE_CDEX[4]= {0.185,0.235,0.285,0.325};
        double Data_RATE_CDEX[4]= {8.2,7.6,7.4,7.9};
        double Data_RE_Err[4]= {0,0,0,0};
        double Data_RE_Rate_Err[4]= {3,3.5,1.5,2};
         
    for(int jjj=0 ; jjj<Data_element ; jjj++)
    {
            RE_DATA[jjj]= Data_RE_CDEX[jjj];
            RE_Rate[jjj]= Data_RATE_CDEX[jjj];
            RE_DATA_Err[jjj]= 0;
            RE_Rate_Err[jjj] = Data_RE_Rate_Err[jjj]*1.64458/0.994;
    }
       
    //============================================================================
    
    //For XENON1T
    /*
        const int Data_element = 8;
    double RE_DATA[Data_element]; double RE_Rate[Data_element]; double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
    double RE_DATA_KL_extracted[Data_element]; double RE_Rate_KL_extracted[Data_element]; double RE_DATA_Err_KL_extracted[Data_element]; double RE_Rate_Err_KL_extracted[Data_element];

        for(int jjj=0 ; jjj<Data_element ; jjj++)
        {
            RE_DATA[jjj]= VrV_le_Xenon1T[jjj][0];//keV
            RE_Rate[jjj]= VrV_le_Xenon1T[jjj][1]*(1e-3);//(count/kg*day*keV)
            RE_DATA_Err[jjj]= 0;
            RE_Rate_Err[jjj] = VrV_le_Xenon1T[jjj][2];
            cout << "RE_DATA[jjj: " << RE_DATA[jjj] << endl;
            cout << "RE_Rate[jjj]: " << RE_Rate[jjj] << endl;
        }
     */
    //===

    //============================================================================
    //For CDMSlite(PhysRevLett.116.071301)
    /*
         const int Data_element = 1;
         double RE_DATA[Data_element]; double RE_Rate[Data_element]; double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
         double RE_DATA_KL_extracted[Data_element]; double RE_Rate_KL_extracted[Data_element]; double RE_DATA_Err_KL_extracted[Data_element]; double RE_Rate_Err_KL_extracted[Data_element];
         
        for(int jjj=0 ; jjj<Data_element ; jjj++)
        {
            RE_DATA[jjj]= 0.17;
            RE_Rate[jjj]= 5.5;
                RE_DATA_Err[jjj]= 0;
            RE_Rate_Err[jjj] = (1.0)*1.64458/0.994;
        }
     */
    //============================================================================

//========For Data
    double qf0, qf1, qf2, qf3, qf4, alpha=1.0;
    qf0 = alpha*0.19816;
    qf1 = alpha*0.05052;
    qf2 = alpha*0.00378;
    qf3 = alpha*0.00192;
    qf4 = alpha*0.0016;

    
    double T_Used[reso_T];double Recoil_X_Used[reso_T];
    
    //for(int kkk=1; kkk<WIMP_mx_Array.size(); kkk++)
    for(int kkk=Index; kkk<Index+1; kkk++)
    {
        double *T_QF_Original_Bef=RecoilX_Event(0,0,WIMP_mx_Array[kkk],0,4,1,File);
        for(int i=0;i<reso_T;i++){T_Used[i]=0;T_Used[i]=T_QF_Original_Bef[i];}
        double *Recoil_X_Bef=RecoilX_Event(1,0,WIMP_mx_Array[kkk],0,4,1,File);
        for(int i=0;i<reso_T;i++){Recoil_X_Used[i]=0;Recoil_X_Used[i]=(Recoil_X_Bef[i]);}
        
        for(int i=0;i<reso_T;i++)
        {
            //cout << "i: " << i << endl;
            //cout << "T_Used[i]: " << T_Used[i] << endl;
            //cout << "Recoil_X_Used[i]: " << Recoil_X_Used[i] << endl;
        }
        
        double Factor[Data_element];
        double pEx = T_Used[0];
        double sig_E; double dEx;
        //double a0_dE = 18/1000.0;
        double a0_dE = 33.4992/1000.0;
        double a1_dE = 13.2145/1000.0;
        //double a_dE = 31.3/100.0;
        //double b_dE = 17.0/100.0;

        for(int jjj=0 ; jjj<Data_element ; jjj++)
        {
            Factor[jjj]=0;
            for(int j=0 ; j<dm_spec_resolution ; j++)
            {
                //Find the theoretical line of the rate of WIMP
                if (Recoil_X_Used[j]>0)
                {
                    sig_E = a0_dE + a1_dE*sqrt(T_Used[j]);
                    //sig_E =  ( b_dE + a_dE/(sqrt(T_Used[j])) )*0.01;
                    //cout << "sig_E: " << sig_E << endl;

                    if(j==0) { dEx = T_Used[0]; } else { dEx = T_Used[j] - pEx; }
                    
                    if( (Recoil_X_Used[j]>0)&&(T_Used[j]>=0.7e-3) )
                    {
                        if(dEx>0)Factor[jjj] = Factor[jjj] + Recoil_X_Used[j]*exp(-pow((T_Used[j]-RE_DATA[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                        /*
                        if( dEx < 0){
                            cout << "j: " << j << endl;
                            cout << "dEx: " << dEx << endl;}
                        if( Recoil_X_Used[j] < 0){
                            cout << "Recoil_X_Used[j]: " << Recoil_X_Used[j] << endl;}
                         */

                    }
                    pEx = T_Used[j];
                }
            }
            //cout << "Factor[jjj]: " << Factor[jjj] << endl;
            //if(kkk==3)cout << "Factor[3][jjj]: " << Factor[3][jjj] << endl;
        }
        /*
        for(int i=0;i<Data_element;i++)
        {
            cout << "i: " << i << endl;
            cout << "RE_DATA[i]: " << RE_DATA[i] << endl;
            cout << "Factor[i]: " << Factor[i] << endl;
        }
         */
        //Find out the scaling factor
        
        //===============================================================
        double Initial_Factor= Factor[0]/(RE_Rate[0]+RE_Rate_Err[0]);
        
        cout << "=======================================" << endl;
            for(int iii=0 ; iii<Data_element ; iii++)
            {
                //cout << "RE_DATA[iii]: " << RE_DATA[iii]  << endl;
                //cout << "Factor[iii]: " << Factor[iii] << endl;
                //cout << "RE_Rate[iii]: " << RE_Rate[iii] << endl;
            }
        cout << "=======================================" << endl;

        int Bin_Possion_Bin=0;
        for(int iii=0 ; iii<Data_element ; iii++)
        {
                double Ratio1 = (Factor[iii]/(RE_Rate[iii]+RE_Rate_Err[iii]));
                if(iii>0 && abs(Ratio1)>Initial_Factor && abs(1/Ratio1)!=0)
                {
                    Bin_Possion_Bin = iii;
                    Initial_Factor = Ratio1;
                }
        }
        cout << "=======================================" << endl;
        cout << "Bin_Possion_Bin: " << Bin_Possion_Bin << endl;
        cout << "1./Initial_Factor: " << 1./Initial_Factor << endl;
        //===============================================================
        for(int jjj=0 ; jjj<Data_element ; jjj++)
        {
            Factor[jjj] = Factor[jjj]*(1./Initial_Factor);
            cout << "RE_DATA[jjj]: " << RE_DATA[jjj] << endl;
            cout << "Factor[jjj]: " << Factor[jjj] << endl;
         }
        double Scaling = sqrt((1./Initial_Factor));
        cout << "CS_Try(1*Scaling,0.5): " << CS_Try(1*Scaling,WIMP_mx_Array[kkk]) << endl;
        cout << "DS_Try(1*Scaling,0.5): " << DS_Try(1e-9*Scaling,WIMP_mx_Array[kkk]) << endl;

        //Cross_section.push_back(CS_Try(1*Scaling,WIMP_mx_Array[kkk]));
            //===============================================================
}
/*
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
leg->AddEntry(gcpkkdX1,Form("m_{#chi}=2.34GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[0],SigmaSI_Scale[0]),"lp");
//leg->AddEntry(gcpkkdX2,Form("m_{#chi}=1GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[1],SigmaSI_Scale[1]),"lp");
//leg->AddEntry(gcpkkdX3,Form("m_{#chi}=2GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[2],SigmaSI_Scale[2]),"lp");
leg->AddEntry(gcpkkdX4,Form("m_{#chi}=3GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[3],SigmaSI_Scale[3]),"lp");
leg->AddEntry(gcpkkdX5,Form("m_{#chi}=4GeV, #sigma^{SI}_{#chi N}=%f #times 10^{-%i}",SigmaSI_Scale_Number[4],SigmaSI_Scale[4]),"lp");

cdexdata->Draw("aplE");
CDEX_DATA->Draw("plsame");
gcpkkdX1->Draw("pl3same");
//gcpkkdX2->Draw("pl3same");
//gcpkkdX3->Draw("pl3same");
gcpkkdX4->Draw("pl3same");
gcpkkdX5->Draw("pl3same");
line1->Draw("lsame");
line2->Draw("lsame");
leg->Draw();
 
  
    //
c1->Print("TEXONO_CDEX.png");
 */
//68% (1sigma) => 90% (2sigma)
}
