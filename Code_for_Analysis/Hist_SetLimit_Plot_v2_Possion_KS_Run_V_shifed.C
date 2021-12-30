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
#include "velocity_distribution_2000_Ave.h"
#include "dsigma_dT2.h"

const int          Simulated_Event_Number = 10;
const double       Max_V                  = 779.;//(km/s)
//Constant
const double       NaI_Density            = 3.67;//3.67(g/cm^3)
const double       NaI_Atomic_Mass        = 22.98*0.5+126*0.5;//
const double       Pb_Density             = 11.29;//3.67(g/cm^3)
const double       Pb_Atomic_Mass         = 207.2;//
const double       Fixed_Length           = 20.;//cm

const double Density_Array[3]           ={2.8 ,Pb_Density    ,NaI_Density};
const double Atomic_Mass_Array[3]={Weighted_Atomic_Number_Cement,Pb_Atomic_Mass,NaI_Atomic_Mass};
const double Number_Density_Array[3]    ={
                                          Density_Array[0]*1./(unified_atomic_mass_g*(Atomic_Mass_Array[0])) ,
                                          Density_Array[1]*1./(unified_atomic_mass_g*(Atomic_Mass_Array[1])) ,
                                          Density_Array[1]*1./(unified_atomic_mass_g*(Atomic_Mass_Array[2])) ,
                                        };

const double Length_Array[3]     ={30*1e2,25,15};//30 meter-water-equivalent, 25cm lead, 15cm NaI

double Mean_free_Path_check(double WIMP_Mass, double Density, double Atomic_Mass, double Sigma_SI)
{
    double MFP_Calculated = 1./((Density)/(unified_atomic_mass_g*(Atomic_Mass))*total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));
    double Collision_Time   = (Density*(MFP_Calculated))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));
    double Collision_Time_2 = (Density*(15.))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));
    /*
    cout << "Atomic_Mass: " << Atomic_Mass << endl;
    cout << "Once: " << (Density*(3))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,1.5e-27,WIMP_Mass,Atomic_Mass)) << endl;;
    cout << "MFP_Calculated: "    << MFP_Calculated << endl;
    cout << "Collision_Time: "    << Collision_Time << endl;
    cout << "Collision_Time_2: "  << Collision_Time_2 << endl;
     */
    return MFP_Calculated;//cm
}

double *Run_Program(double WIMP_Mass, double Density, double Atomic_Mass, double Sigma_SI)//Density, Atomic mass
{
    static double Array[3];
    double Initial_V                      = Max_V;
    double Last_V                         = Max_V;
    double Energy_Loss_Percentage_total   = 0;
    double Collision_Time                 = (Density*(Fixed_Length))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));
    cout << "Collision_Time: " << Collision_Time << endl;
    double Collision_Time_Int             = 0;
    double Energy_Loss_Average            = 0;
    //cout << "Collision_Time: " << Collision_Time << endl;
    while((int)Collision_Time>0 and Energy_DM(WIMP_Mass,Last_V*1e3/3e8)>0.01)
    {
        double *V_aft      = Velocity_Aft_collision_Bent(1,WIMP_Mass,Sigma_SI,Initial_V,Atomic_Mass);
        Initial_V          = V_aft[0];
        double Energy_Loss = Energy_DM(WIMP_Mass,Last_V*1e3/3e8)-Energy_DM(WIMP_Mass,Initial_V*1e3/3e8);//keV
        double Energy_Loss_Percentage = (Energy_Loss)/Energy_DM(WIMP_Mass,Last_V*1e3/3e8);
        //cout << "Energy_Loss: " << Energy_Loss << endl;
        //double Energy_Loss_Percentage = (Energy_Loss)/(Energy_DM(WIMP_Mass,Max_V*1e3/3e8)-0.01);
        //cout << "Energy_Loss_Percentage: " << Energy_Loss_Percentage << endl;
        Energy_Loss_Average          = Energy_Loss_Average + Energy_Loss;
        Energy_Loss_Percentage_total = Energy_Loss_Percentage_total + Energy_Loss_Percentage;
        Last_V             = Initial_V;
        Collision_Time     = Collision_Time - 1;//
        
        Collision_Time_Int = Collision_Time_Int + 1;
    }
    double Energy_Diff = Energy_DM(WIMP_Mass,Max_V*1e3/3e8)-Energy_DM(WIMP_Mass,Last_V*1e3/3e8);
    double Expected_Loss = Energy_DM(WIMP_Mass,Max_V*1e3/3e8)-0.01;
    cout << "Energy_Diff/Fixed_Length: " << Energy_Diff/Fixed_Length << endl;//(dE/dX)
    cout << "(Energy_Diff/Fixed_Length)/(Expected_Loss): " << (Energy_Diff/Fixed_Length)/(Expected_Loss) << endl;//(dE/dX)/(Real_Loss)
    cout << "Collision_Time_Int/Fixed_Length: " << Collision_Time_Int/Fixed_Length << endl;//Collision Time/length
    cout << "Collision_Time_Int: " << Collision_Time_Int << endl;
    cout << "Energy_Diff/Collision_Time_Int: " << Energy_Diff/Collision_Time_Int << endl;//Energy_Loss/Collision Time

    //cout << "Collision_Time_Int: " << Collision_Time_Int  << endl;
    Energy_Loss_Percentage_total = Energy_Loss_Percentage_total/(Collision_Time_Int);
    Energy_Loss_Average          = Energy_Loss_Average/(Collision_Time_Int);
    cout << "Energy_Loss_Average :" << Energy_Loss_Average  << endl;
    //Energy_Loss_Average          = Energy_Loss_Average/(Collision_Time_Int*Energy_Diff);
    Array[0]=Initial_V;Array[1]=Energy_Loss_Percentage_total;Array[2]=Energy_Loss_Average;
    return Array;
}


//Run the program for the individual index and the simulated number of events
void Hist_SetLimit_Plot_v2_Possion_KS_Run_V_shifed()
{
    vector<string> Mass_Point={"0P1","0P09",};

    for(int KKK=0; KKK<Mass_Point.size(); KKK++)
    {
        for(int FILE=0; FILE<27; FILE++)
        {
        TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Flux_CAT/%sGeV/%i.root",Mass_Point[Mass_INT].c_str(),FILE));
        
        TH1F *Flux_HIST_Random;TH1F *Flux_HIST_Aft_Collision_EARTH;
        Flux_HIST_Random=(TH1F*)ROOT_FILE->Get("Flux_HIST_Random");
        Flux_HIST_Aft_Collision_EARTH=(TH1F*)ROOT_FILE->Get("Flux_HIST_Aft_Collision_EARTH");
        }
    }
}
/*
TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.04,"XY");
gStyle->SetTitleFont(62,"XY");
gStyle->SetLegendFont(62);

TGraph * Energy_Aft_Air = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_Air_for_every_cross_section);
Energy_Aft_Air->SetMarkerStyle(20);
Energy_Aft_Air->SetMarkerColor(2);
Energy_Aft_Air->SetMarkerColor(2);
Energy_Aft_Air->GetXaxis()->SetRangeUser(1e-27,1e-26);
Energy_Aft_Air->GetYaxis()->SetRangeUser(0,Energy_DM(Mass_Array[Mass],Max_V*1e3/3e8)+0.01);
Energy_Aft_Air->GetYaxis()->SetTitle("Energy of DM(keV)");
Energy_Aft_Air->GetXaxis()->SetTitle("#sigma_{SI}(cm^2)");
Energy_Aft_Air->Draw("apl");

TGraph * Energy_Aft_NaI = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_NaI_for_every_cross_section);
Energy_Aft_NaI->SetMarkerStyle(20);
Energy_Aft_NaI->SetMarkerColor(3);
Energy_Aft_NaI->SetMarkerColor(3);
Energy_Aft_NaI->GetYaxis()->SetRangeUser(0,20);
Energy_Aft_NaI->Draw("plsame");

TGraph * Energy_Aft_No = new TGraph(Sigma_Array.size(), &Sigma_Array[0], No_ELoss_for_every_cross_section);
Energy_Aft_No->SetMarkerStyle(20);
Energy_Aft_No->SetMarkerColor(4);
Energy_Aft_No->SetMarkerColor(4);
Energy_Aft_No->GetYaxis()->SetRangeUser(0,20);
Energy_Aft_No->Draw("plsame");

TGraph * Energy_Aft_Water = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_Water_for_every_cross_section);
Energy_Aft_Water->SetMarkerStyle(20);
Energy_Aft_Water->SetMarkerColor(5);
Energy_Aft_Water->SetMarkerColor(5);
Energy_Aft_Water->GetYaxis()->SetRangeUser(0,20);
Energy_Aft_Water->Draw("plsame");

TGraph * Energy_Aft_Lead = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_Lead_for_every_cross_section);
Energy_Aft_Lead->SetMarkerStyle(20);
Energy_Aft_Lead->SetMarkerColor(6);
Energy_Aft_Lead->SetMarkerColor(6);
Energy_Aft_Lead->GetYaxis()->SetRangeUser(0,20);
Energy_Aft_Lead->Draw("plsame");

TLegend *leg = new TLegend(0.1,0.1,0.4,0.4);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.04);
leg->SetBorderSize(0);
leg->SetTextFont(22);
leg->AddEntry(Energy_Aft_No,"Max_E","lp");
leg->AddEntry(Energy_Aft_Air,"Air","lp");
leg->AddEntry(Energy_Aft_Water,"Air+30(M.W.E)","lp");
leg->AddEntry(Energy_Aft_Lead,"Air+30(M.W.E)+Lead","lp");
leg->AddEntry(Energy_Aft_NaI,"Air+30(M.W.E)+Lead+NaI","lp");

//c1->SetLogy();
c1->SetLogx();
leg->Draw();
c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/NaI_Proof/0P1GeV.pdf");
*/
