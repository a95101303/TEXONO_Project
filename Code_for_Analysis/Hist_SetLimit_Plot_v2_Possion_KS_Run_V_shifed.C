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

void Hist_SetLimit_Plot_v2_Possion_KS_Run_V_shifed()
{
    const double Total_Bin_Number = 2000;
    vector<string> Mass_Point       ={"0P1","0P09","0P08"};
    vector<double> Mass_Point_Vector={0.1,0.09,0.08};

    double TEXONO_Threshold_E = 0.2;//keV

    for(int KKK=1; KKK<2; KKK++)
    {
        double Min_V = Velocity_DM(Mass_Point_Vector[KKK],TEXONO_Threshold_E);
        int   N_FILE = 0;
        vector<double>Cross_Section;
        vector<double>N_aft_divide_N_Bef_element;
        vector<double>Mean_Energy;
        for(int FILE=1; FILE<28; FILE++)
        {
            //Include the FILE and the HIST
            TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Flux_CAT/%sGeV/Recoil_Spectrum/MD_%i.root",Mass_Point[KKK].c_str(),FILE));
            TFile *FILE_Index = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Flux_CAT/%sGeV/%i.root",Mass_Point[KKK].c_str(),FILE));

            if(ROOT_FILE!=NULL)
            {
                TH1F *Flux_HIST_Random;TH1F *Flux_HIST_Aft_Collision_EARTH;
                ES_HIST_Random             =(TH1F*)ROOT_FILE->Get("ER_Spectrum_Bef");//Energy Spectrum(ES)
                ES_Aft_Collision_EARTH     =(TH1F*)ROOT_FILE->Get("ER_Spectrum_Aft");//
                
                TTree *T1_TREE = (TTree*)FILE_Index->Get("t1");
                Double_t mx,sigma_si;
                T1_TREE->SetBranchAddress("mx",&mx);T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
                T1_TREE->GetEntry(0);
                
                //
                double 
                for(int Bin=0; Bin<Total_Bin_Number; Bin++)
                {
                    double Recoil_Energy = 0.16+0.005*(double)Bin;//keV
                    double ES_Value_Bef  = ES_HIST_Random        ->Eval(Recoil_Energy);
                    double ES_Value_Aft  = ES_Aft_Collision_EARTH->Eval(Recoil_Energy);
                }
                cout << "sigma_si: " << sigma_si << endl;
                //cout << "Event_Aft/Event_Bef: " << Event_Aft/Event_Bef << endl;
                cout << "Mean_Value_Aft: " << Energy_DM(Mass_Point_Vector[KKK],Mean_Value_Aft_V*1e3/3e8) << endl;
                Cross_Section.push_back(sigma_si);
                N_aft_divide_N_Bef_element.push_back(Event_Aft/Event_Bef);
                Mean_Energy.push_back(1e3*Energy_DM(Mass_Point_Vector[KKK],Mean_Value_Aft_V*1e3/3e8));
                N_FILE = N_FILE + 1;
            }
        }
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.04,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);
        
        TGraph *TGraph_Mean_Energy = new TGraph(N_FILE, &Cross_Section[0], &Mean_Energy[0]);
        TGraph_Mean_Energy->SetMarkerStyle(20);
        TGraph_Mean_Energy->SetMarkerColor(2);
        TGraph_Mean_Energy->SetMarkerColor(2);
        TGraph_Mean_Energy->GetXaxis()->SetRangeUser(1e-40,1e-26);
        TGraph_Mean_Energy->GetYaxis()->SetRangeUser(0,1e3);
        TGraph_Mean_Energy->GetYaxis()->SetTitle("<T>(eV)");
        TGraph_Mean_Energy->GetXaxis()->SetTitle("#sigma_{SI}(cm^2)");
        TGraph_Mean_Energy->Draw("apl");

        c1->SetLogx();
        //c1->SetLogy();

    }
}

/*
//Run the program for the individual index and the simulated number of events
void Hist_SetLimit_Plot_v2_Possion_KS_Run_V_shifed()
{
    const double Total_Bin_Number = 2000;
    vector<string> Mass_Point       ={"0P1","0P09","0P08"};
    vector<double> Mass_Point_Vector={0.1,0.09,0.08};

    double TEXONO_Threshold_E = 0.2;//keV

    for(int KKK=1; KKK<2; KKK++)
    {
        double Min_V = Velocity_DM(Mass_Point_Vector[KKK],TEXONO_Threshold_E);
        int   N_FILE = 0;
        vector<double>Cross_Section;
        vector<double>N_aft_divide_N_Bef_element;
        vector<double>Mean_Energy;
        for(int FILE=1; FILE<28; FILE++)
        {
            //Include the FILE and the HIST
            TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Flux_CAT/%sGeV/%i.root",Mass_Point[KKK].c_str(),FILE));
            if(ROOT_FILE!=NULL)
            {
                TH1F *Flux_HIST_Random;TH1F *Flux_HIST_Aft_Collision_EARTH;
                Flux_HIST_Random             =(TH1F*)ROOT_FILE->Get("Flux_HIST_Random");//Flux_HIST_Random(FHR)
                Flux_HIST_Aft_Collision_EARTH=(TH1F*)ROOT_FILE->Get("Flux_HIST_Aft_Collision_EARTH");//Flux_HIST_Aft_Collision_EARTH(FHACE)
                
                TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
                Double_t mx,sigma_si;
                T1_TREE->SetBranchAddress("mx",&mx);T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
                T1_TREE->GetEntry(0);

                double Mean_Value_Aft_V   = Flux_HIST_Aft_Collision_EARTH->GetMean();//(km/s)
                cout << "Mean_Value_Aft_V: " << Mean_Value_Aft_V << endl;
                //
                double TEXONO_Threshold_V = Velocity_DM(Mass_Point_Vector[KKK],TEXONO_Threshold_E);

                //Int_t Maximum_Bin_FHR   = Flux_HIST_Random              ->FindLastBinAbove();
                //Int_t Maximum_Bin_FHACE = Flux_HIST_Aft_Collision_EARTH->FindLastBinAbove();
                double Event_Bef=0;double Event_Aft=0;
                for(int Bin=0; Bin<Total_Bin_Number; Bin++)
                {
                    double Velocity_Bef = Flux_HIST_Random             ->GetXaxis()->GetBinCenter(Bin);//(km/s)
                    double Number_Bef   = Flux_HIST_Random             ->GetBinContent(Bin);//(km/s)
                    double Velocity_Aft = Flux_HIST_Aft_Collision_EARTH->GetXaxis()->GetBinCenter(Bin);//(km/s)
                    double Number_Aft   = Flux_HIST_Aft_Collision_EARTH->GetBinContent(Bin);//(km/s)

                    if(Velocity_Bef>TEXONO_Threshold_V)Event_Bef = Event_Bef + Number_Bef;
                    if(Velocity_Aft>TEXONO_Threshold_V)Event_Aft = Event_Aft + Number_Aft;
                }
                cout << "sigma_si: " << sigma_si << endl;
                //cout << "Event_Aft/Event_Bef: " << Event_Aft/Event_Bef << endl;
                cout << "Mean_Value_Aft: " << Energy_DM(Mass_Point_Vector[KKK],Mean_Value_Aft_V*1e3/3e8) << endl;
                Cross_Section.push_back(sigma_si);
                N_aft_divide_N_Bef_element.push_back(Event_Aft/Event_Bef);
                Mean_Energy.push_back(1e3*Energy_DM(Mass_Point_Vector[KKK],Mean_Value_Aft_V*1e3/3e8));
                N_FILE = N_FILE + 1;
            }
        }
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.04,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);
        
        TGraph *TGraph_Mean_Energy = new TGraph(N_FILE, &Cross_Section[0], &Mean_Energy[0]);
        TGraph_Mean_Energy->SetMarkerStyle(20);
        TGraph_Mean_Energy->SetMarkerColor(2);
        TGraph_Mean_Energy->SetMarkerColor(2);
        TGraph_Mean_Energy->GetXaxis()->SetRangeUser(1e-40,1e-26);
        TGraph_Mean_Energy->GetYaxis()->SetRangeUser(0,1e3);
        TGraph_Mean_Energy->GetYaxis()->SetTitle("<T>(eV)");
        TGraph_Mean_Energy->GetXaxis()->SetTitle("#sigma_{SI}(cm^2)");
        TGraph_Mean_Energy->Draw("apl");

        c1->SetLogx();
        //c1->SetLogy();

    }
}
*/
/*
TGraph *N_aft_divide_N_Bef = new TGraph(N_FILE, &Cross_Section[0], &N_aft_divide_N_Bef_element[0]);
N_aft_divide_N_Bef->SetMarkerStyle(20);
N_aft_divide_N_Bef->SetMarkerColor(2);
N_aft_divide_N_Bef->SetMarkerColor(2);
N_aft_divide_N_Bef->GetXaxis()->SetRangeUser(1e-40,1e-26);
N_aft_divide_N_Bef->GetYaxis()->SetRangeUser(0,1);
N_aft_divide_N_Bef->GetYaxis()->SetTitle("Earth/Vacuum");
N_aft_divide_N_Bef->GetXaxis()->SetTitle("#sigma_{SI}(cm^2)");
N_aft_divide_N_Bef->Draw("apl");
 */

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
