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
void Hist_SetLimit_Plot_v2_Possion_KS_Run_NaI_Lead()
{
    const double       Same_Sigma             = 1e-31;
    const double       Sigma_SI_NaI           = Same_Sigma;
    const double       Sigma_SI_Pb            = Same_Sigma;
    const double       Mx                     = 0.5;//GeV
    const double       Atomic_Mass            = 34.;//GeV
    const double       Velocity               = 700.;//GeV

    const double       Mx_A                   = 0.06;//GeV
    const double       Mx_B                   = 0.1;//GeV
    
    
    TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.04,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);
    
    const int Bin_Number=1000;
    double T_A_array[Bin_Number];double dsigma_dT_A_array[Bin_Number];
    double Max_A = max_recoil_A_keV(Mx_A,Velocity*1e3/3e8,Atomic_Mass);
    double A_Bin_size = Max_A/(double)Bin_Number;
    for(int KKK=0; KKK<Bin_Number; KKK++)
    {
        double T_A = (double)(KKK + 1)*A_Bin_size;
        double dsigma_dT_A = fdsigma_dT_keV(Mx_A,Same_Sigma,Velocity*1e3/3e8,Atomic_Mass,T_A);
        T_A_array[KKK] = T_A;dsigma_dT_A_array[KKK] = dsigma_dT_A;
        cout << "KKK: " << KKK << endl;
        cout << "T_A: " << T_A << endl;
        cout << "dsigma_dT_A: " << (dsigma_dT_A) << endl;
    }
    double Diff_T         = T_A_array[998]-T_A_array[1];
    double Diff_dsigma_dT = TMath::Log10(dsigma_dT_A_array[998]) - TMath::Log10(dsigma_dT_A_array[1]);
    double Ratio          = dsigma_dT_A_array[1]/dsigma_dT_A_array[998];
    cout << "Diff_dsigma_dT: " << Diff_dsigma_dT/Diff_T << endl;
    //cout << "Ratio/Diff_dsigma_dT: " <<  Ratio/Diff_T << endl;
    TGraph * TGraph_dsigma_dT_A = new TGraph(Bin_Number, T_A_array, dsigma_dT_A_array);
    TGraph_dsigma_dT_A->SetMarkerStyle(20);
    TGraph_dsigma_dT_A->SetLineColor(3);
    TGraph_dsigma_dT_A->SetMarkerColor(3);
    TGraph_dsigma_dT_A->GetXaxis()->SetRangeUser(0.,Max_A);
    TGraph_dsigma_dT_A->GetYaxis()->SetRangeUser(1e-28,1e-27);
    TGraph_dsigma_dT_A->GetYaxis()->SetTitle("T[keV]");
    TGraph_dsigma_dT_A->GetXaxis()->SetTitle("dsigma_dT");
    TGraph_dsigma_dT_A->Draw("al");

    /*
    double Max_B = max_recoil_A_keV(Mx_B,Velocity*1e3/3e8,Atomic_Mass);

    TF1 *f_1 = new TF1("f3","fdsigma_dT_keV([0],[1],[2],[3],x)",0,max_recoil_A_keV(Mx_A,Velocity*1e3/3e8,Atomic_Mass));
    f_1->SetParameter(0,Mx_A);f_1->SetParameter(1,Same_Sigma);f_1->SetParameter(2,(Velocity*1e3/3e8));f_1->SetParameter(3,Atomic_Mass);
    f_1->SetLineColor(2);
    
    f_1->GetXaxis()->SetRangeUser(0,max_recoil_A_keV(Mx_A,779*1e3/3e8,Atomic_Mass));
    f_1->GetYaxis()->SetRangeUser(0,1e-20);
    f_1->Draw();

    
    TF1 *f_2 = new TF1("f3","fdsigma_dT_keV([0],[1],[2],[3],x)",0,max_recoil_A_keV(0.1,Velocity*1e3/3e8,Atomic_Mass));
    f_2->SetParameter(0,Mx_B);f_2->SetParameter(1,Same_Sigma);f_2->SetParameter(2,(Velocity*1e3/3e8));f_2->SetParameter(3,Atomic_Mass);
    f_2->SetLineColor(3);
    f_2->Draw("same");
     */
    /*
    TGraph * Energy_Max = new TGraph(Mass_Array_1.size(), &Mass_Array_1[0], &Energy_Array[0]);
    Energy_Max->SetMarkerStyle(20);
    Energy_Max->SetMarkerColor(2);
    Energy_Max->SetMarkerColor(2);
    Energy_Max->GetXaxis()->SetRangeUser(0.05,0.5);
    Energy_Max->GetYaxis()->SetRangeUser(0,2);
    Energy_Max->GetYaxis()->SetTitle("Max Energy(keV)");
    Energy_Max->GetXaxis()->SetTitle("M_{#chi}(GeV)");
    Energy_Max->Draw("apl");
     */
    
    c1->SetLogy();
    //c1->SetLogx();
    c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/NaI_Proof/dsigma_dT_comparison.pdf");

    /*
    vector<double> Cross_Section_Array={8e-28,2e-27};
    //vector<double> Mass_Array={2,0.2};
    int NaI_Event=0;int Pb_Event=0;
    
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);

    TF1 *f3 = new TF1("f3","fdsigma_dT_keV([0],[1],[2],[3],x)",0,max_recoil_A_keV(0.1,779*1e3/3e8,Atomic_Mass));
    f3->SetParameter(0,Mx);f3->SetParameter(1,Same_Sigma);f3->SetParameter(2,(779*1e3/3e8));f3->SetParameter(3,Atomic_Mass);

    TH1F *Try = new TH1F("Try","Try",1000,0,max_recoil_A_keV(0.1,779*1e3/3e8,Atomic_Mass));
    double Energy_Loss_dEdX = MFP_from_DCS_Part(0,779,Same_Sigma,Mx,Atomic_Mass)*1e-3;
    vector<int>    Collision_Time;
    vector<double> Energy_Loss_averaged;
    vector<double> vector_Diff_between_dEdX_and_every;
    const int Event_Number_try=2;
    
    for(int MMM=1; MMM<2; MMM++)
    {
        int Collision_Time_Now = 50000;
        Collision_Time.push_back(Collision_Time_Now);
        
        double Averaged_Energy_Loss=0;
        double Diff_between_dEdX_and_every=0;
        for(int Event=0; Event<Event_Number_try; Event++)
        {
            double Energy_total=0;
            cout << "====================================" << endl;
            cout << "-----------------------------" << endl;
            int Count_Collision=0;
            for(int Collision_Time_Index=0; Collision_Time_Index<(int)Collision_Time_Now; Collision_Time_Index++)
            {
                double Random_Energy= f3->GetRandom();
                //cout << "Energy_Loss_dEdX: " << Energy_Loss_dEdX << endl;
                //cout << "Random_Energy: " << Random_Energy << endl;
                Energy_total = Energy_total + Random_Energy;
                //cout << "Energy_total: " << Energy_total << endl;
                Count_Collision = Count_Collision + 1;
                //cout << "Energy_total/Count_Collision: " << Energy_total/(double)Count_Collision << endl;
                //cout << "(Energy_total-Energy_Loss_dEdX)/Energy_Loss_dEdX: " << abs((Energy_total/(double)Count_Collision-Energy_Loss_dEdX)/Energy_Loss_dEdX) << endl;

            }
            cout << "-----------------------------" << endl;
            Energy_total = Energy_total/Collision_Time_Now;
            Diff_between_dEdX_and_every = Diff_between_dEdX_and_every + abs(Energy_total/Energy_Loss_dEdX);
            cout << "<Energy_total>: " << Energy_total << endl;
            cout << "Energy_Loss_dEdX: " << Energy_Loss_dEdX << endl;
            cout << "(Energy_total/Energy_Loss_dEdX): " << (Energy_total/Energy_Loss_dEdX) << endl;
            cout << "abs(Energy_total/Energy_Loss_dEdX): " << abs(Energy_total/Energy_Loss_dEdX) << endl;
            Averaged_Energy_Loss = Averaged_Energy_Loss + (Energy_total);
            cout << "====================================" << endl;
        }
        Averaged_Energy_Loss = Averaged_Energy_Loss/(double)Event_Number_try;
        cout << "Averaged_Energy_Loss: " << Averaged_Energy_Loss << endl;
        Energy_Loss_averaged.push_back(Averaged_Energy_Loss/Energy_Loss_dEdX);
        vector_Diff_between_dEdX_and_every.push_back((Diff_between_dEdX_and_every)/(double)Event_Number_try);
    }
    
    
    for(int KKK=0; KKK<25-1; KKK++)
    {
        cout << "Collision_Time: " << Collision_Time[KKK] << endl;
        cout << "Energy_Loss_averaged: " << Energy_Loss_averaged[KKK] << endl;
        cout << "vector_Diff_between_dEdX_and_every: " << vector_Diff_between_dEdX_and_every[KKK] << endl;
    }
     */
    
    /*
    for(int KKK=0; KKK<5000; KKK++)
    {
        double Random_Energy= f3->GetRandom();
        cout << "Random_Energy: " << Random_Energy << endl;
        Try->Fill(Random_Energy);
    }
     
    Try->SetMarkerStyle(20);
    Try->SetMarkerColor(2);
    Try->SetMarkerColor(2);
    Try->GetXaxis()->SetRangeUser(0,max_recoil_A_keV(0.1,779*1e3/3e8,34.));
    Try->GetYaxis()->SetRangeUser(0,200);
    Try->Draw("Hist");
     */

    
    /*
    for(int Mass_Index=0; Mass_Index<Mass_Array.size(); Mass_Index++)
{
    cout << "=========Mass_Array[Mass_Index]==========" << Mass_Array[Mass_Index] << endl;
    for(int Cross_Section_IDX=0; Cross_Section_IDX<Cross_Section_Array.size(); Cross_Section_IDX++)
    {
        cout << "=========Cross_Section_Array[Cross_Section_IDX]==========" << Cross_Section_Array[Cross_Section_IDX] << endl;
            double Energy_Ratio_NaI =0;double Energy_Ratio_Pb =0;
            double Energy_Loss_NaI =0; double Energy_Loss_Pb =0;
            for(int Event_Index=0; Event_Index<Simulated_Event_Number; Event_Index++)
            {
                //cout << "=============(NaI)Event: " << Event_Index << endl;
                double *NaI_Parameter = Run_Program(Mass_Array[Mass_Index],NaI_Density,NaI_Atomic_Mass,Cross_Section_Array[Cross_Section_IDX]);
                //cout << "NaI_Energy_Loss_%: "     << NaI_Parameter[1] << endl;
                //cout << "NaI_Energy_Loss_value: " << NaI_Parameter[2] << endl;
                Energy_Ratio_NaI = Energy_Ratio_NaI + NaI_Parameter[1];
                Energy_Loss_NaI = Energy_Loss_NaI + NaI_Parameter[2];

                //if(Final_E_NaI>0.01)NaI_Event = NaI_Event + 1;
            }
             
             
            for(int Event_Index=0; Event_Index<Simulated_Event_Number; Event_Index++)
            {
                //cout << "=============(Pb)Event: " << Event_Index << endl;
                double *Pb_Parameter  = Run_Program(Mass_Array[Mass_Index],Pb_Density,Pb_Atomic_Mass,Cross_Section_Array[Cross_Section_IDX]);
                //cout << "Pb_Energy_Loss_%: "     << Pb_Parameter[1] << endl;
                //cout << "Pb_Energy_Loss_value: " << Pb_Parameter[2] << endl;
                Energy_Ratio_Pb = Energy_Ratio_Pb + Pb_Parameter[1];
                Energy_Loss_Pb  = Energy_Loss_Pb  + Pb_Parameter[2];

                //if(Pb_Parameter[0]>425.){cout << "Big!!Fuck!!" << endl;Pb_Event = Pb_Event + 1;}
            }
            cout << "Energy_Ratio_NaI/Simulated_Event_Number:" << Energy_Ratio_NaI/Simulated_Event_Number << endl;
            cout << "Energy_Ratio_Pb/Simulated_Event_Number:"  << Energy_Ratio_Pb/Simulated_Event_Number << endl;
            cout << "Energy_Loss_NaI/Simulated_Event_Number:"  << Energy_Loss_NaI/Simulated_Event_Number << endl;
            cout << "Energy_Loss_Pb/Simulated_Event_Number:"  << Energy_Loss_Pb/Simulated_Event_Number << endl;

        }
    }
     */
    
    /*
    double Collision_Air_Part=0;
    for(int ATM_Number=0; ATM_Number<19; ATM_Number++)
    {
        double Length_Passed = (atm_table[ATM_Number+1][0]-atm_table[ATM_Number][0])*1e2;//cm
        double Density       = 1e-3*((atm_table[ATM_Number+1][4]+atm_table[ATM_Number][4])*0.5);//g/cm^3
        Collision_Air_Part   = Collision_Air_Part + Length_Passed*Density*1./(unified_atomic_mass_g*(15.)) ;
    }
    
     //Energy Max
    TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.04,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);

    vector<double> Mass_Array_1      = {20.,10,5.,1.,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05};
    vector<double> Energy_Array;
    for(int Mass_Index=0; Mass_Index<Mass_Array_1.size();Mass_Index++){Energy_Array.push_back(Energy_DM(Mass_Array_1[Mass_Index],Max_V*1e3/3e8));}
    
    TGraph * Energy_Max = new TGraph(Mass_Array_1.size(), &Mass_Array_1[0], &Energy_Array[0]);
    Energy_Max->SetMarkerStyle(20);
    Energy_Max->SetMarkerColor(2);
    Energy_Max->SetMarkerColor(2);
    Energy_Max->GetXaxis()->SetRangeUser(0.05,0.5);
    Energy_Max->GetYaxis()->SetRangeUser(0,2);
    Energy_Max->GetYaxis()->SetTitle("Max Energy(keV)");
    Energy_Max->GetXaxis()->SetTitle("M_{#chi}(GeV)");
    Energy_Max->Draw("apl");

    c1->SetLogy();
    c1->SetLogx();
    c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/NaI_Proof/Max_Energy.pdf");
     */
    
    //vector<double> Mass_Array      = {0.1,0.09,0.08,0.07,0.06,0.05};
    //vector<double> Sigma_Array     = {1e-27,2e-27,3e-27,4e-27,5e-27,6e-27,7e-27,8e-27,9e-27,1e-26};
    //vector<double> Mass_Array      = {20,10,5,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2};
    vector<double> Mass_Array      = {10,0.07};
    //vector<double> Mass_Array      = {0.1,0.09,0.08,0.07,0.06,0.05};
    //vector<double> Sigma_Array     = {1e-32,2e-32,3e-32,4e-32,5e-32,6e-32,7e-32,8e-32,9e-32};
    vector<double> Sigma_Array     = {1e-31,2e-31,3e-31,4e-31,5e-31,6e-31,7e-31,8e-31,9e-31};
    //vector<double> Sigma_Array     = {1e-28,2e-28,3e-28,4e-28,5e-28,6e-28,7e-28,8e-28,9e-28};
    //vector<double> Sigma_Array     = {1e-29,2e-29,3e-29,4e-29,5e-29,6e-29,7e-29,8e-29,9e-29};

    //vector<double> Sigma_Array     = {1e-27,2e-27,3e-27,4e-27,5e-27,6e-27,7e-27,8e-27,9e-27};
    //vector<double> Sigma_Array     = {1e-30,2e-30,3e-30,4e-30};

    const int      Material_Number = 3;
    
    double Energy_Loss_per_collision_Array_Water[Mass_Array.size()][Sigma_Array.size()];
    double Energy_Loss_per_collision_Array_Pb[Mass_Array.size()][Sigma_Array.size()];
    double Energy_Loss_per_collision_Array_NaI[Mass_Array.size()][Sigma_Array.size()];

    vector<double> Cross_Section_Set;
    for(int Mass=0; Mass<1; Mass++)//Every Mass I got three values for one Cross-section
    {
    
        double ELoss_Air_for_every_cross_section[Sigma_Array.size()];
        double ELoss_NaI_for_every_cross_section[Sigma_Array.size()];
        double ELoss_Lead_for_every_cross_section[Sigma_Array.size()];
        double ELoss_Water_for_every_cross_section[Sigma_Array.size()];
        double No_ELoss_for_every_cross_section[Sigma_Array.size()];
        
        cout << "============================Mass==============================: "      << Mass_Array[Mass] << endl;
        double  Energy_Max      = Energy_DM(Mass_Array[Mass],Max_V*1e3/3e8);
        double  Final_cross_section =0;
        /*
        for(int Scaling=0; Scaling<10; Scaling++)
        {
            int Scaling_Stop=0;
            if(Scaling_Stop>0)break;
            if(Final_cross_section!=0)break;
         */
            for(int Cross_Section_index=0; Cross_Section_index<Sigma_Array.size(); Cross_Section_index++)
            {
                /*
                if(Final_cross_section!=0)break;
                double CS = Sigma_Array[Cross_Section_index]*TMath::Power(10,Scaling);
                int    Event_Pass=0;
                for(int Event=0; Event<10; Event++)
                {
                    cout << "CS: " << CS << endl;
                    
                    double Collision_Time_Air       = Collision_Air_Part*MFP_from_DCS_Part(2,Max_V,CS,Mass_Array[Mass], 15.);
                    double *V_Aft_Collision_Air     = Velocity_Aft_collision_Bent(Collision_Time_Air,Mass_Array[Mass],CS,Max_V,15.);
                    //cout << "Collision_Time_Air: " << Collision_Time_Air << endl;
                    cout << "Max_V: " << Max_V << endl;
                    double Collision_Time_Water  = 1.0e3*Number_Density_Array[0]*MFP_from_DCS_Part(2,V_Aft_Collision_Air[0],CS,Mass_Array[Mass], Atomic_Mass_Array[0]);
                    double *V_Aft_Collision_Water = Velocity_Aft_collision_Bent(Collision_Time_Water,Mass_Array[Mass],CS,V_Aft_Collision_Air[0],Atomic_Mass_Array[0]);
                    //cout << "Collision_Time_Water: " << Collision_Time_Water << endl;
                    cout << "V_Aft_Collision_Water: " << V_Aft_Collision_Water[0] << endl;

                    double Collision_Time_Lead   = 25*Number_Density_Array[1]*MFP_from_DCS_Part(2,V_Aft_Collision_Water[0],CS,Mass_Array[Mass], Atomic_Mass_Array[1]);
                    double *V_Aft_Collision_Lead  = Velocity_Aft_collision_Bent(Collision_Time_Lead,Mass_Array[Mass],CS,V_Aft_Collision_Water[0],Atomic_Mass_Array[1]);
                    cout << "Collision_Time_Lead: " << Collision_Time_Lead << endl;
                    cout << "V_Aft_Collision_Lead: " << V_Aft_Collision_Lead[0] << endl;
                     
                    double Collision_Time_NaI      = 3.*Number_Density_Array[2]*MFP_from_DCS_Part(2,V_Aft_Collision_Lead[0],CS,Mass_Array[Mass], Atomic_Mass_Array[2]);
                            cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
                    double *V_Aft_Collision_NaI   = Velocity_Aft_collision_Bent(Collision_Time_NaI,Mass_Array[Mass],CS,V_Aft_Collision_Lead[0],Atomic_Mass_Array[2]);
                    cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
                    cout << "V_Aft_Collision_NaI: " << V_Aft_Collision_NaI[0] << endl;
                     
                    double Final_E     = Energy_DM(Mass_Array[Mass],V_Aft_Collision_Lead[0]*1e3/3e8);
                    cout << "Final_E: " << Final_E << endl;
                     if(Final_E>0.16) continue;
                     else
                     {
                         Event_Pass = Event_Pass + 1;
                         cout << "Event_Pass: " << Event_Pass << endl;
                         //cout << "Mass_Array[Mass]: " << Mass_Array[Mass] << endl;
                         //cout << "Sigma_Array[Cross_Section_index]: " << Sigma_Array[Cross_Section_index] << endl;
                     }//else
             }
             if(Event_Pass==10)
             {
                 cout << "SKRRRRRR!!!!!! " << endl;
                 Final_cross_section = CS;
                 Cross_Section_Set.push_back(CS);
                 Scaling_Stop = 1;
                 break;
                 }
                 */
             
                    /*
                double Collision_Time_NaI      = 3.*Number_Density_Array[2]*MFP_from_DCS_Part(2,Max_V,CS,Mass_Array[Mass],Atomic_Mass_Array[2]);
                        cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
                double *V_Aft_Collision_NaI   = Velocity_Aft_collision_Bent(Collision_Time_NaI,Mass_Array[Mass],CS,Max_V,Atomic_Mass_Array[2]);
                cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
                cout << "V_Aft_Collision_NaI: " << V_Aft_Collision_NaI[0] << endl;
                double Final_E     = Energy_DM(Mass_Array[Mass],V_Aft_Collision_NaI[0]*1e3/3e8);
                double Energy_Diff = Energy_Max-Final_E;
                    
                    cout << "Energy_Diff: " << Energy_Diff << endl;
                    //if(Energy_Diff<5.) continue;
                    if(Energy_Diff<0.5) continue;
                    else
                    {
                        Event_Pass = Event_Pass + 1;
                        cout << "Event_Pass: " << Event_Pass << endl;
                        //cout << "Mass_Array[Mass]: " << Mass_Array[Mass] << endl;
                        //cout << "Sigma_Array[Cross_Section_index]: " << Sigma_Array[Cross_Section_index] << endl;
                    }//else
                }
                if(Event_Pass==10)
                {
                    cout << "SKRRRRRR!!!!!! " << endl;
                    Final_cross_section = CS;
                    Cross_Section_Set.push_back(CS);
                    Scaling_Stop = 1;
                    break;
                }
                     */
                /*
                double CS = Sigma_Array[Cross_Section_index];
                double  Energy_Int      = Energy_DM(Mass_Array[Mass],Max_V*1e3/3e8);
                No_ELoss_for_every_cross_section[Cross_Section_index] = Energy_Int;
                cout << "Energy_Int: " << Energy_Int << endl;
                cout << ">>>>>>>>>Sigma_Array[Cross_Section_index]<<<<<<<<<<: "      << Sigma_Array[Cross_Section_index] << endl;
                double Energy_Per_Collision_Air = MFP_from_DCS_Part(0,Max_V,Sigma_Array[Cross_Section_index], Mass_Array[Mass], 15.);
                cout << "Energy_Per_Collision_Air: " << Energy_Per_Collision_Air << endl;
                double Collision_Time_Air       = Collision_Air_Part*MFP_from_DCS_Part(2,Max_V,Sigma_Array[Cross_Section_index], Mass_Array[Mass], 15.);
                cout << "Collision_Time_Air: " << Collision_Time_Air << endl;
                cout << "Total Energy_Air(keV): " << 1e-3*Energy_Per_Collision_Air*Collision_Time_Air << endl;
                double *V_Aft_Collision_Air = Velocity_Aft_collision_Bent(Collision_Time_Air,Mass_Array[Mass],CS,Max_V,15.);

                
                Energy_Int = Energy_Int - 1e-3*Energy_Per_Collision_Air*Collision_Time_Air;
                cout << "After_Air: " << Energy_Int << endl;
                if(Energy_Int>=0)ELoss_Air_for_every_cross_section[Cross_Section_index] = Energy_Int;
                else{ELoss_Air_for_every_cross_section[Cross_Section_index] = 0.;}

                double Velocity_at_Water = Velocity_DM(Mass_Array[Mass],Energy_Int);
                cout << "Velocity_at_Water: " << Velocity_at_Water << endl;
                double Energy_Per_Collision_Water = MFP_from_DCS_Part(0,Velocity_at_Water,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[0]);
                cout << "Energy_Per_Collision_Water: " << Energy_Per_Collision_Water << endl;
                double Collision_Time_Water  = 1.0e3*Number_Density_Array[0]*MFP_from_DCS_Part(2,Velocity_at_Water,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[0]);
                cout << "Collision_Time_Water: " << Collision_Time_Water << endl;
                cout << "Total Energy_Water(keV): " << 1e-3*Energy_Per_Collision_Water*Collision_Time_Water << endl;

                Energy_Int = Energy_Int - 1e-3*Energy_Per_Collision_Water*Collision_Time_Water;
                cout << "After_Water: " << Energy_Int << endl;
                
                if(Energy_Int>0)ELoss_Water_for_every_cross_section[Cross_Section_index] = Energy_Int;
                else{ELoss_Water_for_every_cross_section[Cross_Section_index] = 0;}
                cout << "=======================================" << endl;
                
                double Velocity_at_Lead = Velocity_DM(Mass_Array[Mass],Energy_Int);
                cout << "Velocity_at_Lead: " << Velocity_at_Lead << endl;
                double Energy_Per_Collision_Lead = MFP_from_DCS_Part(0,Velocity_at_Lead,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[1]);
                cout << "Energy_Per_Collision_Lead: " << Energy_Per_Collision_Lead << endl;
                double Collision_Time_Lead       = 25*Number_Density_Array[1]*MFP_from_DCS_Part(2,Velocity_at_Lead,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[1]);
                cout << "Collision_Time_Lead: " << Collision_Time_Lead << endl;
                cout << "Total Energy_Lead(keV): " << 1e-3*Energy_Per_Collision_Lead*Collision_Time_Lead << endl;

                Energy_Int = Energy_Int - 1e-3*Energy_Per_Collision_Lead*Collision_Time_Lead;
                cout << "After_Lead: " << Energy_Int << endl;
                if(Energy_Int>=0)ELoss_Lead_for_every_cross_section[Cross_Section_index] = Energy_Int ;
                else{ELoss_Lead_for_every_cross_section[Cross_Section_index] = 0. ;}
                cout << "=======================================" << endl;


                double Velocity_at_NaI = Velocity_DM(Mass_Array[Mass],Energy_Int);
                cout << "Velocity_at_NaI: " << Velocity_at_NaI << endl;
                double Energy_Per_Collision_NaI = MFP_from_DCS_Part(0,Velocity_at_NaI,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[2]);
                cout << "Energy_Per_Collision_NaI: " << Energy_Per_Collision_NaI << endl;
                double Collision_Time_NaI       = 3.*Number_Density_Array[2]*MFP_from_DCS_Part(2,Velocity_at_NaI,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[2]);
                cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
                cout << "Total Energy_NaI(keV): " << 1e-3*Energy_Per_Collision_NaI*Collision_Time_NaI << endl;
                
                Energy_Int = Energy_Int - 1e-3*Energy_Per_Collision_NaI*Collision_Time_NaI;
                cout << "After_NaI: " << Energy_Int << endl;
                if(Energy_Int>0)ELoss_NaI_for_every_cross_section[Cross_Section_index] = Energy_Int ;
                else{ELoss_NaI_for_every_cross_section[Cross_Section_index] = 0.;}
                cout << "=======================================" << endl;
                 */
                
            }//for(int Cross_Section_index=0; Cross_Section_index<Sigma_Array.size(); Cross_Section_index++)
        //}//for(int Scaling=0; Scaling<6; Scaling++)
        
        /*
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.04,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);
        
        TGraph * Energy_Aft_No = new TGraph(Sigma_Array.size(), &Sigma_Array[0], No_ELoss_for_every_cross_section);
        Energy_Aft_No->SetMarkerStyle(20);
        Energy_Aft_No->SetMarkerColor(4);
        Energy_Aft_No->SetMarkerColor(4);
        Energy_Aft_No->GetXaxis()->SetRangeUser(1e-31,1e-30);
        Energy_Aft_No->GetYaxis()->SetRangeUser(0,Energy_DM(Mass_Array[Mass],Max_V*1e3/3e8)+0.01);
        Energy_Aft_No->GetYaxis()->SetTitle("Energy of DM(keV)");
        Energy_Aft_No->GetXaxis()->SetTitle("#sigma_{SI}(cm^2)");
        Energy_Aft_No->GetYaxis()->SetRangeUser(1e-3,40);
        Energy_Aft_No->Draw("apl");

        TGraph * Energy_Aft_Air = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_Air_for_every_cross_section);
        Energy_Aft_Air->SetMarkerStyle(20);
        Energy_Aft_Air->SetMarkerColor(2);
        Energy_Aft_Air->SetMarkerColor(2);
        //Energy_Aft_Air->Draw("plsame");
        
        TGraph * Energy_Aft_Water = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_Water_for_every_cross_section);
        Energy_Aft_Water->SetMarkerStyle(20);
        Energy_Aft_Water->SetMarkerColor(5);
        Energy_Aft_Water->SetMarkerColor(5);
        Energy_Aft_Water->GetYaxis()->SetRangeUser(0,20);
        //Energy_Aft_Water->Draw("plsame");

        TGraph * Energy_Aft_Lead = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_Lead_for_every_cross_section);
        Energy_Aft_Lead->SetMarkerStyle(20);
        Energy_Aft_Lead->SetMarkerColor(6);
        Energy_Aft_Lead->SetMarkerColor(6);
        Energy_Aft_Lead->GetYaxis()->SetRangeUser(0,20);
        Energy_Aft_Lead->Draw("plsame");

        TGraph * Energy_Aft_NaI = new TGraph(Sigma_Array.size(), &Sigma_Array[0], ELoss_NaI_for_every_cross_section);
        Energy_Aft_NaI->SetMarkerStyle(20);
        Energy_Aft_NaI->SetMarkerColor(3);
        Energy_Aft_NaI->SetMarkerColor(3);
        Energy_Aft_NaI->GetYaxis()->SetRangeUser(0,20);
        Energy_Aft_NaI->Draw("plsame");

        TLegend *leg = new TLegend(0.1,0.1,0.4,0.4);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(Energy_Aft_No,"Max_E","lp");
        //leg->AddEntry(Energy_Aft_Air,"Air","lp");
        //leg->AddEntry(Energy_Aft_Water,"Air+30(M.W.E)","lp");
        leg->AddEntry(Energy_Aft_Lead,"Air+30(M.W.E)+Lead","lp");
        leg->AddEntry(Energy_Aft_NaI,"Air+30(M.W.E)+Lead+NaI","lp");

        c1->SetLogy();
        c1->SetLogx();
        leg->Draw();
        c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/NaI_Proof/10GeV.pdf");
        */
    }
    
    /*
    for(int Cross_Section_Index=0; Cross_Section_Index<Cross_Section_Set.size(); Cross_Section_Index++)
    {
        cout << "Mass_Array: " << Mass_Array[Cross_Section_Index] << endl;
        cout << "Cross_Section: " << Cross_Section_Set[Cross_Section_Index] << endl;
    }

    for(int Cross_Section_Index=0; Cross_Section_Index<Cross_Section_Set.size(); Cross_Section_Index++)
    {
        cout << Mass_Array[Cross_Section_Index] << "," << Cross_Section_Set[Cross_Section_Index] << "," << endl;
    }
     */
    /*
    cout << "================Low-Mass NaI Issue========================: "  << endl;

        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.04,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);

    cout << "================Hello========================: "  << endl;

    vector<double> Sigma_Array_1          = {1e-27,5e-27,1e-26,5e-26,1e-25,5e-25,1e-24,5e-24};
    vector<double> Mass_Array_1           = {0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05};
    string Mass_Point[10]= {"0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05"};

    cout << "================Hello========================: "  << endl;


    for(int Mass_Index=0; Mass_Index<Mass_Array_1.size(); Mass_Index++)
    {
        vector<double> Total_Energy_Loss_NaI;

        for(int Sigma_index=0; Sigma_index<Sigma_Array_1.size(); Sigma_index++)
        {
            double Energy_Per_Collision_NaI = MFP_from_DCS_Part(0,Max_V,Sigma_Array_1[Sigma_index], Mass_Array_1[Mass_Index], Atomic_Mass_Array[2]);
            cout << "Energy_Per_Collision_NaI: " << Energy_Per_Collision_NaI << endl;
            double Collision_Time_NaI       = 10*Number_Density_Array[2]*MFP_from_DCS_Part(2,Max_V,Sigma_Array_1[Sigma_index], Mass_Array_1[Mass_Index], Atomic_Mass_Array[2]);
            cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
            double Total_Energy             = 1e-3*Energy_Per_Collision_NaI*Collision_Time_NaI;//keV
            cout << "Total_Energy: " << Total_Energy << endl;
            if(Total_Energy>=0)Total_Energy_Loss_NaI.push_back(Total_Energy);
            else{Total_Energy_Loss_NaI.push_back(0.);};

        }

        cout << "================Hello========================: "  << endl;

        TGraph * Energy_Loss_NaI = new TGraph(Sigma_Array_1.size(), &Sigma_Array_1[0], &Total_Energy_Loss_NaI[0]);
        Energy_Loss_NaI->SetMarkerStyle(20);
        Energy_Loss_NaI->SetMarkerColor(2);
        Energy_Loss_NaI->SetMarkerColor(2);
        Energy_Loss_NaI->GetXaxis()->SetRangeUser(1e-27,5e-24);
        Energy_Loss_NaI->GetYaxis()->SetRangeUser(0,1e4);
        Energy_Loss_NaI->GetYaxis()->SetTitle("Energy Loss(eV)");
        Energy_Loss_NaI->GetXaxis()->SetTitle("M_{#chi}(GeV)");
        Energy_Loss_NaI->Draw("apl");

        c1->SetLogx();
        c1->SetLogy();
        c1->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/NaI_Proof/NaI_%s.pdf",Mass_Point[Mass_Index].c_str()));

    }
    */
    
        /*
        cout << "================High-Mass NaI Issues========================: "  << endl;

        cout << "================Hello========================: "  << endl;

        double Reference_cross_section        = 1e-31;
        vector<double> Mass_Array_1           = {20,10,8,6,4,2,1,0.9,0.7,0.5,0.3,0.1,0.09,0.07,0.05};
        //vector<double> Mass_Array_1           = {10.,8.,6.,4.,2.,1.};
        //string Mass_Point[10]= {"0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05"};

        cout << "================Hello========================: "  << endl;

        vector<double> Cross_Section_from_NaI;
    
        for(int Mass_Index=0; Mass_Index<Mass_Array_1.size(); Mass_Index++)
        {
            
            cout << "Mass_Array[Mass_Index]: " << Mass_Array_1[Mass_Index] << endl;
            double  Energy_Max      = Energy_DM(Mass_Array_1[Mass_Index],Max_V*1e3/3e8);//keV
            double Energy_Per_Collision_NaI = MFP_from_DCS_Part(0,Max_V,Reference_cross_section, Mass_Array_1[Mass_Index], Atomic_Mass_Array[2]);
            cout << "Energy_Per_Collision_NaI: " << Energy_Per_Collision_NaI << endl;


            double Collision_Time_NaI       = 3.*Number_Density_Array[2]*MFP_from_DCS_Part(2,Max_V,Reference_cross_section, Mass_Array_1[Mass_Index], Atomic_Mass_Array[2]);
            cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
            double Total_Energy             = 1e-3*Energy_Per_Collision_NaI*Collision_Time_NaI;//keV
            //double Scaling                  = 1.0/(Total_Energy);
            //double Scaling                  = 0.1/(Total_Energy);
            cout << "Energy_Max: "     << Energy_Max << endl;
            cout << "(Total_Energy): " << (Total_Energy) << endl;
              double Scaling                  = Energy_Max/(Total_Energy);
            cout << "Scaling: " << Scaling << endl;
            cout << "Final: " << Reference_cross_section*Scaling << endl;
            cout << "Max: " << Energy_DM(Mass_Array_1[Mass_Index],Max_V*1e3/3e8) << endl;
            cout << "================Hello========================: "  << endl;
            Cross_Section_from_NaI.push_back(Reference_cross_section*Scaling);
             

        }
        for(int Mass_Index=0; Mass_Index<Mass_Array_1.size(); Mass_Index++)
        {
            cout << Mass_Array_1[Mass_Index] << "," << Cross_Section_from_NaI[Mass_Index] << "," << endl;
        }
         */
    
        /*
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.04,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);
    
        vector<double> Mass_Array_1  = {20,10,8,6,4,2,1,0.9,0.7,0.5,0.3,0.1,0.09,0.07,0.05};
        vector<double> Cement_Energy_Loss_per_collision_array;
        vector<double> Pb_Energy_Loss_per_collision_array;
        vector<double> NaI_Energy_Loss_per_collision_array;
        double Reference_cross_section        = 1e-31;

        for(int Mass_Index=0; Mass_Index<Mass_Array_1.size(); Mass_Index++)
        {
            double Energy_Per_Collision_Cement = MFP_from_DCS_Part(0,Max_V,Reference_cross_section, Mass_Array_1[Mass_Index], Atomic_Mass_Array[0]);
            cout << "Energy_Per_Collision_Cement: " << Energy_Per_Collision_Cement << endl;
            double Energy_Per_Collision_Pb  = MFP_from_DCS_Part(0,Max_V,Reference_cross_section, Mass_Array_1[Mass_Index], Atomic_Mass_Array[1]);
            cout << "Energy_Per_Collision_Pb: " << Energy_Per_Collision_Pb << endl;
            double Energy_Per_Collision_NaI = MFP_from_DCS_Part(0,Max_V,Reference_cross_section, Mass_Array_1[Mass_Index], Atomic_Mass_Array[2]);
            cout << "Energy_Per_Collision_NaI: " << Energy_Per_Collision_NaI << endl;
            
            Pb_Energy_Loss_per_collision_array.push_back(Energy_Per_Collision_Pb*1E-3);
            NaI_Energy_Loss_per_collision_array.push_back(Energy_Per_Collision_NaI*1E-3);
            Cement_Energy_Loss_per_collision_array.push_back(Energy_Per_Collision_Cement*1E-3);
        }
        
        TGraph *Cement_ELoss_per_collision = new TGraph(Mass_Array_1.size(), &Mass_Array_1[0], &Cement_Energy_Loss_per_collision_array[0]);
        Cement_ELoss_per_collision->SetMarkerStyle(20);
        Cement_ELoss_per_collision->SetMarkerColor(4);
        Cement_ELoss_per_collision->SetMarkerColor(4);
        Cement_ELoss_per_collision->GetXaxis()->SetRangeUser(0.05,20);
        Cement_ELoss_per_collision->GetYaxis()->SetRangeUser(0,2);
        Cement_ELoss_per_collision->GetYaxis()->SetTitle("Energy Loss per collision(keV)");
        Cement_ELoss_per_collision->GetXaxis()->SetTitle("M_{#chi}(GeV)");
        Cement_ELoss_per_collision->Draw("apl");

        TGraph *NaI_ELoss_per_collision = new TGraph(Mass_Array_1.size(), &Mass_Array_1[0], &NaI_Energy_Loss_per_collision_array[0]);
        NaI_ELoss_per_collision->SetMarkerStyle(20);
        NaI_ELoss_per_collision->SetMarkerColor(3);
        NaI_ELoss_per_collision->SetMarkerColor(3);
        NaI_ELoss_per_collision->GetXaxis()->SetRangeUser(0.05,20);
        NaI_ELoss_per_collision->GetYaxis()->SetRangeUser(0,2);
        NaI_ELoss_per_collision->Draw("plsame");

        TGraph *Pb_ELoss_per_collision = new TGraph(Mass_Array_1.size(), &Mass_Array_1[0], &Pb_Energy_Loss_per_collision_array[0]);
        Pb_ELoss_per_collision->SetMarkerStyle(20);
        Pb_ELoss_per_collision->SetMarkerColor(2);
        Pb_ELoss_per_collision->SetMarkerColor(2);
        Pb_ELoss_per_collision->GetXaxis()->SetRangeUser(0.05,0.5);
        Pb_ELoss_per_collision->GetYaxis()->SetRangeUser(0,2);
        Pb_ELoss_per_collision->Draw("plsame");

        TLegend *leg = new TLegend(0.7,0.1,0.9,0.4);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(Cement_ELoss_per_collision,"Cement","lp");
        leg->AddEntry(NaI_ELoss_per_collision,"NaI","lp");
        leg->AddEntry(Pb_ELoss_per_collision ,"Pb","lp");
        leg->Draw();
        
        c1->SetLogx();
        c1->SetLogy();
        c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/NaI_Proof/Energy_Loss_per_collision.pdf");
        */
    
}//End the main function

//double Final_E_NaI    = Energy_DM(WIMP_Mass,NaI_Parameter[0]*1e3/3e8);//keV

    /*
     //Energy Max
    TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.04,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);

    vector<double> Mass_Array_1      = {0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05};
    vector<double> Energy_Array;
    for(int Mass_Index=0; Mass_Index<Mass_Array_1.size();Mass_Index++){Energy_Array.push_back(Energy_DM(Mass_Array_1[Mass_Index],Max_V*1e3/3e8));}
    
    TGraph * Energy_Max = new TGraph(Mass_Array_1.size(), &Mass_Array_1[0], &Energy_Array[0]);
    Energy_Max->SetMarkerStyle(20);
    Energy_Max->SetMarkerColor(2);
    Energy_Max->SetMarkerColor(2);
    Energy_Max->GetXaxis()->SetRangeUser(0.05,0.5);
    Energy_Max->GetYaxis()->SetRangeUser(0,2);
    Energy_Max->GetYaxis()->SetTitle("Max Energy(keV)");
    Energy_Max->GetXaxis()->SetTitle("M_{#chi}(GeV)");
    Energy_Max->Draw("apl");

    c1->SetLogx();
    c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/NaI_Proof/Max_Energy.pdf");
     */
    /*
    TGraph * g_Energy_Loss_per_collision_Array = new TGraph((int)Mass_Array.size(), Energy_H_Check, &TGraph_dsigma_dT_Xe_Ge_Check[0]);
    g_Xe_Ge->SetMarkerStyle(20);
    g_Xe_Ge->SetMarkerColor(2);
    g_Xe_Ge->SetMarkerColor(2);
    g_Xe_Ge->GetYaxis()->SetRangeUser(0,20);
    g_Xe_Ge->GetXaxis()->SetTitle("Recoil Energy(eV)");
    g_Xe_Ge->GetXaxis()->SetTitle("Log10(Ratio of DCS)");
    g_Xe_Ge->Draw("apl");
    
    TGraph * g_Xe_H = new TGraph(7, Energy_H_Check, &TGraph_dsigma_dT_Xe_H_Check[0]);
    g_Xe_H->SetMarkerStyle(20);
    g_Xe_H->SetMarkerColor(3);
    g_Xe_H->Draw("plsame");
     */

/*
for(int Cross_Section_index=0; Cross_Section_index<Sigma_Array.size(); Cross_Section_index++)
{
    cout << "Sigma_Array[Cross_Section_index]: " << Sigma_Array[Cross_Section_index] << endl;
    double NaI_2GeV    = total_Sigma(1,779,Sigma_Array[Cross_Section_index],  2,NaI_Atomic_Mass);
    cout << "NaI_2GeV: " << NaI_2GeV << endl;
    cout << "==================================" << endl;
    double NaI_0P2GeV  = total_Sigma(1,779,Sigma_Array[Cross_Section_index],0.2,NaI_Atomic_Mass);
    cout << "NaI_0P2GeV: " << NaI_0P2GeV << endl;
    cout << "==================================" << endl;
    cout << "Ratio: " << (NaI_0P2GeV)/(NaI_2GeV) << endl;
}
 */
/*
 
double MFP_NaI = Mean_free_Path_check(NaI_Density,NaI_Atomic_Mass,Sigma_SI_NaI);
double MFP_Pb  = Mean_free_Path_check(Pb_Density ,Pb_Atomic_Mass ,Sigma_SI_Pb) ;
cout << "MFP_NaI: " << MFP_NaI << endl;
cout << "MFP_Pb:  " << MFP_Pb  << endl;
 */

    /*
    for(int Material_Index=0; Material_Index<3; Material_Index++)
    {
MFP_Array.push_back(  Mean_free_Path_check(Mass_Array[Mass], Density_Array[Material_Index], Atomic_Mass_Array[Material_Index], Sigma_Array[Cross_Section_index]) );
    }
    cout << "Water_MFP: " << MFP_Array[0] << endl;
    cout << "Pb_MFP: "    << MFP_Array[1] << endl;
    cout << "NaI_MFP: "   << MFP_Array[2] << endl;
    for(int Material_Index=0; Material_Index<3; Material_Index++)
    {
dEdX_Array.push_back(  Number_Density_Array[Material_Index]*MFP_from_DCS_Part(1,Max_V,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[Material_Index]) );
    }
    cout << "Water_dEdX: " << dEdX_Array[0] << endl;
    cout << "Pb_dEdX: "    << dEdX_Array[1] << endl;
    cout << "NaI_dEdX: "   << dEdX_Array[2] << endl;
     */

/*
//Prove that the mean value of energy loss is reasonable
gRandom = new TRandom3(0);
gRandom->SetSeed(0);

TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.04,"XY");
gStyle->SetTitleFont(62,"XY");
gStyle->SetLegendFont(62);

TF1 *f3 = new TF1("f3","fdsigma_dT_keV([0],[1],[2],[3],x)",0,max_recoil_A_keV(0.1,779*1e3/3e8,34.));
f3->SetParameter(0,0.1);f3->SetParameter(1,1e-31);f3->SetParameter(2,(779*1e3/3e8));f3->SetParameter(3,34.);

TH1F *Try = new TH1F("Try","Try",1000,0,max_recoil_A_keV(0.1,779*1e3/3e8,34));
for(int KKK=0; KKK<50000; KKK++)
{
    double Random_Energy= f3->GetRandom();
    cout << "Random_Energy: " << Random_Energy << endl;
    Try->Fill(Random_Energy);
}
Try->SetMarkerStyle(20);
Try->SetMarkerColor(2);
Try->SetMarkerColor(2);
Try->GetXaxis()->SetRangeUser(0,max_recoil_A_keV(0.1,779*1e3/3e8,34.));
Try->GetYaxis()->SetRangeUser(0,200);
Try->Draw("Hist");

double Energy_Loss = MFP_from_DCS_Part(0,779,1e-31, 0.1, 34.)*1e-3;
cout << "Energy_Loss: " << Energy_Loss << endl;
*/
