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
    const double       Same_Sigma             = 4e-28;
    const double       Sigma_SI_NaI           = Same_Sigma;
    const double       Sigma_SI_Pb            = Same_Sigma;

    vector<double> Cross_Section_Array={8e-28,2e-27};

    //vector<double> Mass_Array={2,0.2};
    
    int NaI_Event=0;int Pb_Event=0;
    
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
    double Collision_Air_Part=0;
    for(int ATM_Number=0; ATM_Number<19; ATM_Number++)
    {
        double Length_Passed = (atm_table[ATM_Number+1][0]-atm_table[ATM_Number][0])*1e2;//cm
        double Density       = 1e-3*((atm_table[ATM_Number+1][4]+atm_table[ATM_Number][4])*0.5);//g/cm^3
        //cout << "Length_Passed: " << Length_Passed << endl;
        //cout << "Density: " << Density << endl;
        Collision_Air_Part   = Collision_Air_Part + Length_Passed*Density*1./(unified_atomic_mass_g*(15.)) ;
    }
    
    vector<double> Mass_Array      = {0.1,0.09,0.08,0.07,0.06,0.05};
    vector<double> Sigma_Array     = {1e-27,2e-27,3e-27,4e-27,5e-27,6e-27,7e-27,8e-27,9e-27,1e-26};
    const int      Material_Number = 3;
    
    double Energy_Loss_per_collision_Array_Water[Mass_Array.size()][Sigma_Array.size()];
    double Energy_Loss_per_collision_Array_Pb[Mass_Array.size()][Sigma_Array.size()];
    double Energy_Loss_per_collision_Array_NaI[Mass_Array.size()][Sigma_Array.size()];

    for(int Mass=0; Mass<1; Mass++)//Every Mass I got three values for one Cross-section
    {
    
        double ELoss_Air_for_every_cross_section[Sigma_Array.size()];
        double ELoss_NaI_for_every_cross_section[Sigma_Array.size()];
        double ELoss_Lead_for_every_cross_section[Sigma_Array.size()];
        double ELoss_Water_for_every_cross_section[Sigma_Array.size()];
        double No_ELoss_for_every_cross_section[Sigma_Array.size()];

        cout << "============================Mass==============================: "      << Mass_Array[Mass] << endl;
        for(int Cross_Section_index=0; Cross_Section_index<Sigma_Array.size(); Cross_Section_index++)
        {
            double  Energy_Int      = Energy_DM(Mass_Array[Mass],Max_V*1e3/3e8);
            No_ELoss_for_every_cross_section[Cross_Section_index] = Energy_Int;
            cout << "Energy_Int: " << Energy_Int << endl;
            cout << ">>>>>>>>>Sigma_Array[Cross_Section_index]<<<<<<<<<<: "      << Sigma_Array[Cross_Section_index] << endl;
            double Energy_Per_Collision_Air = MFP_from_DCS_Part(0,Max_V,Sigma_Array[Cross_Section_index], Mass_Array[Mass], 15.);
            cout << "Energy_Per_Collision_Air: " << Energy_Per_Collision_Air << endl;
            double Collision_Time_Air       = Collision_Air_Part*MFP_from_DCS_Part(2,Max_V,Sigma_Array[Cross_Section_index], Mass_Array[Mass], 15.);
            cout << "Collision_Time_Air: " << Collision_Time_Air << endl;
            cout << "Total Energy_Air(keV): " << 1e-3*Energy_Per_Collision_Air*Collision_Time_Air << endl;

            Energy_Int = Energy_Int - 1e-3*Energy_Per_Collision_Air*Collision_Time_Air;
            if(Energy_Int>=0)ELoss_Air_for_every_cross_section[Cross_Section_index] = Energy_Int;
            else{ELoss_Air_for_every_cross_section[Cross_Section_index] = 0.;}
            cout << "=======================================" << endl;

            double Velocity_at_Water = Velocity_DM(Mass_Array[Mass],Energy_Int);
            cout << "Velocity_at_Water: " << Velocity_at_Water << endl;
            double Energy_Per_Collision_Water = MFP_from_DCS_Part(0,Velocity_at_Water,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[0]);
            cout << "Energy_Per_Collision_Water: " << Energy_Per_Collision_Water << endl;
            double Collision_Time_Water  = 1.0e3*Number_Density_Array[0]*MFP_from_DCS_Part(2,Velocity_at_Water,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[0]);
            cout << "Collision_Time_Water: " << Collision_Time_Water << endl;
            cout << "Total Energy_Water(keV): " << 1e-3*Energy_Per_Collision_Water*Collision_Time_Water << endl;
            
            Energy_Int = Energy_Int - 1e-3*Energy_Per_Collision_Air*Collision_Time_Air;
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
            if(Energy_Int>=0)ELoss_Lead_for_every_cross_section[Cross_Section_index] = Energy_Int ;
            else{ELoss_Lead_for_every_cross_section[Cross_Section_index] = 0. ;}
            cout << "=======================================" << endl;


            double Velocity_at_NaI = Velocity_DM(Mass_Array[Mass],Energy_Int);
            cout << "Velocity_at_NaI: " << Velocity_at_NaI << endl;
            double Energy_Per_Collision_NaI = MFP_from_DCS_Part(0,Velocity_at_NaI,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[2]);
            cout << "Energy_Per_Collision_NaI: " << Energy_Per_Collision_NaI << endl;
            double Collision_Time_NaI       = 10*Number_Density_Array[2]*MFP_from_DCS_Part(2,Velocity_at_NaI,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[2]);
            cout << "Collision_Time_NaI: " << Collision_Time_NaI << endl;
            cout << "Total Energy_NaI(keV): " << 1e-3*Energy_Per_Collision_NaI*Collision_Time_NaI << endl;
            
            Energy_Int = Energy_Int - 1e-3*Energy_Per_Collision_NaI*Collision_Time_NaI;
            if(Energy_Int>0)ELoss_NaI_for_every_cross_section[Cross_Section_index] = Energy_Int ;
            else{ELoss_NaI_for_every_cross_section[Cross_Section_index] = 0.;}
            cout << "=======================================" << endl;



            /*
            double First_dEdX = Number_Density_Array[Material_Index];
            double Last_dEdX  = MFP_from_DCS_Part(0,Max_V,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[Material_Index]);
            double dE_dX      = First_dEdX*Last_dEdX;
             */
    //Energy_Loss_per_collision_Array[Mass][Material_Index] = (  MFP_from_DCS_Part(0,Max_V,Sigma_Array[Cross_Section_index], Mass_Array[Mass], Atomic_Mass_Array[Material_Index]) );
            //cout << "Water_ELoss_per_collision "  << Energy_Loss_per_collision_Array[0] << endl;
            //cout << "Pb_ELoss_per_collision: "    << Energy_Loss_per_collision_Array[1] << endl;
            //cout << "NaI_ELoss_per_collision: "   << Energy_Loss_per_collision_Array[2] << endl;
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
        cout << "=====================================================================================: "  << endl;
    }
    
    
    cout << "================Hello========================: "  << endl;

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
}
*/
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
