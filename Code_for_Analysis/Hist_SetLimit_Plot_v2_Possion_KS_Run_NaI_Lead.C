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
//const double       WIMP_Mass              = 0.2;
//Constant
const double       NaI_Density            = 3.67;//3.67(g/cm^3)
const double       NaI_Atomic_Mass        = 22.98*0.5+126*0.5;//
const double       Pb_Density             = 11.29;//3.67(g/cm^3)
const double       Pb_Atomic_Mass         = 207.2;//

const double       Fixed_Length           = 30.;//cm

double Mean_free_Path_check(double WIMP_Mass, double Density, double Atomic_Mass, double Sigma_SI)
{
    double MFP_Calculated = 1./((Density)/(unified_atomic_mass_g*(Atomic_Mass))*total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));
    double Collision_Time   = (Density*(MFP_Calculated))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));
    double Collision_Time_2 = (Density*(15.))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));

    cout << "Atomic_Mass: " << Atomic_Mass << endl;
    cout << "Once: " << (Density*(3))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,1.5e-27,WIMP_Mass,Atomic_Mass)) << endl;;
    cout << "MFP_Calculated: "    << MFP_Calculated << endl;
    cout << "Collision_Time: "    << Collision_Time << endl;
    cout << "Collision_Time_2: "  << Collision_Time_2 << endl;
    return MFP_Calculated;//cm
}

double *Run_Program(double WIMP_Mass, double Density, double Atomic_Mass, double Sigma_SI)//Density, Atomic mass
{
    static double Array[2];
    double Initial_V                = Max_V;
    double Last_V                   = Max_V;
    double Energy_Loss_Percentage_total   = 0;
    double Collision_Time = (Density*(Fixed_Length))/(unified_atomic_mass_g*(Atomic_Mass))*(total_Sigma(1,Max_V,Sigma_SI,WIMP_Mass,Atomic_Mass));
    double Collision_Time_Int       = 0;
    while((int)Collision_Time>0)
    {
        double *V_aft      = Velocity_Aft_collision_Bent(1,WIMP_Mass,Sigma_SI,Initial_V,Atomic_Mass);
        Initial_V          = V_aft[0];
        double Energy_Loss = Energy_DM(WIMP_Mass,Last_V*1e3/3e8)-Energy_DM(WIMP_Mass,Initial_V*1e3/3e8);//keV
        double Energy_Loss_Percentage = (Energy_Loss)/Energy_DM(WIMP_Mass,Last_V*1e3/3e8);
        Energy_Loss_Percentage_total = Energy_Loss_Percentage_total + Energy_Loss_Percentage;
        Last_V             = Initial_V;
        Collision_Time     = Collision_Time - 1;//
        Collision_Time_Int = Collision_Time_Int + 1;
    }
    Energy_Loss_Percentage_total = Energy_Loss_Percentage_total/(Collision_Time_Int);
    Array[0]=Initial_V;Array[1]=Energy_Loss_Percentage_total;
    return Array;
}


//Run the program for the individual index and the simulated number of events
void Hist_SetLimit_Plot_v2_Possion_KS_Run_NaI_Lead()
{
    const double       Same_Sigma             = 1e-27;
    const double       Sigma_SI_NaI           = Same_Sigma;
    const double       Sigma_SI_Pb            = Same_Sigma;

    //vector<double> Mass_Array={2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05};
    vector<double> Mass_Array={2,0.2};

    int NaI_Event=0;int Pb_Event=0;
    
    
    
    for(int Mass_Index=0; Mass_Index<Mass_Array.size(); Mass_Index++)
    {
        double Energy_Ratio_NaI =0;double Energy_Ratio_Pb =0;
        for(int Event_Index=0; Event_Index<Simulated_Event_Number; Event_Index++)
        {
            cout << "=============(NaI)Event: " << Event_Index << endl;
            double *NaI_Parameter = Run_Program(Mass_Array[Mass_Index],NaI_Density,NaI_Atomic_Mass,Sigma_SI_NaI);
            cout << "NaI_Parameter[1]: " << NaI_Parameter[1] << endl;
            Energy_Ratio_NaI = Energy_Ratio_NaI + NaI_Parameter[1];
            //if(Final_E_NaI>0.01)NaI_Event = NaI_Event + 1;
        }
         
         
        for(int Event_Index=0; Event_Index<Simulated_Event_Number; Event_Index++)
        {
            cout << "=============(Pb)Event: " << Event_Index << endl;
            double *Pb_Parameter  = Run_Program(Mass_Array[Mass_Index],Pb_Density,Pb_Atomic_Mass,Sigma_SI_Pb);
            cout << "Pb_Parameter[1]: " << Pb_Parameter[1] << endl;
            Energy_Ratio_Pb = Energy_Ratio_Pb + Pb_Parameter[1];
            //if(Pb_Parameter[0]>425.){cout << "Big!!Fuck!!" << endl;Pb_Event = Pb_Event + 1;}
        }
        cout << "Energy_Ratio_NaI/Simulated_Event_Number:" << Energy_Ratio_NaI/Simulated_Event_Number << endl;
        cout << "Energy_Ratio_Pb/Simulated_Event_Number:"  << Energy_Ratio_Pb/Simulated_Event_Number << endl;
    }
    /*
     
    double MFP_NaI = Mean_free_Path_check(NaI_Density,NaI_Atomic_Mass,Sigma_SI_NaI);
    double MFP_Pb  = Mean_free_Path_check(Pb_Density ,Pb_Atomic_Mass ,Sigma_SI_Pb) ;
    cout << "MFP_NaI: " << MFP_NaI << endl;
    cout << "MFP_Pb:  " << MFP_Pb  << endl;
     */
}

//double Final_E_NaI    = Energy_DM(WIMP_Mass,NaI_Parameter[0]*1e3/3e8);//keV
