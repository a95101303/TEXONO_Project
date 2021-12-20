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

const int          Simulated_Event_Number = 500;
const double       Max_V                  = 779;//(km/s)
const double       WIMP_Mass              = 1;
//Constant
const double       NaI_Density            = 3.67;//3.67(g/cm^3)
const double       NaI_Atomic_Mass        = 22.98*0.5+126*0.5;//
const double       Pb_Density             = 11.29;//3.67(g/cm^3)
const double       Pb_Atomic_Mass         = 207.2;//


double *Run_Program(double Density, double Atomic_Mass, double Sigma_SI)//Density, Atomic mass
{
    static double Array[2];
    double Initial_V      = Max_V;
    double Length_added   = 0;
    int    Collision_Time = 0;
    double Collision_Time = (Density*())/(unified_atomic_mass_g*(APb))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,APb))

    while(Energy_DM(WIMP_Mass,Initial_V*1e3/3e8)>0.01)
    {
        double Length_DM = Length_for_asking_the_collision(0.001,WIMP_Mass,Initial_V,Sigma_SI,Density,Atomic_Mass);//
        double *V_aft    = Velocity_Aft_collision_Bent(1,WIMP_Mass,Sigma_SI,Initial_V,Atomic_Mass);
        Initial_V        = V_aft[0];
        Length_added     = Length_added + Length_DM;//km
    }
    Array[0]=Initial_V;Array[1]=Length_added*1e5;
    return Array;
}


//Run the program for the individual index and the simulated number of events
void Hist_SetLimit_Plot_v2_Possion_KS_Run_NaI_Lead()
{
    const double       Sigma_SI_NaI           = 1e-26;
    const double       Sigma_SI_Pb            = 1e-26;

    int NaI_Event=0;int Pb_Event=0;
    for(int kkk=0; kkk<Simulated_Event_Number; kkk++)
    {
        cout << "=============(NaI)Event: " << kkk << endl;
        double *NaI_Parameter = Run_Program(NaI_Density,NaI_Atomic_Mass,Sigma_SI_NaI);
        double Final_E_NaI    = Energy_DM(WIMP_Mass,NaI_Parameter[0]*1e3/3e8);//keV
        if(Final_E_NaI<0.01)NaI_Event = NaI_Event + 1;
    }
    
    for(int kkk=0; kkk<Simulated_Event_Number; kkk++)
    {
        cout << "=============(Pb)Event: " << kkk << endl;
        double *Pb_Parameter  = Run_Program(Pb_Density,Pb_Atomic_Mass,Sigma_SI_Pb);
        double Final_E_Pb     = Energy_DM(WIMP_Mass,Pb_Parameter[0]*1e3/3e8);//keV
        if(Final_E_Pb<0.01)Pb_Event = Pb_Event + 1;
    }
    cout << "NaI_Event/Simulated_Event_Number:" << NaI_Event/Simulated_Event_Number << endl;
    cout << "Pb_Event/Simulated_Event_Number:"  << Pb_Event/Simulated_Event_Number << endl;
}

