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
#include "TGraph2D.h"

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

#include "VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"
#include "velocity_distribution_2000_Ave.h"
#include "dsigma_dT2.h"

//200eV threshold ==> 1.009keV
//1.009keV        ==> 201.183(km/s)
void Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Material(int Bent_or_not, int Index_Mass, int Simulated_Event_Number, int Index_Sigma, double Sigma_SI)
{
    //Bent_or_not=0 -> Non-scattering case, Bent_or_not=1 -> scattering case
    //Constant
    string Mass_Point[4]={"20","2","0P2","0P05"};
    double WIMP_Mass_Array[4]={20,2,0.2,0.05};//12 for TEXONO
    //Start
    double WIMP_Mass = WIMP_Mass_Array[Index_Mass];
    double DM_mx = WIMP_Mass_Array[Index_Mass];
    cout << "Sigma_SI: " << Sigma_SI << endl;
    cout << "WIMP_Mass: " << WIMP_Mass << endl;

    int Material_Range=0;
    double Density_for_Material[6]={Density_of_Cement,1,11.34,7.86,2.34,8.96};// Cement,H2O,Pb,Fe,B,Cu  g/cm^3
    string ATOM_number_for_Material[6]={"Cement","AH2O","APb","AFe","AB","ACu"};
    double Index_Atom[6]={Weighted_Atomic_Number_Cement,AH2O,APb,AFe,AB,ACu};
    //double Thickness_of_Material[6]={1,1,1,1,1,1};
    //TH1F   *Flux_HIST_Aft_Collision_EARTH = new TH1F("Flux_HIST_Aft_Collision_EARTH","Flux_HIST_Aft_Collision_EARTH",2000,0,791);

    int Event_Number_Check=0;
    double Check_ZERO_COLLISION=0;double Check_ZERO_COLLISION_1=0;

    int kkk=0; int jjj=0; int MMM=0;
    int BorNB=Bent_or_not;
    
    TTree *t1 = new TTree("t1","Information");
    double sigma_si,mx ;int ev,Arrival_air,Arrival_earth,Bent_or_Not;
    double V_Int_A, V_End_A, V_Int_E, V_End_E;
    double Oringal_Length_Air,Path_Length_Air,Oringal_Length_Earth,Path_Length_Earth;
    double Collision_Time_Earth,Collision_Time_Air;
    double Original_Bent_Comparison_Ratio_Air_P,Original_Bent_Comparison_Ratio_Earth_P;
    double Energy_Loss_Percentage_lf;
    
    t1->Branch("sigma_si",&sigma_si,"sigma_si/D");
    t1->Branch("mx",&mx,"mx/D");
    t1->Branch("V_Int_A",&V_Int_A,"V_Int_A/D");
    t1->Branch("V_End_A",&V_End_A,"V_End_A/D");
    t1->Branch("V_Int_E",&V_Int_E,"V_Int_E/D");
    t1->Branch("V_End_E",&V_End_E,"V_End_E/D");

    t1->Branch("ev",&ev,"ev/I");
    t1->Branch("Arrival_air",&Arrival_air,"Arrival_air/I");
    t1->Branch("Arrival_earth",&Arrival_earth,"Arrival_earth/I");
    t1->Branch("Bent_or_Not",&Bent_or_Not,"Bent_or_Not/I");
    
    t1->Branch("Oringal_Length_Air",&Oringal_Length_Air,"Oringal_Length_Air/D");
    t1->Branch("Path_Length_Air",&Path_Length_Air,"Path_Length_Air/D");
    t1->Branch("Oringal_Length_Earth",&Oringal_Length_Earth,"Oringal_Length_Earth/D");
    t1->Branch("Path_Length_Earth",&Path_Length_Earth,"Path_Length_Earth/D");
    
    t1->Branch("Collision_Time_Earth",&Collision_Time_Earth,"Collision_Time_Earth/D");
    t1->Branch("Collision_Time_Air",&Collision_Time_Air,"Collision_Time_Air/D");

    t1->Branch("Energy_Loss_Percentage_lf",&Energy_Loss_Percentage_lf,"Energy_Loss_Percentage_lf/D");


    double Max_V=799.135;
    
    while(jjj<2500)
    //while(jjj<50)
    //while((MMM<500 and Bent_or_not_to_be_Bent==1) or (jjj<500 and Bent_or_not_to_be_Bent==0 ))
    {
        
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;

        sigma_si = Sigma_SI;
        cout << "sigma_si : " << sigma_si << endl;
        mx = WIMP_Mass;
        cout << "mx : " << mx << endl;
        ev = jjj;
        Bent_or_Not = BorNB;
        V_Int_A = Max_V;
                
        double Velocity_X =    1e-50;
        double Velocity_Y =    1e-50;
        double Velocity_Z =   abs(1);
        double V_M[3] = {Velocity_X,Velocity_Y,Velocity_Z};//Velocity_Matrix
                                        
        double Dark_Matter_Energy = Energy_DM(DM_mx,Max_V*1e3/3e8);//KeV
        double Dark_Matter_Velocity = Velocity_DM(DM_mx,Dark_Matter_Energy);//KeV
                 
        //Air_Value
        double *BV =  KS_Collision_Time_Building_Aft_velocity_with_angle(BorNB,Sigma_SI,Max_V,DM_mx,V_M,1,Density_for_Material[Material_Range],Index_Atom[Material_Range]);
        V_End_E = BV[0];Collision_Time_Earth = BV[2];Arrival_earth = BV[3];
        Oringal_Length_Earth = BV[4];Path_Length_Earth = BV[5];
        //Final_Length
        
        jjj = jjj + 1;
        t1->Fill();

    }

    
    
    //================Tree============//

    // fill the tree
    //================================//
       // save the Tree heade; the file will be automatically closed
       // when going out of the function scope
    char fout_name[100];
    if(BorNB==1)sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent/%sGeV/%i_STS_Bent.root",Mass_Point[Index_Mass].c_str(),Index_Sigma));
    if(BorNB==0)sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent/%sGeV/%i_STS_Bent_Comparison.root",Mass_Point[Index_Mass].c_str(),Index_Sigma));
    TFile *fout=new TFile(fout_name,"recreate");
    //Original_Bent_Comparison_Ratio_Earth->Write();
    //Original_Bent_Comparison_Ratio_Air->Write();
    t1->Write();
    fout->Close();
    cout << "fout_name: " << fout_name << endl;

}

