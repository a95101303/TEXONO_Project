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

//200eV threshold ==> 1.009keV
//1.009keV        ==> 201.183(km/s)
void Hist_SetLimit_Plot_v2_Possion_Boundary_Check_SURF()
{
    /*
    const int Event_Number=16;
    double Threshold_of_Detector=0.08; //for CRESST surface
    double WIMP_Mass_Array[Event_Number]={0.5,0.7,0.9,1,2,3,4,5,7,9,11,13,15,17,19,20};//16 for CRESST
     */
    //double WIMP_Mass_Array[Event_Number]={0.5,0.7,0.9,1};//4 for CRESST

    //double WIMP_Mass_Array[Event_Number]={1};//23 for CRESST
    
    const int Event_Number=1;
    double Threshold_of_Detector = 1;//Without attenuating
    //double WIMP_Mass_Array[Event_Number]={2.34,2,3,4,5,7,9,11,13,15,17,19,20};//19 for TEXONO
     double WIMP_Mass_Array[Event_Number]={1e11};
   // double Critical_Point[11] = {779.135000,788.344,777.158,776.366,773.944,768.457,752.241,709.923,621.726,468.272,201.183};
    double Critical_Point[12] = {779.135000,750,700,650,600,550,500,450,400,350,278.036};
    double Sigma_SI_Array[Event_Number];
    //Check
    double Threshold_Criteria_Test=(1000.0*max_recoil_A(2.34,779.135000*1000.0/2.99792458e8, AGe));
    cout << "Threshold_Criteria_Test: " << Threshold_Criteria_Test << endl;
    double Threshold_Criteria_Test_2=(1000.0*max_recoil_A(0.14,779.135000*1000.0/2.99792458e8, 16));//O:0.0174361 Al:0.0104113
    cout << "Threshold_Criteria_Test_2: " << Threshold_Criteria_Test_2 << endl;

for(int WIMP_Event=0; WIMP_Event<Event_Number; WIMP_Event++)
{
   //if(WIMP_Event>=0) continue;
    double WIMP_Mass=WIMP_Mass_Array[WIMP_Event];
    for(int lll=0; lll<7; lll++)
    {
        for(int zzz=10; zzz<11; zzz++)//Z<10
        {
            for(int ppp=0; ppp<1; ppp++ )
            {
                int Count=0;
                cout << "Sigma_SI_Array[WIMP_Event]: " << Sigma_SI_Array[WIMP_Event] << endl;
                if(Sigma_SI_Array[WIMP_Event]>1e-45) continue;
                double Sigma_SI_Default=  ((ppp)+(0.1*zzz)) * (1e-27) * TMath::Power(10,(7-lll));
                cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
                
                int Simulated_Event_Number=1;
                double Velocity[Simulated_Event_Number];double Velocity_Z[Simulated_Event_Number];
                double Collision_Expectation_ATM[Simulated_Event_Number];double Collision_Expectation_EARTH[Simulated_Event_Number];
                double Collision_Expectation_Earth[Simulated_Event_Number];
                cout << "3: " << endl;

                for(int kkk=0 ; kkk<Simulated_Event_Number ; kkk++)
                {
                    cout << "//Event: " << kkk;
                    gRandom = new TRandom3(0);
                    gRandom->SetSeed(0);
                    Velocity[kkk]=778.937;//Max
                    Velocity_Z[kkk]=1;//Minimum
                    double Dark_Matter_Energy = Energy_DM(WIMP_Mass,Velocity[kkk]*1e3/3e8);//KeV
                    double Dark_Matter_Velocity = Velocity_DM(WIMP_Mass,Dark_Matter_Energy);//KeV
                    cout << "Velocity: " << Velocity[kkk] << endl;
                    //cout << "Velocity_Z: " << Velocity_Z[kkk] << endl;
                    cout << "Dark_Matter_Energy: " << Dark_Matter_Energy << endl;
                    cout << "Dark_Matter_Velocity: " << Dark_Matter_Velocity << endl;
                    double Collision_Expectation_velocity_dependent_Earth = 0;
                    double Collision_Expectation_velocity_dependent_EARTH = Collision_Time_EARTH(Sigma_SI_Default,Velocity_Z[kkk],Velocity[kkk],WIMP_Mass);
                    double Collision_Expectation_velocity_dependent_ATM   = Collision_Time_ATM(Sigma_SI_Default, Velocity_Z[kkk],Velocity[kkk],WIMP_Mass);
                    //cout << "COLLISION_TIME_EARTH_Exact: " << Collision_Expectation_velocity_dependent_EARTH << endl;
                    //cout << "COLLISION_TIME_Earth_Rough: " << Collision_Expectation_velocity_dependent_Earth << endl;
                    //cout << "COLLISION_TIME_AIR: " << Collision_Expectation_velocity_dependent_ATM << endl;
                    Collision_Expectation_ATM[kkk] = Collision_Expectation_velocity_dependent_ATM;
                    Collision_Expectation_EARTH[kkk] = Collision_Expectation_velocity_dependent_EARTH;
                    Collision_Expectation_Earth[kkk] = Collision_Expectation_velocity_dependent_Earth;

                }
                cout << "4: " << endl;

                //==============================Third Step==============================
                /*
                TH1F   *Flux_HIST_Aft_Collision_EARTH = new TH1F("Flux_HIST_Aft_Collision_EARTH","Flux_HIST_Aft_Collision_EARTH",2000,0,791);
                TH1F   *Flux_HIST_Aft_Collision_Earth = new TH1F("Flux_HIST_Aft_Collision_Earth","Flux_HIST_Aft_Collision_Earth",2000,0,791);

                double Scaling_Factor=0;double Scaling_Factor_1=0;
                TH1F   *Collision_Time_Hist = new TH1F("Flux_HIST_Aft_Collision","Flux_HIST_Aft_Collision",10,0,10);
                int Event_Number_Check=0;
                TH1F   *Collision_Energy_Lost_E = new TH1F("Collision_Energy_Lost_E","Collision_Energy_Lost_E",100,0,20);
                TH1F   *Collision_Energy_Lost_V = new TH1F("Collision_Energy_Lost_V","Collision_Energy_Lost_V",100,0,20);
                double Check_ZERO_COLLISION=0;double Check_ZERO_COLLISION_1=0;
                 */
                //cout << "Function_of_material_possibility: " << Function_of_material_possibility(2) << endl;
                /*
                double par[3]={1,2,3};double par1[3]={1,2,4};
                double Length = Two_points_difference(par,par1);
                cout << "Length: " << Length << endl;
                double X=0.2; double Y=0.7; double Z=sqrt(1-X*X-Y*Y);
                cout << "X: " << X << endl;cout << "Y: " << Y << endl;cout << "Z: " << Z << endl;
                cout << "Cement_Path_Length:(km) " << Cement_Path_Length(X,Y,Z) << endl;
                cout << "Reactor_Path_Length(): " << Reactor_Path_Length(X,Y,Z) << endl;*/
                
                cout << "5: " << endl;

            for(int jjj=0; jjj<1 ; jjj++)
            {
                double *Collision_Expectation_EARTH_Temp;
                if(jjj==0) Collision_Expectation_EARTH_Temp=Collision_Expectation_EARTH;
                if(jjj==1) Collision_Expectation_EARTH_Temp=Collision_Expectation_Earth;

                for(int kkk=0; kkk<Simulated_Event_Number; kkk++)
                {
                    if(Count!=0)
                    {   cout <<"CONTINUE: " << endl;
                        continue;}
                    cout << "//Event2: " << kkk << endl;;
                    //cout << "Loop2==>EventKIND: " << jjj << endl;
                    //cout << "Loop2=>Event: " << kkk << endl;
                    cout << "Collision_Expectation_ATM[kkk]: " << Collision_Expectation_ATM[kkk] << endl;
                    cout << "Collision_Expectation_EARTH_Temp[kkk]: " << Collision_Expectation_EARTH_Temp[kkk] << endl;

                            //====================Test_Function====================
                    //Two-stage collision
                        cout << "Collision_Expectation_ATM[kkk]: " << Collision_Expectation_ATM[kkk] << endl;
                        cout << "Collision_Expectation_EARTH_Temp[kkk]: " << Collision_Expectation_EARTH_Temp[kkk] << endl;
                        
                        double *V_Aft_Collision_AIR = Velocity_Aft_collision(Collision_Expectation_ATM[kkk],WIMP_Mass,Sigma_SI_Default,Velocity[kkk],2);
                        cout << "V_Aft_AIR: " << V_Aft_Collision_AIR[0] << endl;
                        //cout << "V_Aft_AIR_Energy_Lost_E: " << V_Aft_Collision_AIR[1] << endl;
                        //cout << "V_Aft_AIR_Energy_Lost_V: " << V_Aft_Collision_AIR[2] << endl;

                        double *V_Aft_Collision_EARTH = Velocity_Aft_collision(Collision_Expectation_EARTH_Temp[kkk],WIMP_Mass,Sigma_SI_Default,V_Aft_Collision_AIR[0],3);
                        cout << "V_Aft_EARTH: " << V_Aft_Collision_EARTH[0] << endl;
                        //cout << "V_Aft_EARTH_Energy_Lost_E: " << V_Aft_Collision_EARTH[1] << endl;
                        //cout << "V_Aft_EARTH_Energy_Lost_V: " << V_Aft_Collision_EARTH[2] << endl;
                        double Threshold_Criteria=Energy_DM(WIMP_Mass,V_Aft_Collision_EARTH[0]*1e3/3e8);
                        cout << "Threshold_Criteria: " << Threshold_Criteria << endl;
                  if( Threshold_Criteria > Threshold_of_Detector)
                  {
                      Count = Count +1;
                      cout << "Count: " << Count << endl;
                  }
                 else
                 {
                     continue;
                 }
                    }//Close kkk
                }//Close jjj
                cout << "Count: " << Count << endl;
                cout << "Count/20: " << Count/20 << endl;
                if(Count!=0) {continue;}
                else{
                    cout << "???" << endl;
                    Sigma_SI_Array[WIMP_Event]=Sigma_SI_Default;}
            }//Close ppp
        }//Close zzz
    }// Close lll
}//Close the WIMP
    
    
    for(int WIMP_Event=0; WIMP_Event<Event_Number; WIMP_Event++)
    {
        cout << WIMP_Mass_Array[WIMP_Event] << "," << Sigma_SI_Array[WIMP_Event] << "," << endl;
    }
    

    //E_vs_Collision_Time->Draw("colz");
    //c1->Print("Yes.pdf");
    
    /*
    Collision_Energy_Lost_E->SetLineWidth(1);
    Collision_Energy_Lost_V->SetLineColor(1);
    
    Collision_Energy_Lost_E->SetLineWidth(3);
    Collision_Energy_Lost_V->SetLineColor(2);

    Collision_Energy_Lost_E->Draw();
    Collision_Energy_Lost_V->Draw("histsame");
    c1->Print("34.pdf");*/

/*
    Flux_HIST_Random->SetLineWidth(2);
    Flux_HIST_Random->SetLineColor(2);
    
    Flux_HIST_after_Collision->SetLineWidth(1);
    Flux_HIST_after_Collision->SetLineColor(1);
    
    Flux_HIST_Random->Draw();
    Flux_HIST_after_Collision->Draw("same");
    */
    
    //Collision_Time_Hist->DrawNormalized();
    //Collision_Energy_Lost->Draw();
    //Flux_HIST_Random_Normalized->Draw("same");
    //Collision_Time_Hist->Draw();
    
    //==========================================//
    //==========================================//
}
//cout << "Prove: " << TMath::ACos(Earth_Radius/(Earth_Radius+1)) << endl;
//cout << "Prove1: " << Earth_Radius*TMath::Sin(TMath::ACos(Earth_Radius/(Earth_Radius+1))) << endl;

// (1) Find out the distribution with the random function
// ==> Match between the velocity distribution
// (2) Find out the possion with the sigma(v)*mass density ==> Turn out to be the mean free path
// ==> Earth_radius*2/mean free path to be the expectation value of the possion
// (3) Ask if YES or NO
// ==> If Yes, Find out dsigma/dT distribution ==> Use the SAME function (GetRandom) to get the T we need to substract
// (4) Do all the things again for all 10000 events
//After_Colliding_with_the_possibility
/*
double DM_Energy_Aft_Colliding=0;double DM_Velocity_Aft_Colliding=0;
TF1 *f3 = new TF1("f3","fdsigma_dT_keV([0],[1],[2],[3],x)",0,10);
f3->SetParameter(0,10);f3->SetParameter(1,Sigma_SI_Default);f3->SetParameter(2,(Velocity[kkk]*1e3/3e8));f3->SetParameter(3,AFe);
double Random_Energy= f3->GetRandom();
DM_Energy_Aft_Colliding = (Dark_Matter_Energy - Random_Energy);
DM_Velocity_Aft_Colliding = Velocity_DM(10,DM_Energy_Aft_Colliding);
cout << "Subtracted_Energy: " << Random_Energy << endl;
cout << "DM_Energy_Aft_Colliding: " << DM_Energy_Aft_Colliding << endl;
cout << "DM_Velocity_Aft_Colliding: " << DM_Velocity_Aft_Colliding << endl;
cout << "===============================" << endl;


for(int kkk=0; kkk<1000 ; kkk++)
{
double Collision_Expectionl=Collision_Possibility->GetRandom();
cout << "Collision_Expectionl: " << Collision_Expectionl << endl;
}
cout << "total_Sigma(10,A0): " << total_Sigma(10,AO) << endl;
cout << "total_Sigma(10,AFe): " << total_Sigma(10,AFe) << endl;

cout << "Mean_Free_Path(O):" << Mean_Free_Path(total_Sigma(10,AO)*Mass_Density_of_Earth(atom_mass_O_g))/1e5 << endl;
cout << "Mean_Free_Path(Fe):" << Mean_Free_Path(total_Sigma(10,AFe)*Mass_Density_of_Earth(atom_mass_Fe_g))/1e5 << endl;
cout << "Mean_free_Path(O+Fe+Si): " << 1/(1/(Mean_Free_Path(total_Sigma(10,AO)*Mass_Density_of_Earth(atom_mass_O_g))/1e5)+1/(Mean_Free_Path(total_Sigma(10,AFe)*Mass_Density_of_Earth(atom_mass_Fe_g))/1e5) ) << endl;
cout << "Mean_Free_Path(Fe)[km]:" << Mean_Free_Path(total_Sigma*Mass_Density_of_Earth(atom_mass_Fe_g))/1e5;
cout << "Mass_Density_of_Earth_O: " << Mass_Density_of_Earth(atom_mass_O_g) << endl;
cout << "Mass_Density_of_Earth: " << Mass_Density_of_Earth(atom_mass_Si_g) << endl;
cout << "total_Sigma: " << total_Sigma << endl;
*/

