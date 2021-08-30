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

//Run the program for the individual index and the simulated number of events
void Hist_SetLimit_Plot_v2_Possion_KS_Run_chie(int Index, int Simulated_Event_Number, int Index_Sigma, double Sigma_SI)
{
    cout << "Sigma_SI: " << Sigma_SI << endl;
    //Constant
    //string Mass_Point[12]={"20","19","17","15","13","11","9","7","5","4","3","2P35"};
    //double WIMP_Mass_Array[12]={20,19,17,15,13,11,9,7,5,4,3,2.35};//12 for TEXONO
    
    /*
    string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P11"};
    double WIMP_Mass_Array[10]={0.2,0.19,0.18,0.17,0.16,0.15,0.14,0.13,0.12,0.11};
     */
    /*
    string Mass_Point[19]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10","5","7"};
    double WIMP_Mass_Array[19]={2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,10,5,7};
     */
    string Mass_Point[5]={"3","4","5","10","20"};
    double WIMP_Mass_Array[5]={3,4,5,10,20};

    //Start with the DM Mass
    double WIMP_Mass = WIMP_Mass_Array[Index];
    double DM_mx = WIMP_Mass_Array[Index];

    cout << "WIMP_Mass: " << WIMP_Mass << endl;

    //================DM Halo Model================
    double Vecolity[2000];double Possiblity[2000];double Aft_Vecolity[2000];

    double sum; for(int j=0;j<2000;j++){sum = sum + velo_dist_Ave[j][3];}//Sum all posibilities
    
    for(int j=0;j<2000;j++)
    {
        float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2]);
        Vecolity[j] = v;
        Possiblity[j] = (1/(sum))*velo_dist_Ave[j][3];
    }
    //=================Flux_Hist_for_possibility=================
    TGraph *Flux = new TGraph(2000,Vecolity,Possiblity);
    TH1F   *Flux_HIST = new TH1F("Flux_HIST","Flux_HIST",2000,0,791);
    Flux_HIST->SetLineColor(1);
    Flux_HIST->SetLineWidth(5);
    for(int kkk=0;kkk<2000;kkk++){Flux_HIST->SetBinContent(kkk+1,Possiblity[kkk]);}
    
    //==========================================
    TH1F   *Flux_HIST_Random = new TH1F("Flux_HIST_Random","Flux_HIST_Random",2000,0,791);
    TH1F   *Flux_HIST_after_Collision = new TH1F("Flux_HIST_after_Collision","Flux_HIST_after_Collision",2000,0,791);

    //=================Length of all components=================
    TH1F   *HIST_Cement_Length = new TH1F("HIST_Cement_Length","HIST_Cement_Length",50,0,50);
    TH1F   *HIST_Water_Wall    = new TH1F("HIST_Water_Wall","HIST_Water_Wall",50,0,50);
    TH1F   *HIST_Water_Total   = new TH1F("HIST_Water_Total","HIST_Water_Total",50,0,50);
    TH1F   *Earth_Length       = new TH1F("Earth_Length","Earth_Length",50,0,50);


    double Velocity[Simulated_Event_Number];double Velocity_Z[Simulated_Event_Number];
    double Velocity_X[Simulated_Event_Number];double Velocity_Y[Simulated_Event_Number];
    double Collision_Expectation_ATM[Simulated_Event_Number];double Collision_Expectation_ATM_ER[Simulated_Event_Number];
    double Collision_Expectation_Cement[Simulated_Event_Number];
    double Collision_Expectation_Reactor_Wall[Simulated_Event_Number];double Collision_Expectation_Reactor_Water[Simulated_Event_Number];
    double Collision_Expectation_Shielding[Simulated_Event_Number];
    double Collision_Expectation_EARTH[Simulated_Event_Number];double Collision_Expectation_EARTH_ER[Simulated_Event_Number];
    
    for(int kkk=0 ; kkk<Simulated_Event_Number ; kkk++)
    {
        cout << "===================================" << endl;
        cout << "//Event: " << kkk;
        cout << "WIMP_Mass: " << WIMP_Mass << endl;
        
        //Generate the velocity randomly with a random direction
        gRandom = new TRandom3(0);
        gRandom->SetSeed(0);
        double Random_Velocity = 0;
        Random_Velocity = Flux_HIST->GetRandom();
        //Velocity[kkk] = Random_Velocity;
        //Velocity[kkk] = 232;
        Velocity[kkk] = 776;
        Flux_HIST_Random->Fill(Random_Velocity);//Vacuum case!
        Double_t par[3];
        TRandom *eventGenerator = new TRandom(0);//You can use TRandom(0) or TRandom3(0) to initialize your random function
        eventGenerator->GetSeed();
        eventGenerator->Sphere(par[0],par[1],par[2],Velocity[kkk]);
        
        //X, Y and Z velocity
        
        Velocity_X[kkk] = (par[0]/Velocity[kkk]);
        Velocity_Y[kkk] = (par[1]/Velocity[kkk]);
        Velocity_Z[kkk] = -abs(par[2]/Velocity[kkk]);
         /*
        cout << "Velocity_X[kkk]: " << Velocity_X[kkk] << endl;
        cout << "Velocity_Y[kkk]: " << Velocity_Y[kkk] << endl;
        cout << "Velocity_Z[kkk]: " << Velocity_Z[kkk] << endl;
          */
        /*
        Velocity_X[kkk] = 0;
        Velocity_Y[kkk] = 1;
        Velocity_Z[kkk] = 0;
         */
        //DM Energy and Velocity
        double Dark_Matter_Energy = Energy_DM(DM_mx,Velocity[kkk]*1e3/3e8);//KeV
        double Dark_Matter_Velocity = Velocity_DM(DM_mx,Dark_Matter_Energy);//KeV
        
        cout << "Velocity: " << Velocity[kkk] << endl;

        double Bool_If_Earth_Check= Bool_If_Earth( Velocity_X[kkk], Velocity_Y[kkk],  Velocity_Z[kkk]);
        Earth_Length->Fill(Bool_If_Earth_Check);
        double Cement_Length = KS_Cement_Path_Length( Velocity_X[kkk],  Velocity_Y[kkk],  Velocity_Z[kkk]);
        HIST_Cement_Length->Fill(Cement_Length);
        double Reactor_Length_Total = KS_Reactor_Path_Length(27.5, Velocity_X[kkk],  Velocity_Y[kkk],  Velocity_Z[kkk]);
        HIST_Water_Total->Fill(Reactor_Length_Total);
        double Reactor_Length_Water = KS_Reactor_Path_Length(26.5, Velocity_X[kkk],  Velocity_Y[kkk],  Velocity_Z[kkk]);
        HIST_Water_Wall->Fill(Reactor_Length_Water);
        //cout << "=============Reactor_Length_Water=============: " << Reactor_Length_Water << endl;
double Path_Length_For_Three_Components[4]={Bool_If_Earth_Check,Cement_Length,Reactor_Length_Total,Reactor_Length_Water};//1 for Air check, 2 for the Cement, 3 for the Reactor
        double *A=  KS_Collision_Time_EARTH(Sigma_SI, Velocity_Y[kkk], Velocity_Z[kkk],Velocity[kkk],DM_mx,Path_Length_For_Three_Components);
        Collision_Expectation_EARTH[kkk]=A[1];Collision_Expectation_Cement[kkk]=A[2];
        Collision_Expectation_Reactor_Wall[kkk]=A[3];Collision_Expectation_Reactor_Water[kkk]=A[4];Collision_Expectation_Shielding[kkk]=A[5];Collision_Expectation_EARTH_ER[kkk]=A[6];
        double *B = KS_Collision_Time_ATM(Sigma_SI, Velocity_Y[kkk], Velocity_Z[kkk],Velocity[kkk],DM_mx,Path_Length_For_Three_Components,A[0]);
        Collision_Expectation_ATM[kkk] = B[0];Collision_Expectation_ATM_ER[kkk] = B[1];
    }
    //==============================Third Step==============================
    
    TH1F   *Flux_HIST_Aft_Collision_EARTH = new TH1F("Flux_HIST_Aft_Collision_EARTH","Flux_HIST_Aft_Collision_EARTH",2000,0,791);
    TH1F   *Flux_HIST_Aft_Collision_Earth = new TH1F("Flux_HIST_Aft_Collision_Earth","Flux_HIST_Aft_Collision_Earth",2000,0,791);

    double Scaling_Factor=0;double Scaling_Factor_1=0;
    TH1F   *Collision_Time_Hist_NU = new TH1F("Collision_Time_Hist_NU","Collision_Time_Hist_NU",2000,0,10000);
    TH1F   *Collision_Time_Hist_ER = new TH1F("Collision_Time_Hist_ER","Collision_Time_Hist_ER",2000,0,10000);

    int Event_Number_Check=0;
    TH1F   *Collision_Energy_Lost_E = new TH1F("Collision_Energy_Lost_E","Collision_Energy_Lost_E",100,0,20);
    TH1F   *Collision_Energy_Lost_V = new TH1F("Collision_Energy_Lost_V","Collision_Energy_Lost_V",100,0,20);
    double Check_ZERO_COLLISION=0;double Check_ZERO_COLLISION_1=0;
    //cout << "Function_of_material_possibility: " << Function_of_material_possibility(2) << endl;
    
    /*
    int Earth_Threshold=0;
    int Air_Threshold=0;
    for(int EEE=0; EEE<5; EEE++)//i is earth and 2 is air
    {
        int Check  = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI,779.135,1)[3];
        if(Earth_Threshold<Check)Earth_Threshold=Check;
        int Check1 = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI,779.135,2)[3];
        if(Air_Threshold<Check1)Air_Threshold=Check1;
    }
     
    cout << "Earth_Threshold: " << Earth_Threshold << endl;
    cout << "Air_Threshold: "   << Air_Threshold << endl;
*/

    for(int kkk=0; kkk<Simulated_Event_Number; kkk++)
    {
        
        cout << "===================================" << endl;
        cout << "//Event2: " << kkk << endl;
        cout << "Energy_DM: " << Energy_DM(DM_mx,Velocity[kkk]*1e3/3e8) << endl;
        cout << "***==================================***" << endl;
        cout << "Sigma_SI_Default: " << Sigma_SI << endl;
        cout << "Velocity_X[kkk]: " << Velocity_X[kkk] << endl;
        cout << "Velocity_Y[kkk]: " << Velocity_Y[kkk] << endl;
        cout << "Velocity_Z[kkk]: " << Velocity_Z[kkk] << endl;

        double ECNR  = Collision_Expectation_EARTH[kkk];//EARTH_COLLISION_COUNT_NUCLEAR_RECOIL
        //double ECER  = Collision_Expectation_EARTH_ER[kkk];//EARTH_COLLISION_COUNT_ELECTRON_RECOIL
        double ECER  = ECNR*(1./4.)*(0.5e-3)*(1./4.)*(0.5e-3);//EARTH_COLLISION_COUNT_ELECTRON_RECOIL

        Collision_Time_Hist_NU->Fill(ECNR);
        Collision_Time_Hist_ER->Fill(ECER);

        cout << "ECNR: " << ECNR << endl;
        cout << "ECER: " << ECER << endl;
        
        //double *A = Velocity_Aft_collision_ER(int(ECER),DM_mx,Sigma_SI,Velocity[kkk]);
        
        cout << "LOSS(keV)" << ECER*5/1e3 << endl;
    }
    
    TLegend *leg= new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    
    Flux_HIST_Random->SetLineWidth(1);
    Flux_HIST_Random->SetLineColor(1);
    
    Flux_HIST_Aft_Collision_Earth->SetLineWidth(1);
    Flux_HIST_Aft_Collision_Earth->SetLineColor(2);
    
    Flux_HIST_Aft_Collision_EARTH->SetLineWidth(1);
    Flux_HIST_Aft_Collision_EARTH->SetLineColor(3);
    
    
    
    //================Tree============//
    TTree *t1 = new TTree("t1","Information");
    double sigma_si,mx;int ev;
    t1->Branch("sigma_si",&sigma_si,"sigma_si/D");
    t1->Branch("mx",&mx,"mx/D");
    t1->Branch("ev",&ev,"ev/I");

    // fill the tree
    for (Int_t i=0; i<Simulated_Event_Number; i++) {
        sigma_si = Sigma_SI;
        cout << "sigma_si : " << sigma_si << endl;
        mx = WIMP_Mass;
        cout << "mx : " << mx << endl;
        ev = i;
        cout << " i: " << i << endl;
        t1->Fill();
    }
    //================================//
       // save the Tree heade; the file will be automatically closed
       // when going out of the function scope
    char fout_name[100];
    sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/TEXONO_Flux_CAT/%sGeV/%i.root",Mass_Point[Index].c_str(),Index_Sigma));
    TFile *fout=new TFile(fout_name,"recreate");
    Flux_HIST_Random->Write();
    Flux_HIST_Aft_Collision_Earth->Write();
    Flux_HIST_Aft_Collision_EARTH->Write();
    Collision_Time_Hist_NU->Write();
    Collision_Time_Hist_ER->Write();

    HIST_Cement_Length->Write();
    HIST_Water_Wall->Write();
    HIST_Water_Total->Write();
    Earth_Length->Write();
    
    t1->Write();
    fout->Close();
    cout << "fout_name: " << fout_name << endl;
}

