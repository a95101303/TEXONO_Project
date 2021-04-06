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
void Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_Final(int Index, int Simulated_Event_Number, int Index_Sigma, double Sigma_SI)
{
    cout << "Sigma_SI: " << Sigma_SI << endl;
    //==============================Second step==============================
    
    //Constant
    string Mass_Point[3]={"20","2","0P2"};
    double WIMP_Mass_Array[12]={20,2,0.2};//3 for TEXONO

    //string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P11"};
    //string Mass_Point[19]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10","5","7"};
    //double WIMP_Mass_Array[19]={2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,10,5,7};
    //double WIMP_Mass_Array[10]={0.2,0.19,0.18,0.17,0.16,0.15,0.14,0.13,0.12,0.11};
    //Start
    double WIMP_Mass = WIMP_Mass_Array[Index];
    double DM_mx = WIMP_Mass_Array[Index];
    //Double_t WIMP_Mass_Array[13]={2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08};
    cout << "WIMP_Mass: " << WIMP_Mass << endl;

    double Vecolity[2000];double Possiblity[2000];
    double Aft_Vecolity[2000];

    double sum; for(int j=0;j<2000;j++){sum = sum + velo_dist_Ave[j][3];}
    for(int j=0;j<2000;j++)
    {
        float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2]);
        Vecolity[j] = v;
        Possiblity[j] = (1/(sum))*velo_dist_Ave[j][3];
    }
    TGraph *Flux = new TGraph(2000,Vecolity,Possiblity);
    TH1F   *Flux_HIST = new TH1F("Flux_HIST","Flux_HIST",2000,0,791);
    Flux_HIST->SetLineColor(1);
    Flux_HIST->SetLineWidth(5);

    TH1F   *Flux_HIST_Random = new TH1F("Flux_HIST_Random","Flux_HIST_Random",2000,0,791);
    TH1F   *Flux_HIST_Random_Normalized = new TH1F("Flux_HIST_Random_Normalized","Flux_HIST_Random_Normalized",2000,0,791);
    TH1F   *Flux_HIST_after_Collision = new TH1F("Flux_HIST_after_Collision","Flux_HIST_after_Collision",2000,0,791);
    TH2F   *E_vs_Collision_Time = new TH2F("E_vs_Collision_Time","E_vs_Collision_Time",10000,0,100,10000,0,100);
    
    Flux_HIST_Random_Normalized->SetLineColor(2);

    for(int kkk=0;kkk<2000;kkk++){Flux_HIST->SetBinContent(kkk+1,Possiblity[kkk]);}

    double Velocity[Simulated_Event_Number];double Velocity_Z[Simulated_Event_Number];
    double Velocity_X[Simulated_Event_Number];double Velocity_Y[Simulated_Event_Number];
    double Collision_Expectation_ATM[Simulated_Event_Number];double Collision_Expectation_EARTH[Simulated_Event_Number];
    double Collision_Expectation_Cement[Simulated_Event_Number];
    double Collision_Expectation_Reactor_Wall[Simulated_Event_Number];double Collision_Expectation_Reactor_Water[Simulated_Event_Number];
    double Collision_Expectation_Shielding[Simulated_Event_Number];
    double Collision_Path_Distended[Simulated_Event_Number];
    
    TFile *Bent_File = TFile::Open("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/2_TEXONO_Flux/20GeV/1_STS_Bent.root");
    TH1F *T1_Ratio = (TH1F*)Bent_File->Get("Extended_Length_AIR_Percentage");
         
    for(int kkk=0 ; kkk<Simulated_Event_Number ; kkk++)
    {
        cout << "===================================" << endl;
        cout << "//Event: " << kkk;
        cout << "WIMP_Mass: " << WIMP_Mass << endl;
        gRandom = new TRandom3(0);
        gRandom->SetSeed(0);
        double Random_Velocity = 0;
        Random_Velocity = Flux_HIST->GetRandom();
        
        Collision_Path_Distended[kkk] = T1_Ratio->GetRandom();
        Velocity[kkk] = Random_Velocity;
        Flux_HIST_Random->Fill(Random_Velocity);
        Double_t par[3];
        TRandom *eventGenerator = new TRandom(0);//You can use TRandom(0) or TRandom3(0) to initialize your random function
        eventGenerator->GetSeed();
        eventGenerator->Sphere(par[0],par[1],par[2],Velocity[kkk]);
        
        
        Velocity_X[kkk] = (par[0]/Velocity[kkk]);
        Velocity_Y[kkk] = (par[1]/Velocity[kkk]);
        Velocity_Z[kkk] = (par[2]/Velocity[kkk]);
         
        /*
        Velocity_X[kkk] = 0;
        Velocity_Y[kkk] = -1;
        Velocity_Z[kkk] = sqrt(1-Velocity_Y[kkk]*Velocity_Y[kkk]-Velocity_X[kkk]*Velocity_X[kkk]);
         */

        double Dark_Matter_Energy = Energy_DM(DM_mx,Velocity[kkk]*1e3/3e8);//KeV
        double Dark_Matter_Velocity = Velocity_DM(DM_mx,Dark_Matter_Energy);//KeV
        
        cout << "Velocity: " << Velocity[kkk] << endl;

        double Bool_If_Earth_Check= Bool_If_Earth( Velocity_X[kkk], Velocity_Y[kkk],  Velocity_Z[kkk]);
        double Cement_Length = KS_Cement_Path_Length( Velocity_X[kkk],  Velocity_Y[kkk],  Velocity_Z[kkk]);
        double Reactor_Length_Total = KS_Reactor_Path_Length(27.5, Velocity_X[kkk],  Velocity_Y[kkk],  Velocity_Z[kkk]);
        double Reactor_Length_Water = KS_Reactor_Path_Length(26.5, Velocity_X[kkk],  Velocity_Y[kkk],  Velocity_Z[kkk]);
        //cout << "=============Reactor_Length_Water=============: " << Reactor_Length_Water << endl;
double Path_Length_For_Three_Components[4]={Bool_If_Earth_Check,Cement_Length,Reactor_Length_Total,Reactor_Length_Water};//1 for Air check, 2 for the Cement, 3 for the Reactor
        double *A=KS_Collision_Time_EARTH(Sigma_SI, Velocity_Y[kkk], Velocity_Z[kkk],Velocity[kkk],DM_mx,Path_Length_For_Three_Components);
        Collision_Expectation_EARTH[kkk]=A[1];Collision_Expectation_Cement[kkk]=A[2];
        Collision_Expectation_Reactor_Wall[kkk]=A[3];Collision_Expectation_Reactor_Water[kkk]=A[4];Collision_Expectation_Shielding[kkk]=A[5];
        Collision_Expectation_ATM[kkk] = KS_Collision_Time_ATM(Sigma_SI, Velocity_Y[kkk], Velocity_Z[kkk],Velocity[kkk],DM_mx,Path_Length_For_Three_Components,A[0]);
        
    }
    //==============================Third Step==============================
    
    Bent_File->Close();

    TH1F   *Flux_HIST_Aft_Collision_EARTH = new TH1F("Flux_HIST_Aft_Collision_EARTH","Flux_HIST_Aft_Collision_EARTH",2000,0,791);
    TH1F   *Flux_HIST_Aft_Collision_Earth = new TH1F("Flux_HIST_Aft_Collision_Earth","Flux_HIST_Aft_Collision_Earth",2000,0,791);

    double Scaling_Factor=0;double Scaling_Factor_1=0;
    TH1F   *Collision_Time_Hist = new TH1F("Flux_HIST_Aft_Collision","Flux_HIST_Aft_Collision",10,0,10);
    int Event_Number_Check=0;
    TH1F   *Collision_Energy_Lost_E = new TH1F("Collision_Energy_Lost_E","Collision_Energy_Lost_E",100,0,20);
    TH1F   *Collision_Energy_Lost_V = new TH1F("Collision_Energy_Lost_V","Collision_Energy_Lost_V",100,0,20);
    double Check_ZERO_COLLISION=0;double Check_ZERO_COLLISION_1=0;
    //cout << "Function_of_material_possibility: " << Function_of_material_possibility(2) << endl;
    
    int Earth_Threshold=0;
    int Air_Threshold=0;
    
    TH1F  *HIST_Collision_Expectation_ATM          = new TH1F("HIST_Collision_Expectation_ATM","HIST_Collision_Expectation_ATM",2000,0,2000);
    TH1F  *HIST_Collision_Expectation_ATM_MODIFIED = new TH1F("HIST_Collision_Expectation_ATM_MODIFIED","HIST_Collision_Expectation_ATM_MODIFIED",2000,0,2000);

    for(int EEE=0; EEE<5; EEE++)//i is earth and 2 is air
    {
        int Check  = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI,799.135,1)[3];
        if(Earth_Threshold<Check)Earth_Threshold=Check;
        int Check1 = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI,799.135,2)[3];
        if(Air_Threshold<Check1)Air_Threshold=Check1;
    }
    cout << "Earth_Threshold: " << Earth_Threshold << endl;
    cout << "Air_Threshold: "   << Air_Threshold << endl;


    for(int kkk=0; kkk<Simulated_Event_Number; kkk++)
    {
        
        cout << "===================================" << endl;
        cout << "//Event2: " << kkk << endl;
        
        //cout << "Loop2==>EventKIND: " << jjj << endl;
        //cout << "Loop2=>Event: " << kkk << endl;
        double EC  = Collision_Expectation_EARTH[kkk];//EARTH_COLLISION_COUNT
        double CC  = Collision_Expectation_Cement[kkk];//CEMENT_COLLISION_COUNT
        double RCL = Collision_Expectation_Reactor_Wall[kkk];//REACTOR_COLLISION_WALL_COUNT
        double RCR = Collision_Expectation_Reactor_Water[kkk];//REACTOR_COLLISION_WATER_COUNT
        double AC  = Collision_Expectation_ATM[kkk];//AIR_COLLISION_COUNT
        double ACM = Collision_Expectation_ATM[kkk]*(1+Collision_Path_Distended[kkk]);//AIR_COLLISION_COUNT_MODIFIED
        double CS = Collision_Expectation_Shielding[kkk];
        /*
        cout << "***==================================***" << endl;
        cout << "Sigma_SI_Default: " << Sigma_SI << endl;
        cout << "Velocity_X[kkk]: " << Velocity_X[kkk] << endl;
        cout << "Velocity_Y[kkk]: " << Velocity_Y[kkk] << endl;
        cout << "Velocity_Z[kkk]: " << Velocity_Z[kkk] << endl;

        cout << "Collision_Expectation_EARTH[kkk]: " << Collision_Expectation_EARTH[kkk] << endl;
        cout << "Collision_Expectation_Cement[kkk]: " << Collision_Expectation_Cement[kkk] << endl;
        cout << "Collision_Expectation_Reactor_Wall[kkk]: " << Collision_Expectation_Reactor_Wall[kkk] << endl;
        cout << "Collision_Expectation_Reactor_Water[kkk]: " << Collision_Expectation_Reactor_Water[kkk] << endl;
        cout << "Collision_Expectation_Shielding[kkk]: " << Collision_Expectation_Shielding[kkk] << endl;
         */
        cout << "***==================================***" << endl;
        cout << "Collision_Expectation_ATM[kkk]: " << AC << endl;
        cout << "Collision_Path_Distended[kkk]: " << (1+Collision_Path_Distended[kkk]) << endl;
        cout << "Collision_Expectation_ATM_MODIFIED[kkk]: " << ACM << endl;
        cout << "***==================================***" << endl;
         
        HIST_Collision_Expectation_ATM->Fill(AC);
        HIST_Collision_Expectation_ATM_MODIFIED->Fill(ACM);
        
        if(EC>Earth_Threshold or CC>1e4 or RCL>1e4 or RCR>1e4 or AC>Air_Threshold or CS>1e4){
            cout << "OK GREAT!" << endl;
            Flux_HIST_Aft_Collision_EARTH->Fill(1e-5);}

                //====================Test_Function====================
        //Two-stage collision
        else
        {//1 is earth 2 is air
            cout << "V_Initial: " << Velocity[kkk] << endl;
            
            double *V_Aft_Collision_AIR = Velocity_Aft_collision(ACM,DM_mx,Sigma_SI,Velocity[kkk],2);
            cout << "V_Aft_AIR: " << V_Aft_Collision_AIR[0] << endl;

            double *V_Aft_Collision_EARTH = Velocity_Aft_collision(EC,DM_mx,Sigma_SI,V_Aft_Collision_AIR[0],1);
            //cout << "V_Aft_EARTH: " << V_Aft_Collision_EARTH[0] << endl;
           
            double *V_Aft_Collision_CEMENT = Velocity_Aft_collision(CC,DM_mx,Sigma_SI,V_Aft_Collision_EARTH[0],3);
            //cout << "V_Aft_CEMENT: " << V_Aft_Collision_CEMENT[0] << endl;

            double *V_Aft_Collision_REACTOR_WALL = Velocity_Aft_collision(RCL,DM_mx,Sigma_SI,V_Aft_Collision_CEMENT[0],3);
            //cout << "V_Aft_REACTOR_WALL : " << V_Aft_Collision_REACTOR_WALL[0] << endl;

            double *V_Aft_Collision_REACTOR_WATER = Velocity_Aft_collision(RCR,DM_mx,Sigma_SI,V_Aft_Collision_REACTOR_WALL[0],5);
            //cout << "V_Aft_Collision_REACTOR_WATER : " << V_Aft_Collision_REACTOR_WATER[0] << endl;
            
            double *V_Aft_Collision_SHIELDING = Velocity_Aft_collision(CS,DM_mx,Sigma_SI,V_Aft_Collision_REACTOR_WATER[0],6);
            cout << "V_Aft_Collision_SHIELDING : " << V_Aft_Collision_SHIELDING[0] << endl;

            Flux_HIST_Aft_Collision_EARTH->Fill(V_Aft_Collision_SHIELDING[0]);
        }
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
    for (Int_t i=0; i<2; i++) {
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
    sprintf(fout_name,Form("2_TEXONO_Flux/%sGeV/%i_Bent_Final.root",Mass_Point[Index].c_str(),Index_Sigma));
    TFile *fout=new TFile(fout_name,"recreate");
    cout << "1: " << endl;
    Flux_HIST_Random->Write();
    cout << "2: " << endl;
    Flux_HIST_Aft_Collision_Earth->Write();
    cout << "3: " << endl;
    Flux_HIST_Aft_Collision_EARTH->Write();
    cout << "4: " << endl;
    HIST_Collision_Expectation_ATM->Write();
    cout << "5: " << endl;
    HIST_Collision_Expectation_ATM_MODIFIED->Write();
    cout << "6: " << endl;
    t1->Write();
    cout << "7: " << endl;
    fout->Close();
    cout << "8: " << endl;
    cout << "fout_name: " << fout_name << endl;
}

