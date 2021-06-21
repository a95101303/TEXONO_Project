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
#include "dsigma_dT2_Bent_TEXONO.h"

//200eV threshold ==> 1.009keV
//1.009keV        ==> 201.183(km/s)
void Hist_SetLimit_Plot_v2_Possion_KS_Run_MODIFIED_BentR(int Bent_or_not, int Index_Mass, int Simulated_Event_Number, int Index_Sigma, double Sigma_SI)
{
    cout << "Sigma_SI: " << Sigma_SI << endl;
    //==============================Second step==============================
    
    //Constant
    
    string Mass_Point[5]={"20","10","2","0P2","0P05"};
    double WIMP_Mass_Array[5]={20,10,2,0.2,0.05};//12 for TEXONO

    double WIMP_Mass = WIMP_Mass_Array[Index_Mass];
    double DM_mx = WIMP_Mass_Array[Index_Mass];
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

   // double Velocity[Simulated_Event_Number];double Velocity_Z[Simulated_Event_Number];
   // double Velocity_X[Simulated_Event_Number];double Velocity_Y[Simulated_Event_Number];
    double Collision_Expectation_ATM[Simulated_Event_Number];double Collision_Expectation_EARTH[Simulated_Event_Number];
    double Collision_Expectation_Cement[Simulated_Event_Number];
    double Collision_Expectation_Reactor_Wall[Simulated_Event_Number];double Collision_Expectation_Reactor_Water[Simulated_Event_Number];
    double Collision_Expectation_Shielding[Simulated_Event_Number];

    TH1F   *Flux_HIST_Aft_Collision_EARTH = new TH1F("Flux_HIST_Aft_Collision_EARTH","Flux_HIST_Aft_Collision_EARTH",2000,0,791);
    TH1F   *Flux_HIST_Aft_Collision_Earth = new TH1F("Flux_HIST_Aft_Collision_Earth","Flux_HIST_Aft_Collision_Earth",2000,0,791);
    TH1F   *Flux_HIST_Aft_Collision_Air_I   = new TH1F("Flux_HIST_Aft_Collision_Air_I","Flux_HIST_Aft_Collision_Air_I",2000,0,791);
    TH1F   *Flux_HIST_Aft_Collision_Earth_I = new TH1F("Flux_HIST_Aft_Collision_Earth_I","Flux_HIST_Aft_Collision_Earth",2000,0,791);

    TH1F   *Energy_Loss_Percentage_Hist   = new TH1F("Energy_Loss_Percentage_Hist","Energy_Loss_Percentage_Hist",1e3,0,1);

    double Scaling_Factor=0;double Scaling_Factor_1=0;
    TH1F   *Collision_Time_Hist_Earth = new TH1F("Collision_Time_Hist_Earth","Collision_Time_Hist_Earth",300,0,300);
    TH1F   *Collision_Time_Hist_Air   = new TH1F("Collision_Time_Hist_Air","Collision_Time_Hist_Air",300,0,300);

    int Event_Number_Check=0;
    double Check_ZERO_COLLISION=0;double Check_ZERO_COLLISION_1=0;
    //cout << "Function_of_material_possibility: " << Function_of_material_possibility(2) << endl;

    int kkk=0; int jjj=0; int MMM=0;
    //for(int kkk=0 ; kkk<Simulated_Event_Number ; kkk++)
    //{
    int Bent_or_not_to_be_Bent=Bent_or_not;
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

    
    while(jjj<1)
    //while(jjj<50)
    //while((MMM<500 and Bent_or_not_to_be_Bent==1) or (jjj<500 and Bent_or_not_to_be_Bent==0 ))
    {
        
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;


        sigma_si = Sigma_SI;
        cout << "sigma_si : " << sigma_si << endl;
        mx = WIMP_Mass;
        cout << "mx : " << mx << endl;
        ev = jjj;
        Bent_or_Not = Bent_or_not_to_be_Bent;
        
        gRandom = new TRandom3(0);
        gRandom->SetSeed(0);
        double Random_Velocity = 0;
        Random_Velocity = Flux_HIST->GetRandom();
        V_Int_A = Random_Velocity;
        cout << "Random_Velocity: " << Random_Velocity << endl;
        Double_t par[3];
        TRandom *eventGenerator = new TRandom(0);//You can use TRandom(0) or TRandom3(0) to initialize your random function
        eventGenerator->GetSeed();
        eventGenerator->Sphere(par[0],par[1],par[2],Random_Velocity);
        
        
        double V_X =    (par[0]/Random_Velocity);
        double V_Y =    (par[1]/Random_Velocity);
        double V_Z =    (par[2]/Random_Velocity);
         
        
        /*
        double V_X =    0.0288531;
        double V_Y =    0.303207;
        double V_Z =    -0.952488;
         */
        double DR[3] = {V_X,V_Y,V_Z};

        double *SPF  = Starting_Position();//Starting_Point_Function
        double  SP[3] = {SPF[0],SPF[1],SPF[2]};
        
        double ROP[2] = {SPF[3],SPF[4]};//Restriction_on_Plane
        int ROC1 = SPF[3];int ROC2 = SPF[4];
        if(DR[ROC1]*ROP[1]>0)cout << "Great!" << endl;
        if(DR[ROC1]*ROP[1]<0)
        {
            cout << "Change!" << endl;
            DR[ROC1] = -DR[ROC1];
        }
        /*
        cout << "DR_1[0]: " << DR[0] << endl;
        cout << "DR_1[1]: " << DR[1] << endl;
        cout << "DR_1[2]: " << DR[2] << endl;
         */
        
        //Two times
        Flux_HIST_Random->Fill(Random_Velocity);

        double Dark_Matter_Energy = Energy_DM(DM_mx,Random_Velocity*1e3/3e8);//KeV
        double Dark_Matter_Velocity = Velocity_DM(DM_mx,Dark_Matter_Energy);//KeV

        double *Value = KS_Real_N_With_Angle(Bent_or_not_to_be_Bent, Sigma_SI, Random_Velocity, DM_mx, DR, SP); //Mx(Mass of WIMP),Velocity(km/s) Density(g/cm^3)
        V_End_E=Value[0];Collision_Time_Earth=Value[1];Arrival_earth=Value[2];
        if(Value[2]>0)Flux_HIST_Aft_Collision_EARTH->Fill(Value[0]);
        if(Value[2]>0)Collision_Time_Hist_Earth->Fill(Value[1]);
        //Final_Length

        jjj = jjj + 1;
        t1->Fill();
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;
        cout << "//Event: " << jjj << "//Air: " << kkk << endl;

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

    // fill the tree
    //================================//
       // save the Tree heade; the file will be automatically closed
       // when going out of the function scope
    char fout_name[100];
    if(Bent_or_not_to_be_Bent==1)sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent/%sGeV/%i_STS_Bent.root",Mass_Point[Index_Mass].c_str(),Index_Sigma));
    if(Bent_or_not_to_be_Bent==0)sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_TEXONO_Bent/%sGeV/%i_STS_Bent_Comparison.root",Mass_Point[Index_Mass].c_str(),Index_Sigma));
    TFile *fout=new TFile(fout_name,"recreate");
    Flux_HIST_Random->Write();
    //Flux_HIST_Aft_Collision_Earth->Write();
    Flux_HIST_Aft_Collision_EARTH->Write();
    Collision_Time_Hist_Air->Write();
    Collision_Time_Hist_Earth->Write();
    Flux_HIST_Aft_Collision_Air_I->Write();
    Flux_HIST_Aft_Collision_Earth_I->Write();
    Energy_Loss_Percentage_Hist->Write();
    //Original_Bent_Comparison_Ratio_Earth->Write();
    //Original_Bent_Comparison_Ratio_Air->Write();
    t1->Write();
    fout->Close();
    cout << "fout_name: " << fout_name << endl;

}

