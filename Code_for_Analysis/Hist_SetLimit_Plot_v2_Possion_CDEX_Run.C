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
void Hist_SetLimit_Plot_v2_Possion_CDEX_Run(int Index, int Simulated_Event_Number, int Index_Sigma, double Sigma_SI)
{
    //Constant
    string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    double WIMP_Mass_Array[17]={2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,10};
    //double Minimum_Velocity[17]={30.6725,42.9415,45.3161,48.9865,51.2527,55.2104,60.3555,67.4794,77.7696,95.1836,134.761,141.885,150.592,161.674,174.734,190.169,0}
    //Start
    double WIMP_Mass = WIMP_Mass_Array[Index];
    //Double_t WIMP_Mass_Array[13]={2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08};
    
    int Count=0;

        double Sigma_SI_Default= Sigma_SI;
        cout << "Sigma_SI_Default2: " << Sigma_SI_Default << endl;
        
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
        cout << "2: " << endl;

        for(int kkk=0;kkk<2000;kkk++){Flux_HIST->SetBinContent(kkk+1,Possiblity[kkk]);}

        
        double Velocity[Simulated_Event_Number];double Velocity_Z[Simulated_Event_Number];
        double VX[Simulated_Event_Number];double VY[Simulated_Event_Number];double VZ[Simulated_Event_Number];
        double Collision_Expectation_ATM[Simulated_Event_Number];double Collision_Expectation_EARTH[Simulated_Event_Number];
        double Collision_Expectation_Earth[Simulated_Event_Number];
        /*
        double TRYTRY1 = Collision_Time_Earth(Path_Length_Earth_Rough(-0.8)*1e5,Mean_Free_Path_of_Earth(1,500,Sigma_SI_Default));
        double TRYTRY2 = Collision_Time_EARTH(Sigma_SI_Default,-0.8,500);
        cout << "TRYTRY1: " << TRYTRY1 << endl;cout << "TRYTRY2: " << TRYTRY2 << endl;
        */
        cout << "3: " << endl;

        for(int kkk=0 ; kkk<Simulated_Event_Number ; kkk++)
        {
            cout << "//Event: " << kkk;
            gRandom = new TRandom3(0);
            gRandom->SetSeed(0);
            double Random_Velocity = 0;
            Random_Velocity = Flux_HIST->GetRandom();
            Velocity[kkk] = Random_Velocity;
            Flux_HIST_Random->Fill(Random_Velocity);
            Double_t par[3];
            TRandom *eventGenerator = new TRandom(0);//You can use TRandom(0) or TRandom3(0) to initialize your random function
            eventGenerator->GetSeed();
            eventGenerator->Sphere(par[0],par[1],par[2],Velocity[kkk]);
            VX[kkk] = (par[0]/Velocity[kkk]);
            VY[kkk] = (par[1]/Velocity[kkk]);
            VZ[kkk] = (par[2]/Velocity[kkk]);
            double CJPL_Length=0;//km
            if(VZ[kkk]<0) CJPL_Length=0;
            else if(VZ[kkk]>=0)CJPL_Length = Lcjpl(VX[kkk], VY[kkk], VZ[kkk])/1e5;//km
            cout << "WIMP_Mass: " << WIMP_Mass << endl;
            cout << "CJPL_Length: " << CJPL_Length << endl;//km
            double Dark_Matter_Energy = Energy_DM(WIMP_Mass,Velocity[kkk]*1e3/3e8);//KeV
            double Dark_Matter_Velocity = Velocity_DM(WIMP_Mass,Dark_Matter_Energy);//KeV
            cout << "Velocity: " << Velocity[kkk] << endl;
            double Collision_Expectation_velocity_dependent_Earth = 0;
double Collision_Expectation_velocity_dependent_EARTH
            = CDEX_Collision_Time_EARTH(Sigma_SI_Default,VX[kkk],VY[kkk],VZ[kkk],Velocity[kkk],WIMP_Mass,CJPL_Length);
double Collision_Expectation_velocity_dependent_ATM
            = CDEX_Collision_Time_ATM(Sigma_SI_Default,VZ[kkk],Velocity[kkk],WIMP_Mass,CJPL_Length);
            cout << "COLLISION_TIME_EARTH_Exact: " << Collision_Expectation_velocity_dependent_EARTH << endl;
            cout << "COLLISION_TIME_AIR: " << Collision_Expectation_velocity_dependent_ATM << endl;
            Collision_Expectation_ATM[kkk] = Collision_Expectation_velocity_dependent_ATM;
            Collision_Expectation_EARTH[kkk] = Collision_Expectation_velocity_dependent_EARTH;
            Collision_Expectation_Earth[kkk] = Collision_Expectation_velocity_dependent_Earth;

        }
        for(int kkk=0 ; kkk<2000 ; kkk++){Flux_HIST_Random_Normalized->SetBinContent(kkk+1,(Flux_HIST_Random->GetBinContent(kkk))/(Simulated_Event_Number));};
        cout << "4: " << endl;

        //==============================Third Step==============================
        
        TH1F   *Flux_HIST_Aft_Collision_EARTH = new TH1F("Flux_HIST_Aft_Collision_EARTH","Flux_HIST_Aft_Collision_EARTH",2000,0,791);
        TH1F   *Flux_HIST_Aft_Collision_Earth = new TH1F("Flux_HIST_Aft_Collision_Earth","Flux_HIST_Aft_Collision_Earth",2000,0,791);

        double Scaling_Factor=0;double Scaling_Factor_1=0;
        TH1F   *Collision_Time_Hist = new TH1F("Flux_HIST_Aft_Collision","Flux_HIST_Aft_Collision",10,0,10);
        int Event_Number_Check=0;
        TH1F   *Collision_Energy_Lost_E = new TH1F("Collision_Energy_Lost_E","Collision_Energy_Lost_E",100,0,20);
        TH1F   *Collision_Energy_Lost_V = new TH1F("Collision_Energy_Lost_V","Collision_Energy_Lost_V",100,0,20);
        double Check_ZERO_COLLISION=0;double Check_ZERO_COLLISION_1=0;
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
    //====Check the threshold to 1km/s
    int Earth_Threshold=0;
    for(int EEE=0; EEE<5; EEE++)// 1 is earth 2 is air
    {
        int Check = Velocity_Aft_collision(3000,WIMP_Mass,Sigma_SI_Default,799.135,1)[3];
        if(Earth_Threshold<Check)Earth_Threshold=Check;
    }
    cout << "Earth_Threshold: " << Earth_Threshold << endl;
    
    for(int jjj=0; jjj<1 ; jjj++)
    {
        double *Collision_Expectation_EARTH_Temp;
        if(jjj==0) Collision_Expectation_EARTH_Temp=Collision_Expectation_EARTH;
        if(jjj==1) Collision_Expectation_EARTH_Temp=Collision_Expectation_Earth;

        for(int kkk=0; kkk<Simulated_Event_Number; kkk++)
        {
            cout << "//Event2: " << kkk << endl;;
            //cout << "Loop2==>EventKIND: " << jjj << endl;
            //cout << "Loop2=>Event: " << kkk << endl;
            cout << "Earth_Threshold: " << Earth_Threshold << endl;
            cout << "VX[kkk]: " << VX[kkk] << endl;
            cout << "VY[kkk]: " << VY[kkk] << endl;
            cout << "VZ[kkk]: " << VZ[kkk] << endl;

            cout << "Collision_Expectation_ATM[kkk]: " << Collision_Expectation_ATM[kkk] << endl;
            cout << "Collision_Expectation_EARTH_Temp[kkk]: " << Collision_Expectation_EARTH_Temp[kkk] << endl;

            
            if( (int)Collision_Expectation_EARTH_Temp[kkk]> int(Earth_Threshold) ){
               cout << "OK GREAT!" << endl;
                Flux_HIST_Aft_Collision_EARTH->Fill(1e-5);}
            else if(Collision_Expectation_ATM[kkk]>2e3){
                cout << "OK GREAT!" << endl;
                Flux_HIST_Aft_Collision_EARTH->Fill(1e-5);}
             
                    //====================Test_Function====================
            //Two-stage collision
            else
            {
                //cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
                cout << "Collision_Expectation_ATM[kkk]: " << Collision_Expectation_ATM[kkk] << endl;
                cout << "Collision_Expectation_EARTH_Temp[kkk]: " << Collision_Expectation_EARTH_Temp[kkk] << endl;
                
      double *V_Aft_Collision_AIR = Velocity_Aft_collision(Collision_Expectation_ATM[kkk],WIMP_Mass,Sigma_SI_Default,Velocity[kkk],2);
                //cout << "V_Aft_AIR: " << V_Aft_Collision_AIR[0] << endl;
                //cout << "V_Aft_AIR_Energy_Lost_E: " << V_Aft_Collision_AIR[1] << endl;
                //cout << "V_Aft_AIR_Energy_Lost_V: " << V_Aft_Collision_AIR[2] << endl;

    double *V_Aft_Collision_EARTH = Velocity_Aft_collision(Collision_Expectation_EARTH_Temp[kkk],WIMP_Mass,Sigma_SI_Default,V_Aft_Collision_AIR[0],1);
                cout << "V_Aft_EARTH: " << V_Aft_Collision_EARTH[0] << endl;
                cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
                cout << "E_Aft_Earth: " << Energy_DM(WIMP_Mass,V_Aft_Collision_EARTH[0]*1e3/3e8) << ">=0.8?" << endl;
                if(Energy_DM(WIMP_Mass,V_Aft_Collision_EARTH[0]*1e3/3e8)>0.8) Count = Count+1;
                //cout << "V_Aft_EARTH_Energy_Lost_E: " << V_Aft_Collision_EARTH[1] << endl;
                //cout << "V_Aft_EARTH_Energy_Lost_V: " << V_Aft_Collision_EARTH[2] << endl;

                if(V_Aft_Collision_EARTH[0]>Velocity[kkk]) break;
                Collision_Energy_Lost_E->Fill(V_Aft_Collision_AIR[1]+V_Aft_Collision_EARTH[1]);
                Collision_Energy_Lost_V->Fill(V_Aft_Collision_AIR[2]+V_Aft_Collision_EARTH[2]);
                //Check the consistence
                int Collision_Expectation_ATM_Int = Collision_Expectation_ATM[kkk];
                int Collision_Expectation_EARTH_Int = Collision_Expectation_EARTH_Temp[kkk];
                if( (Collision_Expectation_ATM_Int)+(Collision_Expectation_EARTH_Int)==0 ) Check_ZERO_COLLISION_1 = Check_ZERO_COLLISION_1+1;
                if( (V_Aft_Collision_AIR[1]+V_Aft_Collision_EARTH[1])==0 ) Check_ZERO_COLLISION = Check_ZERO_COLLISION+1;
                
                if(jjj==0)Flux_HIST_Aft_Collision_EARTH->Fill(V_Aft_Collision_EARTH[0]);
                if(jjj==1)Flux_HIST_Aft_Collision_Earth->Fill(V_Aft_Collision_EARTH[0]);
                if(jjj==0)E_vs_Collision_Time->Fill((Collision_Expectation_ATM_Int)+(Collision_Expectation_EARTH_Int),V_Aft_Collision_AIR[1]+V_Aft_Collision_EARTH[1]);
            }
        }
    }
            
    TLegend *leg= new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    
    //cout << "Check_ZERO_COLLISION_1: " << Check_ZERO_COLLISION_1 << endl;
    //cout << "Check_ZERO_COLLISION: " << Check_ZERO_COLLISION << endl;
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
        sigma_si = Sigma_SI_Default;
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
    sprintf(fout_name,Form("1_CDEX_Flux/%sGeV/%i.root",Mass_Point[Index].c_str(),Index_Sigma));
    TFile *fout=new TFile(fout_name,"recreate");
    Flux_HIST_Random->Write();
    Flux_HIST_Aft_Collision_Earth->Write();
    Flux_HIST_Aft_Collision_EARTH->Write();
    t1->Write();
    fout->Close();
    cout << "fout_name: " << fout_name << endl;
    
    delete Flux;
    delete Flux_HIST;
    delete Flux_HIST_Random;
    delete Flux_HIST_Random_Normalized;
    delete Flux_HIST_after_Collision;
    delete E_vs_Collision_Time;
}
