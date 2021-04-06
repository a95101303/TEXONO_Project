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
double WIMP_Mass=0.08;
void Hist_SetLimit_Plot_v2_Possion_CDEX()
{
    /*
    for (int i=0;i<n_point1;i++) {
        for (int j=0;j<lnumber;j++) {
            cout << "Ee[i]" << *(Eem_keV_Array+i) << endl;
            cout << "Eem_keV[i]: " << Eem_keV[i] << endl;
            cout << "dpdE[i][nli]: " << dpdE[i][j] << endl;
              cout << "Prob" << *(dpdE_Array+(i*lnumber)+j) << endl; // energy in eV
          }}
    */
    int Count=0;
for(int ppp=1; ppp<10; ppp++ )
{
    for(int zzz=0; zzz<1; zzz++)//Z<10
    {
        for(int lll=3; lll<4; lll++){
            double Sigma_SI_Default= (ppp+0.2*zzz) * (1e-27) * TMath::Power(10,-lll);
            //double Sigma_SI_Default= (zzz) * (1e-30) * TMath::Power(10,-lll);
            //double Sigma_SI_Default=(1e-29);
            cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
            /*
            TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
            gStyle->SetOptStat(0);
            gStyle->SetTitleSize(0.05,"XY");
            gStyle->SetTitleFont(62,"XY");
            gStyle->SetLegendFont(62);
            */
            //Find the Possion distributions
            /*
             for(int kkk=0 ; kkk<100 ; kkk++)
             {
             {
             cout << "TMath::Poisson(1,0.01*kkk): " << TMath::Poisson(1,0.01*kkk) << endl;
             cout << "0.01*kkk: " << 0.01*kkk << endl;
             }
             }
             TF1 *f1 = new TF1("f1","TMath::Poisson(x,0.001)",0,2);
             f1->Draw();
             */
             //==============================First step==============================
            /*
            cout << "Mean_Free_Path_of_Eath: " << Mean_Free_Path_of_Earth(0,0) << endl;
            double Collision_Expectation = Collision_Time(2*Earth_Radius,Mean_Free_Path_of_Earth(0,0));
            cout << "Collision_Time: " << Collision_Expectation << endl;
            TF1 *f1 = new TF1("f1","TMath::Poisson(x,[0])",0,2);
            f1->SetParameter(0,Collision_Expectation);
            TMath::Poisson(1,Collision_Expectation);
            f1->Draw();
             */
            //==============================Second step==============================
            cout << "1: " << endl;
            
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

            int Simulated_Event_Number=1000;
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
                cout << "VX[kkk]: " << VX[kkk] << endl;
                cout << "VY[kkk]: " << VY[kkk] << endl;
                cout << "VZ[kkk]: " << VZ[kkk] << endl;

                cout << "Collision_Expectation_ATM[kkk]: " << Collision_Expectation_ATM[kkk] << endl;
                cout << "Collision_Expectation_EARTH_Temp[kkk]: " << Collision_Expectation_EARTH_Temp[kkk] << endl;

                
                if(Collision_Expectation_EARTH_Temp[kkk]>200){
                   // cout << "OK GREAT!" << endl;
                    Flux_HIST_Aft_Collision_EARTH->Fill(1e-5);}
                else if(Collision_Expectation_ATM[kkk]>200){
                    //cout << "OK GREAT!" << endl;
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
            
            /*
            Double_t norm = 1;
            Double_t scale = norm/(Flux_HIST_Random->Integral());
            Flux_HIST_Random->Scale(scale);

            Double_t scale1 = norm/(Flux_HIST_Aft_Collision_EARTH->Integral());
            Flux_HIST_Aft_Collision_EARTH->Scale(scale1);
            
            Double_t scale2 = norm/(Flux_HIST_Aft_Collision_Earth->Integral());
            Flux_HIST_Aft_Collision_Earth->Scale(scale2);

            leg->AddEntry(Flux_HIST_Random,"Original_Fluc","l");
            leg->AddEntry(Flux_HIST_Aft_Collision_Earth,"Earth_Average_desntity","l");
            leg->AddEntry(Flux_HIST_Aft_Collision_EARTH,"Earth_non_Average_desntity","l");
            
            Flux_HIST_Random->GetYaxis()->SetRangeUser(0,0.07);
            
            Flux_HIST_Random->Draw("HIST");
            Flux_HIST_Aft_Collision_Earth->Draw("HISTsame");
            Flux_HIST_Aft_Collision_EARTH->Draw("HISTsame");
            leg->Draw();
            c1->Print("Yes.pdf");
            
            char fout_name[100];
            int Sigma_Value=abs(TMath::Log10(Sigma_SI_Default));
             */
            
            char fout_name[100];
            sprintf(fout_name,Form("Sigma_SI_Flux_CDEX_0P08GeV/%iP%ie_M%i.root",ppp,2*zzz,27+lll));
            TFile *fout=new TFile(fout_name,"recreate");
            Flux_HIST_Random->Write();
            Flux_HIST_Aft_Collision_Earth->Write();
            Flux_HIST_Aft_Collision_EARTH->Write();
            
            cout << "Count: " << Count << endl;
            //double PX=0; double PY=0; double PZ=sqrt(1-PX*PX-PY*PY);
        }
    }
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

