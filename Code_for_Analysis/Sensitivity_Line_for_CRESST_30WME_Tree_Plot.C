#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
//#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"

void Sensitivity_Line_for_CRESST_30WME_Tree_Plot()//
{

    string Mass_Point[2]={"2","0P2"};
    double Mass[2]={2,0.2};
    int Bent_or_Not=1;
    string Type[2]={"_Comparison",""};
    vector<double> A;vector<double> B;vector<double> C;vector<double> AC; vector<double> EC;
    vector<double> Energy_Loss_A_all_Array;vector<double> Energy_Loss_S_all_Array;vector<double> Energy_Loss_E_all_Array;

for(int Bent_or_Not=1; Bent_or_Not<2; Bent_or_Not++)
    {
    for(int Mass_INT=1; Mass_INT<2; Mass_INT++)
    {
        
        for(int FILE=31; FILE<50; FILE++)
        {
        string path = Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_CRESST_30MWE_MAT/%sGeV/%i_STS_Bent_Earth.root",Mass_Point[Mass_INT].c_str(),FILE);
            cout << path << endl;
        ifstream fin(path);
            cout << "FILE1: " << FILE << endl;
        if(fin.is_open())
            {
                cout << "FILE2: " << FILE << endl;
        TFile *ROOT_FILE = TFile::Open(path.c_str());
        TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
                
        Double_t mx,sigma_si,V_Int_A,V_End_A,V_End_S,V_End_E;
        Int_t    ev,Arrival_air,Arrival_earth,Bent_or_Not;
    Double_t Oringal_Length_Air,Path_Length_Air,Oringal_Length_Earth,Path_Length_Earth;
        Double_t  Collision_Time_Earth,Collision_Time_Air,Collision_Time_30MWE;
        Double_t Energy_Loss_A,Energy_Loss_S,Energy_Loss_E;
                
        T1_TREE->SetBranchAddress("mx",&mx);
        T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
        T1_TREE->SetBranchAddress("V_Int_A",&V_Int_A);
        T1_TREE->SetBranchAddress("V_End_A",&V_End_A);
        T1_TREE->SetBranchAddress("V_End_S",&V_End_S);
        T1_TREE->SetBranchAddress("V_End_E",&V_End_E);

        T1_TREE->SetBranchAddress("ev",&ev);
        T1_TREE->SetBranchAddress("Arrival_air",&Arrival_air);
        T1_TREE->SetBranchAddress("Arrival_earth",&Arrival_earth);
        T1_TREE->SetBranchAddress("Bent_or_Not",&Bent_or_Not);

        T1_TREE->SetBranchAddress("Oringal_Length_Air",&Oringal_Length_Air);
        T1_TREE->SetBranchAddress("Oringal_Length_Earth",&Oringal_Length_Earth);
        T1_TREE->SetBranchAddress("Path_Length_Air",&Path_Length_Air);
        T1_TREE->SetBranchAddress("Path_Length_Earth",&Path_Length_Earth);
        
        T1_TREE->SetBranchAddress("Collision_Time_Earth",&Collision_Time_Earth);
        T1_TREE->SetBranchAddress("Collision_Time_Air",&Collision_Time_Air);
        T1_TREE->SetBranchAddress("Collision_Time_30MWE",&Collision_Time_30MWE);
        
        T1_TREE->SetBranchAddress("Energy_Loss_A",&Energy_Loss_A);
        T1_TREE->SetBranchAddress("Energy_Loss_S",&Energy_Loss_S);
        T1_TREE->SetBranchAddress("Energy_Loss_E",&Energy_Loss_E);

                TH1F   *H_Energy_Loss_A = new TH1F("H_Energy_Loss_A","H_Energy_Loss_A",1000,0,1);
                TH1F   *H_Energy_Loss_S = new TH1F("H_Energy_Loss_S","H_Energy_Loss_S",1000,0,1);
                TH1F   *H_Energy_Loss_E = new TH1F("H_Energy_Loss_E","H_Energy_Loss_E",1000,0,1);

                TH1F   *H_Length_Earth = new TH1F("H_Length_Air","H_Length_Air ",100,0,1000);
                TH1F   *H_Length_Air = new TH1F("H_Length_Air","H_Length_Air ",100,0,1000);

                TH1F   *Air_Length_per_collision   = new TH1F("Air_Length_per_collision","Air_Length_per_collision ",100,0,1000);
                TH1F   *Earth_Length_per_collision = new TH1F("Earth_Length_per_collision","Earth_Length_per_collision ",100,0,1000);

                double Energy_Loss_A_all=0;double Energy_Loss_S_all=0;double Energy_Loss_E_all=0;
                int kkk=0;
                for(int Entry=0; Entry<T1_TREE->GetEntries(); Entry++)
                {
                    T1_TREE->GetEntry(Entry);
                    double Initial_Energy=Energy_DM(mx,V_Int_A*(1e3/3e8) );
                    double Energy_Aft_A  =Energy_DM(mx,V_End_A*(1e3/3e8) );
                    double Energy_Aft_S  =Energy_DM(mx,V_End_S*(1e3/3e8) );
                    double Energy_Aft_E  =Energy_DM(mx,V_End_E*(1e3/3e8) );
                    

                    if( Energy_Loss_A>0 )
                    {
                        kkk = kkk + 1;
                        Energy_Loss_A_all = Energy_Loss_A_all + (Initial_Energy-Energy_Aft_A)/Initial_Energy;
                        Energy_Loss_S_all = Energy_Loss_S_all + (Energy_Aft_A-Energy_Aft_S)/Initial_Energy;
                        Energy_Loss_E_all = Energy_Loss_E_all + (Energy_Aft_S-Energy_Aft_E)/Initial_Energy;

                        cout << "V_Int_A: " << V_Int_A << endl;
                        cout << "V_End_A: " << V_End_A << endl;
                        cout << "(Initial_Energy-Energy_Aft_A)/Initial_Energy : " << (Initial_Energy-Energy_Aft_A)/Initial_Energy << endl;
                        cout << "Collision_Time_Air: " << Collision_Time_Air << endl;
                        cout << "Path_Length_Air/Collision_Time_Air: " << Path_Length_Air/Collision_Time_Air << endl;
                        cout << "(Energy_Aft_A-Energy_Aft_S)/Initial_Energy: "    << (Energy_Aft_A-Energy_Aft_S)/Initial_Energy << endl;
                        cout << "(Energy_Aft_S-Energy_Aft_E)/Initial_Energy: " << (Energy_Aft_S-Energy_Aft_E)/Initial_Energy << endl;
                        cout << "Path_Length_Earth/Collision_Time_Earth: " << Path_Length_Earth/Collision_Time_Earth << endl;
                        cout << "Collision_Time_Earth: " << Collision_Time_Earth << endl;

                        if(kkk>1)
                        {
                            H_Energy_Loss_A->Fill( (Initial_Energy-Energy_Aft_A)/Initial_Energy );
                            H_Energy_Loss_S->Fill( (Energy_Aft_A-Energy_Aft_S)/Initial_Energy );
                            H_Energy_Loss_E->Fill( (Energy_Aft_S-Energy_Aft_E)/Initial_Energy   );
                            if(Collision_Time_Air>0)Air_Length_per_collision->Fill(Path_Length_Air/Collision_Time_Air);//AC
                            if(Collision_Time_Earth>0)Earth_Length_per_collision->Fill(Path_Length_Earth/Collision_Time_Earth);//EC
                        }
                    }
                }//for(int Entry=0; Entry<T1_TREE->GetEntries(); Entry++)
                Energy_Loss_A_all_Array.push_back(Energy_Loss_A_all/5000.);
                Energy_Loss_S_all_Array.push_back(Energy_Loss_S_all/5000.);
                Energy_Loss_E_all_Array.push_back(Energy_Loss_E_all/5000.);
                A.push_back(H_Energy_Loss_A->GetMean());
                B.push_back(H_Energy_Loss_S->GetMean());
                C.push_back(H_Energy_Loss_E->GetMean());
                AC.push_back(Air_Length_per_collision->GetMean());
                EC.push_back(Earth_Length_per_collision->GetMean());

                char fout_name[100];
    sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/2_CRESST_30MWE_MAT/%sGeV/Hist_%i_STS_Bent_Earth.root",Mass_Point[Mass_INT].c_str(),FILE));
                TFile *fout=new TFile(fout_name,"recreate");
                H_Energy_Loss_A->Write();
                H_Energy_Loss_S->Write();
                H_Energy_Loss_E->Write();
                Air_Length_per_collision->Write();
                Earth_Length_per_collision->Write();
                
                cout << "Yes  " << endl;
                fout->Close();
            }//fin.is_open()
        }//for(int FILE=0; FILE<4; FILE++)
        for(int kkk=0; kkk<20; kkk++)
        {
            cout << "=========================" << endl;
            //cout << "Energy_loss_all_A: " << Energy_Loss_A_all_Array[kkk] << endl;
            //cout << "Energy_loss_all_S: " << Energy_Loss_S_all_Array[kkk] << endl;
            //cout << "Energy_loss_all_E: " << Energy_Loss_E_all_Array[kkk] << endl;
            cout << "A: " << A[kkk] << endl;
            //cout << "AC: " << AC[kkk] << endl;
            cout << "S: " << B[kkk] << endl;
            cout << "E: " << C[kkk] << endl;
            //cout << "EC: " << EC[kkk] << endl;
            cout << "=========================" << endl;
        }
        
    }//for(int Mass_INT=0; Mass_INT<4; Mass_INT++)
    //==============Mass==============//
    }
}
    

