#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
//#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"

void CRESST_Tree_Plot()//
{

    string Mass_Point[3]={"20","2","0P2"};
    double Mass[3]={20,2,0.2};
    int Bent_or_Not=0;
    string Type[2]={"_Comparison",""};

    for(int Mass_INT=0; Mass_INT<3; Mass_INT++)
    {
        for(int FILE=1; FILE<2; FILE++)
        {
            string path = Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/%sGeV/%i_STS_Bent%s.root",Mass_Point[Mass_INT].c_str(),FILE,Type[Bent_or_Not].c_str());
    ifstream fin(path);
        if(fin.is_open())
            {
        TFile *ROOT_FILE = TFile::Open(path.c_str());
        TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
                
        Double_t mx,sigma_si,V_Int_A,V_End_A,V_Int_E,V_End_E;
        Int_t    ev,Arrival_air,Arrival_earth,Bent_or_Not;
    Double_t Oringal_Length_Air,Path_Length_Air,Oringal_Length_Earth,Path_Length_Earth;
        Double_t  Collision_Time_Earth,Collision_Time_Air;
        Double_t Energy_Loss_Percentage_lf;
                
        T1_TREE->SetBranchAddress("mx",&mx);
        T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
        T1_TREE->SetBranchAddress("V_Int_A",&V_Int_A);
        T1_TREE->SetBranchAddress("V_End_A",&V_End_A);
        T1_TREE->SetBranchAddress("V_Int_E",&V_Int_E);
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
    T1_TREE->SetBranchAddress("Energy_Loss_Percentage_lf",&Energy_Loss_Percentage_lf);

                TH1F   *Energy_Diff = new TH1F("Energy_Diff","Energy_Diff",100,0,1);
                TH1F   *Earth_Survive = new TH1F("Earth_Survive","Earth_Survive",3,0,3);

                for(int Entry=0; Entry<T1_TREE->GetEntries(); Entry++)
                {
                    T1_TREE->GetEntry(Entry);
                    Earth_Survive->Fill(Arrival_earth);
                    if(Arrival_earth==1)
                    {
                    Energy_Diff->Fill(Energy_Loss_Percentage_lf);
                    }
                }//for(int Entry=0; Entry<T1_TREE->GetEntries(); Entry++)
                char fout_name[100];
    sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/%sGeV/Hist_%i_STS_Bent%s.root",Mass_Point[Mass_INT].c_str(),FILE,Type[Bent_or_Not].c_str()));
                TFile *fout=new TFile(fout_name,"recreate");
                Energy_Diff->Write();
                Earth_Survive->Write();
                cout << "Yes  " << endl;
                fout->Close();
            }//fin.is_open()
        }//for(int FILE=0; FILE<4; FILE++)
    }//for(int Mass_INT=0; Mass_INT<4; Mass_INT++)
    //==============Mass==============//
}
    

