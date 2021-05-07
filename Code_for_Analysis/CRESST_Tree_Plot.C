#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"

void CRESST_Tree_Plot()//
{

    string path = Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/%i_STS_Bent.root",Mass_Point[Mass_INT].c_str(),FILE);
    ifstream fin(path);
    
    TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/%i_STS_Bent.root",Mass_Point[Mass_INT].c_str(),FILE));
    TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
    Double_t mx,sigma_si;
    T1_TREE->SetBranchAddress("mx",&mx);T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
    T1_TREE->GetEntry(0);
    TH1F   *Flux_HIST_after_Collision = new TH1F("Flux_HIST_after_Collision","Flux_HIST_after_Collision",2000,0,791);
    
    char fout_name[100];
    sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/4_CRESST_Flux/%sGeV/%i_STS_Bent.root",Mass_Point[Index].c_str(),Index_Sigma));
    TFile *fout=new TFile(fout_name,"recreate");
    Flux_HIST_Random->Write();
    fout->Close();
    //==============Mass==============//
}
    

