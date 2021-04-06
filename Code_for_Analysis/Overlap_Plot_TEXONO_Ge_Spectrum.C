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

void Overlap_Plot_TEXONO_Ge_Spectrum()
{
    TCanvas *c3 = new TCanvas("c3");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    const int Number=29;
    double Sigma_SI_Array[Number];
    //Method1 Threshold: 200eV
    double Sigma_SI_With_Threshold_M1[Number];double Sigma_SI_With_Threshold_earth_M1[Number];
    double Ratio_for_Average_M1[Number];double Ratio_for_PREM_M1[Number];

    //Method2 Threshold: 200km/s
    double Sigma_SI_With_Threshold_M2[Number];double Sigma_SI_With_Threshold_earth_M2[Number];
    double Ratio_for_Average_M2[Number];double Ratio_for_PREM_M2[Number];

    //Method1 Threshold: 200eV-300eV
    double Sigma_SI_With_Threshold_M3[Number];double Sigma_SI_With_Threshold_earth_M3[Number];
    double Ratio_for_Average_M3[Number];double Ratio_for_PREM_M3[Number];
    double Sigma_SI_With_Threshold_M3_Error[Number];double Error_X[Number];

    double CPKKD_EXCLUSION[Number];
    double Mass=1;//2.34
    int Take_Plot=0;
    string Type_of_Model="BR"; int Type_of_Model_INT=2;
    cout << "max_recoil_A_EM_keV(): " << max_recoil_A_EM_keV(2.34, 779.135*1000.0/2.99792458e8, AGe) << endl;
        
    int Point_Number=0;
     
    //Sigma_SI_Array[Point_Number]=Multiply_Number*TMath::Power(10,-YYY);
    Sigma_SI_Array[Point_Number] = 1e-36;
    TH1F *Flux_HIST_Random; TH1F *Flux_HIST_Aft_Collision_Earth; TH1F *Flux_HIST_Aft_Collision_EARTH;

    TFile *fin = TFile::Open("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/1GeV/2.root");
    Flux_HIST_Random=(TH1F*)fin->Get("Flux_HIST_Random");
    Flux_HIST_Aft_Collision_EARTH=(TH1F*)fin->Get("Flux_HIST_Aft_Collision_EARTH");

    double T_QF_Original_Bef_Array[reso_T]; double T_QF_Original_Aft_Array[reso_T];
    double Factor1_Original_Bef_Array[reso_T]; double Factor1_Original_Aft_Array[reso_T];

    double *T_QF_Original_Bef=RecoilX_Event(0,Flux_HIST_Random,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        T_QF_Original_Bef_Array[i]=T_QF_Original_Bef[i];}
        cout << "======================================" << endl;
    double *T_QF_Original_Aft=RecoilX_Event(0,Flux_HIST_Aft_Collision_EARTH,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        T_QF_Original_Aft_Array[i]=T_QF_Original_Aft[i];}
        cout << "======================================" << endl;

    double *Factor1_Original_Bef=RecoilX_Event(1,Flux_HIST_Random,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        Factor1_Original_Bef_Array[i]=(Factor1_Original_Bef[i]);}
        cout << "======================================" << endl;

    double *Factor1_Original_Aft=RecoilX_Event(1,Flux_HIST_Aft_Collision_EARTH,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        Factor1_Original_Aft_Array[i]=(Factor1_Original_Aft[i]);}
        cout << "Sigma_SI_Array[Point_Number]: " << Sigma_SI_Array[Point_Number] << endl;
        cout << "======================================" << endl;

    double RecoilX_Event_Original_M1=0; double RecoilX_Event_Aft_EARTH_M1=0;
    double RecoilX_Event_Original_M3=0; double RecoilX_Event_Aft_EARTH_M3=0;
    
        

    double SetOutline_X[1]={0};double SetOutline_Y[1]={0};
    TGraph *SetOutLine = new TGraph(1,SetOutline_X,SetOutline_Y);
    SetOutLine->GetXaxis()->SetRangeUser(0,3);
    SetOutLine->GetYaxis()->SetRangeUser(0,1e+15);

    //Energy recoil Spectrum
    TGraph *ER_Spectrum_Bef = new TGraph(reso_T,T_QF_Original_Bef_Array,Factor1_Original_Bef_Array);
    ER_Spectrum_Bef->GetXaxis()->SetTitle("Energy[keV]");
    ER_Spectrum_Bef->GetYaxis()->SetTitle("Count");
    
    ER_Spectrum_Bef->GetXaxis()->SetLimits(1e-2,1e1);
    ER_Spectrum_Bef->GetYaxis()->SetRangeUser(1e-7,1e7);
     
    ER_Spectrum_Bef->SetLineColor(2);
    ER_Spectrum_Bef->SetLineWidth(5);
    
    ER_Spectrum_Bef->Draw();
    
    char fout_name[100];
    sprintf(fout_name,Form("Thesis_Plot/Spectrum_Ex/%s.root",Type_of_Model.c_str()));
    TFile *fout=new TFile(fout_name,"recreate");
    ER_Spectrum_Bef->Write();

}
    

