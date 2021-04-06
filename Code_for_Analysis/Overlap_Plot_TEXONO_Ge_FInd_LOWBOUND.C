#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"

void Overlap_Plot_TEXONO_Ge_FInd_LOWBOUND()//CDEX:Threshold=160eV, TEXONO:Threshold=200eV
{
    //double WIMP_Mass_Array[16]={2,1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05};//WIMP
    //double WIMP_Mass_Array[10]={0.2,0.19,0.18,0.17,0.16,0.15,0.14,0.13,0.12,0.11};
    //double WIMP_Mass_Array[12]={20,18,16,14,12,10,8,7,6,4,3,2.35};
    double WIMP_Mass_Array[12]={20,19,17,15,13,11,9,7,5,4,3,2.35};
    
    //double WIMP_Mass_Array[11]={1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20};//MIMP
    double Sigma_SI[16];//WIMP
    int Bin_Dominant[16];//WIMP
    string Type_of_Model[4]={"NU","MD","BR","MDMPA"};
    //double Sigma_SI[11];//MIMP
for(int kkk=0; kkk<12;kkk++)
    {
    //double CPKKD_EXCLUSION[Number];
    double Mass = WIMP_Mass_Array[kkk];//GeV
    int Type_of_Model_INT=0;int Data_element_N =4;
    //=======================Recoil Spectrum set for three processes==============================
    double T_QF_Original_Bef_Array[reso_T];double Factor1_Original_Bef_Array[reso_T];
        double Sig=1e-28;
    double *T_QF_Original_Bef=RecoilX_Event(0,0,Mass,Sig,Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){T_QF_Original_Bef_Array[i]=T_QF_Original_Bef[i];}
    double *Factor1_Original_Bef=RecoilX_Event(1,0,Mass,Sig,Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){Factor1_Original_Bef_Array[i]=(Factor1_Original_Bef[i]);}
    //=======================Calculate the Ratio passing through====================//
        
    double *Data_RE_CDEX    =Hist_SetLimit_Plot_v2_Extract_Peak(0);
    double *Data_RATE_CDEX  =Hist_SetLimit_Plot_v2_Extract_Peak(1);
    double *Data_RE_Err     =Hist_SetLimit_Plot_v2_Extract_Peak(2);
    double *Data_RE_Rate_Err=Hist_SetLimit_Plot_v2_Extract_Peak(3);
         
        
        /*
double RE_DATA_Original[Data_element]; double RE_Rate_Original[Data_element]; double RE_DATA_Err_Original[Data_element]; double RE_Rate_Err_Original[Data_element];
        for(int jjj=0 ; jjj<257 ; jjj++)
        {
                RE_DATA_Original[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
                RE_Rate_Original[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
                RE_DATA_Err_Original[jjj]= 0;
                RE_Rate_Err_Original[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*(1.64458/0.994);
            //if(jjj==0) cout << "RE_Rate_Aft[jjj]+RE_Rate_Err_Aft[jjj]=2: " << RE_Rate_Original[jjj]+RE_Rate_Err_Original[jjj] << endl;
        
        }
        double *Factor_TEXONO=cpkkd_calculation_Scaling_Factor(0,T_QF_Original_Bef,Factor1_Original_Bef,RE_DATA_Original,RE_Rate_Original,RE_DATA_Err_Original,RE_Rate_Err_Original,Mass);
*/
        /*
    double *Factor_TEXONO=cpkkd_calculation_Scaling_Factor(0,T_QF_Original_Bef,Factor1_Original_Bef,RE_DATA_1,RE_Rate_1,RE_DATA_Err_1,RE_Rate_Err_1,Mass);
         
        double *Factor_TEXONO=cpkkd_calculation_Scaling_Factor(0,T_QF_Original_Bef,Factor1_Original_Bef,RE_DATA_Original,RE_Rate_Original,RE_DATA_Err_Original,RE_Rate_Err_Original,Mass);
    cout << "1e-36*Factor" << Factor_TEXONO[0]*1e-31 << endl;
        cout << "Factor" << Factor_TEXONO[0]<< endl;
    Sigma_SI[kkk]=Factor_TEXONO[0]*1e-31;
         */
        
        //For CDEX
        /*
        double Data_RE_CDEX[4]= {0.208,0.30,0.35,0.4};//CDEX: 0.208
        double Data_RATE_CDEX[4]= {7.3,7.6,7.4,7.9};
        double Data_RE_Err[4]= {0,0,0,0};
        double Data_RE_Rate_Err[4]= {2.7,3.5,1.5,2};
         */
        /*
        double Data_RE_CDEX[4]= {0.3,0,0,0};//CDEX Brem: 0.265(for threshold 0.25), CDEX Migdal: 0.2, CDEX Brem: 0.206
        double Data_RATE_CDEX[4]= {7.6,0,0,0};
        double Data_RE_Err[4]= {0,0,0,0};
        double Data_RE_Rate_Err[4]= {3.5,0,0,0};
         */
        //Old Data from Dr. HBLi
        /*
        double Data_RE_CDEX[4]= {0.5,0,0,0};//CDEX Brem: 0.265(for threshold 0.25), CDEX Migdal: 0.2, CDEX Brem: 0.206
        double Data_RATE_CDEX[4]= {9,0,0,0};
        double Data_RE_Err[4]= {0,0,0,0};
        double Data_RE_Rate_Err[4]= {3,0,0,0};
         */
    double *Factor_TEXONO=cpkkd_calculation_Scaling_Factor(Sig,0,T_QF_Original_Bef,Factor1_Original_Bef,Data_RE_CDEX,Data_RATE_CDEX,Data_RE_Err,Data_RE_Rate_Err,Mass);


    //=======================Touching Line========================//
    double T[dm_spec_resolution];double Count[dm_spec_resolution];
    for(int kkk=0; kkk<dm_spec_resolution; kkk++)
    {
        T[kkk]    =T_QF_Original_Bef[kkk];
        Count[kkk]=Factor_TEXONO[kkk+2];
    }
        
    Sigma_SI[kkk] = Factor_TEXONO[0] * Sig;
    Bin_Dominant[kkk] = Factor_TEXONO[1];
    cout << "Sigma_SI[kkk] : " << Sigma_SI[kkk]  << endl;
    TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.05,"XY");
    gStyle->SetTitleFont(62,"XY");
    gStyle->SetLegendFont(62);
    //=================Theoretical Line==============//
    TGraph *gcpkkdX1 = new TGraph(dm_spec_resolution,T,Count);
    gcpkkdX1->SetName("gcpkkdX1");
    gcpkkdX1->SetLineColor(4);
    gcpkkdX1->SetName("RS");
    //=================Data Points==============//
    //TGraphErrors *TEXONOData = new TGraphErrors(4,RE_DATA_Original,RE_Rate_Original,RE_DATA_Err_Original,RE_Rate_Err_Original);
    TGraphErrors *TEXONOData = new TGraphErrors(257,Data_RE_CDEX,Data_RATE_CDEX,Data_RE_Err,Data_RE_Rate_Err);
    TEXONOData->SetName("TEXONOData");
    TEXONOData->SetLineColor(2);
    TEXONOData->GetXaxis()->SetLimits(0,2.4);
    TEXONOData->GetYaxis()->SetRangeUser(1,80);
    TEXONOData->SetTitle("The rate of WIMP");
    TEXONOData->GetXaxis()->SetTitle("Edet[keVee]");
    TEXONOData->GetYaxis()->SetTitle("Counts[Evts/day/kg/keV]");
    TEXONOData->SetMarkerStyle(21);
    
    TLegend *leg = new TLegend(0.3,0.5,0.5,0.8);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->AddEntry("",Form("M_{#chi}=%.2fGeV, #sigma^{SI}_{#chi N}:%.2f #times 10^{-%.f}",WIMP_Mass_Array[kkk],CNFNV(0,Sigma_SI[kkk]),CNFNV(1,Sigma_SI[kkk])),"");
    leg->AddEntry("",Form("Dominant Bin: %i",(int)Bin_Dominant[kkk]+1),"");

    TEXONOData->Draw("AP");
    gcpkkdX1->Draw("LPsame");
    leg->Draw();
    c1->Print(Form("LOW_BOUND/High_Mass_Plot/%s_%i_Bin1_TEXONO.pdf",Type_of_Model[Type_of_Model_INT].c_str(),kkk));
        
        char fout_name[300];
        sprintf(fout_name,Form("LOW_BOUND/ROOT/High_Mass/%s_%i_Bin1_TEXONO.root",Type_of_Model[Type_of_Model_INT].c_str(),kkk));
        TFile *fout=new TFile(fout_name,"recreate");
        gcpkkdX1->Write();
        TEXONOData->Write();

    }
    
    
    
    for(int kkk=0; kkk<16;kkk++)
    {
        cout << WIMP_Mass_Array[kkk] << "," << Sigma_SI[kkk] << "," << endl;
    }
    cout << "======================================" << endl;
     for(int kkk=0; kkk<16;kkk++)
     {
         cout << WIMP_Mass_Array[kkk] << "," << Bin_Dominant[kkk] << "," << endl;
     }
    for(int kkk=0; kkk<16;kkk++)
    {
        cout << Bin_Dominant[kkk] << "," ;
    }

}
    

