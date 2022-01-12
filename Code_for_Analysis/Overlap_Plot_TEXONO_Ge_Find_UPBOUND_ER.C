#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
//#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "velocity_distribution_2000_Ave_ER.h"

#include "cpkkd_calculation_New.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/DM_Energy_Loss_Predict.h"


//Unit: Mb/eV from Mukesh, T(keV)
double DM_E(double Mx, double v)//Mx(GeV/c^2),v(km/s)
{
    return 0.5*(Mx*1e6)*(v*1e3/3e8)*(v*1e3/3e8);//keV
}
double Energy_Loss(double a, double b, double Energy_DM_Max, double Scaling)
{
    //cout << "Energy_DM_Max: " << Energy_DM_Max << endl;
    TF1 *fitting_Line = new TF1("fitting_Line","TMath::Power(10,([0]*x+[1])*[2])",0,Energy_DM_Max);
    fitting_Line->SetParameter(0,a);
    fitting_Line->SetParameter(1,b);
    fitting_Line->SetParameter(2,Scaling);

    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);
    double Energy_Loss = fitting_Line->GetRandom();
    //cout << "Energy_Loss: " << Energy_Loss << endl;

    return Energy_Loss;
}
 //Test_fitting_on_TGraph
/*
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()
{
    double x[5] = {1, 2, 3, 4, 5};
    double y[5] = {1.3, 1.4, 1.2, 1.25, 1.35};
    vector<int> AAA;
    TGraph * g = new TGraph((int)AAA.size(), x, y);

    g->SetMarkerStyle(20);
    g->SetMarkerColor(2);
    g->Draw("ap");

    TF1 * f = new TF1("func", "[0]", 0, 6);
    f->SetLineStyle(2);
    f->SetLineColor(8);
    g->Fit(f);
}
 */
/*
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()//Test two fittings
//Test the fitting line
{
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    
    const int Index=0;
    const double Density = 1.8;//g/cm^3
    const double Length  = 1e5;//cm
    //const double Length  = 3e0;//cm
    const double Threshold_keV = 0.1;//keV
    vector<string> File;
    vector<double> WIMP_mx_Array;
    vector<double> Collision_Time_Array;
    vector<double> Energy_Loss_Per_Collision;
    //vector<double> Scaling={1e-18,1e-9};
    vector<double> Scaling={1e-12,1e-9};

    //=========Define the file==========//
    //Xe_c1
    if(Index==0)
    {
    File=
            {"c1_2_0GeV","c1_1_0GeV",
            "c1_0_900GeV","c1_0_800GeV",
            "c1_0_700GeV","c1_0_600GeV",
            "c1_0_500GeV","c1_0_400GeV",
            "c1_0_300GeV","c1_0_200GeV",
            "c1_0_120GeV","c1_0_090GeV",
            "c1_0_070GeV","c1_0_050GeV"};
        WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};
    }
    //Xe_d1
    if(Index==1)
    {
         File=
            {"d1_2_0GeV","d1_1_0GeV",
            "d1_0_900GeV","d1_0_800GeV",
            "d1_0_700GeV","d1_0_600GeV",
            "d1_0_500GeV","d1_0_400GeV",
            "d1_0_300GeV","d1_0_200GeV",
            "d1_0_120GeV","d1_0_090GeV",
            "d1_0_070GeV","d1_0_050GeV"};
        WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};
    }
    //Read the file of DCS for different masses
    
    
    string c1_d1_Xe_Ge_index[4]={"SI_c1_XeData_Vel","SI_d1_XeData_Vel","SI_c1_GeData_Vel","SI_d1_GeData_Vel"};
    //===============electron recoil===================//
    //for(int kkk=0; kkk<File.size(); kkk++)
    for(int kkk=1; kkk<2; kkk++)
    {
        cout << "WIMP_mx_Array: " << WIMP_mx_Array[kkk] << endl;
        TFile *fin   = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str(),File[kkk].c_str()));
        TFile *fin_2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index+2].c_str(),File[kkk].c_str()));

        TH1F*    velocity_TH1F[velocity_N]   ;TH1F*    velocity_TH1F_2[velocity_N];//Mass-based

        for(int LLL=0; LLL<velocity_N; LLL++)
        {
            TH1F *Velocity_TH1F_temp  =(TH1F*)fin  ->Get(velocity_s[LLL].c_str());
            TH1F *Velocity_TH1F_temp_2=(TH1F*)fin_2->Get(velocity_s[LLL].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                cout << "File: "    << velocity_s[LLL].c_str() << endl;
                velocity_TH1F[LLL]    = (Velocity_TH1F_temp);
            }
            if(Velocity_TH1F_temp_2!=NULL)
            {
                cout << "File_2: " << velocity_s[LLL].c_str() << endl;
                velocity_TH1F_2[LLL] = (Velocity_TH1F_temp_2);
            }

        }
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        for(int Applied_Hist=0; Applied_Hist<velocity_N; Applied_Hist++)
        //for(int Applied_Hist=velocity_TH1F.size()-1; Applied_Hist<velocity_TH1F.size(); Applied_Hist++)
        {
            double Total_Cross_Section_Xe=0;double Total_Cross_Section_Ge=0;
            double Mean_Value_Loss_Xe    =0;double Mean_Value_Loss_Ge    =0;
            if((int)velocitykm[Applied_Hist]==300)
            {
                
            cout << "==============================Xe==============================" << endl;
            
                Int_t Minimum_Bin_Xe = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(80);//Min_Recoil_Bin_Xe
                Int_t Maximum_Bin_Xe = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
                for(int LLL=Minimum_Bin_Xe; LLL<Maximum_Bin_Xe+1; LLL++)
                {
                    double DCS_Inidividual = velocity_TH1F[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2);
                    Total_Cross_Section_Xe = Total_Cross_Section_Xe + DCS_Inidividual;
                }
                cout << "Total_Cross_Section_Xe: " << Total_Cross_Section_Xe << endl;
                for(int LLL=Minimum_Bin_Xe; LLL<Maximum_Bin_Xe+1; LLL++)
                {
                    double DCS_Inidividual  = velocity_TH1F[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2);
                    double Recoil_Energy_Xe = velocity_TH1F[Applied_Hist]->GetBinCenter(LLL);//eV
                    Mean_Value_Loss_Xe = Mean_Value_Loss_Xe + Recoil_Energy_Xe*(DCS_Inidividual/Total_Cross_Section_Xe);
                                         
                }
                cout << "velocity_TH1F[Applied_Hist]->GetBinCenter(LLL): " << velocity_TH1F[Applied_Hist]->GetBinCenter(Maximum_Bin_Xe) << endl;
                cout << "Mean_Xe: " << Mean_Value_Loss_Xe << endl;;
                 
                
            cout << "==============================Ge==============================" << endl;
                Int_t Minimum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->GetXaxis()->FindBin(80);//Min_Recoil_Bin_Ge
                Int_t Maximum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
                for(int LLL=Minimum_Bin_Ge; LLL<Maximum_Bin_Ge+1; LLL++)
                {
                    double DCS_Inidividual = velocity_TH1F_2[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2);
                    Total_Cross_Section_Ge = Total_Cross_Section_Ge + DCS_Inidividual;
                    //cout << "Ge_DCS_Inidividual: " << DCS_Inidividual << endl;
                    //double Bin_Width       = velocity_TH1F[Applied_Hist]->GetBinWidth(kkk)*1e-3;
                }
                cout << "velocity_TH1F_2[Applied_Hist]->GetBinCenter(Maximum_Bin_Ge): " << velocity_TH1F_2[Applied_Hist]->GetBinCenter(Maximum_Bin_Ge) << endl;
                cout << "Total_Cross_Section_Ge: " << Total_Cross_Section_Ge << endl;
                for(int LLL=Minimum_Bin_Ge; LLL<Maximum_Bin_Ge+1; LLL++)
                {
                    double DCS_Inidividual  = velocity_TH1F_2[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2);
                    double Recoil_Energy_Ge = velocity_TH1F_2[Applied_Hist]->GetBinCenter(LLL);//eV
                    Mean_Value_Loss_Ge = Mean_Value_Loss_Ge + Recoil_Energy_Ge*(DCS_Inidividual/Total_Cross_Section_Ge);

                }
                cout << "Mean_Ge: " << Mean_Value_Loss_Ge << endl;
                 cour << "Recoil_Energy_Xe[2]: " << Recoil_Energy_Xe[2] <<


            }
        }
        
    }
    
    //==================================================//

}//End_Main
 */
/*
 cout << "================================" << endl;
 cout << "Recoil_Energy_Xe: " << Recoil_Energy_Xe << endl;
 cout << "Mean_Value_Loss_Xe: " << Mean_Value_Loss_Xe << endl;
 cout << "Recoil_Energy_Xe*(DCS_Inidividual/Total_Cross_Section_Xe): " << Recoil_Energy_Xe*(DCS_Inidividual/Total_Cross_Section_Xe) << endl;
 cout << "DCS_Inidividual/Total_Cross_Section_Xe: " << DCS_Inidividual/Total_Cross_Section_Xe << endl;
 cout << "================================" << endl;
*/
 
/*
cout << "================================" << endl;
cout << "Recoil_Energy_Ge: " << Recoil_Energy_Ge << endl;
cout << "Mean_Value_Loss_Ge: " << Mean_Value_Loss_Ge << endl;
cout << "Recoil_Energy_Ge*DCS_Inidividual/Total_Cross_Section_Ge: " << Recoil_Energy_Ge*DCS_Inidividual/Total_Cross_Section_Ge << endl;
cout << "DCS_Inidividual/Total_Cross_Section_Ge: " << DCS_Inidividual/Total_Cross_Section_Ge << endl;
cout << "================================" << endl;
 */

/*
cout << "================================" << endl;
cout << "Recoil_Energy_Ge: " << Recoil_Energy_Ge << endl;
cout << "Mean_Value_Loss_Ge: " << Mean_Value_Loss_Ge << endl;
cout << "Recoil_Energy_Ge*DCS_Inidividual/Total_Cross_Section_Ge: " << Recoil_Energy_Ge*DCS_Inidividual/Total_Cross_Section_Ge << endl;
cout << "DCS_Inidividual/Total_Cross_Section_Ge: " << DCS_Inidividual/Total_Cross_Section_Ge << endl;
cout << "================================" << endl;
 */


const double MeVinverse2_to_GeVinverse2_Plus_UT=1e+4*1e-24;//UT(Unit transform)
const double Energy_H_Point[7]={80,120,150,200,300,400,480};//eV
const double DCS_H_Point[7]={1e-18,1e-19,1e-20,1e-21,1e-22,1e-23,1e-24};//1e-24cm^2/keV,1/MeV^2
//const double Energy_H_Point[11]={17,25,35,50,80,120,150,200,300,400,480};//eV
//const double DCS_H_Point[11]={1e-14,1e-15,1e-16,1e-17,1e-18,1e-19,1e-20,1e-21,1e-22,1e-23,1e-24};//1e-24cm^2/keV,1/MeV^2

/*
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()//Test the fitting
//Test the fitting line
{
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    
    const int Index=0;
    const double Density = 1.8;//g/cm^3
    const double Length  = 1e5;//cm
    //const double Length  = 3e0;//cm
    const double Threshold_keV = 0.1;//keV
    vector<string> File;
    vector<double> WIMP_mx_Array;
    vector<double> Collision_Time_Array;
    vector<double> Energy_Loss_Per_Collision;
    vector<double> Scaling={1e-18,1};

    //=========Define the file==========//
    //Xe_c1
    if(Index==0)
    {
    File=
            {"c1_2_0GeV","c1_1_0GeV",
            "c1_0_900GeV","c1_0_800GeV",
            "c1_0_700GeV","c1_0_600GeV",
            "c1_0_500GeV","c1_0_400GeV",
            "c1_0_300GeV","c1_0_200GeV",
            "c1_0_120GeV","c1_0_090GeV",
            "c1_0_070GeV","c1_0_050GeV"};
        WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};
    }
    //Xe_d1
    if(Index==1)
    {
         File=
            {"d1_2_0GeV","d1_1_0GeV",
            "d1_0_900GeV","d1_0_800GeV",
            "d1_0_700GeV","d1_0_600GeV",
            "d1_0_500GeV","d1_0_400GeV",
            "d1_0_300GeV","d1_0_200GeV",
            "d1_0_120GeV","d1_0_090GeV",
            "d1_0_070GeV","d1_0_050GeV"};
        WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};
    }
    //Read the file of DCS for different masses
    
    
    string c1_d1_Xe_Ge_index[4]={"SI_c1_XeData_Vel","SI_d1_Xeata_Vel","SI_c1_GeData_Vel","SI_d1_GeData_Vel"};
    double UT = 1e-24*1e-12;//Unit change
    double Energy_H_Check[7]={82,120,170,220,320,400,480};//eV
    double DCS_H_Check[7]=
        {
            (1E-2*UT),(1E-3*UT),
            (1E-4*UT),(1E-5*UT),
            (1E-6*UT),(1E-7*UT),
            (1E-8*UT)
        };//1e-24cm^2/keV,1/MeV^2

    
    //===============electron recoil===================//
    //for(int kkk=0; kkk<File.size(); kkk++)
    for(int kkk=1; kkk<2; kkk++)
    {
        
        TFile *fin   = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str(),File[kkk].c_str()));
        TFile *fin_2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index+2].c_str(),File[kkk].c_str()));

        TH1F*    velocity_TH1F[velocity_N];//Mass-based
        TH1F*    velocity_TH1F_2[velocity_N];//Mass-based
        for(int LLL=0; LLL<velocity_N; LLL++)
        {
            TH1F *Velocity_TH1F_temp  =(TH1F*)fin  ->Get(velocity_s[LLL].c_str());
            TH1F *Velocity_TH1F_temp_2=(TH1F*)fin_2->Get(velocity_s[LLL].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                cout << "File: " << velocity_s[LLL].c_str() << endl;
                velocity_TH1F[LLL]=(Velocity_TH1F_temp);
            }
            if(Velocity_TH1F_temp_2!=NULL)
            {
                cout << "File: " << velocity_s[LLL].c_str() << endl;
                velocity_TH1F_2[LLL]=(Velocity_TH1F_temp_2);
            }
        }
        
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        //for(int Applied_Hist=0; Applied_Hist<velocity_N; Applied_Hist++)
        for(int Applied_Hist=5; Applied_Hist<6; Applied_Hist++)
        {
            vector<double> TGraph_Recoil_Energy_Xe;vector<double> TGraph_Recoil_Energy_Ge;
            vector<double> TGraph_dsigma_dT_Xe          ;vector<double> TGraph_dsigma_dT_Ge;
            vector<double> TGraph_dsigma_dT_Xe_Check    ;vector<double> TGraph_dsigma_dT_Ge_Check;
            vector<double> TGraph_dsigma_dT_Xe_Ge_Check  ;vector<double> TGraph_dsigma_dT_Xe_H_Check  ;
            double Total_dsigma_dT_Xe=0;double Total_dsigma_dT_Ge=0;
            //Xe
            Int_t Minimum_Bin_Xe = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(80);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin_Xe = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Lowest_Y_Xe   = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Minimum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Highest_Y_Xe  = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Maximum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Length_Xe_X   = velocity_TH1F[Applied_Hist]->GetBinCenter(Maximum_Bin_Xe);//eV
            double Slope_Xe      = (Highest_Y_Xe - Lowest_Y_Xe)/(Length_Xe_X);
            //Ge
            Int_t Minimum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->GetXaxis()->FindBin(80);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Lowest_Y_Ge    = TMath::Log10(velocity_TH1F_2[Applied_Hist]->GetBinContent(Minimum_Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2));
            double Highest_Y_Ge   = TMath::Log10(velocity_TH1F_2[Applied_Hist]->GetBinContent(Maximum_Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2));
            double Length_Ge_X    = velocity_TH1F[Applied_Hist]->GetBinCenter(Maximum_Bin_Ge);//eV
            double Slope_Ge       = (Highest_Y_Ge - Lowest_Y_Ge)/(Length_Ge_X);

            double Length_H_X     = (480.-80.);
            double Slope_H        = (TMath::Log10(DCS_H_Point[6]) - TMath::Log10(DCS_H_Point[0]))/(Length_H_X);
            
            cout << "Lowest_Y_Xe : " << Lowest_Y_Xe  << endl;
            cout << "Lowest_Y_Ge : " << Lowest_Y_Ge  << endl;

            cout << "Slope_Xe: " << Slope_Xe << endl;
            cout << "Slope_Ge: " << Slope_Ge << endl;
            cout << "Slope_H: "  << Slope_H  << endl;
            
        }
        
    }
    
    //==================================================//

}//End_Main
*/


double Linear_Line(double x, double a, double b)
{
    return a*x+b;
}

void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()//Test the fitting
//Test the fitting line
{
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    string velocitykm_s[velocity_N];for(int kkk=0; kkk<velocity_N; kkk++){velocitykm_s[kkk]=std::to_string((int)velocitykm[kkk]);}//(km/s)}
    
    const int Index=0;
    const double Density = 1.8;//g/cm^3
    const double Length  = 1e5;//cm
    //const double Length  = 3e0;//cm
    const double Threshold_keV = 0.1;//keV
    vector<string> File;
    vector<double> WIMP_mx_Array;
    vector<double> Collision_Time_Array;
    vector<double> Energy_Loss_Per_Collision;
    vector<double> Scaling={1e-18,1e-9};

    const int    Mass_N=14;
    //=========Define the file==========//
    //Xe_c1
    if(Index==0)
    {
    File=
            {"c1_2_0GeV","c1_1_0GeV",
            "c1_0_900GeV","c1_0_800GeV",
            "c1_0_700GeV","c1_0_600GeV",
            "c1_0_500GeV","c1_0_400GeV",
            "c1_0_300GeV","c1_0_200GeV",
            "c1_0_120GeV","c1_0_090GeV",
            "c1_0_070GeV","c1_0_050GeV"};
        WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};
    }
    //Xe_d1
    if(Index==1)
    {
         File=
            {"d1_2_0GeV","d1_1_0GeV",
            "d1_0_900GeV","d1_0_800GeV",
            "d1_0_700GeV","d1_0_600GeV",
            "d1_0_500GeV","d1_0_400GeV",
            "d1_0_300GeV","d1_0_200GeV",
            "d1_0_120GeV","d1_0_090GeV",
            "d1_0_070GeV","d1_0_050GeV"};
        WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};
    }
    //Read the file of DCS for different masses
    string Mass_s[Mass_N]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P12","0P09","0P07","0P05"};

    
    string c1_d1_Xe_Ge_index[4]={"SI_c1_XeData_Vel","SI_d1_Xeata_Vel","SI_c1_GeData_Vel","SI_d1_GeData_Vel"};
    double UT = 1e-24*1e-12;//Unit change
    double Energy_H_Check[7]={82,120,170,220,320,400,480};//eV
    double DCS_H_Check[7]=
        {
            (1E-2*UT),(1E-3*UT),
            (1E-4*UT),(1E-5*UT),
            (1E-6*UT),(1E-7*UT),
            (1E-8*UT)
        };//1e-24cm^2/keV,1/MeV^2

    double Energy_H[12]={12,20,28,38,55,82,120,170,220,320,400,480};//eV
    double DCS_H[12]=
        {
            TMath::Log10(2.5E+2*UT),TMath::Log10(1E+2*UT),
            TMath::Log10(1E+1*UT),TMath::Log10(1E+0*UT),
            TMath::Log10(1E-1*UT),TMath::Log10(1E-2*UT),
            TMath::Log10(1E-3*UT),TMath::Log10(1E-4*UT),
            TMath::Log10(1E-5*UT),TMath::Log10(1E-6*UT),
            TMath::Log10(1E-7*UT),TMath::Log10(1E-8*UT)
        };//1e-24cm^2/keV,1/MeV^2

    //===============electron recoil===================//
    
    double Collision_Air_Part=0;
    for(int ATM_Number=0; ATM_Number<19; ATM_Number++)
    {
        double Length_Passed = (atm_table[ATM_Number+1][0]-atm_table[ATM_Number][0])*1e2;//cm
        double Density       = 1e-3*((atm_table[ATM_Number+1][4]+atm_table[ATM_Number][4])*0.5);//g/cm^3
        Collision_Air_Part   = Collision_Air_Part + Length_Passed*Density*1./(unified_atomic_mass_g*(15.)) ;
    }
    const double       NaI_Density            = 3.67;//3.67(g/cm^3)
    const double       NaI_Atomic_Mass        = 22.98*0.5+126*0.5;//
    const double       Pb_Density             = 11.29;//3.67(g/cm^3)
    const double       Pb_Atomic_Mass         = 207.2;//
    const double       Fixed_Length           = 20.;//cm
    
    const double Density_Array[3]           ={2.8 ,Pb_Density    ,NaI_Density};
    const double Atomic_Mass_Array[3]={Weighted_Atomic_Number_Cement,Pb_Atomic_Mass,NaI_Atomic_Mass};
    const double Number_Density_Array[3]    ={
                                              Density_Array[0]*1./(unified_atomic_mass_g*(Atomic_Mass_Array[0])) ,
                                              Density_Array[1]*1./(unified_atomic_mass_g*(Atomic_Mass_Array[1])) ,
                                              Density_Array[1]*1./(unified_atomic_mass_g*(Atomic_Mass_Array[2])) ,
                                            };
    
    double   Total_Cross_Section_Mass_dependent_HH[Mass_N];
    double   Total_Cross_Section_Mass_dependent_Xe[Mass_N];
    double   Total_Cross_Section_Mass_dependent_Ge[Mass_N];
    double   Ratio_Xe[Mass_N][velocity_N-1];
    double   Ratio_Ge[Mass_N][velocity_N-1];
    double   Ratio_averaged[Mass_N][velocity_N-1];
    
    double  Total_Cross_Section_Ge_array[velocity_N];
    double  Total_Cross_Section_HH_at_1GeV_300kms        = 1.87548e-36;
    double  Total_Cross_Section_times_velocity_square    = Total_Cross_Section_HH_at_1GeV_300kms*300.*300.;
    
    //double A = Total_Cross_Section_HH_at_1GeV_300kms*Constant_scaling[1][5];
    
    for(int kkk=0; kkk<File.size()-1; kkk++)
    //for(int kkk=10; kkk<11; kkk++)
    //for(int kkk=7; kkk<8; kkk++)
    //for(int kkk=8; kkk<9; kkk++)//0.12GeV
    {
        
        TFile *fin   = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str()  ,File[kkk].c_str()));
        TFile *fin_2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index+2].c_str(),File[kkk].c_str()));
        TH1F*    velocity_TH1F[velocity_N];//Mass-based
        TH1F*    velocity_TH1F_2[velocity_N];//Mass-based
        int      Check_Coherent[velocity_N];
        
        for(int LLL=0; LLL<velocity_N; LLL++)
        {
            cout << "velocity_s[LLL].c_str(): " << velocity_s[LLL].c_str() << endl;
            TH1F *Velocity_TH1F_temp  =(TH1F*)fin  ->Get(velocity_s[LLL].c_str());
            TH1F *Velocity_TH1F_temp_2=(TH1F*)fin_2->Get(velocity_s[LLL].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                cout << "File_1: " << velocity_s[LLL].c_str() << endl;
                velocity_TH1F[LLL]=(Velocity_TH1F_temp);
            }
            if(Velocity_TH1F_temp_2!=NULL)
            {
                cout << "File_2: " << velocity_s[LLL].c_str() << endl;
                velocity_TH1F_2[LLL]=(Velocity_TH1F_temp_2);
            }
            cout << "=================" << endl;
            if(Velocity_TH1F_temp!=NULL and Velocity_TH1F_temp_2!=NULL)
            {
                cout << "LLL_1: " << LLL << endl;
                Check_Coherent[LLL]=1;
            }
            else
            {
                cout << "LLL_0: " << LLL << endl;
                Check_Coherent[LLL]=0;
            }
            cout << "=================" << endl;
        }
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        //for(int Applied_Hist=0; Applied_Hist<velocity_N; Applied_Hist++)
        double Power_Final_array[velocity_N];
        double Ratio1_array[velocity_N];
        double Ratio2_array[velocity_N];
        
        double Total_Cross_Section_Xe_array[velocity_N];
        double Total_Cross_Section_Ge_array[velocity_N];
        
        double Energy_Loss_Xe_array[velocity_N];
        double Energy_Loss_Ge_array[velocity_N];
        //Total_Cross_Section_Xe_array[0]=0;Total_Cross_Section_Xe_array[1]=0;
       // for(int Applied_Hist=0; Applied_Hist<velocity_N; Applied_Hist++)
        for(int Applied_Hist=0; Applied_Hist<velocity_N; Applied_Hist++)
        //for(int Applied_Hist=0; Applied_Hist<6; Applied_Hist++)
        //for(int Applied_Hist=3; Applied_Hist<velocity_N; Applied_Hist++)
        {
            Energy_Loss_Xe_array[Applied_Hist]=0;
            Energy_Loss_Ge_array[Applied_Hist]=0;
            cout << "================================================" << endl;
            if(Check_Coherent[Applied_Hist]!=1)
            {
                Total_Cross_Section_Xe_array[Applied_Hist]=0;
                Total_Cross_Section_Ge_array[Applied_Hist]=0;
                continue;
            }
            cout << "Applied_Hist: " << Applied_Hist << endl;
            cout << "Check_Coherent[Applied_Hist]: " << Check_Coherent[Applied_Hist] << endl;
            
            //Find the relation between the total cross sections and the velocity
            cout << "velocitykm: " << velocitykm[Applied_Hist] << endl;
            cout << "velocity_d[Applied_Hist]" << velocity_d[Applied_Hist] << endl;
            Int_t Minimum_Bin_Xe = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin_Xe = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Xe_Y1 = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Minimum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Xe_Y2 = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Maximum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Xe_X1 = velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(Minimum_Bin_Xe);
            double Xe_X2 = velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(Maximum_Bin_Xe);
            double A_Slope_Xe = (Xe_Y2-Xe_Y1)/(Xe_X2-Xe_X1);
            cout << "A_Slope_Xe: " << A_Slope_Xe << endl;
            double B_Xe       =  Xe_Y1 - Xe_X1*A_Slope_Xe;
            
            Int_t Maximum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Ge_Y1 = Xe_Y1;
            double Ge_Y2 = TMath::Log10(velocity_TH1F_2[Applied_Hist]->GetBinContent(Maximum_Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2));
            double Ge_X1 = Xe_X1;
            double Ge_X2 = velocity_TH1F_2[Applied_Hist]->GetXaxis()->GetBinCenter(Maximum_Bin_Ge);
            double A_Slope_Ge = (Ge_Y2-Ge_Y1)/(Ge_X2-Ge_X1);
            cout << "A_Slope_Ge: " << A_Slope_Ge << endl;
            double B_Ge       =  Ge_Y1 - Ge_X1*A_Slope_Ge;

            //cout << "TMath::Power(10, Linear_Line(12.,A_Slope_HH,B_HH)): " << TMath::Power(10, Linear_Line(12.,A_Slope_HH,B_HH)) << endl;
            double Total_Xe_DCS=0;double Total_Ge_DCS=0;double Total_HH_DCS=0;
            double Bin = (480.-12.)/10000.;
            for(int KKK=0; KKK<10000; KKK++)
            {
                double DCS_Xe     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Xe,B_Xe));
                double DCS_Ge     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Ge,B_Ge));
                Total_Xe_DCS = Total_Xe_DCS + DCS_Xe;
                Total_Ge_DCS = Total_Ge_DCS + DCS_Ge;
            }
            double Atomic_Number_Array[3]={1.00784,28.,54.};
            double Atomic_Mass_Array[3]  ={1.00784,72.64,131.293};

            Total_Cross_Section_Xe_array[Applied_Hist]=Total_Xe_DCS;
            Total_Cross_Section_Ge_array[Applied_Hist]=Total_Ge_DCS;
            cout << "velocitykm: " << velocitykm[Applied_Hist] << endl;
            cout << "Total_Cross_Section_Xe_array[Applied_Hist]: " << Total_Cross_Section_Xe_array[Applied_Hist] << endl;
            cout << "Total_Cross_Section_Ge_array[Applied_Hist]: " << Total_Cross_Section_Ge_array[Applied_Hist] << endl;
            //cout << "Total_Cross_Section_Xe_array[Applied_Hist]/A^3: " << Total_Cross_Section_Xe_array[Applied_Hist]/(131.293*131.293*131.293) << endl;
            
            Total_Cross_Section_Mass_dependent_Xe[kkk]=Total_Xe_DCS;
            Total_Cross_Section_Mass_dependent_Ge[kkk]=Total_Ge_DCS;
             
            
            //Calculating the collision Time and mean value of energy loss
            /*
            vector<double> A_array               = {28.0855,15.99};
            vector<double> Z_Array               = {14.,8.};
            double Total_Cross_Section_const     = 8e-39;
            //First case is for CDEX
            double CDEX_Threshold            = 160;//eV
            double Collision_Time            = 1.4e3*1e2*( 2.33*1./(unified_atomic_mass_g*(A_array[0])) )*Total_Cross_Section_const*TMath::Power(A_array[0],3);
            double Energy_Loss_per_collision = 10*sqrt(Z_Array[0]);
            double Total_Energy              = Collision_Time*Energy_Loss_per_collision;
            cout << "Total_Energy: " << Total_Energy << endl;
            double Scaling                   = CDEX_Threshold/(Total_Energy);
            cout << "CS_Try: " << CS_Try((Scaling),1) << endl;
            */
            /*
            Int_t Minimum_Bin_Xe = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin_Xe = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Xe_Y1 = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Minimum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Xe_Y2 = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Maximum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Xe_X1 = velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(Minimum_Bin_Xe);
            double Xe_X2 = velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(Maximum_Bin_Xe);
            double A_Slope_Xe = (Xe_Y2-Xe_Y1)/(Xe_X2-Xe_X1);
            cout << "A_Slope_Xe: " << A_Slope_Xe << endl;
            double B_Xe       =  Xe_Y1 - Xe_X1*A_Slope_Xe;
            cout << "B_Xe: " << B_Xe << endl;

            Int_t Maximum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Ge_Y1 = Xe_Y1;
            double Ge_Y2 = TMath::Log10(velocity_TH1F_2[Applied_Hist]->GetBinContent(Maximum_Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2));
            double Ge_X1 = Xe_X1;
            double Ge_X2 = velocity_TH1F_2[Applied_Hist]->GetXaxis()->GetBinCenter(Maximum_Bin_Ge);
            double A_Slope_Ge = (Ge_Y2-Ge_Y1)/(Ge_X2-Ge_X1);
            cout << "A_Slope_Ge: " << A_Slope_Ge << endl;
            double B_Ge       =  Ge_Y1 - Ge_X1*A_Slope_Ge;
            cout << "B_Ge: " << B_Ge << endl;
            //cout << "TMath::Power(10, Linear_Line(12.,A_Slope_HH,B_HH)): " << TMath::Power(10, Linear_Line(12.,A_Slope_HH,B_HH)) << endl;
            double Total_Xe_DCS=0;double Total_Ge_DCS=0;
            double Bin = (480.-12.)/10000;
            
            for(int KKK=0; KKK<10000; KKK++)
            {
                double DCS_Xe     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Xe,B_Xe));
                double DCS_Ge     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Ge,B_Ge));
                Total_Xe_DCS = Total_Xe_DCS + DCS_Xe;
                Total_Ge_DCS = Total_Ge_DCS + DCS_Ge;
            }
            cout << "Total_Xe_DCS: " << Total_Xe_DCS << endl;
            cout << "Total_Ge_DCS: " << Total_Ge_DCS << endl;
            double Energy_Loss_Xe_DCS=0;double Energy_Loss_Ge_DCS=0;double CCC=0;
            for(int KKK=0; KKK<10000; KKK++)
            {
                double DCS_Xe     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Xe,B_Xe));
                double DCS_Ge     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Ge,B_Ge));
                
                CCC = CCC + (12.+Bin*(KKK+1.))*DCS_Xe/Total_Xe_DCS;
                Energy_Loss_Ge_DCS = Energy_Loss_Ge_DCS + (12.+Bin*(KKK+1.))*DCS_Ge/Total_Ge_DCS;
            }
            cout << "BBB: " << CCC << endl;
            cout << "Energy_Loss_Xe_DCS: " << CCC << endl;
            cout << "Energy_Loss_Ge_DCS: " << Energy_Loss_Ge_DCS << endl;

            Energy_Loss_Xe_array[Applied_Hist]=CCC;
            Energy_Loss_Ge_array[Applied_Hist]=Energy_Loss_Ge_DCS;
             */
             //Find the mean value of energy loss and its relationship with Z
            /*
            cout << "velocitykm: " << velocitykm[Applied_Hist] << endl;
            cout << "velocity_d[Applied_Hist]" << velocity_d[Applied_Hist] << endl;
            Int_t Minimum_Bin_Xe = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin_Xe = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Xe_Y1 = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Minimum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Xe_Y2 = TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(Maximum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            double Xe_X1 = velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(Minimum_Bin_Xe);
            double Xe_X2 = velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(Maximum_Bin_Xe);
            double A_Slope_Xe = (Xe_Y2-Xe_Y1)/(Xe_X2-Xe_X1);
            cout << "A_Slope_Xe: " << A_Slope_Xe << endl;
            double B_Xe       =  Xe_Y1 - Xe_X1*A_Slope_Xe;
            
            Int_t Maximum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Ge_Y1 = Xe_Y1;
            double Ge_Y2 = TMath::Log10(velocity_TH1F_2[Applied_Hist]->GetBinContent(Maximum_Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2));
            double Ge_X1 = Xe_X1;
            double Ge_X2 = velocity_TH1F_2[Applied_Hist]->GetXaxis()->GetBinCenter(Maximum_Bin_Ge);
            double A_Slope_Ge = (Ge_Y2-Ge_Y1)/(Ge_X2-Ge_X1);
            cout << "A_Slope_Ge: " << A_Slope_Ge << endl;
            double B_Ge       =  Ge_Y1 - Ge_X1*A_Slope_Ge;

            double HH_Y1 = TMath::Log10(DCS_H_Check[0]);
            double HH_Y2 = TMath::Log10(DCS_H_Check[6]);
            double HH_X1 = Energy_H_Check[0];
            double HH_X2 = Energy_H_Check[6];
            double A_Slope_HH = (HH_Y2-HH_Y1)/(HH_X2-HH_X1);
            cout << "A_Slope_HH: " << A_Slope_HH << endl;
            double B_HH       =  HH_Y1 - HH_X1*A_Slope_HH;

            //cout << "TMath::Power(10, Linear_Line(12.,A_Slope_HH,B_HH)): " << TMath::Power(10, Linear_Line(12.,A_Slope_HH,B_HH)) << endl;
            double Total_Xe_DCS;double Total_Ge_DCS;double Total_HH_DCS;
            double Bin = (480.-12.)/10000;
            for(int KKK=0; KKK<10000; KKK++)
            {
                double DCS_Xe     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Xe,B_Xe));
                double DCS_Ge     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Ge,B_Ge));
                Total_Xe_DCS = Total_Xe_DCS + DCS_Xe;
                Total_Ge_DCS = Total_Ge_DCS + DCS_Ge;
            }
            for(int KKK=0; KKK<12; KKK++)
            {
                double DCS_HH     =  TMath::Power(10,DCS_H[KKK]);
                Total_HH_DCS = Total_HH_DCS + DCS_HH;
            }
            double Energy_Loss_Xe_DCS;double Energy_Loss_Ge_DCS;double Energy_Loss_HH_DCS;
            double CCC;
            for(int KKK=0; KKK<10000; KKK++)
            {
                double DCS_Xe     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Xe,B_Xe));
                double DCS_Ge     =  TMath::Power(10, Linear_Line(12.+Bin*(KKK+1.),A_Slope_Ge,B_Ge));
                
                CCC = CCC + (12.+Bin*(KKK+1.))*DCS_Xe/Total_Xe_DCS;
                //cout << "(12.+Bin*(KKK+1.))*DCS_Xe/Total_Xe_DCS: " << (12.+Bin*(KKK+1.))*DCS_Xe/Total_Xe_DCS << endl;
                //cout << "Energy_Loss_Xe_DCS: " << CCC << endl;
                
                Energy_Loss_Ge_DCS = Energy_Loss_Ge_DCS + (12.+Bin*(KKK+1.))*DCS_Ge/Total_Ge_DCS;
            }
            cout << "BBB: " << CCC << endl;
            for(int KKK=0; KKK<12; KKK++)
            {
                double DCS_HH     =  TMath::Power(10,DCS_H[KKK]);
                Energy_Loss_HH_DCS = Energy_Loss_HH_DCS + Energy_H[KKK]*DCS_HH/Total_HH_DCS;
            }
            cout << "Energy_Loss_Xe_DCS: " << CCC << endl;
            cout << "Energy_Loss_Ge_DCS: " << Energy_Loss_Ge_DCS << endl;
            cout << "Energy_Loss_HH_DCS: " << Energy_Loss_HH_DCS << endl;

            double Test_variable_Xe=54.;
            double Test_variable_Ge=28.;
            double Test_variable_HH=1.;

            double Power_Final=0;
            double Diff_1;double Diff_2;double Ratio_1;double Ratio_2;
            double Start_Point_Diff=100.;
            for(int KKK = 0; KKK< 1000; KKK++)
            {
                double Power_Constant = 0.001*(double)KKK;
                double Xe_Energy_Loss_const =CCC/TMath::Power(Test_variable_Xe,Power_Constant);
                double Ge_Energy_Loss_const =Energy_Loss_Ge_DCS/TMath::Power(Test_variable_Ge,Power_Constant);
                double HH_Energy_Loss_const =Energy_Loss_HH_DCS/TMath::Power(Test_variable_HH,Power_Constant);
                //cout << "Power_Constant: " << Power_Constant << endl;
                //cout << "Xe_Energy_Loss_const: " << Xe_Energy_Loss_const << endl;
                //cout << "Ge_Energy_Loss_const: " << Ge_Energy_Loss_const << endl;
                //cout << "HH_Energy_Loss_const: " << HH_Energy_Loss_const << endl;
                Diff_1 = abs(Xe_Energy_Loss_const-HH_Energy_Loss_const);
                Diff_2 = abs(Ge_Energy_Loss_const-HH_Energy_Loss_const);
                if( Diff_1+Diff_2<Start_Point_Diff )
                {
                    Power_Final = Power_Constant;
                    Ratio_1 = Xe_Energy_Loss_const/HH_Energy_Loss_const;
                    Ratio_2 = Ge_Energy_Loss_const/HH_Energy_Loss_const;
                    Start_Point_Diff = Diff_1+Diff_2;
                }
            }
            Power_Final_array[Applied_Hist] = Power_Final;
            Ratio1_array[Applied_Hist]      = Ratio_1;
            Ratio2_array[Applied_Hist]      = Ratio_2;
            */
            
            //double Test_variable_Xe=131.;
            //double Test_variable_Ge=72.;
            //double Test_variable_HH=1.;
             /*
            double Test_variable_Xe=54.;
            double Test_variable_Ge=28.;
            double Test_variable_HH=1.;

            double Atomic_Mass_Array[3]={1.00784,28.,54.};
            //double Atomic_Mass_Array[3]={1.00784,72.64,131.293};
            double Energy_Loss_DCS_0P5Power[3];double Energy_Loss_DCS_1Power[3];double Energy_Loss_DCS_2Power[3];
            
            Energy_Loss_DCS_0P5Power[0]=Energy_Loss_Xe_DCS/TMath::Power(Test_variable_Xe,0.5);
            Energy_Loss_DCS_0P5Power[1]=Energy_Loss_Ge_DCS/TMath::Power(Test_variable_Ge,0.5);
            Energy_Loss_DCS_0P5Power[2]=Energy_Loss_HH_DCS/TMath::Power(Test_variable_HH,0.5);

            Energy_Loss_DCS_1Power[0]  =Energy_Loss_Xe_DCS/TMath::Power(Test_variable_Xe,1);
            Energy_Loss_DCS_1Power[1]  =Energy_Loss_Ge_DCS/TMath::Power(Test_variable_Ge,1);
            Energy_Loss_DCS_1Power[2]  =Energy_Loss_HH_DCS/TMath::Power(Test_variable_HH,1);

            Energy_Loss_DCS_2Power[0]  =Energy_Loss_Xe_DCS/TMath::Power(Test_variable_Xe,2);
            Energy_Loss_DCS_2Power[1]  =Energy_Loss_Ge_DCS/TMath::Power(Test_variable_Ge,2);
            Energy_Loss_DCS_2Power[2]  =Energy_Loss_HH_DCS/TMath::Power(Test_variable_HH,2);

            TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
            gStyle->SetOptStat(0);
            gStyle->SetTitleSize(0.04,"XY");
            gStyle->SetTitleFont(62,"XY");
            gStyle->SetLegendFont(62);
            
            TGraph *ELoss_over_0P5Power_of_Z = new TGraph(3,Atomic_Mass_Array,Energy_Loss_DCS_0P5Power);//Total Cross Section(TCS), Electron Number(EN)
            ELoss_over_0P5Power_of_Z->SetMarkerStyle(20);
            ELoss_over_0P5Power_of_Z->SetMarkerColor(2);
            ELoss_over_0P5Power_of_Z->SetLineColor(2);
            ELoss_over_0P5Power_of_Z->SetLineWidth(5);
            ELoss_over_0P5Power_of_Z->GetXaxis()->SetTitle("Atomic Number(Z)");
            ELoss_over_0P5Power_of_Z->GetXaxis()->SetRangeUser(0,500);
            ELoss_over_0P5Power_of_Z->GetYaxis()->SetRangeUser(1e-2,15);
            ELoss_over_0P5Power_of_Z->Draw("apl");

            TGraph *ELoss_over_1Power_of_Z = new TGraph(3,Atomic_Mass_Array,Energy_Loss_DCS_1Power);//Total Cross Section(TCS), Electron Number(EN)
            ELoss_over_1Power_of_Z->SetMarkerStyle(20);
            ELoss_over_1Power_of_Z->SetMarkerColor(3);
            ELoss_over_1Power_of_Z->SetLineColor(3);
            ELoss_over_1Power_of_Z->SetLineWidth(5);
            ELoss_over_1Power_of_Z->Draw("plsame");
        
            TGraph *ELoss_over_2Power_of_Z = new TGraph(3,Atomic_Mass_Array,Energy_Loss_DCS_2Power);//Total Cross Section(TCS), Electron Number(EN)
            ELoss_over_2Power_of_Z->SetMarkerStyle(20);
            ELoss_over_2Power_of_Z->SetMarkerColor(4);
            ELoss_over_2Power_of_Z->SetLineColor(4);
            ELoss_over_2Power_of_Z->SetLineWidth(5);
            ELoss_over_2Power_of_Z->Draw("plsame");

            TLegend *leg = new TLegend(0.7,0.1,0.9,0.3);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            leg->AddEntry(ELoss_over_0P5Power_of_Z,"<E>/Z^{#frac{1}{2}}","lP");
            leg->AddEntry(ELoss_over_1Power_of_Z,"<E>/Z^{1}","lP");
            leg->AddEntry(ELoss_over_2Power_of_Z,"<E>/Z^{2}","lP");
            leg->Draw();

            c1->SetLogy();
            c1->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/ELoss_over_Power_of_Z_%s.pdf",velocitykm_s[Applied_Hist].c_str()));
            */
            
            
             //Get the histogram of DCS and the mean value of energy loss
            /*
            vector<double> TGraph_Recoil_Energy_Xe;vector<double> TGraph_Recoil_Energy_Ge;
            vector<double> TGraph_dsigma_dT_Xe          ;vector<double> TGraph_dsigma_dT_Ge;
            vector<double> TGraph_dsigma_dT_Xe_Check    ;vector<double> TGraph_dsigma_dT_Ge_Check;
            vector<double> TGraph_dsigma_dT_Xe_Ge_Check  ;vector<double> TGraph_dsigma_dT_Xe_H_Check  ;
            double Total_dsigma_dT_Xe=0;double Total_dsigma_dT_Ge=0;
            //Xe
            Int_t Minimum_Bin_Xe = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin_Xe = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Lowest_Y_Xe    = velocity_TH1F[Applied_Hist]->GetBinContent(Minimum_Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2);
            for(int LLL=Minimum_Bin_Xe; LLL<Maximum_Bin_Xe+1; LLL++)
            {
                Total_dsigma_dT_Xe = Total_dsigma_dT_Xe + velocity_TH1F[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2) ;
            }
            double Avergae_Total_Xe = 0;
            for(int LLL=Minimum_Bin_Xe; LLL<Maximum_Bin_Xe+1; LLL++)
            {
                double ER  = velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(LLL);
                double DCS = velocity_TH1F[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2);
                Avergae_Total_Xe = Avergae_Total_Xe + ER*DCS/Total_dsigma_dT_Xe;
                TGraph_Recoil_Energy_Xe.push_back( ER );
                //TGraph_dsigma_dT_Xe.push_back( TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2)/Total_dsigma_dT_Xe)  );
                TGraph_dsigma_dT_Xe.push_back( TMath::Log10(DCS)  );

            }
            for(int MMM=0; MMM<7; MMM++)
            {
                Int_t Bin_Xe = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Energy_H_Point[MMM]);
                TGraph_dsigma_dT_Xe_Check.push_back(velocity_TH1F[Applied_Hist]->GetBinContent(Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2));
            }
            
            //Ge
            
            Int_t Minimum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->GetXaxis()->FindBin(80);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin_Ge = velocity_TH1F_2[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            double Lowest_Y_Ge    = velocity_TH1F_2[Applied_Hist]->GetBinContent(Minimum_Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2);
            
            double Ratio_Used    = Lowest_Y_Xe/Lowest_Y_Ge;
            cout << "Ratio_Used: " <<  Ratio_Used << endl;
            for(int LLL=Minimum_Bin_Ge; LLL<Maximum_Bin_Ge+1; LLL++)
            {
                Total_dsigma_dT_Ge = Total_dsigma_dT_Ge + velocity_TH1F_2[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2) ;
            }
            double Avergae_Total_Ge = 0;
            for(int LLL=Minimum_Bin_Ge; LLL<Maximum_Bin_Ge+1; LLL++)
            {
                double ER  = velocity_TH1F_2[Applied_Hist]->GetXaxis()->GetBinCenter(LLL);
                double DCS = velocity_TH1F_2[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2);
                Avergae_Total_Ge = Avergae_Total_Ge + ER*DCS/Total_dsigma_dT_Ge;
                TGraph_Recoil_Energy_Ge.push_back( ER );
                TGraph_dsigma_dT_Ge.push_back( TMath::Log10(DCS)  );
            }
            for(int MMM=0; MMM<7; MMM++)
            {
                Int_t Bin_Ge = velocity_TH1F_2[Applied_Hist]->GetXaxis()->FindBin(Energy_H_Point[MMM]);
                cout << "Bin_Ge: " << Bin_Ge << endl;
                TGraph_dsigma_dT_Ge_Check.push_back(velocity_TH1F_2[Applied_Hist]->GetBinContent(Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2));
            }
            for(int MMM=0; MMM<7; MMM++)
            {
                TGraph_dsigma_dT_Xe_Ge_Check.push_back(TMath::Log10(TGraph_dsigma_dT_Xe_Check[MMM]/TGraph_dsigma_dT_Ge_Check[MMM]));
            }
            
            for(int MMM=0; MMM<7; MMM++)
            {
                cout << "TGraph_dsigma_dT_Xe_Check[MMM]: " << TGraph_dsigma_dT_Xe_Check[MMM] << endl;
                cout << "DCS_H[MMM]: " << DCS_H[MMM] << endl;
                cout << "TGraph_dsigma_dT_Xe_Check[MMM]/DCS_H[MMM]: " << TGraph_dsigma_dT_Xe_Check[MMM]/DCS_H[MMM] << endl;
                cout << "TMath::Log10(TGraph_dsigma_dT_Xe_Check[MMM]/DCS_H[MMM]): " << TMath::Log10(TGraph_dsigma_dT_Xe_Check[MMM]/DCS_H[MMM]) << endl;
                TGraph_dsigma_dT_Xe_H_Check.push_back(TMath::Log10(TGraph_dsigma_dT_Xe_Check[MMM]/DCS_H[MMM]));
            }

            for(int MMM=0; MMM<7; MMM++)
            {
                cout << "DCS_H[MMM]: " << DCS_H[MMM] << endl;
                cout << "TGraph_dsigma_dT_Xe_Ge_Check: " << TGraph_dsigma_dT_Xe_Ge_Check[MMM] << endl;
                cout << "TGraph_dsigma_dT_Xe_H_Check: " << TGraph_dsigma_dT_Xe_H_Check[MMM] << endl;
            }
            
            cout << "Total_dsigma_dT_Xe: " << Total_dsigma_dT_Xe << endl;
            cout << "Total_dsigma_dT_Ge "  << Total_dsigma_dT_Ge << endl;
            */
            /*
            TGraph * g_Xe = new TGraph((int)TGraph_Recoil_Energy_Xe.size(), &TGraph_Recoil_Energy_Xe[0], &TGraph_dsigma_dT_Xe[0]);
            g_Xe->SetMarkerStyle(20);
            g_Xe->SetMarkerColor(2);
            g_Xe->SetMarkerColor(2);
            g_Xe->GetXaxis()->SetRangeUser(0,500);
            g_Xe->GetYaxis()->SetRangeUser(-45,-30);
            g_Xe->Draw("apl");

            TGraph * g_Ge = new TGraph((int)TGraph_Recoil_Energy_Ge.size(), &TGraph_Recoil_Energy_Ge[0], &TGraph_dsigma_dT_Ge[0]);
            g_Ge->SetMarkerStyle(20);
            g_Ge->SetMarkerColor(3);
            g_Ge->Draw("plsame");

            TGraph * g_H = new TGraph(12, Energy_H, DCS_H);
            g_H->SetMarkerStyle(20);
            g_H->SetMarkerColor(4);
            g_H->SetLineWidth(5);
            g_H->Draw("plsame");

            TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            leg->AddEntry(g_Xe,"Xe","lP");
            leg->AddEntry(g_Ge,"Ge","lP");
            leg->AddEntry(g_H,"H","lP");
            leg->Draw();
             */
            /*
            TGraph * g_Xe_Ge = new TGraph(7, Energy_H_Point, &TGraph_dsigma_dT_Xe_Ge_Check[0]);
            g_Xe_Ge->SetMarkerStyle(20);
            g_Xe_Ge->SetMarkerColor(2);
            g_Xe_Ge->SetMarkerColor(2);
            g_Xe_Ge->GetYaxis()->SetRangeUser(0,20);
            g_Xe_Ge->GetXaxis()->SetTitle("Recoil Energy(eV)");
            g_Xe_Ge->GetYaxis()->SetTitle("Log10(Ratio of DCS)");
            g_Xe_Ge->Draw("apl");

            TGraph * g_Xe_H = new TGraph(7, Energy_H_Point, &TGraph_dsigma_dT_Xe_H_Check[0]);
            g_Xe_H->SetMarkerStyle(20);
            g_Xe_H->SetMarkerColor(3);
            g_Xe_H->Draw("plsame");
             */
                    
            
            //Predict the total cross sections of Silicon and oxygen
            /*
            TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
            gStyle->SetOptStat(0);
            gStyle->SetTitleSize(0.04,"XY");
            gStyle->SetTitleFont(62,"XY");
            gStyle->SetLegendFont(62);

            double Total_Cross_Section_Xe_Check=0;
            double Total_Cross_Section_Ge_Check=0;
            double Total_Cross_Section_HH_Check=0;

            double Atomic_Mass_Array[3]={1.00784,72.64,131.293};
            double TCS_over_EN_Element[3];double TCS_over_EN_square_Element[3];double TCS_over_EN_third_Element[3];double TCS_over_EN_fourth_Element[3];
            double TCS_over_AM_Element[3];double TCS_over_AM_square_Element[3];double TCS_over_AM_third_Element[3];double TCS_over_AM_fourth_Element[3];//Atomic_Mass(AM)
            for(int Energy_Index=0; Energy_Index<7; Energy_Index++)
            {
                Int_t Bin_Xe = velocity_TH1F[Applied_Hist]    ->GetXaxis()->FindBin(Energy_H_Check[Energy_Index]);//Min_Recoil_Bin_Xe
                Int_t Bin_Ge = velocity_TH1F_2[Applied_Hist]  ->GetXaxis()->FindBin(Energy_H_Check[Energy_Index]);//Min_Recoil_Bin_Xe
                
                Total_Cross_Section_Xe_Check = Total_Cross_Section_Xe_Check + velocity_TH1F[Applied_Hist]  ->GetBinContent(Bin_Xe)*1e-15*TMath::Power(Scaling[Index],2);
                Total_Cross_Section_Ge_Check = Total_Cross_Section_Ge_Check + velocity_TH1F_2[Applied_Hist]->GetBinContent(Bin_Ge)*1e-15*TMath::Power(Scaling[Index],2);
                Total_Cross_Section_HH_Check = Total_Cross_Section_HH_Check + DCS_H_Check[Energy_Index];
            }
             
            cout << "Total_Cross_Section_Xe_Check: " << Total_Cross_Section_Xe_Check << endl;
            cout << "Total_Cross_Section_Ge_Check: " << Total_Cross_Section_Ge_Check << endl;
            cout << "Total_Cross_Section_HH_Check: " << Total_Cross_Section_HH_Check << endl;
            
            
            TCS_over_EN_Element[0]=Total_Cross_Section_HH_Check/1.;
            TCS_over_EN_Element[1]=Total_Cross_Section_Ge_Check/28.;
            TCS_over_EN_Element[2]=Total_Cross_Section_Xe_Check/54.;

            TCS_over_EN_square_Element[0]=Total_Cross_Section_HH_Check/(1.*1.);
            TCS_over_EN_square_Element[1]=Total_Cross_Section_Ge_Check/(28.*28.);
            TCS_over_EN_square_Element[2]=Total_Cross_Section_Xe_Check/(54.*54.);

            TCS_over_EN_third_Element[0]=Total_Cross_Section_HH_Check/(1.*1.*1.);
            TCS_over_EN_third_Element[1]=Total_Cross_Section_Ge_Check/(28.*28.*28.);
            TCS_over_EN_third_Element[2]=Total_Cross_Section_Xe_Check/(54.*54.*54.);

            TCS_over_EN_fourth_Element[0]=Total_Cross_Section_HH_Check/(1.*1.*1.*1.);
            TCS_over_EN_fourth_Element[1]=Total_Cross_Section_Ge_Check/(28.*28.*28.*28.);
            TCS_over_EN_fourth_Element[2]=Total_Cross_Section_Xe_Check/(54.*54.*54.*54.);

            TCS_over_AM_Element[0]=Total_Cross_Section_HH_Check/(1.00784);
            TCS_over_AM_Element[1]=Total_Cross_Section_Ge_Check/(72.64);
            TCS_over_AM_Element[2]=Total_Cross_Section_Xe_Check/(131.293);

            TCS_over_AM_square_Element[0]=Total_Cross_Section_HH_Check/(1.00784*1.00784);
            TCS_over_AM_square_Element[1]=Total_Cross_Section_Ge_Check/(72.64*72.64);
            TCS_over_AM_square_Element[2]=Total_Cross_Section_Xe_Check/(131.293*131.293);

            TCS_over_AM_third_Element[0]=Total_Cross_Section_HH_Check/(1.00784*1.00784*1.00784);
            TCS_over_AM_third_Element[1]=Total_Cross_Section_Ge_Check/(72.64*72.64*72.64);
            TCS_over_AM_third_Element[2]=Total_Cross_Section_Xe_Check/(131.293*131.293*131.293);
             
            cout << "TCS_over_AM_third_Element[0]: " << TCS_over_AM_third_Element[0] << endl;
            cout << "TCS_over_AM_third_Element[1]: " << TCS_over_AM_third_Element[1] << endl;
            cout << "TCS_over_AM_third_Element[2]: " << TCS_over_AM_third_Element[2] << endl;
            */
            /*
            TGraph *TCS_over_AM_third = new TGraph(3,Atomic_Mass_Array,TCS_over_AM_third_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_AM_third->SetMarkerStyle(20);
            TCS_over_AM_third->SetMarkerColor(4);
            TCS_over_AM_third->SetLineWidth(5);
            TCS_over_AM_third->GetXaxis()->SetTitle("Atomic Mass");
            TCS_over_AM_third->GetYaxis()->SetTitle("Total Cross Section/Atomic Mass");
            TCS_over_AM_third->Draw("apl");
            
            c1->SetLogy();
            c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/TCS_over_AM_third.pdf");
             */
            /*
            TGraph *TCS_over_AM_square = new TGraph(3,Atomic_Mass_Array,TCS_over_AM_square_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_AM_square->SetMarkerStyle(20);
            TCS_over_AM_square->SetMarkerColor(4);
            TCS_over_AM_square->SetLineWidth(5);
            TCS_over_AM_square->GetXaxis()->SetTitle("Atomic Mass");
            TCS_over_AM_square->GetYaxis()->SetTitle("Total Cross Section/Atomic Mass");
            TCS_over_AM_square->Draw("apl");
            
            c1->SetLogy();
            c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/TCS_over_AM_square.pdf");
             */
            
            /*
            TGraph *TCS_over_AM = new TGraph(3,Atomic_Mass_Array,TCS_over_AM_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_AM->SetMarkerStyle(20);
            TCS_over_AM->SetMarkerColor(4);
            TCS_over_AM->SetLineWidth(5);
            TCS_over_AM->GetXaxis()->SetTitle("Atomic Mass");
            TCS_over_AM->GetYaxis()->SetTitle("Total Cross Section/Atomic Mass");
            TCS_over_AM->Draw("apl");
            
            c1->SetLogy();
            c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/TCS_over_AM.pdf");
             */
            /*
            TGraph *TCS_over_EN_fourth = new TGraph(3,Atomic_Mass_Array,TCS_over_EN_fourth_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_EN_fourth->SetMarkerStyle(20);
            TCS_over_EN_fourth->SetMarkerColor(4);
            TCS_over_EN_fourth->SetLineWidth(5);
            TCS_over_EN_fourth->GetXaxis()->SetTitle("Atomic Mass");
            TCS_over_EN_fourth->GetYaxis()->SetTitle("Total Cross Section/electron number^4");
            TCS_over_EN_fourth->Draw("apl");
            
            c1->SetLogy();
            c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/TCS_over_EN_fourth.pdf");
             */
            
            /*
            TGraph *TCS_over_AM = new TGraph(3,Atomic_Mass_Array,TCS_over_AM_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_AM->SetMarkerStyle(20);
            TCS_over_AM->SetMarkerColor(2);
            TCS_over_AM->SetLineColor(2);
            TCS_over_AM->SetLineWidth(5);
            TCS_over_AM->GetXaxis()->SetTitle("Atomic Mass(AM)");
            TCS_over_AM->GetXaxis()->SetRangeUser(0,500);
            TCS_over_AM->GetYaxis()->SetRangeUser(1E-40,1E-33);
            TCS_over_AM->Draw("apl");

            TGraph *TCS_over_AM_square = new TGraph(3,Atomic_Mass_Array,TCS_over_AM_square_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_AM_square->SetMarkerStyle(20);
            TCS_over_AM_square->SetMarkerColor(3);
            TCS_over_AM_square->SetLineColor(3);
            TCS_over_AM_square->SetLineWidth(5);
            TCS_over_AM_square->Draw("plsame");
        
            TGraph *TCS_over_AM_third = new TGraph(3,Atomic_Mass_Array,TCS_over_AM_third_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_AM_third->SetMarkerStyle(20);
            TCS_over_AM_third->SetMarkerColor(4);
            TCS_over_AM_third->SetLineColor(4);
            TCS_over_AM_third->SetLineWidth(5);
            TCS_over_AM_third->Draw("plsame");


            TLegend *leg = new TLegend(0.1,0.7,0.4,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            leg->AddEntry(TCS_over_AM,"#sigma_{total}/AM","lP");
            leg->AddEntry(TCS_over_AM_square,"#sigma_{total}/(AM^2)","lP");
            leg->AddEntry(TCS_over_AM_third,"#sigma_{total}/(AM^3)","lP");
            leg->Draw();

            c1->SetLogy();
            c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/All_TCS_over_AM.pdf");
             */
            /*
            TGraph *TCS_over_EN = new TGraph(3,Atomic_Mass_Array,TCS_over_EN_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_EN->SetMarkerStyle(20);
            TCS_over_EN->SetMarkerColor(2);
            TCS_over_EN->SetLineWidth(5);
            TCS_over_EN->SetLineColor(2);
            TCS_over_EN->GetXaxis()->SetTitle("Atomic Mass");
            TCS_over_EN->GetXaxis()->SetRangeUser(0,500);
            TCS_over_EN->GetYaxis()->SetRangeUser(1E-40,1E-33);
            TCS_over_EN->Draw("alp");
            
            TGraph *TCS_over_EN_square = new TGraph(3,Atomic_Mass_Array,TCS_over_EN_square_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_EN_square->SetMarkerStyle(20);
            TCS_over_EN_square->SetMarkerColor(3);
            TCS_over_EN_square->SetLineColor(3);
            TCS_over_EN_square->SetLineWidth(5);
            TCS_over_EN_square->Draw("lpsame");

            TGraph *TCS_over_EN_third = new TGraph(3,Atomic_Mass_Array,TCS_over_EN_third_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_EN_third->SetMarkerStyle(20);
            TCS_over_EN_third->SetMarkerColor(4);
            TCS_over_EN_third->SetLineColor(4);
            TCS_over_EN_third->SetLineWidth(5);
            TCS_over_EN_third->Draw("lpsame");

            TGraph *TCS_over_EN_fourth = new TGraph(3,Atomic_Mass_Array,TCS_over_EN_fourth_Element);//Total Cross Section(TCS), Electron Number(EN)
            TCS_over_EN_fourth->SetMarkerStyle(20);
            TCS_over_EN_fourth->SetMarkerColor(6);
            TCS_over_EN_fourth->SetLineColor(6);
            TCS_over_EN_fourth->SetLineWidth(5);
            TCS_over_EN_fourth->Draw("lpsame");

            TLegend *leg = new TLegend(0.1,0.7,0.4,0.9);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);
            leg->AddEntry(TCS_over_EN,"#sigma_{total}/EN","lP");
            leg->AddEntry(TCS_over_EN_square,"#sigma_{total}/(EN^2)","lP");
            leg->AddEntry(TCS_over_EN_third,"#sigma_{total}/(EN^3)","lP");
            leg->AddEntry(TCS_over_EN_fourth,"#sigma_{total}/(EN^4)","lP");
            leg->Draw();

            c1->SetLogy();
            c1->Print("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/All_TCS_over_EN.pdf");
             */
            
            
        }//for(int Applied_Hist=5; Applied_Hist<6; Applied_Hist++)
        /*
        double Energy_Loss_divide_V_square_Xe[velocity_N];double Energy_Loss_divide_V_square_Ge[velocity_N];
        for(int VN=0; VN<velocity_N; VN++)
        {
            cout << "velocity: " << velocitykm[VN] << endl;
            cout << "Energy_Loss_Xe_array[Applied_Hist]: " << Energy_Loss_Xe_array[VN] << endl;
            cout << "Energy_Loss_Ge_array[Applied_Hist]: " << Energy_Loss_Ge_array[VN] << endl;
            cout << "Energy_Loss_Xe_array[Applied_Hist]/velocity: " << Energy_Loss_Xe_array[VN]/velocitykm[VN] << endl;
            cout << "Energy_Loss_Ge_array[Applied_Hist]/velocity: " << Energy_Loss_Ge_array[VN]/velocitykm[VN] << endl;
            cout << "Energy_Loss_Xe_array[Applied_Hist]/velocity^2 " << Energy_Loss_Xe_array[VN]/(velocitykm[VN]*velocitykm[VN]) << endl;
            cout << "Energy_Loss_Ge_array[Applied_Hist]/velocity^2 " << Energy_Loss_Ge_array[VN]/(velocitykm[VN]*velocitykm[VN]) << endl;
            Energy_Loss_divide_V_square_Xe[VN]=Energy_Loss_Xe_array[VN];
            Energy_Loss_divide_V_square_Ge[VN]=Energy_Loss_Ge_array[VN];
        }
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.03,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);

        TGraph * TG_Energy_Loss_divide_V_square_Xe = new TGraph(velocity_N,velocitykm,Energy_Loss_divide_V_square_Xe);
        TG_Energy_Loss_divide_V_square_Xe->SetMarkerStyle(20);
        TG_Energy_Loss_divide_V_square_Xe->SetMarkerColor(2);
        TG_Energy_Loss_divide_V_square_Xe->SetMarkerColor(2);
        TG_Energy_Loss_divide_V_square_Xe->GetYaxis()->SetRangeUser(1,1e4);
        TG_Energy_Loss_divide_V_square_Xe->GetXaxis()->SetTitle("V(km/s)");
        TG_Energy_Loss_divide_V_square_Xe->GetYaxis()->SetTitle("<E>");
        TG_Energy_Loss_divide_V_square_Xe->Draw("apl");

        TGraph * TG_Energy_Loss_divide_V_square_Ge = new TGraph(velocity_N,velocitykm,Energy_Loss_divide_V_square_Ge);
        TG_Energy_Loss_divide_V_square_Ge->SetMarkerStyle(20);
        TG_Energy_Loss_divide_V_square_Ge->SetMarkerColor(3);
        TG_Energy_Loss_divide_V_square_Ge->SetMarkerColor(3);
        TG_Energy_Loss_divide_V_square_Ge->Draw("plsame");

        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(TG_Energy_Loss_divide_V_square_Xe,"Xe","lP");
        leg->AddEntry(TG_Energy_Loss_divide_V_square_Ge,"Ge","lP");
        leg->Draw();

        c1->SetLogy();
        c1->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/Energy_Loss_%sGeV.pdf",Mass_s[kkk].c_str()));
        */
        /*
        double Energy_Loss_divide_V_square_Xe[velocity_N];double Energy_Loss_divide_V_square_Ge[velocity_N];
        for(int VN=0; VN<velocity_N; VN++)
        {
            cout << "velocity: " << velocitykm[VN] << endl;
            cout << "Energy_Loss_Xe_array[Applied_Hist]: " << Energy_Loss_Xe_array[VN] << endl;
            cout << "Energy_Loss_Ge_array[Applied_Hist]: " << Energy_Loss_Ge_array[VN] << endl;
            cout << "Energy_Loss_Xe_array[Applied_Hist]/velocity: " << Energy_Loss_Xe_array[VN]/velocitykm[VN] << endl;
            cout << "Energy_Loss_Ge_array[Applied_Hist]/velocity: " << Energy_Loss_Ge_array[VN]/velocitykm[VN] << endl;
            cout << "Energy_Loss_Xe_array[Applied_Hist]/velocity^2 " << Energy_Loss_Xe_array[VN]/(velocitykm[VN]*velocitykm[VN]) << endl;
            cout << "Energy_Loss_Ge_array[Applied_Hist]/velocity^2 " << Energy_Loss_Ge_array[VN]/(velocitykm[VN]*velocitykm[VN]) << endl;
            Energy_Loss_divide_V_square_Xe[VN]=Energy_Loss_Xe_array[VN]/(velocitykm[VN]*velocitykm[VN]);
            Energy_Loss_divide_V_square_Ge[VN]=Energy_Loss_Ge_array[VN]/(velocitykm[VN]*velocitykm[VN]);
        }
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.03,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);

        TGraph * TG_Energy_Loss_divide_V_square_Xe = new TGraph(velocity_N,velocitykm,Energy_Loss_divide_V_square_Xe);
        TG_Energy_Loss_divide_V_square_Xe->SetMarkerStyle(20);
        TG_Energy_Loss_divide_V_square_Xe->SetMarkerColor(2);
        TG_Energy_Loss_divide_V_square_Xe->SetMarkerColor(2);
        TG_Energy_Loss_divide_V_square_Xe->GetYaxis()->SetRangeUser(1e-6,1e3);
        TG_Energy_Loss_divide_V_square_Xe->GetXaxis()->SetTitle("V(km/s)");
        TG_Energy_Loss_divide_V_square_Xe->GetYaxis()->SetTitle("<E>/V^2");
        TG_Energy_Loss_divide_V_square_Xe->Draw("apl");

        TGraph * TG_Energy_Loss_divide_V_square_Ge = new TGraph(velocity_N,velocitykm,Energy_Loss_divide_V_square_Ge);
        TG_Energy_Loss_divide_V_square_Ge->SetMarkerStyle(20);
        TG_Energy_Loss_divide_V_square_Ge->SetMarkerColor(3);
        TG_Energy_Loss_divide_V_square_Ge->SetMarkerColor(3);
        TG_Energy_Loss_divide_V_square_Ge->Draw("plsame");

        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(TG_Energy_Loss_divide_V_square_Xe,"Xe","lP");
        leg->AddEntry(TG_Energy_Loss_divide_V_square_Ge,"Ge","lP");
        leg->Draw();

        c1->SetLogy();
        c1->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/Energy_Loss_divide_V_square_%sGeV.pdf",Mass_s[kkk].c_str()));
         */
        /*
        double Total_Cross_Section_Xe_Times_V_square[velocity_N];
        double Total_Cross_Section_Ge_Times_V_square[velocity_N];
        double Total_Cross_Section_Xe_Times_V[velocity_N];
        double Total_Cross_Section_Ge_Times_V[velocity_N];
        double Total_Cross_Section_HH_Times_V_square[velocity_N];
        double VReverse[velocity_N];
        for(int VN=0; VN<velocity_N; VN++)
        {
            //cout << "velocitykm: " << velocitykm[VN] << endl;
            //cout << "1./(velocitykm[VN]*velocitykm[VN]): " << 1./(velocitykm[VN]*velocitykm[VN]) << endl;
            VReverse[VN] = 1./(velocitykm[VN]*velocitykm[VN]);
            //cout << "VReverse[VN]: " << VReverse[VN]  << endl;
            //cout << "velocity_s[VN]: " << velocity_s[VN].c_str() << endl;
            //cout << "Total_Cross_Section_Xe_array: " << Total_Cross_Section_Xe_array[VN] << endl;
            //cout << "Total_Cross_Section_Ge_array: " << Total_Cross_Section_Ge_array[VN] << endl;
            //cout << "Total_Cross_Section_Xe_array[VN]/A^3: " << Total_Cross_Section_Xe_array[VN]/(131.293*131.293*131.293) << endl;
            //cout << "Total_Cross_Section_Ge_array[VN]/A^3: " << Total_Cross_Section_Ge_array[VN]/(72.64*72.64*72.64) << endl;
            cout << "Total_Cross_Section_Xe_scaled: " << Total_Cross_Section_Xe_array[VN]*(velocitykm[VN]) << endl;
            cout << "Total_Cross_Section_Ge_scaled: " << Total_Cross_Section_Ge_array[VN]*(velocitykm[VN]) << endl;
            Total_Cross_Section_Xe_Times_V_square[VN]=Total_Cross_Section_Xe_array[VN]*(velocitykm[VN]*velocitykm[VN]);
            Total_Cross_Section_Ge_Times_V_square[VN]=Total_Cross_Section_Ge_array[VN]*(velocitykm[VN]*velocitykm[VN]);
            Total_Cross_Section_Xe_Times_V[VN]=Total_Cross_Section_Xe_array[VN]*(velocitykm[VN]);
            Total_Cross_Section_Ge_Times_V[VN]=Total_Cross_Section_Ge_array[VN]*(velocitykm[VN]);
        }
        for(int VN=0; VN<velocity_N-1; VN++)
        {
            cout << "=============================" << endl;
            cout << "Total_Cross_Section_Xe_Times_V_square[VN+1]/Total_Cross_Section_Xe_Times_V_square[VN]: " << Total_Cross_Section_Xe_Times_V_square[VN+1]/Total_Cross_Section_Xe_Times_V_square[VN] << endl;
            cout << "Total_Cross_Section_Ge_Times_V_square[VN+1]/Total_Cross_Section_Ge_Times_V_square[VN]: " << Total_Cross_Section_Ge_Times_V_square[VN+1]/Total_Cross_Section_Ge_Times_V_square[VN] << endl;
            cout << "=============================" << endl;
            Ratio_Xe[kkk][VN] = Total_Cross_Section_Xe_Times_V_square[VN+1]/Total_Cross_Section_Xe_Times_V_square[VN];
            Ratio_Ge[kkk][VN] = Total_Cross_Section_Ge_Times_V_square[VN+1]/Total_Cross_Section_Ge_Times_V_square[VN];
        }
         */
         
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.03,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);
        TGraph * TG_Total_Cross_Section_Xe_Times_V_square = new TGraph(velocity_N,velocitykm,Total_Cross_Section_Xe_array);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerStyle(20);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerColor(2);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerColor(2);
        TG_Total_Cross_Section_Xe_Times_V_square->GetYaxis()->SetRangeUser(1e-35,1e-25);
        TG_Total_Cross_Section_Xe_Times_V_square->GetXaxis()->SetTitle("V(km/s)");
        TG_Total_Cross_Section_Xe_Times_V_square->GetYaxis()->SetTitle("#sigma_{SI}");
        TG_Total_Cross_Section_Xe_Times_V_square->Draw("apl");

        TGraph * TG_Total_Cross_Section_Ge_Times_V_square = new TGraph(velocity_N,velocitykm,Total_Cross_Section_Ge_array);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerStyle(20);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerColor(3);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerColor(3);
        TG_Total_Cross_Section_Ge_Times_V_square->Draw("plsame");

        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(TG_Total_Cross_Section_Xe_Times_V_square,"Xe","lP");
        leg->AddEntry(TG_Total_Cross_Section_Ge_Times_V_square,"Ge","lP");
        leg->Draw();

        c1->SetLogy();
        c1->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/Total_Cross_Section_%sGeV.pdf",Mass_s[kkk].c_str()));
          
        /*
        TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
        gStyle->SetOptStat(0);
        gStyle->SetTitleSize(0.03,"XY");
        gStyle->SetTitleFont(62,"XY");
        gStyle->SetLegendFont(62);
        TGraph * TG_Total_Cross_Section_Xe_Times_V_square = new TGraph(velocity_N,velocitykm,Total_Cross_Section_Xe_Times_V);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerStyle(20);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerColor(2);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerColor(2);
        TG_Total_Cross_Section_Xe_Times_V_square->GetYaxis()->SetRangeUser(1e-30,1e-25);
        TG_Total_Cross_Section_Xe_Times_V_square->GetXaxis()->SetTitle("V(km/s)");
        TG_Total_Cross_Section_Xe_Times_V_square->GetYaxis()->SetTitle("#sigma_{SI} #times V(km/s)");
        TG_Total_Cross_Section_Xe_Times_V_square->Draw("apl");

        TGraph * TG_Total_Cross_Section_Ge_Times_V_square = new TGraph(velocity_N,velocitykm,Total_Cross_Section_Ge_Times_V);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerStyle(20);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerColor(3);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerColor(3);
        TG_Total_Cross_Section_Ge_Times_V_square->Draw("plsame");

        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(TG_Total_Cross_Section_Xe_Times_V_square,"Xe","lP");
        leg->AddEntry(TG_Total_Cross_Section_Ge_Times_V_square,"Ge","lP");
        leg->Draw();

        c1->SetLogy();
        c1->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/Total_Cross_Section_Times_V_%sGeV.pdf",Mass_s[kkk].c_str()));
         */
        /*
        TGraph * TG_Total_Cross_Section_Xe_Times_V_square = new TGraph(velocity_N,velocitykm,Total_Cross_Section_Xe_Times_V_square);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerStyle(20);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerColor(2);
        TG_Total_Cross_Section_Xe_Times_V_square->SetMarkerColor(2);
        TG_Total_Cross_Section_Xe_Times_V_square->GetYaxis()->SetRangeUser(1e-27,1e-22);
        TG_Total_Cross_Section_Xe_Times_V_square->GetXaxis()->SetTitle("V(km/s)");
        TG_Total_Cross_Section_Xe_Times_V_square->GetYaxis()->SetTitle("#sigma_{SI} #times V(km/s)");
        TG_Total_Cross_Section_Xe_Times_V_square->Draw("apl");

        TGraph * TG_Total_Cross_Section_Ge_Times_V_square = new TGraph(velocity_N,velocitykm,Total_Cross_Section_Ge_Times_V_square);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerStyle(20);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerColor(3);
        TG_Total_Cross_Section_Ge_Times_V_square->SetMarkerColor(3);
        TG_Total_Cross_Section_Ge_Times_V_square->Draw("plsame");

        TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry(TG_Total_Cross_Section_Xe_Times_V_square,"Xe","lP");
        leg->AddEntry(TG_Total_Cross_Section_Ge_Times_V_square,"Ge","lP");
        leg->Draw();

        c1->SetLogy();
        c1->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Predict_DCS_ER/Total_Cross_Section_Times_V_%sGeV.pdf",Mass_s[kkk].c_str()));
         */
    }//for(int kkk=1; kkk<2; kkk++)
    /*
    for(int Plot_1=0; Plot_1<Mass_N; Plot_1++)
    {
        cout << "WIMP_mx_Array: " << WIMP_mx_Array[Plot_1] << endl;
        cout << "Total_Cross_Section_Mass_dependent_Xe[Plot_1]: " << Total_Cross_Section_Mass_dependent_Xe[Plot_1] << endl;
        cout << "Total_Cross_Section_Mass_dependent_Ge[Plot_1]: " << Total_Cross_Section_Mass_dependent_Ge[Plot_1] << endl;
        cout << "Total_Cross_Section_Mass_dependent_Xe[Plot_1]/Mx: " << Total_Cross_Section_Mass_dependent_Xe[Plot_1]/WIMP_mx_Array[Plot_1] << endl;
        cout << "Total_Cross_Section_Mass_dependent_Ge[Plot_1]/Mx: " << Total_Cross_Section_Mass_dependent_Ge[Plot_1]/WIMP_mx_Array[Plot_1] << endl;
        cout << "Total_Cross_Section_Mass_dependent_Xe[Plot_1]/(Mx*Mx): " << Total_Cross_Section_Mass_dependent_Xe[Plot_1]/(WIMP_mx_Array[Plot_1]*WIMP_mx_Array[Plot_1]) << endl;
        cout << "Total_Cross_Section_Mass_dependent_Ge[Plot_1]/(Mx*Mx): " << Total_Cross_Section_Mass_dependent_Ge[Plot_1]/(WIMP_mx_Array[Plot_1]*WIMP_mx_Array[Plot_1]) << endl;
    }
     */
    /*
    for(int kkk=0; kkk<Mass_N; kkk++)
    {
        cout << "WIMP_mx_Array: " << WIMP_mx_Array[kkk] << endl;
        for(int VN=0; VN<velocity_N-1; VN++)
        {
            double Averaged_value=0;
            cout << "Ratio_Xe[kkk][VN]: " << Ratio_Xe[kkk][VN] << endl;
            cout << "Ratio_Ge[kkk][VN]: " << Ratio_Ge[kkk][VN] << endl;
            Averaged_value = (Ratio_Xe[kkk][VN]+Ratio_Ge[kkk][VN])*0.5;
            if(Averaged_value>0 and Averaged_value<100)Ratio_averaged[kkk][VN] = Averaged_value;
            else{Ratio_averaged[kkk][VN] = 0.;}
        }
    }
    for(int kkk=0; kkk<Mass_N; kkk++)
    {
        cout << "WIMP_mx_Array: " << WIMP_mx_Array[kkk] << endl;
        cout << "" << endl;
        for(int VN=0; VN<velocity_N-1; VN++)
        {
            cout << Ratio_averaged[kkk][VN] << ",";
        }
        cout << "" << endl;
    }
     */
    //==================================================//

}//End_Main
 

/*
TGraph * g_Xe_Ge = new TGraph(7, Energy_H_Point, &TGraph_dsigma_dT_Xe_Ge_Check[0]);
g_Xe_Ge->SetMarkerStyle(20);
g_Xe_Ge->SetMarkerColor(2);
g_Xe_Ge->SetMarkerColor(2);
g_Xe_Ge->GetYaxis()->SetRangeUser(0,20);
g_Xe_Ge->GetXaxis()->SetTitle("Recoil Energy(eV)");
g_Xe_Ge->GetXaxis()->SetTitle("Log10(Ratio of DCS)");
g_Xe_Ge->Draw("apl");

TGraph * g_Xe_H = new TGraph(7, Energy_H_Point, &TGraph_dsigma_dT_Xe_H_Check[0]);
g_Xe_H->SetMarkerStyle(20);
g_Xe_H->SetMarkerColor(3);
g_Xe_H->Draw("plsame");
 */

/*
 TGraph * g_Xe = new TGraph((int)TGraph_Recoil_Energy_Xe.size(), &TGraph_Recoil_Energy_Xe[0], &TGraph_dsigma_dT_Xe[0]);
 g_Xe->SetMarkerStyle(20);
 g_Xe->SetMarkerColor(2);
 g_Xe->SetMarkerColor(2);
 g_Xe->GetXaxis()->SetRangeUser(0,500);
 g_Xe->GetYaxis()->SetRangeUser(-60,-30);
 g_Xe->Draw("apl");

 TGraph * g_Ge = new TGraph((int)TGraph_Recoil_Energy_Ge.size(), &TGraph_Recoil_Energy_Ge[0], &TGraph_dsigma_dT_Ge[0]);
 g_Ge->SetMarkerStyle(20);
 g_Ge->SetMarkerColor(3);
 g_Ge->Draw("plsame");

 TGraph * g_H = new TGraph(11, Energy_H_Point, DCS_H);
 g_H->SetMarkerStyle(20);
 g_H->SetMarkerColor(4);
 g_H->SetLineWidth(5);
 g_H->Draw("plsame");

 TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
 leg->SetFillColor(0);
 leg->SetFillStyle(0);
 leg->SetTextSize(0.04);
 leg->SetBorderSize(0);
 leg->SetTextFont(22);
 leg->AddEntry(g_Xe,"Xe","lP");
 leg->AddEntry(g_Ge,"Ge","lP");
 leg->AddEntry(g_H,"H","lP");
 leg->Draw();
  
*/
 


/*
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()
{
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    
    const int Index      = 1;
     //CDEX
    const double Density = 1.8;//g/cm^3
    const double Length  = 1e5;//cm
     
     //TEXONO
    
    const double Density = 1;//g/cm^3
    const double Length  = 3e3;//cm
     

    //const double Length  = 3e0;//cm
    const double Threshold_keV = 0.01;//keV
    vector<string> File;
    vector<double> WIMP_mx_Array;
    vector<double> Scaling={1e-18,1e-9};
    //=========Define the file==========//
    //Xe_c1
    if(Index==0)
    {
    File=
            {"c1_20_0GeV",
            "c1_10_0GeV","c1_5_0GeV",
            "c1_4_0GeV","c1_3_0GeV",
            "c1_2_0GeV","c1_1_0GeV",
            "c1_0_900GeV","c1_0_800GeV",
            "c1_0_700GeV","c1_0_600GeV",
            "c1_0_500GeV","c1_0_400GeV",
            "c1_0_300GeV","c1_0_200GeV",
            "c1_0_120GeV","c1_0_090GeV",
            "c1_0_070GeV","c1_0_050GeV",
            "c1_0_040GeV","c1_0_030GeV",
            "c1_0_020GeV","c1_0_010GeV",
            };
        WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05,0.04,0.03,0.02,0.01};
    }
    //Xe_d1
    if(Index==1)
    {
         File=
            {
            "d1_20_0GeV",
            "d1_10_0GeV","d1_5_0GeV",
            "d1_4_0GeV","d1_3_0GeV",
            "d1_2_0GeV","d1_1_0GeV",
            "d1_0_900GeV","d1_0_800GeV",
            "d1_0_700GeV","d1_0_600GeV",
            "d1_0_500GeV","d1_0_400GeV",
            "d1_0_300GeV","d1_0_200GeV",
            "d1_0_120GeV","d1_0_090GeV",
            "d1_0_080GeV","d1_0_070GeV",
            "d1_0_050GeV",
            "d1_0_040GeV","d1_0_030GeV",
            "d1_0_020GeV","d1_0_010GeV",
            };
        WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.08,0.07,0.05,0.04,0.03,0.02,0.01};
    }
    //Read the file of DCS for different masses
    double Cross_Section_Set[File.size()];
    double Collision_Time_Array[File.size()];
    double Energy_Loss_Per_Collision[File.size()];
    double Energy_Loss_Per_Count_Only[File.size()];
    //int File_Start=File.size()-3;
    int File_Start=0;

    string c1_d1_Xe_Ge_index[4]={"SI_c1_XeData_Vel","SI_d1_XeData_Vel","SI_c1_GeData_Vel","SI_d1_GeData_Vel"};
    //===============electron recoil===================//
    //for(int FIle_Index=File_Start; FIle_Index<File.size(); FIle_Index++)
    for(int FIle_Index=6; FIle_Index<7; FIle_Index++)
    {
        cout << "File[FIle_Index].c_str(): " << File[FIle_Index].c_str() << endl;
        TFile *fin = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str(),File[FIle_Index].c_str()));
        vector<TH1F*> velocity_TH1F;//Mass-based
        vector<double>   velocity_Used;//km/s
        velocity_Used.push_back(0);
        for(int LLL=0; LLL<velocity_N; LLL++)
        {
            TH1F *Velocity_TH1F_temp=(TH1F*)fin->Get(velocity_s[LLL].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                //cout << "File: " << velocity_s[kkk].c_str() << endl;
                //cout << "velocitykm: " << velocitykm[kkk] << endl;
                velocity_TH1F.push_back(Velocity_TH1F_temp);
                velocity_Used.push_back(velocitykm[LLL]);
            }
        }
        vector<double> Fitting_a;vector<double> Fitting_b;
        for(int Applied_Hist=0; Applied_Hist<velocity_TH1F.size(); Applied_Hist++)
        //for(int Applied_Hist=0; Applied_Hist<3; Applied_Hist++)
        {
            vector<double> TGraph_Recoil_Energy;vector<double> TGraph_dsigma_dT;
            
            Int_t Minimum_Bin = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            for(int LLL=Minimum_Bin; LLL<Maximum_Bin+1; LLL++)
            {
                TGraph_Recoil_Energy.push_back( velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(LLL) );//keV recoil energy
                TGraph_dsigma_dT.push_back( TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2)) );//Real dsigma_dT
            }
            TGraph * g = new TGraph((int)TGraph_Recoil_Energy.size(), &TGraph_Recoil_Energy[0], &TGraph_dsigma_dT[0]);
            TF1 *fitting_Line = new TF1("fitting_Line","[0]*x+[1]",12,velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter( Maximum_Bin-1 ));//Fitting
            g->Fit(fitting_Line,"R");
            Double_t par[9];
            fitting_Line->GetParameters(&par[0]);
            Fitting_a.push_back(par[0]);
            Fitting_b.push_back(par[1]);
            cout << "par[0]: " << par[0] << endl;
            cout << "par[1]: " << par[1] << endl;
        }
        vector<TF1*> TF1_Fitting_Line;


        //=======================================Start interacting=======================================//
        
        int Event_Number=0;
        double collision_Time_Total=0;
        double Energy_Loss_eV_collision=0;
        int    Applied_Hist=0;

        while(Event_Number<20)
        {
            //cout << "=====Eve====: " << Event_Number <<endl; ;
            double Velocity_DM_Temp = 779;//Initial Energy and iteration =>km/s
            double Energy_DM_Temp   = DM_E(WIMP_mx_Array[FIle_Index],779);//Initial Energy and iteration =>keV
            int collision_Time=0;
            double Energy_Loss_Inidivual=0;
            //for(int kkk=0; kkk<1; kkk++)
            //cout << "=====================================================================" << endl;
            while(Energy_DM_Temp>0.01)
            {
                //Confirm the TH1F
                gRandom = new TRandom3(0);
                gRandom->SetSeed(0);
                
                //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;
                //cout << "Energy_DM_Temp: " << Energy_DM_Temp << endl;
                //Find the file
                Applied_Hist=0;
                
                if(Velocity_DM_Temp>velocity_Used[velocity_Used.size()-1])
                {
                    Applied_Hist=(velocity_TH1F.size())-1;
                }
                if(Applied_Hist==0)
                {
                    for(int kkk=0; kkk<velocity_Used.size(); kkk++)
                    {
                        if(Velocity_DM_Temp>velocity_Used[kkk] and Velocity_DM_Temp<velocity_Used[kkk+1])
                        {
                            Applied_Hist=kkk;
                        }
                    }
                }
                //Set the range
                double Energy_Loss_eV = Energy_Loss(Fitting_a[Applied_Hist],Fitting_b[Applied_Hist],Energy_DM_Temp*1e3,1);
                double Energy_Loss_keV = Energy_Loss_eV*1e-3;//Energy_Loss(keV)
                //cout << "---------------------------" << endl;
                //cout << "Energy_DM_Temp(eV): " << Energy_DM_Temp*1e3 << endl;

                //cout << "Energy_Loss(eV): " << Energy_Loss_keV*1e3 << endl;
                Energy_DM_Temp   = Energy_DM_Temp - Energy_Loss_keV;
                Velocity_DM_Temp = sqrt(2*Energy_DM_Temp/(WIMP_mx_Array[FIle_Index]*1e6))*(3e8/1e3);//km/s
                //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;//km/s
                Energy_Loss_Inidivual = Energy_Loss_Inidivual + Energy_Loss_keV;
                collision_Time = collision_Time + 1;
                Applied_Hist =0;
            }
            //cout << "collision_Time: " << collision_Time << endl;
            //if(collision_Time>0)cout << "Energy_Loss_Inidivual/collision_Time: " << Energy_Loss_Inidivual/collision_Time << endl;
            Energy_Loss_eV_collision = Energy_Loss_eV_collision + Energy_Loss_Inidivual/collision_Time;
            collision_Time_Total = collision_Time_Total + collision_Time;
            Event_Number = Event_Number + 1;
        }
        double collision_Time_Ave = collision_Time_Total / Event_Number;
        Collision_Time_Array[FIle_Index] = collision_Time_Ave;//Count for every DM
        //cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        double Energy_Loss_Real = Energy_DM(WIMP_mx_Array[FIle_Index],779*(1e3/3e8))-Threshold_keV;
        //cout << "Energy_Loss_Real: " << Energy_Loss_Real << endl;
        //cout << "Energy_Loss_eV_collision/Event_Number: " << Energy_Loss_eV_collision/Event_Number << endl;
        Energy_Loss_Per_Count_Only[FIle_Index] = Energy_Loss_eV_collision/Event_Number;
        Energy_Loss_eV_collision = (1./Energy_Loss_Real) * (Energy_Loss_eV_collision/Event_Number) ;
        Energy_Loss_Per_Collision[FIle_Index] = (Energy_Loss_eV_collision);
        //cout << "Threshold: " << sqrt(2*0.16/(WIMP_mx_Array[kkk]*1e6))*(3e8/1e3);
        //cout << "DM_E(WIMP_mx_Array[kkk],779) : " << DM_E(WIMP_mx_Array[kkk],779) << endl;;
        //cout << "Energy_Loss_eV_collision: " << Energy_Loss_eV_collision << endl;
        //==================================================//
        
        //double Target_Total_Cross_Section = collision_Time_Ave*(unified_atomic_mass_g*(ASi))/(Density*Length);
        Applied_Hist=0;
        double Velocity_DM_Now=779;//km/s
        if(Velocity_DM_Now>velocity_Used[velocity_Used.size()-1])Applied_Hist=(velocity_TH1F.size())-1;
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        if(Applied_Hist==0)
        {
            for(int kkk=0; kkk<velocity_Used.size(); kkk++)
            {
                if(Velocity_DM_Now>velocity_Used[kkk] and Velocity_DM_Now<velocity_Used[kkk+1]){Applied_Hist=kkk;}
            }
        }
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        double v_c = Velocity_DM_Now*1e3/3e8;//c
        double Max_Recoil = DM_E(WIMP_mx_Array[FIle_Index],Velocity_DM_Now);//Max_Recoil(keV) given by mass
        Int_t Maximum_Bin = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Max_Recoil*1e3);//Max_Recoil_Bin
        //cout << "Maximum_Bin: " << Maximum_Bin << endl;
        int LastBin = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
        //cout << "LastBin: " << LastBin << endl;

        
        double Total_Count; double Area_Integral;
        double X_Range     =  Max_Recoil*1e3;//eV
        double X_Range_Bin = (Max_Recoil*1e3)*2e-3;//eV
        double T           =  0;
        for(int kkk=1; kkk<499; kkk++)
        {
            //cout << kkk << endl;
            double T_Temp    = T_Temp + X_Range_Bin;
            double dT        = (T_Temp - T)*1e-3;//keV
            //cout << "Max_Recoil: " << Max_Recoil << endl;
            //cout << "dT: " << dT << endl;
            double T_Central = (T + T_Temp)*0.5;
            //cout << "T_Central: " << T_Central << endl;
            TF1 *fitting_Line = new TF1("fitting_Line","TMath::Power(10,([0]*x+[1]))",0,X_Range);
            //cout << "X_Range: " << X_Range << endl;
            fitting_Line->SetParameter(0,Fitting_a[Applied_Hist]);
            fitting_Line->SetParameter(1,Fitting_b[Applied_Hist]);
            //cout << "velocity_Used[Applied_Hist+1]: " << velocity_Used[Applied_Hist+1] << endl;
            //cout << "Applied_Hist: " << Applied_Hist << endl;
            //cout << "T_Central: " << T_Central << endl;
            //cout << "fitting_Line->Eval(T_Central): " << fitting_Line->Eval(T_Central) << endl;
            double dsigma_dT  = fitting_Line->Eval(T_Central);
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            //cout << "1e-15*TMath::Power(Scaling[Index],2): " << 1e-15*TMath::Power(Scaling[Index],2) << endl;
            //cout << "dsigma_dT: " << dsigma_dT << endl;
            //cout << "dT: " << dT << endl;
            
            Total_Count = Total_Count + Length*(Density/((unified_atomic_mass_g*(ASi))))*(dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
            Area_Integral = Area_Integral + (dsigma_dT)*(dT);
            //cout << "Total_Count: " << Total_Count << endl;
            T = T_Temp;
        }
        //cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        cout << "Area_Integral: " << Area_Integral << endl;
        cout << "Total_Count: " << Total_Count << endl;
        
        double Scaling = (collision_Time_Ave)/(Total_Count);
        //cout << "Scaling: " << Scaling << endl;
        double Scaling_Factor = sqrt(Scaling);
        
        cout << "CS_Try(1*Scaling_Factor,WIMP_mx_Array[kkk]): " << CS_Try(1,WIMP_mx_Array[FIle_Index]) << endl;
        cout << "DS_Try(1e-9*Scaling_Factor,WIMP_mx_Array[kkk]): " << DS_Try(1e-9,WIMP_mx_Array[FIle_Index]) << endl;
        if(Index==0)Cross_Section_Set[FIle_Index]= ( CS_Try(1*Scaling_Factor,WIMP_mx_Array[FIle_Index]) );
        if(Index==1)Cross_Section_Set[FIle_Index]= ( DS_Try(1e-9*Scaling_Factor,WIMP_mx_Array[FIle_Index]) );

        //cout << "WIMP_mx_Array[kkk]: " << WIMP_mx_Array[kkk] << endl;
    }
    
    cout << "========================================================" << endl;
    cout << "File_Number: " << File.size() << endl;
    cout << "WIMP_mx_Array_Number: " << WIMP_mx_Array.size() << endl;
    for(int kkk=0; kkk< File.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Cross_Section_Set[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< File.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Collision_Time_Array[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< File.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Energy_Loss_Per_Collision[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< File.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Energy_Loss_Per_Count_Only[kkk] << "," << endl;
    }

    //==================================================//

}//End_Main
*/

/*
for(int Applied_Hist=0; Applied_Hist<velocity_TH1F.size(); Applied_Hist++)
{
    //cout << "velocity_Used[Applied_Hist+1]: " << Max_X_Array[Applied_Hist+1] << endl;
    //cout << "Max_X_Array[Applied_Hist]: " << Max_X_Array[Applied_Hist] << endl;
    TF1 *fitting_Line = new TF1("fitting_Line","TMath::Power(10,[0]*x+[1])",0,Energy_DM_Max);
    fitting_Line->SetParameter(0,Fitting_a[Applied_Hist]);
    fitting_Line->SetParameter(1,Fitting_b[Applied_Hist]);
    TF1_Fitting_Line.push_back(fitting_Line);
}
 */

double dE_dX_ER_from_Paper(double E_0, double E_d, double Length)//Arxiv: 1905.06348, E_0(keV), E_d(keV), Length(km)
{
    return (E_0-E_d)/(Length);
}

double Max_Recoil_elastic(double mx, double Velocity)//Mass(GeV/c^2),velocity(km/s)
{
    return TMath::Power(10,9)*(mx*Velocity*1e3/3e8)*(mx*Velocity*1e3/3e8)/(2*28);//eV
}
double q_max_from_paper(double mx, double Velocity)//Calculated with q2
{
    double reduce_mass_DM = mx*1e9*1*1e9/((mx*1e9)+(1*1e9));//(eV/c^2)
    double v_c            = Velocity*1e3/3e8;//c
    //cout << "q_max: " << (4*reduce_mass_DM*reduce_mass_DM*v_c*v_c)/(2.0*1e9*28.0) << endl;
    return sqrt(4*reduce_mass_DM*reduce_mass_DM*v_c*v_c);//(eV/c)
}
double v_min(double mx, double q){return q/(2*mx);}

double Energy_Transfer_ER_atomic_scattering(double mx, double Velocity)//mx(GeV/c^2),v(km/s)
{
    double v_c      = Velocity*1e3/3e8;
    double momentum = mx*v_c;//(GeV/c)
    double Energy   = momentum*momentum/(2.0*28.0);//  (GeV/c)^2/(GeV/c^2)=GeV
    
    return Energy*1e9;
}

const double Length_for_exp[3]={0.107,2,1.78};//(km), for MINOS (107 m underground),SNOLAB (2000 m underground), and DAMIC-M in Modane (1780 m underground)
const double Cross_Section_for_Exp[3]={9e-22,3e-23,3.5e-23};//(km), for MINOS (107 m underground),SNOLAB (2000 m underground), and DAMIC-M in Modane (1780 m underground)

 
 //CDEX
const double Density = 1.8;//g/cm^3
const double Length  = 1e5;//cm
 
   // XENON10/XENON100
/*
const double Density = 1.8;//g/cm^3
const double Length  = 1.4e5;//cm
*/
/* //TEXONO
const double Density = 1;//g/cm^3
const double Length  = 3e3;//cm
*/

/* //MINOS
const double Density = 1.8;//g/cm^3
const double Length  = 0.107e5;//cm
*/

/*
const double Density = 11.29;//g/cm^3
const double Length  = 50.;//cm
*/

/*
const double Density = 1.8;//g/cm^3
const double Length  = 1.07e4;//cm
*/

/*
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()//Upper boundary
{
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    
    const int Index=0;
    
    //const double Length  = 3e0;//cm
    const double Threshold_keV = 0.2;//keV
    vector<string> File;
    vector<double> WIMP_mx_Array;
    vector<double> Collision_Time_Array;
    vector<double> dE_dX_1_from_mukesh_array;
    vector<double> dE_dX_2_from_Paper_array;

    vector<double> Energy_Loss_Per_Collision;
    vector<double> Scaling={1e-18,1e-9};
    //=========Define the file==========//
    //Xe_c1
    if(Index==0)
    {
    File=
            {"c1_20_0GeV",
            "c1_10_0GeV","c1_5_0GeV",
            "c1_4_0GeV","c1_3_0GeV",
            "c1_2_0GeV","c1_1_0GeV",
            "c1_0_900GeV","c1_0_800GeV",
            "c1_0_700GeV","c1_0_600GeV",
            "c1_0_500GeV","c1_0_400GeV",
            "c1_0_300GeV","c1_0_200GeV",
            "c1_0_120GeV","c1_0_090GeV",
            "c1_0_070GeV","c1_0_050GeV",
            "c1_0_040GeV",
            "c1_0_020GeV","c1_0_010GeV",
            };
        WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05,0.04,0.02,0.01};
    }
    //Xe_d1
    if(Index==1)
    {
         File=
            {
            "d1_20_0GeV",
            "d1_10_0GeV","d1_5_0GeV",
            "d1_4_0GeV","d1_3_0GeV",
            "d1_2_0GeV","d1_1_0GeV",
            "d1_0_900GeV","d1_0_800GeV",
            "d1_0_700GeV","d1_0_600GeV",
            "d1_0_500GeV","d1_0_400GeV",
            "d1_0_300GeV","d1_0_200GeV",
            "d1_0_120GeV","d1_0_090GeV",
            "d1_0_080GeV","d1_0_070GeV",
            "d1_0_050GeV",
            "d1_0_040GeV",
            "d1_0_020GeV","d1_0_010GeV",
            };
        WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.08,0.07,0.05,0.04,0.02,0.01};
    }
    //Read the file of DCS for different masses
    vector<double> Cross_Section_Set;

    
    string c1_d1_Xe_Ge_index[4]={"SI_c1_XeData_Vel","SI_d1_XeData_Vel","SI_c1_GeData_Vel","SI_d1_GeData_Vel"};
    //===============electron recoil===================//
    //for(int Index_File=0; Index_File<File.size(); Index_File++)
    //for(int Index_File=6; Index_File<7; Index_File++) //For 1.0GeV
    for(int Index_File=15; Index_File<16; Index_File++) //For 0.1GeV
    //for(int Index_File=File.size()-1; Index_File<File.size(); Index_File++) //For 0.01GeV
    {
        TFile *fin  = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str()  ,File[Index_File].c_str()));

        vector<TH1F*> velocity_TH1F   ;//Mass-based
        vector<double>   velocity_Used;//km/s
        velocity_Used.push_back(0);
        for(int kkk=0; kkk<velocity_N; kkk++)
        {
            TH1F *Velocity_TH1F_temp=(TH1F*)fin->Get(velocity_s[kkk].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                cout << "File: " << velocity_s[kkk].c_str() << endl;
                //cout << "velocitykm: " << velocitykm[kkk] << endl;
                velocity_TH1F.push_back(Velocity_TH1F_temp);
                velocity_Used.push_back(velocitykm[kkk]);
            }
        }
        
        int    Applied_Hist=0;
        double Velocity_DM_Now=779;//km/s
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        double v_c = Velocity_DM_Now*1e3/3e8;//c
        double Max_Recoil = DM_E(WIMP_mx_Array[Index_File],Velocity_DM_Now);//Max_Recoil(keV) given by mass
        Int_t Maximum_Bin = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Max_Recoil*1e3);//Max_Recoil_Bin
        //cout << "Maximum_Bin: " << Maximum_Bin << endl;
        int LastBin = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
        //cout << "LastBin: " << LastBin << endl;

        //
        vector<double> Fitting_a;vector<double> Fitting_b;
        Fitting_a.clear();Fitting_b.clear();
        for(int Applied_Hist_1=velocity_TH1F.size()-1; Applied_Hist_1<velocity_TH1F.size(); Applied_Hist_1++)
        //for(int Applied_Hist=0; Applied_Hist<3; Applied_Hist++)
        {
            vector<double> TGraph_Recoil_Energy;vector<double> TGraph_dsigma_dT;
            
            Int_t Minimum_Bin = velocity_TH1F[Applied_Hist_1]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin = velocity_TH1F[Applied_Hist_1]->FindLastBinAbove();//Max_Recoil_Bin
            for(int LLL=Minimum_Bin; LLL<Maximum_Bin+1; LLL++)
            {
                TGraph_Recoil_Energy.push_back( velocity_TH1F[Applied_Hist_1]->GetXaxis()->GetBinCenter(LLL) );//keV recoil energy
                TGraph_dsigma_dT.push_back( TMath::Log10(velocity_TH1F[Applied_Hist_1]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2)) );//Real dsigma_dT
            }
            TGraph * g = new TGraph((int)TGraph_Recoil_Energy.size(), &TGraph_Recoil_Energy[0], &TGraph_dsigma_dT[0]);
            TF1 *fitting_Line = new TF1("fitting_Line","[0]*x+[1]",12,velocity_TH1F[Applied_Hist_1]->GetXaxis()->GetBinCenter( Maximum_Bin-1 ));//Fitting
            g->Fit(fitting_Line,"R");
            Double_t par[9];
            fitting_Line->GetParameters(&par[0]);
            Fitting_a.push_back(par[0]);
            Fitting_b.push_back(par[1]);
            cout << "par[0]: " << par[0] << endl;
            cout << "par[1]: " << par[1] << endl;
        }

        
        //=======================================Start interacting=======================================//
        
        int Event_Number=0;
        double collision_Time_Total=0;
        double Energy_Loss_eV_collision=0;
        

        while(Event_Number<50)
        {
            //cout << "=====Event_Number: " << Event_Number << endl;
            double Velocity_DM_Temp = 784;//Initial Energy and iteration =>km/s
            double Energy_DM_Temp   = DM_E(WIMP_mx_Array[Index_File],784);//Initial Energy and iteration =>keV
            int collision_Time=0;
            double Energy_Loss_Inidivual=0;
            //for(int kkk=0; kkk<1; kkk++)
            //while(Energy_DM_Temp>0.01)
            while(Energy_DM_Temp>0.0011)
            {
                //Confirm the TH1F
                gRandom = new TRandom3(0);
                gRandom->SetSeed(0);
                
                if(Velocity_DM_Temp>velocity_Used[velocity_Used.size()-1])
                {
                    Applied_Hist=(velocity_TH1F.size())-1;
                }
                if(Applied_Hist==0)
                {
                    for(int kkk=0; kkk<velocity_Used.size(); kkk++)
                    {
                        if(Velocity_DM_Temp>velocity_Used[kkk] and Velocity_DM_Temp<velocity_Used[kkk+1])
                        {
                            Applied_Hist=kkk;
                        }
                    }
                }
                //Set the range
                int Maximum_Bin_Loss = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Energy_DM_Temp*1e3);//Max_Recoil_Bin
                //cout << "Applied_Hist: " << Applied_Hist << endl;
                int LastBin_Loss = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_from_Hist
                if(Maximum_Bin_Loss>LastBin_Loss)Maximum_Bin_Loss=LastBin_Loss;
                double  Max_X        = velocity_TH1F[Applied_Hist]->GetBinCenter(Maximum_Bin_Loss);
                velocity_TH1F[Applied_Hist]->GetXaxis()->SetRange(12,Max_X);
                double Energy_Loss_eV = 0;//Energy_Loss(eV)
                
                while(Energy_Loss_eV<12)
                {Energy_Loss_eV = velocity_TH1F[Applied_Hist]->GetRandom();}
                                
                double Energy_Loss_keV = Energy_Loss_eV*1e-3;//Energy_Loss(keV)
                //cout << "Energy_Loss_(eV): " << Energy_Loss_eV << endl;
                Energy_Loss_Inidivual = Energy_Loss_Inidivual + Energy_Loss_eV;
                //cout << "Energy_Loss_Inidivual: " << Energy_Loss_Inidivual << endl;
                //Energy_DM_Temp   = Energy_DM_Temp - Energy_Loss_keV;
                //Energy_DM_Temp   = Energy_DM_Temp - 1.9*1e-5;
                Energy_DM_Temp   = Energy_DM_Temp - Energy_Loss_eV*1e-3;

                //cout << "Energy_DM_Temp(keV): " << Energy_DM_Temp << endl;
                Velocity_DM_Temp = sqrt(2*Energy_DM_Temp/(WIMP_mx_Array[Index_File]*1e6))*(3e8/1e3);//km/s
                //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;//km/s
                collision_Time = collision_Time + 1;
                //cout << "Average_Energy_Loss(eV): " << Energy_Loss_Inidivual/collision_Time << endl;
                
                //cout << "collision_Time: " << collision_Time << endl;
                Applied_Hist =0;
            }
            //cout << "collision_Time: " << collision_Time << endl;
            //if(collision_Time>0)cout << "Energy_Loss_Inidivual/collision_Time: " << Energy_Loss_Inidivual*1e3/collision_Time << endl;
            Energy_Loss_eV_collision = Energy_Loss_eV_collision + Energy_Loss_Inidivual*1e3/collision_Time;
            collision_Time_Total = collision_Time_Total + collision_Time;
            Event_Number = Event_Number + 1;
        }
        double collision_Time_Ave = collision_Time_Total / Event_Number;
        //double collision_Time_Ave = DM_E(WIMP_mx_Array[Index_File],784) / ( 4e-04);

        Collision_Time_Array.push_back(collision_Time_Ave);//Count for every DM
        //cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        double Energy_Loss_Real = Energy_DM(WIMP_mx_Array[Index_File],779*(1e3/3e8))-Threshold_keV;
        //cout << "Energy_Loss_Real: " << Energy_Loss_Real << endl;
        //Energy_Loss_eV_collision = (1./Energy_Loss_Real) * (Energy_Loss_eV_collision/Event_Number) ;
        Energy_Loss_Per_Collision.push_back(Energy_Loss_eV_collision/Event_Number);
        //cout << "Threshold: " << sqrt(2*0.16/(WIMP_mx_Array[kkk]*1e6))*(3e8/1e3);
        //cout << "DM_E(WIMP_mx_Array[kkk],779) : " << DM_E(WIMP_mx_Array[kkk],779) << endl;;
        //cout << "Energy_Loss_eV_collision: " << Energy_Loss_eV_collision << endl;
        //==================================================//
        
        //double Target_Total_Cross_Section = collision_Time_Ave*(unified_atomic_mass_g*(ASi))/(Density*Length);
        
        if(Velocity_DM_Now>velocity_Used[velocity_Used.size()-1])Applied_Hist=(velocity_TH1F.size())-1;
        if(Applied_Hist==0)
        {
            for(int kkk=0; kkk<velocity_Used.size(); kkk++)
            {
                if(Velocity_DM_Now>velocity_Used[kkk] and Velocity_DM_Now<velocity_Used[kkk+1]){Applied_Hist=kkk;}
            }
        }
        
        for(int kkk=0; kkk<3; kkk++)
        {
            cout << "c_1(Cross_Section_for_Exp[kkk],WIMP_mx_Array[Index_File]): " << c_1(Cross_Section_for_Exp[kkk],WIMP_mx_Array[Index_File]) << endl;;
        }

        double Total_Count=0;
        const double constant_for_total_count = Length*(Density/((unified_atomic_mass_g*(ASi))));
        double dE_dX_1_from_mukesh = 0;//N_atom_1kg_Ge_Electron(Number density)

        for(int kkk=1; kkk<Maximum_Bin-1; kkk++)
        {
            double dsigma_dT  = velocity_TH1F[Applied_Hist]->GetBinContent(kkk)*1e-15*TMath::Power(Scaling[Index],2);
            double dT         = velocity_TH1F[Applied_Hist]->GetBinWidth(kkk)*1e-3;
            double T          = velocity_TH1F[Applied_Hist]->GetBinCenter(kkk)*1e-3;

            if(velocity_TH1F[Applied_Hist]->GetBinCenter(kkk)>12)Total_Count = Total_Count + (dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
            if(velocity_TH1F[Applied_Hist]->GetBinCenter(kkk)>12)dE_dX_1_from_mukesh = dE_dX_1_from_mukesh + 1e-3*(dsigma_dT)*(dT)*T;//
            //cout << "Total_Count: " << Total_Count << endl;
        }
        dE_dX_1_from_mukesh = TMath::Power(10,12)*dE_dX_1_from_mukesh*N_atom_1kg_Ge_Electron;
        
        double X_Range         =  Max_Recoil*1e3;//eV
        double X_Pre_Range     =  DM_E(WIMP_mx_Array[Index_File],650)*1e3;//eV

        double X_Range_Bin = (Max_Recoil*1e3)*2e-3;//eV
        double T           =  0;
        double Total_Count_1= 0;
        double T_Temp       = 0;double dT      = 0;double T_Central = 0;
        //cout << "Fitting_a[0]: " << Fitting_a[0] << endl;
        //cout << "Fitting_b[0]: " << Fitting_b[0] << endl;
        for(int kkk=1; kkk<499; kkk++)
        {
            //cout << kkk << endl;
            T_Temp    = T_Temp + X_Range_Bin;
            dT        = (T_Temp - T)*1e-3;//keV
            //cout << "Max_Recoil: " << Max_Recoil << endl;
            T_Central = (T + T_Temp)*0.5;//eV
            TF1 *fitting_Line = new TF1("fitting_Line","TMath::Power(10,([0]*x+[1]))",0,X_Range);
            fitting_Line->SetParameter(0,Fitting_a[0]);
            fitting_Line->SetParameter(1,Fitting_b[0]);
            double dsigma_dT  = fitting_Line->Eval(T_Central);
            //cout << "T_Central: " << T_Central <<endl;
            //cout << "dsigma_dT: " << dsigma_dT << endl;
            if(T_Central>X_Pre_Range and T_Central<X_Range) Total_Count_1 = Total_Count_1 + (dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
            T = T_Temp;
        }
        cout << "Total_Count: " << Total_Count  << endl;
        //cout << "Total_Count_1: " << Total_Count_1  << endl;
        Total_Count = Total_Count + Total_Count_1;
        //cout << "X_Range: " << X_Range << endl;
        Total_Count = Total_Count*constant_for_total_count;
        //cout << "velocity_TH1F[Applied_Hist]->GetMean(): " << velocity_TH1F[Applied_Hist]->GetMean() << endl;
        //cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        //cout << "Total_Count: " << Total_Count << endl;
        double Scaling = (collision_Time_Ave)/(Total_Count);
        //cout << "WIMP_mx_Array[Index_File]: " << WIMP_mx_Array[Index_File] << endl;
        //cout << "Scaling: " << Scaling << endl;
        double c1_GeV2 = sqrt(Scaling);
        double d1      = sqrt(Scaling);
        
        cout << "CS_Try(1*c1_GeV2,WIMP_mx_Array[Index_File]): " << TMath::Power(10,12)*CS_Try(1,WIMP_mx_Array[Index_File]) << endl;
        cout << "DS_Try(1e-9*d1,WIMP_mx_Array[Index_File]): " << DS_Try(1e-9,WIMP_mx_Array[Index_File]) << endl;
        if(Index==0)Cross_Section_Set.push_back( CS_Try(1*c1_GeV2,WIMP_mx_Array[Index_File]) );
        if(Index==1)Cross_Section_Set.push_back( DS_Try(1e-9*d1,WIMP_mx_Array[Index_File]) );
                
        double E_0_Paper = DM_E(WIMP_mx_Array[Index_File],784);
        double E_d_Paper = 1.1e-3;
        double dE_dX_2_from_Paper  = dE_dX_ER_from_Paper(E_0_Paper,E_d_Paper,Length/(1e5));//dE/dX from 1905.06348-2.pdf
        
        dE_dX_1_from_mukesh_array.push_back(dE_dX_1_from_mukesh);
        dE_dX_2_from_Paper_array.push_back(dE_dX_2_from_Paper);
        //cout << "WIMP_mx_Array[kkk]: " << WIMP_mx_Array[kkk] << endl;
    }
    //From the paper of others'
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << dE_dX_1_from_mukesh_array[kkk] << "," << endl;
        cout << WIMP_mx_Array[kkk] << "," << dE_dX_2_from_Paper_array[kkk] << "," << endl;
    }
    
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << dE_dX_1_from_mukesh_array[kkk] << "," << endl;
        cout << WIMP_mx_Array[kkk] << "," << dE_dX_2_from_Paper_array[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Cross_Section_Set[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Collision_Time_Array[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Energy_Loss_Per_Collision[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;

    //==================================================//
     
}//End_Main
*/
/*
int CCC=0;
while(CCC==0)
{
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);

    TF1 *Atomic_scattering = new TF1("Atomic_scattering","1",0,q_max_from_paper(WIMP_mx_Array[Index_File],Velocity_DM_Temp));
    double q_chosen        = Atomic_scattering->GetRandom();//eV/c
    double e_chosen        = q_chosen*q_chosen/(2.0*1e9*28.0);//eV
    double v_min_chosen    = v_min(WIMP_mx_Array[Index_File]*1e9,q_chosen)*3e8/1e3;
    if(Velocity_DM_Temp>v_min_chosen)
    {
        CCC=1;
        Energy_Loss_eV = e_chosen;
        cout << "q_chosen: " << q_chosen << endl;
        cout << "e_chosen: " << e_chosen << endl;
    }
}
 */

//cout << "Applied_Hist: " << Applied_Hist << endl;
//cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
//double dsigma_dT  = velocity_TH1F[Applied_Hist]->GetBinContent(kkk)*1e-15*TMath::Power(Scaling[Index],2)*TMath::Power(c_1(Cross_Section_for_Exp[0],WIMP_mx_Array[Index_File]),2);
//cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
//cout << "1e-15*TMath::Power(Scaling[Index],2): " << 1e-15*TMath::Power(Scaling[Index],2) << endl;
//cout << "dsigma_dT: " << dsigma_dT << endl;
//cout << "dT: " << dT << endl;
//cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
//cout << "velocity_TH1F[Applied_Hist]->GetBinCenter(kkk): " << velocity_TH1F[Applied_Hist]->GetBinCenter(kkk) << endl;


/*
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()
{
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    
    const int Index=0;
    
    //const double Length  = 3e0;//cm
    const double Threshold_keV = 0.2;//keV
    vector<string> File;
    vector<double> WIMP_mx_Array;
    vector<double> Collision_Time_Array;
    vector<double> dE_dX_1_from_mukesh_array;
    vector<double> dE_dX_2_from_Paper_array;

    vector<double> Energy_Loss_Per_Collision;
    vector<double> Scaling={1e-18,1e-9};
    //=========Define the file==========//
    //Xe_c1
    if(Index==0)
    {
    File=
            {"c1_20_0GeV",
            "c1_10_0GeV","c1_5_0GeV",
            "c1_4_0GeV","c1_3_0GeV",
            "c1_2_0GeV","c1_1_0GeV",
            "c1_0_900GeV","c1_0_800GeV",
            "c1_0_700GeV","c1_0_600GeV",
            "c1_0_500GeV","c1_0_400GeV",
            "c1_0_300GeV","c1_0_200GeV",
            "c1_0_120GeV","c1_0_090GeV",
            "c1_0_070GeV","c1_0_050GeV",
            "c1_0_040GeV",
            "c1_0_020GeV","c1_0_010GeV",
            };
        WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05,0.04,0.02,0.01};
    }
    //Xe_d1
    if(Index==1)
    {
         File=
            {
            "d1_20_0GeV",
            "d1_10_0GeV","d1_5_0GeV",
            "d1_4_0GeV","d1_3_0GeV",
            "d1_2_0GeV","d1_1_0GeV",
            "d1_0_900GeV","d1_0_800GeV",
            "d1_0_700GeV","d1_0_600GeV",
            "d1_0_500GeV","d1_0_400GeV",
            "d1_0_300GeV","d1_0_200GeV",
            "d1_0_120GeV","d1_0_090GeV",
            "d1_0_080GeV","d1_0_070GeV",
            "d1_0_050GeV",
            "d1_0_040GeV",
            "d1_0_020GeV","d1_0_010GeV",
            };
        WIMP_mx_Array ={20.0,10.0,5.0,4.0,3.0,2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.08,0.07,0.05,0.04,0.02,0.01};
    }
    //Read the file of DCS for different masses
    vector<double> Cross_Section_Set;

    
    string c1_d1_Xe_Ge_index[4]={"SI_c1_XeData_Vel","SI_d1_XeData_Vel","SI_c1_GeData_Vel","SI_d1_GeData_Vel"};
    //===============electron recoil===================//
    //for(int Index_File=0; Index_File<File.size(); Index_File++)
    //for(int Index_File=6; Index_File<7; Index_File++) //For 1.0GeV
    for(int Index_File=15; Index_File<16; Index_File++) //For 0.1GeV
    //for(int Index_File=File.size()-1; Index_File<File.size(); Index_File++) //For 0.01GeV
    {
        TFile *fin = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str(),File[Index_File].c_str()));
        vector<TH1F*> velocity_TH1F;//Mass-based
        vector<double>   velocity_Used;//km/s
        velocity_Used.push_back(0);
        for(int kkk=0; kkk<velocity_N; kkk++)
        {
            TH1F *Velocity_TH1F_temp=(TH1F*)fin->Get(velocity_s[kkk].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                cout << "File: " << velocity_s[kkk].c_str() << endl;
                //cout << "velocitykm: " << velocitykm[kkk] << endl;
                velocity_TH1F.push_back(Velocity_TH1F_temp);
                velocity_Used.push_back(velocitykm[kkk]);
            }
        }
        
        int    Applied_Hist=0;
        double Velocity_DM_Now=779;//km/s
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        double v_c = Velocity_DM_Now*1e3/3e8;//c
        double Max_Recoil = DM_E(WIMP_mx_Array[Index_File],Velocity_DM_Now);//Max_Recoil(keV) given by mass
        Int_t Maximum_Bin = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Max_Recoil*1e3);//Max_Recoil_Bin
        //cout << "Maximum_Bin: " << Maximum_Bin << endl;
        int LastBin = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
        //cout << "LastBin: " << LastBin << endl;

        //
        vector<double> Fitting_a;vector<double> Fitting_b;
        Fitting_a.clear();Fitting_b.clear();
        for(int Applied_Hist_1=velocity_TH1F.size()-1; Applied_Hist_1<velocity_TH1F.size(); Applied_Hist_1++)
        //for(int Applied_Hist=0; Applied_Hist<3; Applied_Hist++)
        {
            vector<double> TGraph_Recoil_Energy;vector<double> TGraph_dsigma_dT;
            
            Int_t Minimum_Bin = velocity_TH1F[Applied_Hist_1]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin = velocity_TH1F[Applied_Hist_1]->FindLastBinAbove();//Max_Recoil_Bin
            for(int LLL=Minimum_Bin; LLL<Maximum_Bin+1; LLL++)
            {
                TGraph_Recoil_Energy.push_back( velocity_TH1F[Applied_Hist_1]->GetXaxis()->GetBinCenter(LLL) );//keV recoil energy
                TGraph_dsigma_dT.push_back( TMath::Log10(velocity_TH1F[Applied_Hist_1]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2)) );//Real dsigma_dT
            }
            TGraph * g = new TGraph((int)TGraph_Recoil_Energy.size(), &TGraph_Recoil_Energy[0], &TGraph_dsigma_dT[0]);
            TF1 *fitting_Line = new TF1("fitting_Line","[0]*x+[1]",12,velocity_TH1F[Applied_Hist_1]->GetXaxis()->GetBinCenter( Maximum_Bin-1 ));//Fitting
            g->Fit(fitting_Line,"R");
            Double_t par[9];
            fitting_Line->GetParameters(&par[0]);
            Fitting_a.push_back(par[0]);
            Fitting_b.push_back(par[1]);
            cout << "par[0]: " << par[0] << endl;
            cout << "par[1]: " << par[1] << endl;
        }

        
        //=======================================Start interacting=======================================//
        
        int Event_Number=0;
        double collision_Time_Total=0;
        double Energy_Loss_eV_collision=0;
        while(Event_Number<10)
        {
            //cout << "=====Event_Number: " << Event_Number << endl;
            double Velocity_DM_Temp = 779;//Initial Energy and iteration =>km/s
            double Energy_DM_Temp   = DM_E(WIMP_mx_Array[Index_File],779);//Initial Energy and iteration =>keV
            int collision_Time=0;
            double Energy_Loss_Inidivual=0;
            //for(int kkk=0; kkk<1; kkk++)
            //while(Energy_DM_Temp>0.01)
            while(Energy_DM_Temp>0.01)
            {
                //Confirm the TH1F
                gRandom = new TRandom3(0);
                gRandom->SetSeed(0);
                
                if(Velocity_DM_Temp>velocity_Used[velocity_Used.size()-1])
                {
                    Applied_Hist=(velocity_TH1F.size())-1;
                }
                if(Applied_Hist==0)
                {
                    for(int kkk=0; kkk<velocity_Used.size(); kkk++)
                    {
                        if(Velocity_DM_Temp>velocity_Used[kkk] and Velocity_DM_Temp<velocity_Used[kkk+1])
                        {
                            Applied_Hist=kkk;
                        }
                    }
                }
                //Set the range
                int Maximum_Bin_Loss = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Energy_DM_Temp*1e3);//Max_Recoil_Bin
                //cout << "Applied_Hist: " << Applied_Hist << endl;
                int LastBin_Loss = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_from_Hist
                if(Maximum_Bin_Loss>LastBin_Loss)Maximum_Bin_Loss=LastBin_Loss;
                double  Max_X        = velocity_TH1F[Applied_Hist]->GetBinCenter(Maximum_Bin_Loss);
                velocity_TH1F[Applied_Hist]->GetXaxis()->SetRange(12,Max_X);
                double Energy_Loss_eV = 0;//Energy_Loss(eV)
                while(Energy_Loss_eV<12)
                {Energy_Loss_eV = velocity_TH1F[Applied_Hist]->GetRandom();}
                double Energy_Loss_keV = Energy_Loss_eV*1e-3;//Energy_Loss(keV)
                //cout << "Energy_Loss_(keV): " << Energy_Loss_eV*1e-3 << endl;
                Energy_DM_Temp   = Energy_DM_Temp - Energy_Loss_keV;
                //Energy_DM_Temp   = Energy_DM_Temp - 0.55*1e-3;
                //cout << "Energy_DM_Temp(keV): " << Energy_DM_Temp << endl;
                Velocity_DM_Temp = sqrt(2*Energy_DM_Temp/(WIMP_mx_Array[Index_File]*1e6))*(3e8/1e3);//km/s
                //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;//km/s
                Energy_Loss_Inidivual = Energy_Loss_Inidivual + Energy_Loss_keV;
                collision_Time = collision_Time + 1;
                //cout << "collision_Time: " << collision_Time << endl;
                Applied_Hist =0;
            }
            //cout << "collision_Time: " << collision_Time << endl;
            //if(collision_Time>0)cout << "Energy_Loss_Inidivual/collision_Time: " << Energy_Loss_Inidivual*1e3/collision_Time << endl;
            Energy_Loss_eV_collision = Energy_Loss_eV_collision + Energy_Loss_Inidivual*1e3/collision_Time;
            collision_Time_Total = collision_Time_Total + collision_Time;
            Event_Number = Event_Number + 1;
        }
        double collision_Time_Ave = collision_Time_Total / Event_Number;
        Collision_Time_Array.push_back(collision_Time_Ave);//Count for every DM
        //cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        double Energy_Loss_Real = Energy_DM(WIMP_mx_Array[Index_File],779*(1e3/3e8))-Threshold_keV;
        //cout << "Energy_Loss_Real: " << Energy_Loss_Real << endl;
        //Energy_Loss_eV_collision = (1./Energy_Loss_Real) * (Energy_Loss_eV_collision/Event_Number) ;
        Energy_Loss_Per_Collision.push_back(Energy_Loss_eV_collision/Event_Number);
        //cout << "Threshold: " << sqrt(2*0.16/(WIMP_mx_Array[kkk]*1e6))*(3e8/1e3);
        //cout << "DM_E(WIMP_mx_Array[kkk],779) : " << DM_E(WIMP_mx_Array[kkk],779) << endl;;
        //cout << "Energy_Loss_eV_collision: " << Energy_Loss_eV_collision << endl;
        //==================================================//
        
        //double Target_Total_Cross_Section = collision_Time_Ave*(unified_atomic_mass_g*(ASi))/(Density*Length);
        
        if(Velocity_DM_Now>velocity_Used[velocity_Used.size()-1])Applied_Hist=(velocity_TH1F.size())-1;
        if(Applied_Hist==0)
        {
            for(int kkk=0; kkk<velocity_Used.size(); kkk++)
            {
                if(Velocity_DM_Now>velocity_Used[kkk] and Velocity_DM_Now<velocity_Used[kkk+1]){Applied_Hist=kkk;}
            }
        }
        
        for(int kkk=0; kkk<3; kkk++)
        {
            cout << "c_1(Cross_Section_for_Exp[kkk],WIMP_mx_Array[Index_File]): " << c_1(Cross_Section_for_Exp[kkk],WIMP_mx_Array[Index_File]) << endl;;
        }

        double Total_Count=0;
        const double constant_for_total_count = Length*(Density/((unified_atomic_mass_g*(ASi))));
        for(int kkk=1; kkk<Maximum_Bin-1; kkk++)
        {
            //cout << "Applied_Hist: " << Applied_Hist << endl;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            double dsigma_dT  = velocity_TH1F[Applied_Hist]->GetBinContent(kkk)*1e-15*TMath::Power(Scaling[Index],2);
            double dT         = velocity_TH1F[Applied_Hist]->GetBinWidth(kkk)*1e-3;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            //cout << "1e-15*TMath::Power(Scaling[Index],2): " << 1e-15*TMath::Power(Scaling[Index],2) << endl;
            //cout << "dsigma_dT: " << dsigma_dT << endl;
            //cout << "dT: " << dT << endl;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinCenter(kkk): " << velocity_TH1F[Applied_Hist]->GetBinCenter(kkk) << endl;

            if(velocity_TH1F[Applied_Hist]->GetBinCenter(kkk)>12)Total_Count = Total_Count + (dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
            //cout << "Total_Count: " << Total_Count << endl;
        }
        double X_Range         =  Max_Recoil*1e3;//eV
        double X_Pre_Range     =  DM_E(WIMP_mx_Array[Index_File],650)*1e3;//eV

        double X_Range_Bin = (Max_Recoil*1e3)*2e-3;//eV
        double T           =  0;
        double Total_Count_1= 0;
        double T_Temp       = 0;double dT      = 0;double T_Central = 0;
        //cout << "Fitting_a[0]: " << Fitting_a[0] << endl;
        //cout << "Fitting_b[0]: " << Fitting_b[0] << endl;
        for(int kkk=1; kkk<499; kkk++)
        {
            //cout << kkk << endl;
            T_Temp    = T_Temp + X_Range_Bin;
            dT        = (T_Temp - T)*1e-3;//keV
            //cout << "Max_Recoil: " << Max_Recoil << endl;
            T_Central = (T + T_Temp)*0.5;//eV
            TF1 *fitting_Line = new TF1("fitting_Line","TMath::Power(10,([0]*x+[1]))",0,X_Range);
            fitting_Line->SetParameter(0,Fitting_a[0]);
            fitting_Line->SetParameter(1,Fitting_b[0]);
            double dsigma_dT  = fitting_Line->Eval(T_Central);
            //cout << "T_Central: " << T_Central <<endl;
            //cout << "dsigma_dT: " << dsigma_dT << endl;
            if(T_Central>X_Pre_Range and T_Central<X_Range) Total_Count_1 = Total_Count_1 + (dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
            T = T_Temp;
        }
        //cout << "Total_Count: " << Total_Count  << endl;
        //cout << "Total_Count_1: " << Total_Count_1  << endl;
        Total_Count = Total_Count + Total_Count_1;
        //cout << "X_Range: " << X_Range << endl;
        Total_Count = Total_Count*constant_for_total_count;
        //cout << "velocity_TH1F[Applied_Hist]->GetMean(): " << velocity_TH1F[Applied_Hist]->GetMean() << endl;
        //cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        //cout << "Total_Count: " << Total_Count << endl;
        
        double Scaling = (collision_Time_Ave)/(Total_Count);
        //cout << "WIMP_mx_Array[Index_File]: " << WIMP_mx_Array[Index_File] << endl;
        //cout << "Scaling: " << Scaling << endl;
        double c1_GeV2 = sqrt(Scaling);
        double d1      = sqrt(Scaling);
        cout << "c1_GeV2: " << c1_GeV2 << endl;
        cout << "CS_Try(1*c1_GeV2,WIMP_mx_Array[Index_File]): " << CS_Try(1*c1_GeV2,WIMP_mx_Array[Index_File]) << endl;
        cout << "DS_Try(1e-9*d1,WIMP_mx_Array[Index_File]): " << DS_Try(1e-9*d1,WIMP_mx_Array[Index_File]) << endl;
        if(Index==0)Cross_Section_Set.push_back( CS_Try(1*c1_GeV2,WIMP_mx_Array[Index_File]) );
        if(Index==1)Cross_Section_Set.push_back( DS_Try(1e-9*d1,WIMP_mx_Array[Index_File]) );
        
        
        double E_0_Paper = DM_E(WIMP_mx_Array[Index_File],779);
        double E_d_Paper = 1.1e-3;
        double dE_dX_1_from_mukesh = 0;
        double dE_dX_2_from_Paper  = dE_dX_ER_from_Paper(E_0_Paper,E_d_Paper,Length/(1e5));
        dE_dX_1_from_mukesh_array.push_back(dE_dX_1_from_mukesh);
        dE_dX_2_from_Paper_array.push_back(dE_dX_1_from_mukesh);
        //cout << "WIMP_mx_Array[kkk]: " << WIMP_mx_Array[kkk] << endl;
    }
    
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << dE_dX_1_from_mukesh_array[kkk] << "," << endl;
        cout << WIMP_mx_Array[kkk] << "," << dE_dX_2_from_Paper_array[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Cross_Section_Set[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Collision_Time_Array[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Energy_Loss_Per_Collision[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;

    //==================================================//

}//End_Main
*/
     


