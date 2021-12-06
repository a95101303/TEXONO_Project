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
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()//Test the fitting
//Test the fitting line
{
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    
    const int Index=1;
    const double Density = 1.8;//g/cm^3
    const double Length  = 1e5;//cm
    //const double Length  = 3e0;//cm
    const double Threshold_keV = 0.1;//keV
    vector<string> File;
    vector<double> WIMP_mx_Array;
    vector<double> Collision_Time_Array;
    vector<double> Energy_Loss_Per_Collision;
    vector<double> Scaling={1e-18,1e-9};

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
    for(int kkk=0; kkk<1; kkk++)
    {
        
        TFile *fin = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str(),File[kkk].c_str()));
        vector<TH1F*>    velocity_TH1F;//Mass-based
        for(int LLL=0; LLL<velocity_N; LLL++)
        {
            TH1F *Velocity_TH1F_temp=(TH1F*)fin->Get(velocity_s[LLL].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                cout << "File: " << velocity_s[LLL].c_str() << endl;
                velocity_TH1F.push_back(Velocity_TH1F_temp);
            }
        }
        
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        //for(int Applied_Hist=0; Applied_Hist<velocity_N; Applied_Hist++)
        for(int Applied_Hist=velocity_TH1F.size()-1; Applied_Hist<velocity_TH1F.size(); Applied_Hist++)
        {
            vector<double> TGraph_Recoil_Energy;vector<double> TGraph_dsigma_dT;

            Int_t Minimum_Bin = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(12);//Min_Recoil_Bin_Xe
            Int_t Maximum_Bin = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
            for(int LLL=Minimum_Bin; LLL<Maximum_Bin+1; LLL++)
            {
                TGraph_Recoil_Energy.push_back( velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter(LLL) );
                TGraph_dsigma_dT.push_back( TMath::Log10(velocity_TH1F[Applied_Hist]->GetBinContent(LLL)*1e-15*TMath::Power(Scaling[Index],2) ) );
            }
            cout << "Minimum_Bin: " << Minimum_Bin << endl;
            cout << "Maximum_Bin: " << Maximum_Bin << endl;
            
            TGraph * g = new TGraph((int)TGraph_Recoil_Energy.size(), &TGraph_Recoil_Energy[0], &TGraph_dsigma_dT[0]);
            TF1 *fitting_Line = new TF1("fitting_Line","[0]*x+[1]",12,velocity_TH1F[Applied_Hist]->GetXaxis()->GetBinCenter( Maximum_Bin-1 ));
            g->Fit(fitting_Line,"R");
            Double_t par[9];
            fitting_Line->GetParameters(&par[0]);
            cout << "par[0](a): " << par[0] << endl;cout << "par[1](b): " << par[1] << endl;
            g->SetMarkerStyle(20);
            g->SetMarkerColor(2);
            g->Draw("ap");
            
        }
        
    }
    
    //==================================================//

}//End_Main
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
            //cout << "Applied_Hist: " << Applied_Hist << endl;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            //double dsigma_dT  = velocity_TH1F[Applied_Hist]->GetBinContent(kkk)*1e-15*TMath::Power(Scaling[Index],2)*TMath::Power(c_1(Cross_Section_for_Exp[0],WIMP_mx_Array[Index_File]),2);
            double dsigma_dT  = velocity_TH1F[Applied_Hist]->GetBinContent(kkk)*1e-15*TMath::Power(Scaling[Index],2);
            double dT         = velocity_TH1F[Applied_Hist]->GetBinWidth(kkk)*1e-3;
            double T          = velocity_TH1F[Applied_Hist]->GetBinCenter(kkk)*1e-3;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            //cout << "1e-15*TMath::Power(Scaling[Index],2): " << 1e-15*TMath::Power(Scaling[Index],2) << endl;
            //cout << "dsigma_dT: " << dsigma_dT << endl;
            //cout << "dT: " << dT << endl;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinCenter(kkk): " << velocity_TH1F[Applied_Hist]->GetBinCenter(kkk) << endl;

            if(velocity_TH1F[Applied_Hist]->GetBinCenter(kkk)>12)Total_Count = Total_Count + (dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
            if(velocity_TH1F[Applied_Hist]->GetBinCenter(kkk)>12)dE_dX_1_from_mukesh = dE_dX_1_from_mukesh + (dsigma_dT)*(dT)*T;//
            //cout << "Total_Count: " << Total_Count << endl;
        }
        dE_dX_1_from_mukesh = dE_dX_1_from_mukesh*N_atom_1kg_Ge_Electron;
        
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
        
        cout << "CS_Try(1*c1_GeV2,WIMP_mx_Array[Index_File]): " << CS_Try(1*c1_GeV2,WIMP_mx_Array[Index_File]) << endl;
        cout << "DS_Try(1e-9*d1,WIMP_mx_Array[Index_File]): " << DS_Try(1e-9*d1,WIMP_mx_Array[Index_File]) << endl;
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
    /*
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
     */
}//End_Main


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
     


