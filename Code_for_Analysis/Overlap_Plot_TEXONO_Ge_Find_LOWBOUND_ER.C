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
    return 0.5*(Mx*1e6)*(v*1e3/3e8)*(v*1e3/3e8);
}
void Overlap_Plot_TEXONO_Ge_Find_LOWBOUND_ER()
{
    const double WIMP_mx=1;//GeV
    const int    velocity_N=13;
    const double dr_mukesh_c_to_kms = 1e-3*(3e8/1e3);
    string velocity_s[velocity_N]={"01667","03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};
    double velocity_d[velocity_N]={0.1667,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};
    double velocitykm[velocity_N];//Filled in the next line
    for(int kkk=0; kkk<velocity_N; kkk++){velocitykm[kkk]=velocity_d[kkk]*dr_mukesh_c_to_kms;}//(km/s)}
    
    
    //Read the file of DCS for different masses
    TFile *fin = TFile::Open("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/SI_c1_XeData_Vel/c1_1_0GeV/DCS.root");
    vector<TH1F*> velocity_TH1F;//Mass-based
    vector<double>   velocity_Used;//km/s
    velocity_Used.push_back(0);
    for(int kkk=0; kkk<velocity_N; kkk++)
    {
        TH1F *Velocity_TH1F_temp=(TH1F*)fin->Get(velocity_s[kkk].c_str());
        if(Velocity_TH1F_temp!=NULL)
        {
            cout << "File: " << velocity_s[kkk].c_str() << endl;
            cout << "velocitykm: " << velocitykm[kkk] << endl;
            velocity_TH1F.push_back(Velocity_TH1F_temp);
            velocity_Used.push_back(velocitykm[kkk]);
        }
    }
    
    int    Applied_Hist=0;
    double Velocity_DM_Now=799;//km/s
    if(Velocity_DM_Now>velocity_Used[velocity_Used.size()-1])Applied_Hist=(velocity_TH1F.size())-1;
    if(Applied_Hist==0)
    {
        for(int kkk=0; kkk<velocity_Used.size(); kkk++)
        {
            if(Velocity_DM_Now>velocity_Used[kkk] and Velocity_DM_Now<velocity_Used[kkk+1]){Applied_Hist=kkk;}
        }
    }
    cout << "Applied_Hist: " << Applied_Hist << endl;
    double v_c = Velocity_DM_Now*1e3/3e8;//c
    double Max_Recoil = DM_E(WIMP_mx,Velocity_DM_Now);//Max_Recoil(keV) given by mass
    Int_t Maximum_Bin = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Max_Recoil*1e3);//Max_Recoil_Bin
    cout << "Maximum_Bin: " << Maximum_Bin << endl;
    int LastBin = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
    cout << "LastBin: " << LastBin << endl;

    double Total_Cross_Section=0;
    double c1_GeV2 = 1;
    double Scaling = 3e3;
    double c1_eV2 = c1_GeV2*1e-18*Scaling;

    for(int kkk=0; kkk<Maximum_Bin; kkk++)
    {
        double dsigma_dT  = velocity_TH1F[Applied_Hist]->GetBinContent(kkk)*1e-15*c1_eV2*c1_eV2;
        double dT         = velocity_TH1F[Applied_Hist]->GetBinWidth(kkk)*1e-3;
        //cout << "dsigma_dT: " << dsigma_dT << endl;
        //cout << "dT: " << dT << endl;
        Total_Cross_Section = Total_Cross_Section + (dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
        //cout << "Total_Cross_Section: " << Total_Cross_Section << endl;
    }
    cout << "Total_Cross_Section: " << Total_Cross_Section << endl;
    cout << "CS_Try: " << CS_Try(c1_GeV2*Scaling,0.5);
    int N = (1.8*1*1e5*Total_Cross_Section)/(unified_atomic_mass_g*(28));
    cout << "Count: " << N << endl;
    //=======================================Start interacting=======================================//
    
    int Event_Number=0;
    double collision_Time_Total=0;
    while(Event_Number<100)
    {
        cout << "Event_Number: " << Event_Number << endl;
        double Velocity_DM_Temp = 779;//Initial Energy and iteration =>km/s
        double Energy_DM_Temp   = DM_E(WIMP_mx,779);//Initial Energy and iteration =>keV
        int collision_Time=0;
        //for(int kkk=0; kkk<1; kkk++)
        while(Energy_DM_Temp>0.16)
        {
            //Confirm the TH1F
            gRandom = new TRandom3(0);
            gRandom->SetSeed(0);
            
            //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;
            //cout << "Energy_DM_Temp: " << Energy_DM_Temp << endl;
            //Find the file
            if(Velocity_DM_Temp>velocity_Used[velocity_Used.size()-1])Applied_Hist=(velocity_TH1F.size())-1;
            if(Applied_Hist==0)
            {
                for(int kkk=0; kkk<velocity_Used.size(); kkk++)
                {
                    if(Velocity_DM_Temp>velocity_Used[kkk] and Velocity_DM_Temp<velocity_Used[kkk+1]){Applied_Hist=kkk;}
                }
            }
            //Set the range
            int Maximum_Bin_Loss = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Energy_DM_Temp*1e3);//Max_Recoil_Bin
            int LastBin_Loss = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_from_Hist
            if(Maximum_Bin_Loss>LastBin_Loss)Maximum_Bin_Loss=LastBin_Loss;
            double  Max_X        = velocity_TH1F[Applied_Hist]->GetBinCenter(Maximum_Bin_Loss);
            velocity_TH1F[Applied_Hist]->GetXaxis()->SetRange(0,Max_X);
            double Energy_Loss_eV = velocity_TH1F[Applied_Hist]->GetRandom();//Energy_Loss(eV)
            double Energy_Loss_keV = Energy_Loss_eV*1e-3;//Energy_Loss(keV)
            //cout << "Energy_Loss_(keV): " << Energy_Loss_keV << endl;
            Energy_DM_Temp   = Energy_DM_Temp - Energy_Loss_keV;
            //cout << "Energy_DM_Temp(keV): " << Energy_DM_Temp << endl;
            Velocity_DM_Temp = sqrt(2*Energy_DM_Temp/(WIMP_mx*1e6))*(3e8/1e3);//km/s
            //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;//km/s
            collision_Time = collision_Time + 1;
        }
        collision_Time_Total = collision_Time_Total + collision_Time;
        Event_Number = Event_Number + 1;
    }
    collision_Time_Total = collision_Time_Total / Event_Number;
    cout << "collision_Time_Ave: " << collision_Time_Total << endl;

}//End_Main
    
     

