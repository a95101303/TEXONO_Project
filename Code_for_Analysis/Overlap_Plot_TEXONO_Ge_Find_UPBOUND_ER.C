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
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_ER()
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
    vector<double> Cross_Section_Set;
    
    
    string c1_d1_Xe_Ge_index[4]={"SI_c1_XeData_Vel","SI_d1_XeData_Vel","SI_c1_GeData_Vel","SI_d1_GeData_Vel"};
    //===============electron recoil===================//
    for(int kkk=0; kkk<File.size(); kkk++)
    //for(int kkk=2; kkk<3; kkk++)
    {
        TFile *fin = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%s/%s/DCS.root",c1_d1_Xe_Ge_index[Index].c_str(),File[kkk].c_str()));
        vector<TH1F*> velocity_TH1F;//Mass-based
        vector<double>   velocity_Used;//km/s
        velocity_Used.push_back(0);
        for(int kkk=0; kkk<velocity_N; kkk++)
        {
            TH1F *Velocity_TH1F_temp=(TH1F*)fin->Get(velocity_s[kkk].c_str());
            if(Velocity_TH1F_temp!=NULL)
            {
                //cout << "File: " << velocity_s[kkk].c_str() << endl;
                //cout << "velocitykm: " << velocitykm[kkk] << endl;
                velocity_TH1F.push_back(Velocity_TH1F_temp);
                velocity_Used.push_back(velocitykm[kkk]);
            }
        }
        
        int    Applied_Hist=0;
        double Velocity_DM_Now=779;//km/s
        if(Velocity_DM_Now>velocity_Used[velocity_Used.size()-1])Applied_Hist=(velocity_TH1F.size())-1;
        if(Applied_Hist==0)
        {
            for(int kkk=0; kkk<velocity_Used.size(); kkk++)
            {
                if(Velocity_DM_Now>velocity_Used[kkk] and Velocity_DM_Now<velocity_Used[kkk+1]){Applied_Hist=kkk;}
            }
        }
        //cout << "Applied_Hist: " << Applied_Hist << endl;
        double v_c = Velocity_DM_Now*1e3/3e8;//c
        double Max_Recoil = DM_E(WIMP_mx_Array[kkk],Velocity_DM_Now);//Max_Recoil(keV) given by mass
        Int_t Maximum_Bin = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Max_Recoil*1e3);//Max_Recoil_Bin
        cout << "Maximum_Bin: " << Maximum_Bin << endl;
        int LastBin = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_Bin
        cout << "LastBin: " << LastBin << endl;

        /*
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
        cout << "CS_Try: " << CS_Try(c1_GeV2*Scaling,1);
        int N = (1.8*1*1e5*Total_Cross_Section)/(unified_atomic_mass_g*(28));
        cout << "Count: " << N << endl;
         */
        //=======================================Start interacting=======================================//
        
        int Event_Number=0;
        double collision_Time_Total=0;
        double Energy_Loss_eV_collision=0;
        while(Event_Number<50)
        {
            //cout << "=====Event_Number: " << Event_Number << endl;
            double Velocity_DM_Temp = 779;//Initial Energy and iteration =>km/s
            double Energy_DM_Temp   = DM_E(WIMP_mx_Array[kkk],779);//Initial Energy and iteration =>keV
            int collision_Time=0;
            double Energy_Loss_Inidivual=0;
            //for(int kkk=0; kkk<1; kkk++)
            while(Energy_DM_Temp>0.01)
            {
                //Confirm the TH1F
                gRandom = new TRandom3(0);
                gRandom->SetSeed(0);
                
                //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;
                //cout << "Energy_DM_Temp: " << Energy_DM_Temp << endl;
                //Find the file
                if(Velocity_DM_Temp>velocity_Used[velocity_Used.size()-1])
                {
                  //  cout << "vel>650: " << endl;
                  //  cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;
                  //  cout << "velocity_Used[velocity_Used.size()-1]: " << velocity_Used[velocity_Used.size()-1] << endl;
                    Applied_Hist=(velocity_TH1F.size())-1;
                  //  cout << "Applied_Hist: " << Applied_Hist << endl;
                }
                if(Applied_Hist==0)
                {
                    for(int kkk=0; kkk<velocity_Used.size(); kkk++)
                    {
                        if(Velocity_DM_Temp>velocity_Used[kkk] and Velocity_DM_Temp<velocity_Used[kkk+1])
                        {
                            //cout << "velocity_Used[kkk]: " << velocity_Used[kkk] << endl;
                            //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;
                            //cout << "velocity_Used[kkk+1]: " << velocity_Used[kkk+1] << endl;
                            Applied_Hist=kkk;
                            //cout << "Applied_Hist: " << Applied_Hist << endl;
                        }
                    }
                }
                //Set the range
                int Maximum_Bin_Loss = velocity_TH1F[Applied_Hist]->GetXaxis()->FindBin(Energy_DM_Temp*1e3);//Max_Recoil_Bin
                //cout << "Applied_Hist: " << Applied_Hist << endl;
                int LastBin_Loss = velocity_TH1F[Applied_Hist]->FindLastBinAbove();//Max_Recoil_from_Hist
                if(Maximum_Bin_Loss>LastBin_Loss)Maximum_Bin_Loss=LastBin_Loss;
                double  Max_X        = velocity_TH1F[Applied_Hist]->GetBinCenter(Maximum_Bin_Loss);
                velocity_TH1F[Applied_Hist]->GetXaxis()->SetRange(0,Max_X);
                double Energy_Loss_eV = velocity_TH1F[Applied_Hist]->GetRandom();//Energy_Loss(eV)
                double Energy_Loss_keV = Energy_Loss_eV*1e-3;//Energy_Loss(keV)
                //cout << "Energy_Loss_(keV): " << Energy_Loss_keV << endl;
                Energy_DM_Temp   = Energy_DM_Temp - Energy_Loss_keV;
                //cout << "Energy_DM_Temp(keV): " << Energy_DM_Temp << endl;
                Velocity_DM_Temp = sqrt(2*Energy_DM_Temp/(WIMP_mx_Array[kkk]*1e6))*(3e8/1e3);//km/s
                //cout << "Velocity_DM_Temp: " << Velocity_DM_Temp << endl;//km/s
                Energy_Loss_Inidivual = Energy_Loss_Inidivual + Energy_Loss_keV;
                collision_Time = collision_Time + 1;
                Applied_Hist =0;
            }
            //if(collision_Time>0)cout << "Energy_Loss_Inidivual/collision_Time: " << Energy_Loss_Inidivual/collision_Time << endl;
            Energy_Loss_eV_collision = Energy_Loss_eV_collision + Energy_Loss_Inidivual/collision_Time;
            collision_Time_Total = collision_Time_Total + collision_Time;
            Event_Number = Event_Number + 1;
        }
        double collision_Time_Ave = collision_Time_Total / Event_Number;
        Collision_Time_Array.push_back(collision_Time_Ave);//Count for every DM
        cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        double Energy_Loss_Real = Energy_DM(WIMP_mx_Array[kkk],779*(1e3/3e8))-Threshold_keV;
        
        Energy_Loss_eV_collision = (1./Energy_Loss_Real) * (Energy_Loss_eV_collision/Event_Number) ;
        Energy_Loss_Per_Collision.push_back(Energy_Loss_eV_collision);
        //cout << "Threshold: " << sqrt(2*0.16/(WIMP_mx_Array[kkk]*1e6))*(3e8/1e3);
        //cout << "DM_E(WIMP_mx_Array[kkk],779) : " << DM_E(WIMP_mx_Array[kkk],779) << endl;;
        //cout << "Energy_Loss_eV_collision: " << Energy_Loss_eV_collision << endl;
        //==================================================//
        
        //double Target_Total_Cross_Section = collision_Time_Ave*(unified_atomic_mass_g*(ASi))/(Density*Length);
        
        double Total_Count;
        
        for(int kkk=1; kkk<Maximum_Bin-1; kkk++)
        {
            //cout << "kkk:: " << kkk << endl;
            double dsigma_dT  = velocity_TH1F[Applied_Hist]->GetBinContent(kkk)*1e-15*TMath::Power(Scaling[Index],2);
            double dT         = velocity_TH1F[Applied_Hist]->GetBinWidth(kkk)*1e-3;
            //cout << "velocity_TH1F[Applied_Hist]->GetBinContent(kkk): " << velocity_TH1F[Applied_Hist]->GetBinContent(kkk) << endl;
            //cout << "1e-15*TMath::Power(Scaling[Index],2): " << 1e-15*TMath::Power(Scaling[Index],2) << endl;
            //cout << "dsigma_dT: " << dsigma_dT << endl;
            //cout << "dT: " << dT << endl;
            Total_Count = Total_Count + Length*(Density/((unified_atomic_mass_g*(ASi))))*(dsigma_dT)*(dT);//(cm^2/keV)*(keV)=cm^2
            //cout << "Total_Count: " << Total_Count << endl;
        }
        //cout << "velocity_TH1F[Applied_Hist]->GetMean(): " << velocity_TH1F[Applied_Hist]->GetMean() << endl;
        //cout << "collision_Time_Ave: " << collision_Time_Ave << endl;
        //cout << "Total_Count: " << Total_Count << endl;
        
        double Scaling = (collision_Time_Ave)/(Total_Count);
        //cout << "Scaling: " << Scaling << endl;
        double c1_GeV2 = sqrt(Scaling);
        double d1      = sqrt(Scaling);
        if(Index==0)Cross_Section_Set.push_back( CS_Try(1*c1_GeV2,WIMP_mx_Array[kkk]) );
        if(Index==1)Cross_Section_Set.push_back( DS_Try(1e-9*d1,WIMP_mx_Array[kkk]) );

        //cout << "WIMP_mx_Array[kkk]: " << WIMP_mx_Array[kkk] << endl;
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
    
     


