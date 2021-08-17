#include "DataJune_ShortRange/DataJune_0_06GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_07GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_08GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_09GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_10GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_12GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_20GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_25GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_30GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_40GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_50GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_60GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_70GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_80GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_0_90GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_1_0GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_1_3GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_1_5GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_2_0GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_2_5GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_3_0GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_3_5GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_4_0GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_5_0GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_10_0GeV_c1_dcs.h"
#include "DataJune_ShortRange/DataJune_20_0GeV_c1_dcs.h"

//These parameters are for calculating the sigma_e(Unit transformation!)
const double hbar_eV_s = 6.582119569E-16;//(eV*s)
const double Light_Speed = 2.99792458E8;//(m/s)
const double GeV_for_transform = 1E9;//(eV)
const double hbarC_divide_GeV_m  = hbar_eV_s*Light_Speed/(GeV_for_transform);//(m)
const double hbarC_divide_GeV_cm = 1e2*hbarC_divide_GeV_m;//(cm)
const double hbarC_divide_GeV_cm_Square = pow(hbarC_divide_GeV_cm,2);//(cm^2)

//Other parameters that should be used
const double PI = 3.141592653589793238;
const double AGe = 72.64;
const double avogadro_number = 6.02214129e+23;
const double density = 0.3; //GeV cm^{-3} for (DM)
const double Averaged_Velocity=232.;//(km/s)
const double Electron_Mass_MeV=0.511;//MeV
const double Me_times_alpha=3.7;//keV/c^2

double RMe(double mx)//Reduce_Mass_e(RMe)
{
    //cout << "RMe: " << 1e-3*Electron_Mass_MeV*1e-3*1000*mx/((1e-3*Electron_Mass_MeV)+1e-3*1000*mx) << endl;
    return 1e-3*Electron_Mass_MeV*1e-3*1000*mx/((1e-3*Electron_Mass_MeV)+1e-3*1000*mx);//GeV
}

double c_1(double CS, double mx)//Cross-Section(CS) to c1
{
    return sqrt( ( CS*PI ) / ( RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square )  );
}
double d_1(double CS, double mx)//Cross-Section(CS) to d1
{
    return sqrt( (CS*PI*pow(Me_times_alpha,4)) /( RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square*1e24 ) );
}
double electron_number(double A = AGe)
{
    return (1000.0*avogadro_number)/(A);//kg^{-1}
}

double dsigma_dT_keV_ER(double CS, double Month, double T, double mx, double A)//dsigma_dT_keV for Electron Recoil, Month1: June, Month2: Dec
{
    double rate_factor = pow(d_1(CS,mx),2)*1E-18*1E+3*86400*3E+10*(density*electron_number(A))/(mx);

    if(Month==1)return rate_factor*LongRangeDcs_June(T*1E+3 , mx )*1e-29;
    if(Month==2)return rate_factor*LongRangeDcs_Dec(T*1E+3 , mx )*1e-29;
}

double Max_Recoil_A_keV_ER(double V, double mx)//Electron Recoil Max recoil energy, mx(GeV/c^2), velocity(V)
{
    double Beta = (V*1e3)/(3E+8);//Beta(c)
    return 0.5*mx*(1e6)*(Beta)*(Beta);//Energy of DM basically
}

double Total_Sigma_ER(double Month, double CS, double V, double mx, double A)//Total Sigma for Electron Recoil
{
    int reso_T=1000;double T[reso_T];double total_Sigma=0;
     double WIMP_max_T = Max_Recoil_A_keV_ER(V,mx); //keV
    //======
    for(int i=0;i<reso_T;i++)
    {
        T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
    }
    //======
    double dEx=0;
    double pEx = T[0];
    for(int i=0;i<reso_T;i++)
    {
        if(i==0) { dEx = T[0]; }
        else { dEx = T[i] - pEx; }
        total_Sigma = total_Sigma + (dsigma_dT_keV_ER(CS, 1, T[i], mx, A)*dEx);
        pEx = T[i];
    }
    return total_Sigma;
}

/*
double
{
    TF1 *fa2 = new TF1("fa2","dsigma_dT_keV_ER([0],x,[1])",95.0*1E-3,95.0*1E-3+5);
    fa2->SetParameter(0,2);
    fa2->SetParameter(1,mx);
}
*/

double CS_Try(double c_1, double mx)
{
    return c_1*c_1*RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square/(PI);
}
double DS_Try(double d_1, double mx)
{
    return d_1*d_1*RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square*1e24/(PI*pow(Me_times_alpha,4));
}

