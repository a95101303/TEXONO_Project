const int DataBin = 255;

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_06GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_07GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_08GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_09GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_10GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_12GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_20GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_25GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_30GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_40GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_50GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_60GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_70GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_80GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_90GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_1_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_1_3GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_1_5GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_2_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_2_5GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_3_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_3_5GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_4_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_5_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_10_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_20_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_LR_DCS_June.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_LR_DCS_Dec.h"


//These parameters are for calculating the sigma_e(Unit transformation!)
const double hbar_eV_s = 6.582119569E-16;//(eV*s)
const double Light_Speed = 2.99792458E8;//(m/s)
const double GeV_for_transform = 1E9;//(eV)
const double hbarC_divide_GeV_m  = hbar_eV_s*Light_Speed/(GeV_for_transform);//(m)
const double hbarC_divide_GeV_cm = 1e2*hbarC_divide_GeV_m;//(cm)
const double hbarC_divide_GeV_cm_Square = pow(hbarC_divide_GeV_cm,2);//(cm^2)

//Other parameters that should be used
const double avogadro_number = 6.02214129e+23;
const double density = 0.3; //GeV cm^{-3} for (DM)
const double Averaged_Velocity=232.;//(km/s)
const double Electron_Mass_MeV=0.511;//MeV
const double Me_times_alpha=3.7;//keV/c^2

const double V_to_C = 1e3/3e8;//km/s to beta
const double MeV_to_GeV = 1E-3;//
const double GeV_to_keV = 1E+6;//

double Final_Energy(int Times, double Initial_Energy)
{
    double A=Initial_Energy; double B;
    for(int kkk=0; kkk<Times ; kkk++)
    {
        B = A - (A*0.5);
        A = B;
    }
    return A;
}

double RMe(double mx)//mx(GeV/c^2),Reduce_Mass_e(RMe)
{
    //cout << "RMe: " << 1e-3*Electron_Mass_MeV*1e-3*1000*mx/((1e-3*Electron_Mass_MeV)+1e-3*1000*mx) << endl;
    return 1e-3*Electron_Mass_MeV*mx/((1e-3*Electron_Mass_MeV)+mx);//GeV/c^2
}
double max_recoil_A_for_ER_keV(double velocity, double mx)//mx(GeV/c^2),Velocity(km/s)
{
  double max_recoil_A_0 = 0.5*mx*1e6*(velocity*V_to_C)*(velocity*V_to_C); //(keV)
    cout << "max_recoil_A_0: " << max_recoil_A_0 << endl;
  return max_recoil_A_0;
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

/*
double dsigma_dT_keV_ER(double CS, double Month, double T, double mx, double A)//dsigma_dT_keV for Electron Recoil, Month1: June, Month2: Dec
{
    double rate_factor = pow(d_1(CS,mx),2)*1E-18*1E+3*86400*3E+10*(density*electron_number(A))/(mx);
    if(Month==1)return rate_factor*LongRangeDcs_June(T*1E+3 , mx )*1e-29;
    if(Month==2)return rate_factor*LongRangeDcs_Dec(T*1E+3 , mx )*1e-29;
}
*/
double Max_Recoil_A_keV_ER(double V, double mx)//Electron Recoil Max recoil energy, mx(GeV/c^2), velocity(V)
{
    double Beta = (V*1e3)/(3E+8);//Beta(c)
    return 0.5*mx*(1e6)*(Beta)*(Beta);//Energy of DM basically
    //return 0.5*RMe(mx)*(1e6)*(Beta)*(Beta);//Energy of DM basically

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
//================From paper 1905.06348v2===============//
/*
double Max_Recoil_A_keV_ER_Free(double mx, double V_Initial)//Electron Recoil Free electron, Max recoil energy, mx(GeV/c^2), velocity(V)
{
    //cout << "V_Initial: " << V_Initial << endl;
    double V_Final = (V_Initial)*(2*mx/(mx+Electron_Mass_MeV*MeV_to_GeV));
    //cout << "V_Final" << V_Final << endl;
    double Beta    = (V_Final*1e3)/(3E+8);//Beta(c)
    //cout << "Beta: " << Beta << endl;
    //cout << "0.5*Electron_Mass_MeV*1e3*(Beta)*(Beta): " << 0.5*Electron_Mass_MeV*1e3*(Beta)*(Beta) << endl;
    return 0.5*Electron_Mass_MeV*1e3*(Beta)*(Beta);//Energy of DM basically
}
 */

const double Alpha_FS = 1./137.;//Fine Sturcture constant
const double q_ref    = Alpha_FS*0.511*1e3;//keV/c^2

double Momentum_Transfer(double T)//T(keV)
{
    return sqrt(2*Electron_Mass_MeV*1e3*T);//(keV/c)
}
double F_DM(double T, double MA_Plus)//T(keV), MA_Plus(Mass of Dark Photon), ER_Square
{
    cout << "T: " << T << endl;
    double  q = Momentum_Transfer(T);
    cout << "q: " << q << endl;
    return (q_ref*q_ref+MA_Plus*MA_Plus)/(q*q+MA_Plus*MA_Plus);//No Unit
}
double dsigma_dT_keV_ERSS(double Sigma_SI, double mx, double V)//ER_Square
{
    //return Sigma_SI*(1./(4*RMe(mx)*RMe(mx)*V*V_to_C*V*V_to_C))*F_DM(mx,V)*F_DM(mx,V);
    return Sigma_SI*(1./(4*RMe(mx)*1e6*RMe(mx)*1e6*V*V_to_C*V*V_to_C));//(cm^2*c^2/keV^2)
}
double dsigma_dT_keV_ER(double Sigma_SI, double mx, double V, double T, double M_DP)//Sigma_SI(cm^2), mx(GeV/c^2), V(km/s), q(momentum from above), M_DP(keV/c^2)[Dark Photon Mass], ER
{
    //return Sigma_SI*(1./(4*RMe(mx)*RMe(mx)*V*V_to_C*V*V_to_C))*F_DM(mx,V)*F_DM(mx,V);
    return 2*Electron_Mass_MeV*1e3*(dsigma_dT_keV_ERSS(Sigma_SI,mx,V)*F_DM(T, M_DP)*F_DM(T, M_DP));//
}


double Total_Sigma_ER(double CS, double V, double mx)//Total Sigma for Electron Recoil, V(km/s)
{
    int reso_T=1000;double T[reso_T];double total_Sigma=0;
     double WIMP_max_T   = max_recoil_A_for_ER_keV(V,mx); //keV
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
        //cout << "T[i]: " << T[i] << endl;
        if(i==0) { dEx = T[0]; }
        else { dEx = T[i] - pEx; }
        total_Sigma = total_Sigma + (dsigma_dT_keV_ER(CS, mx, V, T[i], 1e6)*dEx);
        pEx = T[i];
    }
    cout << "mx: " << mx << endl;
    cout << "total_Sigma: " << total_Sigma << endl;
    return total_Sigma;
}

//=====================================================//
//Find the dsigma_dT
