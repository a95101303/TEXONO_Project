const int DataBin = 255;

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_LR_DCS_June.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_LR_DCS_Dec.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_SR_DCS_June.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_SR_DCS_Dec.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/constant.h"



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

double RMe(double mx)//Reduce_Mass_e(RMe)
{
    //cout << "RMe: " << 1e-3*Electron_Mass_MeV*1e-3*1000*mx/((1e-3*Electron_Mass_MeV)+1e-3*1000*mx) << endl;
    return 1e-3*Electron_Mass_MeV*1e-3*1000*mx/((1e-3*Electron_Mass_MeV)+1e-3*1000*mx);//GeV
}
double max_recoil_A_for_ER_keV(double velocity, double mx)//mx(GeV/c^2),Velocity(km/s)
{
  double max_recoil_A_0 = 0.5*RMe(mx)*1e6*(velocity*V_to_C)*(velocity*V_to_C); //(keV)
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
    //return 0.5*RMe(mx)*(1e6)*(Beta)*(Beta);//Energy of DM basically

}

double Total_Sigma_ER(double Month, double CS, double V, double mx, double A)//Total Sigma for Electron Recoil, V(km/s)
{
    int reso_T=1000;double T[reso_T];double total_Sigma=0;
     double WIMP_max_T   = Max_Recoil_A_keV_ER(V,mx); //keV
     double WIMP_max_T_2 = max_recoil_A_for_ER_keV(V,mx);//keV
    //cout << "WIMP_max_T: " << WIMP_max_T << endl;
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
        total_Sigma = total_Sigma + (dsigma_dT_keV_ER(CS, 1, T[i], mx, A)*dEx);
        //cout << "dsigma_dT_keV_ER(CS, 1, T[i], mx, A): " << dsigma_dT_keV_ER(CS, 1, T[i], mx, A) << endl;
        pEx = T[i];
    }
    //cout << "total_Sigma: " << total_Sigma << endl;
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


//mx: mass of dark matter, unit: GeV
//velocity: natural unit
//sigma_SI: in unit of cm^2
//T: recoil nuclear energy, in unit of MeV
double max_recoil_A(double mx, double velocity, double atomic_mass)
{
    double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));

    double r_A = 4.0*pow(reduce_mass_A,2.0)/(mx*1000.0*unified_atomic_mass_MeV*atomic_mass);
    
    double max_recoil_A_0 = 0.5*mx*1000.0*velocity*velocity*r_A; //MeV
    
    return max_recoil_A_0;
}

double max_recoil_A_keV(double mx, double velocity, double atomic_mass)
{ return 1000.0*max_recoil_A(mx, velocity, atomic_mass); }
//

///////////////////////////////////////////////////////////////////////////////
// no unit
// T: MeV
double F2(double atomic_mass, double T)
{
    double s_fm = 0.9;//Thickness, fm
    double s = s_fm*fm_MeV1; //unit MeV^-1
    
    double q = sqrt(2.0*unified_atomic_mass_MeV*atomic_mass*T); //unit MeV
    double rn2_fm = pow(((1.23*pow(atomic_mass,(1.0/3.0)))-0.6),2.0)+(7.0/3.0)*pow((0.52*TMath::Pi()),2.0)-(5.0*s_fm*s_fm); //unit fm^2
    double rn = sqrt(rn2_fm)*fm_MeV1; //unit MeV^-1
    double form_factor = 3.0*((sin(q*rn)-(q*rn)*cos(q*rn))/pow((q*rn),3.0))*exp(-1.0*pow((q*s),2.0)/2.0);
    
    return form_factor*form_factor;
}

//cm^2/MeV
//mx: GeV
//sigma_SI: cm^2
//velocity: ratio of speed of light
//T: MeV
double fdsigma_dT(double mx, double sigma_SI, double velocity, double atomic_mass, double T)
{
    double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));
    double reduce_mass_n = mx*1000.0*unified_atomic_mass_MeV/((mx*1000.0)+unified_atomic_mass_MeV);
    
    double sigma_chiN0 = sigma_SI*pow(reduce_mass_A,2.0)*pow(atomic_mass,2.0)/pow(reduce_mass_n,2.0);
    double dsigma_dq2 = sigma_chiN0*F2(atomic_mass, T)/(4.0*pow(reduce_mass_A,2.0)*pow(velocity,2.0));
    double dsigma_dT_0;
    
    
    if(max_recoil_A(mx, velocity, atomic_mass)>T)
    {
        dsigma_dT_0 = 2.0*unified_atomic_mass_MeV*atomic_mass*dsigma_dq2;
    }
    else {
        dsigma_dT_0 = 0.0; }
     
    
    return dsigma_dT_0;
}
// cm^2/keV
//mx: GeV
//sigma_SI: cm^2
//velocity: ratio of speed of light
//T: keV
double fdsigma_dT_keV(double mx=10, double sigma_SI=1e-40, double velocity=0, double atomic_mass=0, double T=0)
{
    cout << "(1.0/1000.0)*fdsigma_dT(mx, sigma_SI, velocity, atomic_mass, (T/1000.0)): " << (1.0/1000.0)*fdsigma_dT(mx, sigma_SI, velocity, atomic_mass, (T/1000.0)) << endl;
    cout << "T: " << T << endl;
    return (1.0/1000.0)*fdsigma_dT(mx, sigma_SI, velocity, atomic_mass, (T/1000.0));
}

double fdsigma_dT_keV_ER(double mx=10, double sigma_SI=1e-40, double velocity=0, double atomic_mass=0, double T=0)
{
    return (1/atomic_mass)*(1/atomic_mass)*(1./4.)*(0.5e-3)*(1./4.)*(0.5e-3)*(1.0/1000.0)*fdsigma_dT(mx, sigma_SI, velocity, atomic_mass, (T/1000.0));
}

