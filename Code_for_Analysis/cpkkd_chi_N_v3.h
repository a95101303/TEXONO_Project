#include "math.h"

#ifndef INCLUDE_CONSTANT_H
#define INCLUDE_CONSTANT_H
#include "constant.h"
#endif

#ifndef INCLUDE_DM_PARAMETERS_H
#define INCLUDE_DM_PARAMETERS_H
//#include "dm_parameters.h"
#endif



double cpkkd_chi_N(double mx, // GeV
                   double T, // keV
                   double time, // days
                   double phase, // days, should default to 152.5 for HALO
                   double atomic_mass
                  )
{
const double sig_p0 = 0.0001; // the final dm-spec will be in 10^-40 cm^2 unit

double Energy = T/1e6;
double r, E0, cpd90, R0, sig_N;
double q, s, rrms, rn, qrn, F;
double sig_N_band;
double red_mGe, red_mp;
double first_factor, second_factor;
double k1k0;
double mt, mp;
double ve, vmin;

k1k0 = 1.0/0.9965;
mp = Mp/1000.0; //GeV
//mt= 67.56;
//mt = atom_mass_Ge_MeV/1000.0; //GeV
mt = atomic_mass*unified_atomic_mass_MeV/1000.0; //GeV


    // compare with proton   
    sig_N_band = pow(atomic_mass,2)*sig_p0*pow(((mx*mt)/(mx+mt)),2)/pow(((mx*mp)/(mx+mp)),2);
    // form factor
    q = 6.92*0.001*(sqrt(atomic_mass))*(sqrt(Energy*1.0e+6));
    s = 1.0;
    rrms = pow(atomic_mass,(1.0/3.0));
    rn = pow((5.0/3.0*pow(rrms,2)-5.0*pow(s,2)),0.5);
    qrn = q*rn;

//    if(Energy>0.01)
    { F = 3*(sin(qrn)/pow(qrn,2)-cos(qrn)/qrn)/qrn*exp(-0.5*pow((q*s),2)); }
//    else
//    { F = cos(qrn)*exp(-0.5*pow((q*s),2)); }

    sig_N = sig_N_band*pow(F,2);
    //
    r = 4.0*mx*mt/pow((mx+mt),2);
///////////
    red_mGe= mx*mt/(mx+mt);
    red_mp = mx*mp/(mx+mp);

    E0 = 0.5*mx*(pow((v0/3.0e5),2));

    R0 = pow(atomic_mass,2)*sig_p0*pow(((mx*mt)/(mx+mt)),2)/pow(((mx*mp)/(mx+mp)),2)*5.47*rohx*v0/(mx*mt);

    vmin = v0*(sqrt(Energy/(r*E0)));

    //ve = 232.0;
    ve = 232.0 + 15.0*cos(2*TMath::Pi()*(time-phase)/365.25);

    first_factor = (R0*sqrt(TMath::Pi())*v0)/(E0*r*4.0*ve)*(TMath::Erf((vmin+ve)/v0) - TMath::Erf((vmin-ve)/v0));
    second_factor = k1k0*(first_factor - (R0/(r*E0)*exp(-1*(vesc*vesc)/(v0*v0))));

    if(second_factor>0) { return second_factor*1.0e-6*pow(F,2); } 
                   else { return 0.0; }

}
