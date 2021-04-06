#include "math.h"
#ifndef INCLUDE_CONSTANT_H
#define INCLUDE_CONSTANT_H
#include "constant.h"
#endif

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

double min_v_at_max_E(double mx, double T, double atomic_mass)
{
  double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));
  double r_A = 4.0*pow(reduce_mass_A,2.0)/(mx*1000.0*unified_atomic_mass_MeV*atomic_mass);

  double velocity0 = sqrt(2.0*T/(mx*1000.0*r_A));

  return velocity0;
}


double fdsigma_dT(double mx, double sigma_SI, double velocity, double atomic_mass, double T)
{
  double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));
  double reduce_mass_n = mx*1000.0*unified_atomic_mass_MeV/((mx*1000.0)+unified_atomic_mass_MeV);
  double r_A = 4.0*pow(reduce_mass_A,2.0)/(mx*1000.0*unified_atomic_mass_MeV*atomic_mass);

  double s_fm = 0.9;//Thickness, fm
  double s = s_fm*fm_MeV1; //unit MeV^-1

  double q = sqrt(2.0*unified_atomic_mass_MeV*atomic_mass*T); //unit MeV
  double rn2_fm = pow(((1.23*pow(atomic_mass,(1.0/3.0)))-0.6),2.0)+(7.0/3.0)*pow((0.52*TMath::Pi()),2.0)-(5.0*s_fm*s_fm); //unit fm^2
  double rn = sqrt(rn2_fm)*fm_MeV1; //unit MeV^-1
  double form_factor = 3.0*((sin(q*rn)-(q*rn)*cos(q*rn))/pow((q*rn),3.0))*exp(-1.0*pow((q*s),2.0)/2.0);

  double sigma_chiN0 = sigma_SI*4.0*pow(reduce_mass_A,2.0)*pow(atomic_mass,2.0)/pow(reduce_mass_n,2.0);
  double dsigma_dq2 = sigma_chiN0*pow(form_factor,2.0)/(4.0*pow(reduce_mass_A,2.0)*pow(velocity,2.0));
  double dsigma_dT_0 = 2.0*unified_atomic_mass_MeV*atomic_mass*dsigma_dq2;

  return dsigma_dT_0;
}

double insert(double x, double x0, double y0, double x1, double y1)
{
  double y  = y1 - (x1-x)*(y1-y0)/(x1-x0);
  return y;
}
