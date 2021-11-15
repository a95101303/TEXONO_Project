const int data_bin = 156;
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_03333V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_05000V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_06667V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_08333V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_10000V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_11670V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_13330V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_15000V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_16670V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_18330V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_20000V.h"
#include "SI_c1_XeData_Vel/c1_1_0GeV/DM_1_0GeV_21670V.h"
//#include "c1_1_0GeV/DM_1_0GeV_25460V.h"



double sub_MB_dist(double dm_vel)
{
  double MB_factor;
  double V0    = 220E+5; //cm sec^{-1}
  double Ve    = 232E+5; //cm sec^{-1} 
  return MB_factor = dm_vel*exp(-pow((dm_vel-Ve)/V0,2));

 }

double sub_MB_dist_Test(double dm_vel, double dm_vel_0)
{
  double MB_factor;
  double V0    = 220E+5; //cm sec^{-1}
  double Ve    = 232E+5; //cm sec^{-1}
  return MB_factor = dm_vel*exp(-pow((dm_vel-Ve)/V0,2))*(dm_vel-dm_vel_0);

 }

double add_MB_dist(double dm_vel)
{
  double MB_factor;
  double V0    = 220E+5; //cm sec^{-1}
  double Ve    = 232E+5; //cm sec^{-1} 
  return MB_factor = dm_vel*exp(-pow((dm_vel+Ve)/V0,2));

 }


void dm_dcsXe_1GeV_avDCS_Chih()
{
  
   
  const double PI = 3.141592653589793238;
  const double AGe = 131.293;
  const double avogadro_number = 6.02214129e+23; 
  const double N_atom_1kg_Ge = (1000.0*avogadro_number)/AGe;
  
  double density = 0.4; //GeV cm^{-3}
  double electron_number = N_atom_1kg_Ge; //kg^{-1}
  
  double dm_mass = 1.0; //GeV
 
  const double percm_GeV = 1.97326971780039025e-14;
  float  me = 0.510998928E-3;
  float red_mass = ((me*dm_mass)/(me+dm_mass));
  //float alpha = (1E-37*TMath::Pi())/(pow(red_mass,2)*pow(percm_GeV,2));

  double V0    = 220E+5; //cm sec^{-1} 
  double Vesc  = 544E+5; //cm sec^{-1} 
  double Ve    = 232E+5; //cm sec^{-1} 
  double c1 = 1e-18;
  
  
  const double rate_factor = pow(c1,2)*1E-18*1E+3*86400;
  const int bin = 13; 
  double vel_dist[bin] = {
    1.667e-04,
    3.333e-04,
    5.000e-04,
    6.667e-04,
    8.333e-04,
    1.000e-03,
    1.167e-03,
    1.333e-03,
    1.500e-03,
    1.667e-03,
    1.833e-03,
    2.000e-03,
    2.167e-03,
    //2.546e-03
  }; 

  double norm_factor =0;
  
  for(int ll = 0; ll<bin;ll++)
    {
      vel_dist[ll] =  vel_dist[ll]*3E+10; //cm sec^{-1}
      if(ll==0)
	{
	  norm_factor += vel_dist[ll]*exp(-pow((vel_dist[ll]-Ve)/V0,2))*vel_dist[ll];
	}
      else
	{
	  norm_factor += vel_dist[ll]*exp(-pow((vel_dist[ll]-Ve)/V0,2))*(vel_dist[ll]-vel_dist[ll-1]);  
	}    
}
  norm_factor = 1.0/norm_factor;
  printf(" Value   %e \n", norm_factor);  

  
  int number_electron =0;
  double temp =0;
  
  double rate[data_bin], energy[data_bin];
  for(int rr = 0; rr < data_bin-1;  rr++)
    {
      energy[rr] = DM_1_0GeV_10000V[rr][0]*1E-3;

      rate[rr] =  (rate_factor*norm_factor*                  
		   ( vel_dist[1]*  DM_1_0GeV_03333V[rr][1]*sub_MB_dist(vel_dist[1])*(vel_dist[1]-vel_dist[0])   +
		     vel_dist[2]*  DM_1_0GeV_06667V[rr][1]*sub_MB_dist(vel_dist[2])*(vel_dist[2]-vel_dist[1])   +
		     vel_dist[3]*  DM_1_0GeV_06667V[rr][1]*sub_MB_dist(vel_dist[3])*(vel_dist[3]-vel_dist[2])   +
		     vel_dist[4]*  DM_1_0GeV_08333V[rr][1]*sub_MB_dist(vel_dist[4])*(vel_dist[4]-vel_dist[3])   +
		     vel_dist[5]*  DM_1_0GeV_10000V[rr][1]*sub_MB_dist(vel_dist[5])*(vel_dist[5]-vel_dist[4])   +
		     vel_dist[6]*  DM_1_0GeV_11670V[rr][1]*sub_MB_dist(vel_dist[6])*(vel_dist[6]-vel_dist[5])   +
		     vel_dist[7]*  DM_1_0GeV_13330V[rr][1]*sub_MB_dist(vel_dist[7])*(vel_dist[7]-vel_dist[6])   +
		     vel_dist[8]*  DM_1_0GeV_15000V[rr][1]*sub_MB_dist(vel_dist[8])*(vel_dist[8]-vel_dist[7])   +
		     vel_dist[9]*  DM_1_0GeV_16670V[rr][1]*sub_MB_dist(vel_dist[9])*(vel_dist[9]-vel_dist[8])   +
		     vel_dist[10]* DM_1_0GeV_18330V[rr][1]*sub_MB_dist(vel_dist[10])*(vel_dist[10]-vel_dist[9]) +
		     vel_dist[11]* DM_1_0GeV_20000V[rr][1]*sub_MB_dist(vel_dist[11])*(vel_dist[11]-vel_dist[10]) +
		     vel_dist[12]* DM_1_0GeV_21670V[rr][1]*sub_MB_dist(vel_dist[12])*(vel_dist[12]-vel_dist[11])
		     //vel_dist[13]* DM_1_0GeV_25460V[rr][1]*sub_MB_dist(vel_dist[13])*(vel_dist[13]-vel_dist[12]) 
		     ));

    }
    cout << "norm_factor: " << norm_factor << endl;

    for(int kkk=0; kkk<12; kkk++)
    {
        cout << "sub_MB_dist_Test(vel_dist[x]): "  << sub_MB_dist_Test(vel_dist[kkk+1],vel_dist[kkk])*norm_factor << endl;
    }

  
  TCanvas *plot = new TCanvas("plot","",800,600);
 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFrameLineWidth(2);	
  gStyle->SetLabelOffset(0.00,"Y");
  gStyle->SetLabelOffset(-0.01,"X");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetNdivisions(510,"Y");
  gStyle->SetLabelSize(0.04, "X");
  gStyle->SetNdivisions(510,"X");
  gStyle->SetTitleSize(0.045, "X" );
  gStyle->SetTitleSize(0.05, "Y" );
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetLabelFont(22,"X");
  gStyle->SetLabelFont(22,"Y");
  gStyle->SetTitleOffset(0.85,"X");
  gStyle->SetTitleOffset(0.95,"Y");


  gPad->SetLogy(1);
  gPad->SetLogx(1);
 

 
  
    TH2F *frame = new TH2F("frame","",10,1E-2,10,10,1E-30,1E-21);
  frame->GetYaxis()->SetTitle("#frac{d<#sigma v>}{dT} cm^{3} keV^{-1} Day^{-1}");
  frame->GetXaxis()->SetTitle(" T (keV)");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  //frame->GetXaxis()->SetMoreLogLabels();
  frame->Draw();

  
  TGraph *gr = new TGraph(data_bin-1,energy, rate);
  gr->SetLineColor(1);
  gr->SetLineWidth(10);
  gr->GetXaxis()->SetRangeUser(1e-2,2.8);
  gr->Draw("Al");
  
}









