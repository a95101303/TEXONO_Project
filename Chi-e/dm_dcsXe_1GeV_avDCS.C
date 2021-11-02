const int data_bin = 156;
#include "Xe_c1_1_0GeV/DM_1_0GeV_03333V.h"
#include "Xe_c1_1_0GeV/DM_1_0GeV_05000V.h"
#include "Xe_c1_1_0GeV/DM_1_0GeV_06667V.h"
#include "Xe_c1_1_0GeV/DM_1_0GeV_08333V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_10000V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_11670V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_13330V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_15000V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_16670V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_18330V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_20000V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_21670V.h" 
#include "Xe_c1_1_0GeV/DM_1_0GeV_25460V.h"



double sub_MB_dist(double dm_vel)
{
  double MB_factor;
  double V0    = 220E+5; //cm sec^{-1}
  double Ve    = 232E+5; //cm sec^{-1} 
  return MB_factor = dm_vel*exp(-pow((dm_vel-Ve)/V0,2));

 }

double add_MB_dist(double dm_vel)
{
  double MB_factor;
  double V0    = 220E+5; //cm sec^{-1}
  double Ve    = 232E+5; //cm sec^{-1} 
  return MB_factor = dm_vel*exp(-pow((dm_vel+Ve)/V0,2));

 }






void dm_dcsXe_1GeV_avDCS()
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
  float alpha = (1E-37*TMath::Pi())/(pow(red_mass,2)*pow(percm_GeV,2));



  

  double V0    = 220E+5; //cm sec^{-1} 
  double Vesc  = 544E+5; //cm sec^{-1} 
  double Ve    = 232E+5; //cm sec^{-1} 
  double c1 = 1E-18;
  
  double k = pow((TMath::Pi()*pow(V0,2)),3.0/2.0)*(TMath::Erf(Vesc/V0) -(2.0/sqrt(TMath::Pi()))*(Vesc/V0)*exp(-1*pow(Vesc/V0,2)));

  double kf = TMath::Pi()*pow(V0,2)/(k*Ve);
  
  printf("%e   \n ", kf); 
  
  const double rate_factor = pow(c1,2)*1E-18*1E+3*86400.0*alpha;
  const int bin = 14; 
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
    2.546e-03
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
		     vel_dist[12]* DM_1_0GeV_21670V[rr][1]*sub_MB_dist(vel_dist[12])*(vel_dist[12]-vel_dist[11]) +
		     vel_dist[13]* DM_1_0GeV_25460V[rr][1]*sub_MB_dist(vel_dist[13])*(vel_dist[13]-vel_dist[12]) 
		    ));
      

      //rate[rr] = log10(rate[rr]); 
      
     
      //printf("%f  %e   \n",energy[rr], rate[rr] );
    }



 double ratef[data_bin], energyf[data_bin];
  for(int rr = 0; rr < data_bin-1;  rr++)
    {
      energyf[rr] = DM_1_0GeV_10000V[rr][0]*1E-3;

      ratef[rr] =  (rate_factor*kf*(                  
				    (DM_1_0GeV_03333V[rr][1]*sub_MB_dist(vel_dist[1])*(vel_dist[1]-vel_dist[0])   +
				     DM_1_0GeV_06667V[rr][1]*sub_MB_dist(vel_dist[2])*(vel_dist[2]-vel_dist[1])   +
				     DM_1_0GeV_06667V[rr][1]*sub_MB_dist(vel_dist[3])*(vel_dist[3]-vel_dist[2])   +
				     DM_1_0GeV_08333V[rr][1]*sub_MB_dist(vel_dist[4])*(vel_dist[4]-vel_dist[3])   +
				     DM_1_0GeV_10000V[rr][1]*sub_MB_dist(vel_dist[5])*(vel_dist[5]-vel_dist[4])   +	  
				     DM_1_0GeV_11670V[rr][1]*sub_MB_dist(vel_dist[6])*(vel_dist[6]-vel_dist[5])   +	  
				     DM_1_0GeV_13330V[rr][1]*sub_MB_dist(vel_dist[7])*(vel_dist[7]-vel_dist[6])   +	  
				     DM_1_0GeV_15000V[rr][1]*sub_MB_dist(vel_dist[8])*(vel_dist[8]-vel_dist[7])   +	  
				     DM_1_0GeV_16670V[rr][1]*sub_MB_dist(vel_dist[9])*(vel_dist[9]-vel_dist[8])   +
				     DM_1_0GeV_18330V[rr][1]*sub_MB_dist(vel_dist[10])*(vel_dist[10]-vel_dist[9]) +
				     DM_1_0GeV_20000V[rr][1]*sub_MB_dist(vel_dist[11])*(vel_dist[11]-vel_dist[10]) +
				     DM_1_0GeV_21670V[rr][1]*sub_MB_dist(vel_dist[12])*(vel_dist[12]-vel_dist[11]) +
				     DM_1_0GeV_25460V[rr][1]*sub_MB_dist(vel_dist[13])*(vel_dist[13]-vel_dist[12]) 
				     )
				    -( DM_1_0GeV_03333V[rr][1]*add_MB_dist(vel_dist[1])*(vel_dist[1]-vel_dist[0])   +
				       DM_1_0GeV_06667V[rr][1]*add_MB_dist(vel_dist[2])*(vel_dist[2]-vel_dist[1])   +
				       DM_1_0GeV_06667V[rr][1]*add_MB_dist(vel_dist[3])*(vel_dist[3]-vel_dist[2])   +
				       DM_1_0GeV_08333V[rr][1]*add_MB_dist(vel_dist[4])*(vel_dist[4]-vel_dist[3])   +
				       DM_1_0GeV_10000V[rr][1]*add_MB_dist(vel_dist[5])*(vel_dist[5]-vel_dist[4])   	  
				       )
				    -(
				      DM_1_0GeV_11670V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[6]-vel_dist[5])   +	  
				      DM_1_0GeV_13330V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[7]-vel_dist[6])   +	  
				      DM_1_0GeV_15000V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[8]-vel_dist[7])   +	  
				      DM_1_0GeV_16670V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[9]-vel_dist[8])   +
				      DM_1_0GeV_18330V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[10]-vel_dist[9])  +
				      DM_1_0GeV_20000V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[11]-vel_dist[10]) +
				      DM_1_0GeV_21670V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[12]-vel_dist[11]) +
				      DM_1_0GeV_25460V[rr][1]*exp(-pow(Vesc/V0,2))*(vel_dist[13]-vel_dist[12]) 
				      )
						      )
		    );
      

      //rate[rr] = log10(rate[rr]); 
      
     
      //printf("%f  %e   \n",energy[rr], rate[rr] );
    }











  
  TCanvas *plot = new TCanvas("plot","",800,600);
  TPad *pad1 = new TPad("pad1", "",0.0,0.3,1.0,1.0);

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

  
  pad1->SetFillColor(0);
  pad1->Draw();
  pad1->cd();
  
  
  pad1->SetFillColor(0);
  pad1->SetBorderSize(2);
  pad1->SetRightMargin(0.03);
  pad1->SetTopMargin(0.03);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.1);

  gPad->SetLogy(1);
  gPad->SetLogx(1);
 

  /*
  TCanvas *plot = new TCanvas("plot");
  gPad->SetLogy(1);
  gPad->SetLogx(1);
  gStyle->SetOptStat(0);
  */
  
  TH2F *frame = new TH2F("frame","",10,1E-2,1.8,10,1E-30,1E-21);
  frame->GetYaxis()->SetTitle("#frac{d<#sigma v>}{dT} cm^{3} keV^{-1} Day^{-1}");
  frame->GetXaxis()->SetTitle(" T (keV)");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  //frame->GetXaxis()->SetMoreLogLabels();
  frame->Draw();

  
  TGraph *gr = new TGraph(data_bin-1,energy, rate);
  gr->SetLineColor(1);
  gr->SetLineWidth(10);
  gr->Draw("al");

  TGraph *grf = new TGraph(data_bin-1,energyf, ratef);
  grf->SetLineColor(4);
  grf->SetLineWidth(3);
  grf->Draw("l");




  
  TGraph *grBM19 = new TGraph("Xe1904_shortRange_1GeV.txt");
  grBM19->SetLineColor(2);
  grBM19->SetLineWidth(3);
  grBM19->Draw("l");
  
//TLatex *  tex = new TLatex(0.01864007,-24.25754,"Our RRPA");
//tex->SetTextFont(22);
//tex->SetLineWidth(2);
//tex->Draw();
//tex = new TLatex(0.01944321,-24.9254,"BM Roberts arXiv1904.");
//tex->SetTextColor(2);
//tex->SetTextFont(22);
//tex->SetLineWidth(2);
//tex->Draw();
//

  
 
  //FILE *Etafile = fopen("data_c1_dcs_1_000GeV.txt","r");
  FILE *Etafile = fopen("data_c1_Xedcs_1_000GeV.txt","r");
  //const int DataPt = 155;
  const int DataPt = 183;
  
  double EtaRate[DataPt], EtaEnergy[DataPt];
  for(int rr = 0; rr < DataPt;  rr++)
    {
      fscanf(Etafile,"%lf  %lf \n", &EtaEnergy[rr], &EtaRate[rr]);
      EtaEnergy[rr] = EtaEnergy[rr]*1E-3;
      
      EtaRate[rr] =  (rate_factor*EtaRate[rr]*3E+10);
    }


  TGraph *grEta = new TGraph(DataPt, EtaEnergy, EtaRate);
  grEta->SetLineColor(6);
  grEta->SetLineWidth(3);
  grEta->SetMarkerStyle(20);
  grEta->Draw("l");

  TLegend *leg = new TLegend(0.17,0.24,0.45,0.43);
  leg ->AddEntry(grBM19,"BM Roberts ","l");
  leg ->AddEntry(gr,"Old Discrete averaged","l");
  leg ->AddEntry(grf,"New Discrete averaged","l");
  
  leg ->AddEntry(grEta,"With #eta-fun","l");
  leg ->SetBorderSize(0);
  leg ->SetTextFont(22);
  leg ->SetTextSize(0.04);
  leg ->Draw();
   



  plot->cd();
      
  TPad *pad2 = new TPad("pad2", "",0.0,0.0,1.0,0.3);
  pad2->Draw();
  pad2->cd();

  pad2->SetFillColor(0);
  pad2->SetBorderSize(2);
  pad2->SetRightMargin(0.03);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.15);
  pad2->SetLeftMargin(0.1);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetFrameLineWidth(2);	
  gStyle->SetLabelOffset(0.010,"Y");
  gStyle->SetLabelSize(0.09, "Y");
  gStyle->SetNdivisions(507,"Y");
  gStyle->SetTitleSize(0.12, "Y" );
  gStyle->SetTitleOffset(0.25,"Y");
  
  gStyle->SetLabelOffset(-0.01,"X");
  gStyle->SetLabelSize(0.09, "X");
  gStyle->SetNdivisions(510,"X");
  gStyle->SetTitleSize(0.1, "X" );
  gStyle->SetTitleOffset(0.6,"X");
  gStyle->SetTitleFont(22,"X");
  gStyle->SetTitleFont(22,"Y");
  gStyle->SetLabelFont(22,"X");
  gStyle->SetLabelFont(22,"Y");
  gStyle->SetTickLength(0.1,"X");

  gPad->SetLogx(1);
  /*
 TCanvas *plot2 = new TCanvas("plot2");
  gPad->SetLogx(1);
  */
  
  TH2F *frame2 = new TH2F("frame2","",10,1E-2,1.8,10,0.1,2.3);
  frame2->GetYaxis()->SetTitle("Ratio");
  frame2->GetXaxis()->SetTitle(" T (keV)");
  frame2->GetXaxis()->CenterTitle();
  frame2->GetYaxis()->CenterTitle();
  //frame->GetXaxis()->SetMoreLogLabels();
  frame2->Draw();


  double ratio1[data_bin];
  double ratio2[data_bin];
  for(int rr = 0; rr < data_bin-1;  rr++)
    {
      energy[rr] = DM_1_0GeV_10000V[rr][0]*1E-3;
      //ratio[rr] = pow(10,grBM19->Eval(energy[rr]))/pow(10, rate[rr]);
      ratio1[rr] = grf->Eval(energy[rr])/grEta-> Eval(energy[rr]);
      ratio2[rr] = grBM19->Eval(energy[rr])/gr-> Eval(energy[rr]);
      //printf("%f  %.3f %.3f    \n",energy[rr], grBM19->Eval(energy[rr]), rate[rr] );
    }

   TGraph *gratio = new TGraph(data_bin-1,energy, ratio1);
  gratio->SetLineColor(1);
  gratio->SetLineWidth(3);
  gratio->SetMarkerStyle(20);
  gratio->Draw("pl");
  
  TGraph *grat2 = new TGraph(data_bin-1,energy, ratio2);
  grat2->SetLineColor(2);
  grat2->SetMarkerColor(2);
  grat2->SetLineWidth(3);
  grat2->SetMarkerStyle(20);
  grat2->Draw("pl");
  

  TLatex *tex = new TLatex(0.201039,1.872379,"Disc-New / #eta_{fun}");
   tex->SetTextFont(22);
   tex->SetTextSize(0.1275362);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.02648188,1.902387,"BMRobert / Disc-Old");
   tex->SetTextColor(2);
   tex->SetTextFont(22);
   tex->SetTextSize(0.1275362);
   tex->SetLineWidth(2);
   tex->Draw();
}









