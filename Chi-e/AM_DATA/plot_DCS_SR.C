const int DataBin = 255;
#include "Header_SR_DCS_June.h"
#include "Header_SR_DCS_Dec.h"
void plot_DCS_SR()
{

  double massDM[26] =  {0.06, 0.07,  0.08, 0.09,  0.10, 0.12, 0.20,
			0.25, 0.30,  0.40, 0.50,  0.60, 0.70, 0.80,
			0.90, 1.00,  1.30, 1.50,  2.00, 2.50, 3.00,
			3.50, 4.00,  5.00, 10.00, 20.00};


  const int NumPt = 1000;
  
  double Energy[NumPt];
  double JuneRate[NumPt];
  double DecRate[NumPt];

  double Ratio[NumPt];


  const double PI = 3.141592653589793238;
  const double AGe = 72.64;
  const double avogadro_number = 6.02214129e+23; 
  const double N_atom_1kg_Ge = (1000.0*avogadro_number)/AGe;
  
  const double density = 0.3; //GeV cm^{-3}
  const double electron_number = N_atom_1kg_Ge; //kg^{-1}
  
  double rate_factor;
 
  double c1 = 1E-18;
  double mx; //GeV

  mx = massDM[10];
 
  printf("mx = %f\n",mx);
  
  for(int rr = 0; rr<NumPt;rr++)
    {
      rate_factor = pow(c1,2)*1E-18*1E+3*86400*3E+10*(density*electron_number)/(mx);
      
      Energy[rr] = 95.0E-3 + 5.0*(double)rr/(double)NumPt;
      JuneRate[rr] = rate_factor*ShortRangeDcs_June(Energy[rr]*1E+3 , mx);
      DecRate[rr]  = rate_factor*ShortRangeDcs_Dec(Energy[rr]*1E+3  , mx);
      Ratio[rr]  = DecRate[rr]-JuneRate[rr];

      
      //printf("%f  %e \n", Energy[rr], JuneRate[rr]);
    }


  TCanvas *plot = new TCanvas("plot");
  gPad->SetLogy(1);
  gPad->SetLogx(1);

      
      
  TGraph *grJ = new TGraph(NumPt , Energy , JuneRate );
  grJ->SetLineWidth(3);
  grJ->SetLineColor(kRed);
  grJ->Draw("al");

     
  TGraph *grD = new TGraph(NumPt , Energy , DecRate );
  grD->SetLineWidth(3);
  grD->SetLineColor(kBlue);
  grD->Draw("l");



  TCanvas *plot1 = new TCanvas("plot1");
  //gPad->SetLogy(1);
  gPad->SetLogx(1);
    
  TGraph *grF = new TGraph(NumPt , Energy , Ratio);
  grF->SetLineWidth(3);
  grF->SetLineColor(2);
  grF->Draw("al");


     


  
  
}
