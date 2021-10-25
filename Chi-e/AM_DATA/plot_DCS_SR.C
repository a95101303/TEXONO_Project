#include "Header_SL_Parameter.h"

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
 
  
    float SE_to_NT_P = 1.64458;//68%->90%
    const double Lower_Bound = 49.432*SE_to_NT_P;

    double Sigma_e_Dec[26];double Sigma_e_Jun[26];

for(int M_Idx=10; M_Idx<11; M_Idx++)
{
      double mx = massDM[M_Idx] ;
    cout << "mx: " << mx << endl;
    cout << "=================================================" << endl;
    int DecRate_Check=0;int JuneRate_Check=0;
    
        for(int Sig_Idx=1; Sig_Idx<2; Sig_Idx++)
        {
            double c1 = 9*1e-20;
            cout << "c1: " << CS_Try(c1,mx) << endl;

          for(int rr = 0; rr<NumPt;rr++)
            {
              rate_factor = pow(c1,2)*1E-18*1E+3*86400*3E+10*(density*electron_number)/(mx);
              
              Energy[rr] = 95.0E-3 + 5.0*(double)rr/(double)NumPt;
              JuneRate[rr] = rate_factor*ShortRangeDcs_June(Energy[rr]*1E+3 , mx);
              DecRate[rr]  = rate_factor*ShortRangeDcs_Dec(Energy[rr]*1E+3  , mx);
              Ratio[rr]  = DecRate[rr]-JuneRate[rr];

              if(Energy[rr]==0.20)
              {
                  cout << "=================================================" << endl;
                  cout << "DecRate_Final[rr]: " << DecRate[rr] << endl;
                  cout << "=================================================" << endl;

                  if(JuneRate_Check==0)
                  {
                      double scale_June = Lower_Bound/JuneRate[rr];
                      double scale_June_c1 = sqrt(scale_June);
                      cout << "mx: " << mx << endl;
                      cout << "JuneRate_Final[rr]: " << JuneRate[rr]*scale_June << endl;
                      cout << "June_Cross-section_sigma_e: " << CS_Try(c1*scale_June_c1,mx) << endl;
                      Sigma_e_Jun[M_Idx] = CS_Try(c1*scale_June_c1,mx);
                      JuneRate_Check = JuneRate_Check + 1;
                  }
                  if(DecRate_Check==0)
                  {
                      double scale_Dec = Lower_Bound/DecRate[rr];
                      double scale_Dec_c1 = sqrt(scale_Dec);
                      cout << "mx: " << mx << endl;
                      cout << "DecRate_Final[rr]: " << DecRate[rr]*scale_Dec << endl;
                      cout << "Dec_Cross-section_sigma_e: " << CS_Try(c1*scale_Dec_c1,mx) << endl;
                      Sigma_e_Dec[M_Idx] = CS_Try(c1*scale_Dec_c1,mx);
                      DecRate_Check = DecRate_Check + 1;
                  }
              }

              //printf("%f  %e \n", Energy[rr], JuneRate[rr]);
            }
        }
    
}
    /*
    for(int M_Idx=0; M_Idx<26; M_Idx++)
    {
        cout << "mx: " << massDM[M_Idx] << endl;
        cout << "Sigma_e_Jun[M_Idx]: " << Sigma_e_Jun[M_Idx] << endl;
        cout << "Sigma_e_Dec[M_Idx]: " << Sigma_e_Dec[M_Idx] << endl;
     }
     
    for(int M_Idx=0; M_Idx<26; M_Idx++){cout << massDM[M_Idx] << "," << Sigma_e_Jun[M_Idx] << "," << endl;}
    for(int M_Idx=0; M_Idx<26; M_Idx++){cout << massDM[M_Idx] << "," << Sigma_e_Dec[M_Idx] << "," << endl;}
     */

  TCanvas *plot = new TCanvas("plot");
  gPad->SetLogy(1);
  gPad->SetLogx(1);

      
      
  TGraph *grJ = new TGraph(NumPt , Energy , JuneRate );
  grJ->SetLineWidth(3);
  grJ->SetLineColor(kRed);
  grJ->Draw("al");

    /*
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
     */

     


  
  
}
