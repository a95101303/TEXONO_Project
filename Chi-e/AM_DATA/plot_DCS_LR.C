#include "Header_SL_Parameter.h"
#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"

void plot_DCS_LR()
{
    
  double massDM[26] =  {0.06, 0.07,  0.08, 0.09,  0.10, 0.12, 0.20,
			0.25, 0.30,  0.40, 0.50,  0.60, 0.70, 0.80,
			0.90, 1.00,  1.30, 1.50,  2.00, 2.50, 3.00,
			3.50, 4.00,  5.00, 10.00, 20.00};
     
    //double massDM[2] =  {2,0.2};

  const int NumPt = 1000;
  
  double Energy[NumPt];
  double JuneRate[NumPt];
  double DecRate[NumPt];

  double Ratio[NumPt];

  double mx; //GeV
  double rate_factor;
  
    float SE_to_NT_P = 1.64458;//68%->90%
    //mx = 2;//10: 0.5GeV
  const double Lower_Bound = 49.432*SE_to_NT_P;
  printf("mx = %f\n",mx);
    
    //cout << "DS_Try(4.89e-11,0.5): " << DS_Try(4.89e-11,0.5) << endl;
    double Sigma_e_Dec[26];double Sigma_e_Jun[26];
    
    for(int M_Idx=25; M_Idx<26; M_Idx++)
    {
          double mx = massDM[M_Idx] ;
        cout << "mx: " << mx << endl;
        cout << "=================================================" << endl;
        int DecRate_Check=0;int JuneRate_Check=0;
        
            for(int Sig_Idx=1; Sig_Idx<2; Sig_Idx++)
            {
                  double d1 = 1E-9;
                  for(int rr = 0; rr<NumPt;rr++)
                    {
                          rate_factor = pow(d1,2)*1E-18*1E+3*86400*3E+10*(density*electron_number(AGe))/(mx);
                          
                          Energy[rr] = 95.0E-3 + 5.0*(double)rr/(double)NumPt;
                          JuneRate[rr] = rate_factor*LongRangeDcs_June(Energy[rr]*1E+3 , mx );
                          DecRate[rr]  = rate_factor*LongRangeDcs_Dec(Energy[rr]*1E+3  , mx );

                          Ratio[rr]  = JuneRate[rr]-DecRate[rr];
                        if(Energy[rr]==0.20)
                        {
                            cout << "=================================================" << endl;
                            cout << "DecRate_Final[rr]: " << DecRate[rr] << endl;
                            cout << "=================================================" << endl;

                            if(JuneRate_Check==0)
                            {
                                double scale_June = Lower_Bound/JuneRate[rr];
                                double scale_June_d1 = sqrt(scale_June);
                                cout << "mx: " << mx << endl;
                                cout << "JuneRate_Final[rr]: " << JuneRate[rr]*scale_June << endl;
                                cout << "June_Cross-section_sigma_e: " << DS_Try(d1*scale_June_d1,mx) << endl;
                                Sigma_e_Jun[M_Idx] = DS_Try(d1*scale_June_d1,mx);
                                JuneRate_Check = JuneRate_Check + 1;
                            }
                            if(DecRate_Check==0)
                            {
                                double scale_Dec = Lower_Bound/DecRate[rr];
                                double scale_Dec_d1 = sqrt(scale_Dec);
                                cout << "mx: " << mx << endl;
                                cout << "DecRate_Final[rr]: " << DecRate[rr]*scale_Dec << endl;
                                cout << "Dec_Cross-section_sigma_e: " << DS_Try(d1*scale_Dec_d1,mx) << endl;
                                Sigma_e_Dec[M_Idx] = DS_Try(d1*scale_Dec_d1,mx);
                                DecRate_Check = DecRate_Check + 1;
                            }
                        }
                      //printf("%f  %e \n", Energy[rr], JuneRate[rr]);
                    }
            }
        cout << "=================================================" << endl;
    }
    
    for(int M_Idx=0; M_Idx<26; M_Idx++)
    {
        cout << "mx: " << mx << endl;
        cout << "Sigma_e_Jun[M_Idx]: " << Sigma_e_Jun[M_Idx] << endl;
        cout << "Sigma_e_Dec[M_Idx]: " << Sigma_e_Dec[M_Idx] << endl;
     }
     
    for(int M_Idx=0; M_Idx<26; M_Idx++){cout << massDM[M_Idx] << "," << Sigma_e_Jun[M_Idx] << "," << endl;}
    for(int M_Idx=0; M_Idx<26; M_Idx++){cout << massDM[M_Idx] << "," << Sigma_e_Dec[M_Idx] << "," << endl;}
     
    
    TCanvas *plot = new TCanvas("plot");
    gPad->SetLogy(1);
    //gPad->SetLogx(1);

    /*
    double X[1]={0};double Y[1]={0};
    TGraph *SetFrame = new TGraph(1,X,Y);

    SetFrame->GetXaxis()->SetLimits(1e-1,1e2);
    SetFrame->GetYaxis()->SetRangeUser(1e-38,1e-26);
    SetFrame->GetYaxis()->SetTitle("#sigma_{e} and #sigma_{N}");
    SetFrame->GetXaxis()->SetTitle("ER energy[kev]");
    SetFrame->Draw("ALP");
    
    cout << "START: " << endl;
    TF1 *fa2 = new TF1("fa2","dsigma_dT_keV_ER([0],[1],x,[2],[3])",0,20);
    fa2->SetParameter(0,1.72334e-31);
    fa2->SetParameter(1,1);
    fa2->SetParameter(2,3);
    fa2->SetParameter(3,AGe);
    
    
    TF1 *fa3 = new TF1("fa3","fdsigma_dT_keV([0],[1],[2],[3],x)",1e-1,1e1);
    fa3->SetParameter(0,3);
    fa3->SetParameter(1,1e-31);
    fa3->SetParameter(2,779.*V_to_C);
    fa3->SetParameter(3,AGe);
    fa3->Draw("Lsame");
     */
    
    //fa2->Draw("Lsame");

    //cout << "Total_Sigma_ER: " << Total_Sigma_ER(1,300,10,34) << endl;
    /*
    cout << "c1: " << c_1(9e-42,0.5) << endl;
    cout << "d1: " << d_1(1.72334e-31,0.5) << endl;
    cout << "CS_Try(): " << CS_Try(5.279*1e-4,0.5) << endl;
    cout << "DS_Try(): " << DS_Try(1e-9,0.5) << endl;
     */

    
  TCanvas *plot1 = new TCanvas("plot");
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
   
   
   /*
   TCanvas *plot1 = new TCanvas("plot1");
   //gPad->SetLogy(1);
   gPad->SetLogx(1);
     
    
   TGraph *grF = new TGraph(NumPt , Energy , Ratio);
   grF->SetLineWidth(3);
   grF->SetLineColor(2);
   grF->Draw("al");
   */
    
     


  
  
}
