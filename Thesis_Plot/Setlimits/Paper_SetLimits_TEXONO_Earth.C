#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "TCutG.h"
#include "TCut.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/cdms_si_allow68.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/superCDMS2014_band.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/cogent2013.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CMB.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/XQC.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/cresstII.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/cresst_surf.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/cresst_surf_Ours.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_surf_Ours.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Ours.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_real_Ours.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_real_Migdal_Lower.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_real_Migdal_Upper.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Lower.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Upper.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Brem_Lower.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Brem_Upper.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Lower_All_Bins.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Upper_Bend.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Brem_Upper_Bend.h"

void Paper_SetLimits_TEXONO_Earth()
{
    TLegend *leg = new TLegend(0.15,0.14,0.4,0.35);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);

  int Number_of_Candidates=20;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetFrameLineWidth(2);	
  //gStyle->SetLabelOffset(0.007,"Y");
  //gStyle->SetLabelSize(0.04, "Y");
  //gStyle->SetNdivisions(510,"Y");
  //gStyle->SetLabelSize(0.04, "X");
  //gStyle->SetNdivisions(510,"X");
  gStyle->SetTitleSize(0.055, "X" );
  gStyle->SetTitleSize(0.055, "Y" );
  //gStyle->SetTitleFont(42,"X");
  //gStyle->SetTitleFont(42,"Y");
  //gStyle->SetLabelFont(42,"X");
  //gStyle->SetLabelFont(42,"Y");
  //gStyle->SetTitleOffset(1.0,"X");
  //gStyle->SetTitleOffset(1.5,"Y");

  TCanvas *plot = new TCanvas("plot","",800,600);
  
  TPad *pad1 = new TPad("pad1", "",0.0,0.0,1.0,1.0);

  pad1->Draw();
  pad1->cd();
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad1->SetRightMargin(0.02);
  pad1->SetTopMargin(0.02);
  pad1->SetBottomMargin(0.12);
  pad1->SetLeftMargin(0.13);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameBorderMode(0);
  pad1->SetTickx(1);
  pad1->SetTicky(1);

//  TH2F *frame = new TH2F("frame","",100,3.0,200,100,5e-43,1e-39);
  TH2F *frame = new TH2F("frame","",100,0.05,2,100,1e-32,1e-26);
  //frame->GetXaxis()->SetTitle("WIMP Mass (GeV/c^{2})");
  //frame->GetYaxis()->SetTitle("SI Corss section (cm^{2})");
  //frame->GetXaxis()->SetTitle("m_{#chi} (GeV/c^{2})");
//frame->GetYaxis()->SetTitle("#sigma_{SI} (cm^{2})");

  //frame->GetXaxis()->CenterTitle();
  //frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetMoreLogLabels();
  frame->Draw();
    frame->GetXaxis()->SetTitle("M_{#chi}[GeV]");
    frame->GetYaxis()->SetTitle("#sigma_{SI}[cm^{2}]");
  gPad->SetLogx();
  gPad->SetLogy();

    
    TFile *f=new TFile("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/Set_Limits_Plot_for_ALL.root");
    


//TEXONO
    TGraph *TEXONO_Down=(TGraph*)f->Get("TEXONO_With_Subtracting_the_KL_shells");
    TEXONO_Down->SetLineWidth(8);
    TEXONO_Down->SetName("TEXONO_Real");
    TEXONO_Down->SetFillColor(kOrange-3);
    TEXONO_Down->SetLineColor(kOrange+10);
    TEXONO_Down->SetLineWidth(5);
    TEXONO_Down->SetFillStyle(3144);
    TEXONO_Down->Draw("lfsame");

    
    double x, y0;

     double TEXONO_x[Number_of_Candidates], TEXONO_y[Number_of_Candidates];
    for(int i=0;i<Number_of_Candidates;i++)
      {
        TEXONO_Down->GetPoint(i,x,y0);
        TEXONO_x[i] = x;
        TEXONO_y[i] = y0;
          cout << "TEXONO_x[i]: " << TEXONO_x[i] << endl;
          cout << "TEXONO_y[i]: " << TEXONO_y[i] << endl;

      }

    TGraph *texono_real = new TGraph(); texono_real->SetName("texono_real");
    texono_real->SetPoint(0, TEXONO_x[0],TEXONO_y[0]);

    cout << "TEXONO_x[0]: " << TEXONO_x[0] << endl;
    cout << "TEXONO_y[0]: " << TEXONO_y[0] << endl;

    for(int i=0;i<12;i++)
    {
        cout << "TEXONO_real_chiN[i][0]: " << TEXONO_real_chiN[i][0] << endl;
        cout << "TEXONO_real_chiN[i][1]: " << TEXONO_real_chiN[i][1] << endl;

        texono_real->SetPoint((i+1), TEXONO_real_chiN[i][0],TEXONO_real_chiN[i][1]);
    }
    texono_real->SetPoint(13, TEXONO_x[18],TEXONO_y[18]);
    cout << "TEXONO_x[20]: " << TEXONO_x[18] << endl;
    cout << "TEXONO_y[20]: " << TEXONO_y[18] << endl;

    texono_real->SetName("texono_real");
    texono_real->SetFillColor(kOrange-3);
    texono_real->SetFillStyle(3144);
    texono_real->SetLineColor(kOrange+10);
    texono_real->SetLineWidth(5);
    texono_real->Draw("lfsame");
    
    TLatex *tex11 = new TLatex(2,1e-39,"TEXONO(NU)");
    tex11->SetTextFont(42);
    tex11->SetTextSize(0.04);
    tex11->SetLineWidth(2);
    tex11->SetTextColor(kOrange-3);
    tex11->SetTextAngle(-25);
    tex11->Draw();

     
//
         
     //=====MAT
     TGraph *texono_real_Migdal_Bend = new TGraph(); texono_real_Migdal_Bend->SetName("texono_real_Migdal_Bend");

     for(int i=0;i<4;i++)
     {
         cout << "TEXONO_real_Migdal_Upper_Bend[i][0]: " << TEXONO_real_Migdal_Upper_Bend[i][0] << endl;
         cout << "TEXONO_real_Migdal_Upper_Bend[i][1]: " << TEXONO_real_Migdal_Upper_Bend[i][1] << endl;

         texono_real_Migdal_Bend->SetPoint((i), TEXONO_real_Migdal_Upper_Bend[i][0],TEXONO_real_Migdal_Upper_Bend[i][1]);
     }
     

     texono_real_Migdal_Bend->SetLineColor(2);
     texono_real_Migdal_Bend->SetLineWidth(10);
     texono_real_Migdal_Bend->Draw("lsame");

    TLatex *tex20;
    tex20 = new TLatex(0.05,9e-28,"MD: TEXONO(MAT:Real)");
    tex20->SetTextFont(42);
    tex20->SetTextSize(0.04);
    tex20->SetLineWidth(2);
    tex20->SetTextColor(2);
    tex20->SetTextAngle(0);
    tex20->Draw("same");


        //=====CAT
    TGraph *texono_real_Migdal = new TGraph(); texono_real_Migdal->SetName("texono_real_MD");

    for(int i=0;i<16;i++)
    {
        cout << "TEXONO_real_Migdal_Lower[i][0]: " << TEXONO_real_Migdal_Lower[i][0] << endl;
        cout << "TEXONO_real_Migdal_Lower[i][1]: " << TEXONO_real_Migdal_Lower[i][1] << endl;

        texono_real_Migdal->SetPoint((i), TEXONO_real_Migdal_Lower[i][0],TEXONO_real_Migdal_Lower[i][1]);
    }
    
    for(int i=0;i<16;i++)
    {
        cout << "CDEX_real_Migdal_Upper[i][0]: " << TEXONO_real_Migdal_Upper[i][0] << endl;
        cout << "CDEX_real_Migdal_Upper[i][1]: " << TEXONO_real_Migdal_Upper[i][1] << endl;

        texono_real_Migdal->SetPoint((i+16), TEXONO_real_Migdal_Upper[i][0],TEXONO_real_Migdal_Upper[i][1]);
    }
    
     texono_real_Migdal->SetPoint((32), TEXONO_real_Migdal_Lower[0][0],TEXONO_real_Migdal_Lower[0][1]);

    texono_real_Migdal->SetFillColor(kCyan-4);
    texono_real_Migdal->SetLineColor(kOrange+10);
    texono_real_Migdal->SetLineWidth(5);
    texono_real_Migdal->Draw("fsame");
    texono_real_Migdal->SetFillStyle(3144);

    TLatex *tex6;
    tex6 = new TLatex(0.06,2e-33,"MD: TEXONO(CAT:Tradition)");
    tex6->SetTextFont(42);
    tex6->SetTextSize(0.04);
    tex6->SetLineWidth(2);
    tex6->SetTextColor(kCyan-4);
    tex6->SetTextAngle(0);
    tex6->Draw();

     //=====MAT
     TGraph *texono_real_brem_Bend = new TGraph(); texono_real_brem_Bend->SetName("texono_real_brem_Bend");

     for(int i=0;i<4;i++)
     {
         cout << "TEXONO_real_Brem_Upper_Bend[i][0]: " << TEXONO_real_Brem_Upper_Bend[i][0] << endl;
         cout << "TEXONO_real_Brem_Upper_Bend[i][1]: " << TEXONO_real_Brem_Upper_Bend[i][1] << endl;

         texono_real_brem_Bend->SetPoint((i), TEXONO_real_Brem_Upper_Bend[i][0],TEXONO_real_Brem_Upper_Bend[i][1]);
     }
     

     texono_real_brem_Bend->SetLineColor(2);
     texono_real_brem_Bend->SetLineWidth(10);
     texono_real_brem_Bend->Draw("lsame");

    TLatex *tex21;
    tex21 = new TLatex(0.3,1e-29,"Brem: TEXONO(MAT:Real)");
    tex21->SetTextFont(42);
    tex21->SetTextSize(0.04);
    tex21->SetLineWidth(2);
    tex21->SetTextColor(2);
    tex21->SetTextAngle(0);
    tex21->Draw("same");

    //=====CAT
    TGraph *texono_real_brem = new TGraph(); texono_real_brem->SetName("texono_real_BR");

    for(int i=0;i<19;i++)
    {
        cout << "TEXONO_real_Brem_Lower[i][0]: " << TEXONO_real_Brem_Lower[i][0] << endl;
        cout << "TEXONO_real_Brem_Lower[i][1]: " << TEXONO_real_Brem_Lower[i][1] << endl;

        texono_real_brem->SetPoint((i), TEXONO_real_Brem_Lower[i][0],TEXONO_real_Brem_Lower[i][1]);
    }
    
    for(int i=0;i<19;i++)
    {
        cout << "TEXONO_real_Brem_Lower[i][0]: " << TEXONO_real_Brem_Upper[i][0] << endl;
        cout << "TEXONO_real_Brem_Lower[i][1]: " << TEXONO_real_Brem_Upper[i][1] << endl;

        texono_real_brem->SetPoint((i+19), TEXONO_real_Brem_Upper[i][0],TEXONO_real_Brem_Upper[i][1]);
    }

     
    texono_real_brem->SetFillColor(kGreen+2);
    texono_real_brem->SetLineColor(kOrange+10);
    texono_real_brem->SetLineWidth(5);
    texono_real_brem->Draw("fsame");
    texono_real_brem->SetFillStyle(3144);
   // texono_real_brem->SetMarkerStyle(8);
   // texono_real_brem->SetMarkerSize(0.5);
   // texono_real_brem->SetMarkerColor(2);

    TLatex *tex7;
    tex7 = new TLatex(0.1,1e-29,"Brem: TEXONO(CAT:Tradition)");
    tex7->SetTextFont(42);
    tex7->SetTextSize(0.04);
    tex7->SetLineWidth(2);
    tex7->SetTextColor(kGreen+2);
    tex7->SetTextAngle(-15);
    tex7->Draw();
    

     
    leg->Draw();
      plot->Print("TEXONO_Earth_Bend.pdf");
}
