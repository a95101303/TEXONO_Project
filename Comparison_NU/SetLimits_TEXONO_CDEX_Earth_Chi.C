#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "TCutG.h"
#include "TCut.h"
#include "cdms_si_allow68.h"
#include "superCDMS2014_band.h"
#include "cogent2013.h"

#include "CMB.h"
#include "XQC.h"
#include "cresstII.h"
#include "cresst_surf.h"
#include "cresst_surf_Ours.h"
#include "TEXONO_surf_Ours.h"
#include "TEXONO_real_Ours.h"
#include "CDEX_real_Ours.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_ER_Lower_c1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_ER_Lower_d1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_ER_Upper_c1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_ER_Upper_d1.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_ER_Lower_c1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_ER_Lower_d1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_ER_Upper_c1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_ER_Upper_d1.h"


void SetLimits_TEXONO_CDEX_Earth_Chi()
{
    TLegend *leg = new TLegend(0.15,0.14,0.4,0.35);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);

    TLegend *leg1 = new TLegend(0.55,0.14,0.75,0.35);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.04);
    leg1->SetBorderSize(0);
    leg1->SetTextFont(22);

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
  TH2F *frame = new TH2F("frame","",100,0.05,15,100,4e-41,1e-15);
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

//=================================================================================================
    /*
    cout << "===========================================================================================" << endl;
  TGraph *Graph_TEXONO_ER_Lower_c1 = new TGraph(); Graph_TEXONO_ER_Lower_c1->SetName("Graph_TEXONO_ER_Lower_c1");
  for(int i=0;i<13;i++)
  {
      Graph_TEXONO_ER_Lower_c1->SetPoint((i), TEXONO_ER_Lower_c1[i][0], TEXONO_ER_Lower_c1[i][1]);
      cout << "TEXONO_ER_Lower_c1.h[i][0]: " << TEXONO_ER_Lower_c1[i][0] <<  endl;
      cout << "TEXONO_ER_Lower_c1.h[i][1]: " << TEXONO_ER_Lower_c1[i][1] <<  endl;
  }
    Graph_TEXONO_ER_Lower_c1->SetFillColor(kAzure+6);
    //Graph_TEXONO_ER_Lower_c1->SetFillStyle(3114);
    Graph_TEXONO_ER_Lower_c1->SetLineStyle(1);
    Graph_TEXONO_ER_Lower_c1->SetLineColor(2);
    Graph_TEXONO_ER_Lower_c1->Draw("L");

    TGraph *Graph_TEXONO_ER_Lower_d1 = new TGraph(); Graph_TEXONO_ER_Lower_d1->SetName("Graph_TEXONO_ER_Lower_d1");
    for(int i=0;i<13;i++)
    {
        Graph_TEXONO_ER_Lower_d1->SetPoint((i), TEXONO_ER_Lower_d1[i][0], TEXONO_ER_Lower_d1[i][1]);
        cout << "TEXONO_ER_Lower_d1.h[i][0]: " << TEXONO_ER_Lower_d1[i][0] <<  endl;
        cout << "TEXONO_ER_Lower_d1.h[i][1]: " << TEXONO_ER_Lower_d1[i][1] <<  endl;
    }
      Graph_TEXONO_ER_Lower_d1->SetFillColor(kAzure+6);
      Graph_TEXONO_ER_Lower_d1->SetLineStyle(1);
      Graph_TEXONO_ER_Lower_d1->SetLineColor(3);
      //Graph_TEXONO_ER_Lower_c1->SetFillStyle(3114);
      Graph_TEXONO_ER_Lower_d1->Draw("Lsame");

    TGraph *Graph_TEXONO_ER_Upper_c1 = new TGraph(); Graph_TEXONO_ER_Upper_c1->SetName("Graph_TEXONO_ER_Upper_c1");
    for(int i=0;i<13;i++)
    {
        Graph_TEXONO_ER_Upper_c1->SetPoint((i), TEXONO_ER_Upper_c1[i][0], TEXONO_ER_Upper_c1[i][1]);
        cout << "TEXONO_ER_Upper_c1[i][0][i][0]: " << TEXONO_ER_Upper_c1[i][0] <<  endl;
        cout << "TEXONO_ER_Upper_c1[i][1]: " << TEXONO_ER_Upper_c1[i][1] <<  endl;
    }
      Graph_TEXONO_ER_Upper_c1->SetFillColor(kAzure+6);
      Graph_TEXONO_ER_Upper_c1->SetLineStyle(5);
      Graph_TEXONO_ER_Upper_c1->SetLineColor(2);
      //Graph_TEXONO_ER_Lower_c1->SetFillStyle(3114);
      Graph_TEXONO_ER_Upper_c1->Draw("Lsame");

    
      TGraph *Graph_TEXONO_ER_Upper_d1 = new TGraph(); Graph_TEXONO_ER_Upper_d1->SetName("Graph_TEXONO_ER_Upper_d1");
      for(int i=0;i<13;i++)
      {
          Graph_TEXONO_ER_Upper_d1->SetPoint((i), TEXONO_ER_Upper_d1[i][0], TEXONO_ER_Upper_d1[i][1]);
          cout << "TEXONO_ER_Upper_d1[i][0]: " << TEXONO_ER_Upper_d1[i][0] <<  endl;
          cout << "TEXONO_ER_Upper_d1[i][1]: " << TEXONO_ER_Upper_d1[i][1] <<  endl;
      }
        Graph_TEXONO_ER_Upper_d1->SetFillColor(kAzure+6);
        Graph_TEXONO_ER_Upper_d1->SetLineColor(3);
        Graph_TEXONO_ER_Upper_d1->SetLineStyle(5);
        //Graph_TEXONO_ER_Lower_c1->SetFillStyle(3114);
        Graph_TEXONO_ER_Upper_d1->Draw("Lsame");
    cout << "===========================================================================================" << endl;
    */
      cout << "===========================================================================================" << endl;
  TGraph *Graph_CDEX_ER_Lower_c1 = new TGraph(); Graph_CDEX_ER_Lower_c1->SetName("Graph_CDEX_ER_Lower_c1");
  for(int i=0;i<13;i++)
  {
      Graph_CDEX_ER_Lower_c1->SetPoint((i), CDEX_ER_Lower_c1[i][0], CDEX_ER_Lower_c1[i][1]);
      cout << "CDEX_ER_Lower_c1.h[i][0]: " << CDEX_ER_Lower_c1[i][0] <<  endl;
      cout << "CDEX_ER_Lower_c1.h[i][1]: " << CDEX_ER_Lower_c1[i][1] <<  endl;
  }
    Graph_CDEX_ER_Lower_c1->SetFillColor(kAzure+6);
    //Graph_CDEX_ER_Lower_c1->SetFillStyle(3114);
    Graph_CDEX_ER_Lower_c1->SetLineStyle(1);
    Graph_CDEX_ER_Lower_c1->SetLineColor(2);
    Graph_CDEX_ER_Lower_c1->Draw("L");

    TGraph *Graph_CDEX_ER_Lower_d1 = new TGraph(); Graph_CDEX_ER_Lower_d1->SetName("Graph_CDEX_ER_Lower_d1");
    for(int i=0;i<13;i++)
    {
        Graph_CDEX_ER_Lower_d1->SetPoint((i), CDEX_ER_Lower_d1[i][0], CDEX_ER_Lower_d1[i][1]);
        cout << "CDEX_ER_Lower_d1.h[i][0]: " << CDEX_ER_Lower_d1[i][0] <<  endl;
        cout << "CDEX_ER_Lower_d1.h[i][1]: " << CDEX_ER_Lower_d1[i][1] <<  endl;
    }
      Graph_CDEX_ER_Lower_d1->SetFillColor(kAzure+6);
      Graph_CDEX_ER_Lower_d1->SetLineStyle(1);
      Graph_CDEX_ER_Lower_d1->SetLineColor(3);
      //Graph_CDEX_ER_Lower_c1->SetFillStyle(3114);
      Graph_CDEX_ER_Lower_d1->Draw("Lsame");

    TGraph *Graph_CDEX_ER_Upper_c1 = new TGraph(); Graph_CDEX_ER_Upper_c1->SetName("Graph_CDEX_ER_Upper_c1");
    for(int i=0;i<13;i++)
    {
        Graph_CDEX_ER_Upper_c1->SetPoint((i), CDEX_ER_Upper_c1[i][0], CDEX_ER_Upper_c1[i][1]);
        cout << "CDEX_ER_Upper_c1[i][0][i][0]: " << CDEX_ER_Upper_c1[i][0] <<  endl;
        cout << "CDEX_ER_Upper_c1[i][1]: " << CDEX_ER_Upper_c1[i][1] <<  endl;
    }
      Graph_CDEX_ER_Upper_c1->SetFillColor(kAzure+6);
      Graph_CDEX_ER_Upper_c1->SetLineStyle(5);
      Graph_CDEX_ER_Upper_c1->SetLineColor(2);
      //Graph_CDEX_ER_Lower_c1->SetFillStyle(3114);
      Graph_CDEX_ER_Upper_c1->Draw("Lsame");

    
      TGraph *Graph_CDEX_ER_Upper_d1 = new TGraph(); Graph_CDEX_ER_Upper_d1->SetName("Graph_CDEX_ER_Upper_d1");
      for(int i=0;i<13;i++)
      {
          Graph_CDEX_ER_Upper_d1->SetPoint((i), CDEX_ER_Upper_d1[i][0], CDEX_ER_Upper_d1[i][1]);
          cout << "CDEX_ER_Upper_d1[i][0]: " << CDEX_ER_Upper_d1[i][0] <<  endl;
          cout << "CDEX_ER_Upper_d1[i][1]: " << CDEX_ER_Upper_d1[i][1] <<  endl;
      }
        Graph_CDEX_ER_Upper_d1->SetFillColor(kAzure+6);
        Graph_CDEX_ER_Upper_d1->SetLineColor(3);
        Graph_CDEX_ER_Upper_d1->SetLineStyle(5);
        //Graph_CDEX_ER_Lower_c1->SetFillStyle(3114);
        Graph_CDEX_ER_Upper_d1->Draw("Lsame");
    cout << "===========================================================================================" << endl;

    
    /*
    TLatex *CMB_Text = new TLatex(0.055,2.5e-26,"CMB");
    CMB_Text->SetTextFont(42);
    CMB_Text->SetTextSize(0.04);
    CMB_Text->SetLineWidth(2);
    CMB_Text->SetTextColor(1);
    CMB_Text->SetTextAngle(5);
     */
    //CMB_Text->Draw();
    
    //leg->AddEntry(gCMB,"CMB","f");
    /*
    leg->AddEntry(gXQC,"XQC","f");
    leg->AddEntry(gcresst_surf,"CRESST(2017) Surface","f");
    leg->AddEntry(gcresstII,"CRESST II","f");
    leg->AddEntry(g_damah_90,"DAMA2009","f");
    leg->AddEntry(g_cogent2013,"CoGeNT(2013)","f");
    leg->Draw();
    leg1->Draw();
     */
    //plot->Print("Both_Earth_Chi_TEXONO.pdf");

  plot->Print("Both_Earth_Chi_CDEX.pdf");
}
