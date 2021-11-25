#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "TCutG.h"
#include "TCut.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison_ER/TEXONO_ER_Lower_d1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison_ER/TEXONO_ER_Upper_d1.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison_ER/CDEX_ER_Lower_d1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison_ER/CDEX_ER_Upper_d1.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison_ER/CDMSlite_ER_Lower_d1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison_ER/XENON_ER_Upper_d1.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison_ER/XENON1T_ER_Lower_d1.h"


void SetLimits_d1()
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
  TH2F *frame = new TH2F("frame","",100,1e-2,10,100,1e-35,1e-21);
  //frame->GetXaxis()->SetTitle("WIMP Mass (GeV/c^{2})");
  //frame->GetYaxis()->SetTitle("SI Corss section (cm^{2})");
  //frame->GetXaxis()->SetTitle("m_{#chi} (GeV/c^{2})");
//frame->GetYaxis()->SetTitle("#sigma_{SI} (cm^{2})");

  //frame->GetXaxis()->CenterTitle();
  //frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetMoreLogLabels();
  frame->Draw();
    frame->GetXaxis()->SetTitle("M_{#chi}[GeV]");
    frame->GetYaxis()->SetTitle("#sigma_{e}[cm^{2}]");
  gPad->SetLogx();
  gPad->SetLogy();
    //=========================================CDMSlite=======================================
    TGraph *Graph_CDMSlite_ER_Lower_d1 = new TGraph(); Graph_CDMSlite_ER_Lower_d1->SetName("Graph_CDMSlite_ER_Lower_d1");
    int Row_CDMSlite = sizeof(CDMSlite_ER_Lower_d1)/sizeof(CDMSlite_ER_Lower_d1[0]);
    cout << "Row_CDMSlite: " << Row_CDMSlite << endl;
    for(int i=0;i<Row_CDMSlite;i++)
    {
        Graph_CDMSlite_ER_Lower_d1->SetPoint((i), CDMSlite_ER_Lower_d1[i][0], CDMSlite_ER_Lower_d1[i][1]);
        cout << "CDMSlite_ER_Lower_d1[i][0]: " << CDMSlite_ER_Lower_d1[i][0] <<  endl;
        cout << "CDMSlite_ER_Lower_d1[i][1]: " << CDMSlite_ER_Lower_d1[i][1] <<  endl;
    }
      Graph_CDMSlite_ER_Lower_d1->SetLineStyle(1);
      Graph_CDMSlite_ER_Lower_d1->SetLineColor(2);
      Graph_CDMSlite_ER_Lower_d1->SetLineWidth(8);
      Graph_CDMSlite_ER_Lower_d1->Draw("L");

    TLatex *Line_CDMSlite_Index = new TLatex(3,2e-32,"CDMSlite");
    Line_CDMSlite_Index->SetTextFont(42);
    Line_CDMSlite_Index->SetTextSize(0.04);
    Line_CDMSlite_Index->SetLineWidth(2);
    Line_CDMSlite_Index->SetTextColor(2);
    Line_CDMSlite_Index->SetTextAngle(5);
    Line_CDMSlite_Index->Draw("Lsame");

    //=========================================XENON1T(12e-)=======================================
    TGraph *Graph_XENON1T_12e_ER_Lower_d1 = new TGraph(); Graph_XENON1T_12e_ER_Lower_d1->SetName("Graph_XENON1T_12e_ER_Lower_d1");
    int Row_XENON1T_12e = sizeof(XENON1T_ER_Lower_d1_12e)/sizeof(XENON1T_ER_Lower_d1_12e[0]);
    cout << "Row_XENON1T_12e: " << Row_XENON1T_12e << endl;
    for(int i=0;i<Row_XENON1T_12e;i++)
    {
        Graph_XENON1T_12e_ER_Lower_d1->SetPoint((i), XENON1T_ER_Lower_d1_12e[i][0], XENON1T_ER_Lower_d1_12e[i][1]);
        cout << "XENON1T_ER_Lower_d1_12e[i][0]: " << XENON1T_ER_Lower_d1_12e[i][0] <<  endl;
        cout << "XENON1T_ER_Lower_d1_12e[i][1]: " << XENON1T_ER_Lower_d1_12e[i][1] <<  endl;
    }
      Graph_XENON1T_12e_ER_Lower_d1->SetLineStyle(1);
      Graph_XENON1T_12e_ER_Lower_d1->SetLineColor(7);
      Graph_XENON1T_12e_ER_Lower_d1->SetLineWidth(8);
      Graph_XENON1T_12e_ER_Lower_d1->Draw("Lsame");

    TLatex *Line_XENON1T_12e_Index = new TLatex(1,2e-34,"XENON1T(>12e-)");
    Line_XENON1T_12e_Index->SetTextFont(42);
    Line_XENON1T_12e_Index->SetTextSize(0.04);
    Line_XENON1T_12e_Index->SetLineWidth(2);
    Line_XENON1T_12e_Index->SetTextColor(7);
    Line_XENON1T_12e_Index->SetTextAngle(5);
    Line_XENON1T_12e_Index->Draw("Lsame");

    //=========================================XENON1T(12e-)=======================================
    TGraph *Graph_XENON1T_ER_Lower_d1 = new TGraph(); Graph_XENON1T_ER_Lower_d1->SetName("Graph_XENON1T_ER_Lower_d1");
    int Row_XENON1T = sizeof(XENON1T_ER_Lower_d1)/sizeof(XENON1T_ER_Lower_d1[0]);
    cout << "Row_XENON1T: " << Row_XENON1T << endl;
    for(int i=0;i<Row_XENON1T;i++)
    {
        Graph_XENON1T_ER_Lower_d1->SetPoint((i), XENON1T_ER_Lower_d1[i][0], XENON1T_ER_Lower_d1[i][1]);
        cout << "XENON1T_ER_Lower_d1[i][0]: " << XENON1T_ER_Lower_d1[i][0] <<  endl;
        cout << "XENON1T_ER_Lower_d1[i][1]: " << XENON1T_ER_Lower_d1[i][1] <<  endl;
    }
      Graph_XENON1T_ER_Lower_d1->SetLineStyle(1);
      Graph_XENON1T_ER_Lower_d1->SetLineColor(6);
      Graph_XENON1T_ER_Lower_d1->SetLineWidth(8);
      Graph_XENON1T_ER_Lower_d1->Draw("Lsame");

    TLatex *Line_XENON1T_Index = new TLatex(3,2e-35,"XENON1T");
    Line_XENON1T_Index->SetTextFont(42);
    Line_XENON1T_Index->SetTextSize(0.04);
    Line_XENON1T_Index->SetLineWidth(2);
    Line_XENON1T_Index->SetTextColor(6);
    Line_XENON1T_Index->SetTextAngle(5);
    Line_XENON1T_Index->Draw("Lsame");

  //===========================================CDEX_Lower=============================================
  TGraph *Graph_CDEX_ER_Lower_d1 = new TGraph(); Graph_CDEX_ER_Lower_d1->SetName("Graph_CDEX_ER_Lower_d1");
  int Row_CDEX_Lower = sizeof(CDEX_ER_Lower_d1)/sizeof(CDEX_ER_Lower_d1[0]);
  for(int i=0;i<Row_CDEX_Lower;i++)
  {
      Graph_CDEX_ER_Lower_d1->SetPoint((i), CDEX_ER_Lower_d1[i][0], CDEX_ER_Lower_d1[i][1]);
      cout << "CDEX_ER_Lower_d1.h[i][0]: " << CDEX_ER_Lower_d1[i][0] <<  endl;
      cout << "CDEX_ER_Lower_d1.h[i][1]: " << CDEX_ER_Lower_d1[i][1] <<  endl;
  }
    Graph_CDEX_ER_Lower_d1->SetLineStyle(1);
    Graph_CDEX_ER_Lower_d1->SetLineColor(3);
    Graph_CDEX_ER_Lower_d1->SetLineWidth(8);
    Graph_CDEX_ER_Lower_d1->Draw("Lsame");

    //===========================================CDEX_Upper=============================================
    TGraph *Graph_CDEX_ER_Upper_d1 = new TGraph(); Graph_CDEX_ER_Upper_d1->SetName("Graph_CDEX_ER_Upper_d1");
    int Row_CDEX_Upper = sizeof(CDEX_ER_Upper_d1)/sizeof(CDEX_ER_Upper_d1[0]);
    for(int i=0;i<Row_CDEX_Upper;i++)
    {
        Graph_CDEX_ER_Upper_d1->SetPoint((i), CDEX_ER_Upper_d1[i][0], CDEX_ER_Upper_d1[i][1]);
        cout << "CDEX_ER_Upper_d1.h[i][0]: " << CDEX_ER_Upper_d1[i][0] <<  endl;
        cout << "CDEX_ER_Upper_d1.h[i][1]: " << CDEX_ER_Upper_d1[i][1] <<  endl;
    }
      Graph_CDEX_ER_Upper_d1->SetLineStyle(2);
      Graph_CDEX_ER_Upper_d1->SetLineColor(3);
      Graph_CDEX_ER_Upper_d1->SetLineWidth(8);
      Graph_CDEX_ER_Upper_d1->Draw("Lsame");

    //===========================================CDEX_Boundaries=============================================
    TLine *Line_CDEX = new TLine(CDEX_ER_Lower_d1[Row_CDEX_Lower-1][0],CDEX_ER_Lower_d1[Row_CDEX_Lower-1][1],CDEX_ER_Upper_d1[Row_CDEX_Upper-1][0],CDEX_ER_Upper_d1[Row_CDEX_Upper-1][1]);
    Line_CDEX->SetLineColor(3);
    Line_CDEX->SetLineWidth(6);
    Line_CDEX->SetLineStyle(10);
    Line_CDEX->Draw("Lsame");
    TLatex *Line_CDEX_Index = new TLatex(3,2e-31,"CDEX-1B");
    Line_CDEX_Index->SetTextFont(42);
    Line_CDEX_Index->SetTextSize(0.04);
    Line_CDEX_Index->SetLineWidth(2);
    Line_CDEX_Index->SetTextColor(3);
    Line_CDEX_Index->SetTextAngle(5);
    Line_CDEX_Index->Draw("Lsame");

    
      //===========================================TEXONO_Lower=============================================
    TGraph *Graph_TEXONO_ER_Lower_d1 = new TGraph(); Graph_TEXONO_ER_Lower_d1->SetName("Graph_TEXONO_ER_Lower_d1");
    int Row_TEXONO_Lower = sizeof(TEXONO_ER_Lower_d1)/sizeof(TEXONO_ER_Lower_d1[0]);
    for(int i=0;i<Row_TEXONO_Lower;i++)
    {
        Graph_TEXONO_ER_Lower_d1->SetPoint((i), TEXONO_ER_Lower_d1[i][0], TEXONO_ER_Lower_d1[i][1]);
        cout << "TEXONO_ER_Lower_d1.h[i][0]: " << TEXONO_ER_Lower_d1[i][0] <<  endl;
        cout << "TEXONO_ER_Lower_d1.h[i][1]: " << TEXONO_ER_Lower_d1[i][1] <<  endl;
    }
      Graph_TEXONO_ER_Lower_d1->SetLineStyle(1);
      Graph_TEXONO_ER_Lower_d1->SetLineColor(4);
      Graph_TEXONO_ER_Lower_d1->SetLineWidth(8);
      Graph_TEXONO_ER_Lower_d1->Draw("Lsame");

    //===========================================TEXONO_Upper=============================================
    TGraph *Graph_TEXONO_ER_Upper_d1 = new TGraph(); Graph_TEXONO_ER_Upper_d1->SetName("Graph_TEXONO_ER_Upper_d1");
    int Row_TEXONO_Upper = sizeof(TEXONO_ER_Upper_d1)/sizeof(TEXONO_ER_Upper_d1[0]);
    for(int i=0;i<Row_TEXONO_Upper;i++)
    {
        Graph_TEXONO_ER_Upper_d1->SetPoint((i), TEXONO_ER_Upper_d1[i][0], TEXONO_ER_Upper_d1[i][1]);
        cout << "TEXONO_ER_Upper_d1.h[i][0]: " << TEXONO_ER_Upper_d1[i][0] <<  endl;
        cout << "TEXONO_ER_Upper_d1.h[i][1]: " << TEXONO_ER_Upper_d1[i][1] <<  endl;
    }
      Graph_TEXONO_ER_Upper_d1->SetLineStyle(2);
      Graph_TEXONO_ER_Upper_d1->SetLineColor(4);
      Graph_TEXONO_ER_Upper_d1->SetLineWidth(8);
      Graph_TEXONO_ER_Upper_d1->Draw("Lsame");

    //===========================================TEXONO_Boundaries=============================================
    TLine *Line_TEXONO = new TLine(TEXONO_ER_Lower_d1[Row_TEXONO_Lower-1][0],TEXONO_ER_Lower_d1[Row_TEXONO_Lower-1][1],TEXONO_ER_Upper_d1[Row_TEXONO_Upper-1][0],TEXONO_ER_Upper_d1[Row_TEXONO_Upper-1][1]);
    Line_TEXONO->SetLineColor(4);
    Line_TEXONO->SetLineWidth(6);
    Line_TEXONO->SetLineStyle(10);
    Line_TEXONO->Draw("Lsame");

    TLatex *Line_TEXONO_Index = new TLatex(3,2e-30,"TEXONO");
    Line_TEXONO_Index->SetTextFont(42);
    Line_TEXONO_Index->SetTextSize(0.04);
    Line_TEXONO_Index->SetLineWidth(2);
    Line_TEXONO_Index->SetTextColor(4);
    Line_TEXONO_Index->SetTextAngle(5);
    Line_TEXONO_Index->Draw("Lsame");

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

  plot->Print("Chie_d1.pdf");
}
