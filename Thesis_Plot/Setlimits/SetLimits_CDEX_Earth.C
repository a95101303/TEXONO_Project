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

void SetLimits_CDEX_Earth()
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
  TH2F *frame = new TH2F("frame","",100,0.05,15,100,4e-47,1e-22);
  //frame->GetXaxis()->SetTitle("WIMP Mass (GeV/c^{2})");
  //frame->GetYaxis()->SetTitle("SI Corss section (cm^{2})");
  //frame->GetXaxis()->SetTitle("m_{#chi} (GeV/c^{2})");
//frame->GetYaxis()->SetTitle("#sigma_{SI} (cm^{2})");

  //frame->GetXaxis()->CenterTitle();
  //frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetMoreLogLabels();
  frame->Draw();
  gPad->SetLogx();
  gPad->SetLogy();
    frame->GetXaxis()->SetTitle("M_{#chi}");
    frame->GetYaxis()->SetTitle("#sigma_{SI}");

//=================================================================================================
  TGraph *gCMB = new TGraph(); gCMB->SetName("gCMB"); 
  for(int i=0;i<8;i++)
  { gCMB->SetPoint((i), pow(10.0,CMB_chiN[i][0]), pow(10.0,CMB_chiN[i][1]));
      cout << "pow(10.0,CMB_chiN[i][0]): " << pow(10.0,CMB_chiN[i][0]) <<  endl;
      cout << "pow(10.0,CMB_chiN[i][1])): " << pow(10.0,CMB_chiN[i][1]) << endl;
  }
    gCMB->SetFillColor(kAzure+6);
    gCMB->SetFillStyle(3114);
    gCMB->Draw("f");

    TLatex *CMB_Text = new TLatex(0.055,2.5e-26,"CMB");
    CMB_Text->SetTextFont(42);
    CMB_Text->SetTextSize(0.04);
    CMB_Text->SetLineWidth(2);
    CMB_Text->SetTextColor(1);
    CMB_Text->SetTextAngle(5);
    //CMB_Text->Draw();
//=================================================================================================
  TGraph *gXQC = new TGraph(); gXQC->SetName("gXQC");
  for(int i=0;i<15;i++)
  { gXQC->SetPoint((i), pow(10.0,XQC_chiN[i][0]), pow(10.0,XQC_chiN[i][1]));
    cout << "pow(10.0,XQC_chiN[i][0]):  " << pow(10.0,XQC_chiN[i][0]) << endl;
    cout << "pow(10.0,XQC_chiN[i][1]):  " << pow(10.0,XQC_chiN[i][1]) << endl;
  }
    gXQC->SetFillColor(kYellow-7);
    gXQC->SetFillStyle(3114);
    gXQC->Draw("f");

    TLatex *XQC_Text = new TLatex(0.15,6e-25,"XQC");
    XQC_Text->SetTextFont(42);
    XQC_Text->SetTextSize(0.04);
    XQC_Text->SetLineWidth(2);
    XQC_Text->SetTextColor(1);
    XQC_Text->SetTextAngle(-40);
    //XQC_Text->Draw();
//==================================================================================================
    TGraph *gcresst_surf = new TGraph(); gcresst_surf->SetName("gcresst_surf");
    for(int i=0;i<34;i++)
    { gcresst_surf->SetPoint((i), pow(10.0,cresst_surf_chiN[i][0]), pow(10.0,cresst_surf_chiN[i][1])); }
      gcresst_surf->SetPoint((34), pow(10.0,cresst_surf_chiN[0][0]), pow(10.0,cresst_surf_chiN[0][1]));
      gcresst_surf->SetName("CDEX-1a");
      gcresst_surf->SetFillColor(kViolet-4);
      gcresst_surf->SetFillStyle(3144);
      gcresst_surf->Draw("f");

      TLatex *CRESST_Surface_text = new TLatex(0.15,3e-28,"CRESST(2017) Surface");
      CRESST_Surface_text->SetTextFont(42);
      CRESST_Surface_text->SetTextSize(0.04);
      CRESST_Surface_text->SetLineWidth(2);
      CRESST_Surface_text->SetTextColor(1);
      CRESST_Surface_text->SetTextAngle(-2);
      //CRESST_Surface_text->Draw();
//==================================================================================================
    TGraph *gcresstII = new TGraph(); gcresstII->SetName("gcresstII");
    for(int i=0;i<37;i++)
    { gcresstII->SetPoint((i), pow(10.0,cresstII_chiN[i][0]), pow(10.0,cresstII_chiN[i][1])); }
      gcresstII->SetPoint((37), pow(10.0,cresstII_chiN[0][0]), pow(10.0,cresstII_chiN[0][1]));
      gcresstII->SetName("CDEX-1a");
      gcresstII->SetFillColor(kGreen-4);
      gcresstII->SetFillStyle(3144);
      gcresstII->Draw("f");

      TLatex *CRESST_2_text = new TLatex(0.58,1e-36,"CRESST II");
      CRESST_2_text->SetTextFont(42);
      CRESST_2_text->SetTextSize(0.04);
      CRESST_2_text->SetLineWidth(2);
      CRESST_2_text->SetTextColor(1);
      CRESST_2_text->SetTextAngle(-30);
      //CRESST_2_text->Draw();
//==================================================================================================
    TGraph *gcresst_surf_Ours = new TGraph(); gcresst_surf_Ours->SetName("gcresst_surf_Ours");
    for(int i=0;i<16;i++)
    { gcresst_surf_Ours->SetPoint((i), cresst_surf_chiN_Ours[i][0], cresst_surf_chiN_Ours[i][1]); }
    gcresst_surf_Ours->SetLineColor(2);
    gcresst_surf_Ours->SetLineWidth(8);
    gcresst_surf_Ours->SetLineStyle(4);

    //gcresst_surf_Ours->Draw("l");
//============================
     

    

     
    
    /*
    TGraph *texono_surf = new TGraph(); texono_surf->SetName("texono_surf");
    for(int i=0;i<12;i++)
    { texono_surf->SetPoint((i), TEXONO_surf_chiN[i][0],TEXONO_surf_chiN[i][1]); }
    texono_surf->SetLineColor(6);
    texono_surf->SetLineWidth(8);
    texono_surf->SetLineStyle(5);
    texono_surf->Draw("l");
    
    
    
    TLatex *tex1 = new TLatex(3,1e-27,"TEXONO(Surface)");
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04);
    tex1->SetLineWidth(2);
    tex1->SetTextColor(6);
    tex1->Draw();
    */
    
    TFile *f=new TFile("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/Set_Limits_Plot_for_ALL.root");
    
    //CDEX-1a
    TGraph *CDEX_1a=(TGraph*)f->Get("CDEX-1a");
    CDEX_1a->SetName("CDEX-1a");
    CDEX_1a->SetFillColor(kYellow+0);
    CDEX_1a->SetLineWidth(5);
    CDEX_1a->SetLineColor(kOrange+7);
    CDEX_1a->SetFillStyle(3144);
    CDEX_1a->Draw("lfsame");

    double x1, y1;

     double CDEX_x[Number_of_Candidates], CDEX_y[Number_of_Candidates];
    for(int i=0;i<Number_of_Candidates;i++)
      {
        CDEX_1a->GetPoint(i,x1,y1);
        CDEX_x[i] = x1;
        CDEX_y[i] = y1;
          cout << "CDEX_x[i]: " << CDEX_x[i] << endl;
          cout << "CDEX_y[i]: " << CDEX_y[i] << endl;

      }

    TGraph *cdex_real = new TGraph(); cdex_real->SetName("cdex_real");
    cdex_real->SetPoint(0, CDEX_x[0],CDEX_y[0]);

    cout << "CDEX_x[0]: " << CDEX_x[0] << endl;
    cout << "CDEX_y[0]: " << CDEX_y[0] << endl;

    for(int i=0;i<12;i++)
    {
        cout << "CDEX_real_chiN[i][0]: " << CDEX_real_chiN[i][0] << endl;
        cout << "CDEX_real_chiN[i][1]: " << CDEX_real_chiN[i][1] << endl;

        cdex_real->SetPoint((i+1), CDEX_real_chiN[i][0],CDEX_real_chiN[i][1]);
    }
    cdex_real->SetPoint(13, CDEX_x[18],CDEX_y[18]);
    cout << "CDEX_x[20]: " << CDEX_x[18] << endl;
    cout << "CDEX_y[20]: " << CDEX_y[18] << endl;

    cdex_real->SetName("cdex_real");
    cdex_real->SetFillColor(kYellow+0);
    cdex_real->SetFillStyle(3144);
    cdex_real->SetLineColor(kOrange+10);
    cdex_real->SetLineWidth(5);

    cdex_real->Draw("lfsame");

    TLatex *tex10 = new TLatex(2,1e-40,"CDEX-1b(NU)");
    tex10->SetTextFont(42);
    tex10->SetTextSize(0.04);
    tex10->SetLineWidth(2);
    tex10->SetTextColor(kCyan-3);
    tex10->SetTextAngle(-20);
    tex10->Draw();


     
    TGraph *cdex_real_Migdal = new TGraph(); cdex_real_Migdal->SetName("cdex_real_MD");

    for(int i=0;i<17;i++)
    {
        cout << "CDEX_real_Migdal_Lower[i][0]: " << CDEX_real_Migdal_Lower[i][0] << endl;
        cout << "CDEX_real_Migdal_Lower[i][1]: " << CDEX_real_Migdal_Lower[i][1] << endl;

        cdex_real_Migdal->SetPoint((i), CDEX_real_Migdal_Lower[i][0],CDEX_real_Migdal_Lower[i][1]);
    }

    for(int i=0;i<17;i++)
    {
        cout << "CDEX_real_Migdal_Upper[i][0]: " << CDEX_real_Migdal_Upper[i][0] << endl;
        cout << "CDEX_real_Migdal_Upper[i][1]: " << CDEX_real_Migdal_Upper[i][1] << endl;

        cdex_real_Migdal->SetPoint((i+16), CDEX_real_Migdal_Upper[i][0],CDEX_real_Migdal_Upper[i][1]);
    }

    cdex_real_Migdal->SetPoint((32), CDEX_real_Migdal_Lower[0][0],CDEX_real_Migdal_Lower[0][1]);

    cdex_real_Migdal->SetFillColor(kCyan+4);
    cdex_real_Migdal->SetLineColor(kOrange+10);
    cdex_real_Migdal->SetLineWidth(5);
    cdex_real_Migdal->Draw("lfsame");
    cdex_real_Migdal->SetFillStyle(3144);

    TLatex *tex5;
    tex5 = new TLatex(0.05,1e-33,"CDEX-1b(Migdal Effect)");
    tex5->SetTextFont(42);
    tex5->SetTextSize(0.04);
    tex5->SetLineWidth(2);
    tex5->SetTextColor(kCyan+4);
    tex5->SetTextAngle(-15);
    tex5->Draw();
    
     
    
    TGraph *texono_real_Migdal = new TGraph(); texono_real_Migdal->SetName("texono_real_MD");

    for(int i=0;i<15;i++)
    {
        cout << "TEXONO_real_Migdal_Lower[i][0]: " << TEXONO_real_Migdal_Lower[i][0] << endl;
        cout << "TEXONO_real_Migdal_Lower[i][1]: " << TEXONO_real_Migdal_Lower[i][1] << endl;

        texono_real_Migdal->SetPoint((i), TEXONO_real_Migdal_Lower[i][0],TEXONO_real_Migdal_Lower[i][1]);
    }
    
    for(int i=0;i<15;i++)
    {
        cout << "CDEX_real_Migdal_Upper[i][0]: " << TEXONO_real_Migdal_Upper[i][0] << endl;
        cout << "CDEX_real_Migdal_Upper[i][1]: " << TEXONO_real_Migdal_Upper[i][1] << endl;

        texono_real_Migdal->SetPoint((i+15), TEXONO_real_Migdal_Upper[i][0],TEXONO_real_Migdal_Upper[i][1]);
    }
    
     
    texono_real_Migdal->SetFillColor(kCyan-4);
    texono_real_Migdal->SetLineColor(1);
    texono_real_Migdal->SetLineWidth(5);
    //texono_real_Migdal->Draw("fsame");
    texono_real_Migdal->SetFillStyle(3144);

    TLatex *tex6;
    tex6 = new TLatex(0.05,5e-35,"TEXONO(Migdal Effect)");
    tex6->SetTextFont(42);
    tex6->SetTextSize(0.04);
    tex6->SetLineWidth(2);
    tex6->SetTextColor(kCyan-4);
    //tex6->Draw();

    /*
    TGraph *texono_real_Migdal_all_bins = new TGraph(); texono_real_Migdal_all_bins->SetName("texono_real_MD_all_Bins");

    for(int i=0;i<15;i++)
    {
        cout << "TEXONO_real_Migdal_Lower_All_Bins[i][0]: " << TEXONO_real_Migdal_Lower_All_Bins[i][0] << endl;
        cout << "TEXONO_real_Migdal_Lower_All_Bins[i][1]: " << TEXONO_real_Migdal_Lower_All_Bins[i][1] << endl;

        texono_real_Migdal_all_bins->SetPoint((i), TEXONO_real_Migdal_Lower_All_Bins[i][0],TEXONO_real_Migdal_Lower_All_Bins[i][1]);
    }

    texono_real_Migdal_all_bins->SetFillColor(kCyan-4);
    texono_real_Migdal_all_bins->SetLineColor(1);
    texono_real_Migdal_all_bins->SetLineWidth(5);
    texono_real_Migdal_all_bins->Draw("Lsame");
    texono_real_Migdal_all_bins->SetFillStyle(3144);

    */
    
    TGraph *texono_real_brem = new TGraph(); texono_real_brem->SetName("texono_real_BR");

    for(int i=0;i<18;i++)
    {
        cout << "TEXONO_real_Brem_Lower[i][0]: " << TEXONO_real_Brem_Lower[i][0] << endl;
        cout << "TEXONO_real_Brem_Lower[i][1]: " << TEXONO_real_Brem_Lower[i][1] << endl;

        texono_real_brem->SetPoint((i), TEXONO_real_Brem_Lower[i][0],TEXONO_real_Brem_Lower[i][1]);
    }
    
    for(int i=0;i<18;i++)
    {
        cout << "TEXONO_real_Brem_Lower[i][0]: " << TEXONO_real_Brem_Upper[i][0] << endl;
        cout << "TEXONO_real_Brem_Lower[i][1]: " << TEXONO_real_Brem_Upper[i][1] << endl;

        texono_real_brem->SetPoint((i+18), TEXONO_real_Brem_Upper[i][0],TEXONO_real_Brem_Upper[i][1]);
    }

     
    texono_real_brem->SetFillColor(kGreen+2);
    texono_real_brem->SetLineColor(1);
    texono_real_brem->SetLineWidth(5);
    //texono_real_brem->Draw("fsame");
    texono_real_brem->SetFillStyle(3144);
   // texono_real_brem->SetMarkerStyle(8);
   // texono_real_brem->SetMarkerSize(0.5);
   // texono_real_brem->SetMarkerColor(2);

    TLatex *tex7;
    tex7 = new TLatex(0.05,1e-29,"TEXONO(Brem)");
    tex7->SetTextFont(42);
    tex7->SetTextSize(0.04);
    tex7->SetLineWidth(2);
    tex7->SetTextColor(kGreen+2);
    //tex7->Draw();
    
    /*
    TGraph *texono_real_Migdal = new TGraph(); texono_real_Migdal->SetName("texono_real");

    for(int i=0;i<16;i++)
    {
        cout << "TEXONO_real_Migdal_Lower[i][0]: " << TEXONO_real_Migdal_Lower[i][0] << endl;
        cout << "TEXONO_real_Migdal_Lower[i][1]: " << TEXONO_real_Migdal_Lower[i][1] << endl;

        texono_real_Migdal->SetPoint((i+1), TEXONO_real_Migdal_Lower[i][0],TEXONO_real_Migdal_Lower[i][1]);
    }
     */


       // DAMA legend
       //CoGent AM legend
       TLatex *tex;
       TPave *pave = new TPave(3.2,5.6e-42,3.8,8e-42,4,"br");
       pave->SetFillColor(kBlue-4);
       pave->SetLineColor(kBlue-4);
       pave->SetLineWidth(0);
       pave->SetFillStyle(3001);
       pave->SetShadowColor(0);
       //pave->Draw();

       pave = new TPave(3.35,6.15e-42,3.635,7.44e-42,4,"br");
       pave->SetFillColor(kBlue);
       pave->SetLineColor(kBlue);
       pave->SetFillStyle(1001);
       pave->SetBorderSize(1);
       pave->SetShadowColor(0);
       //pave->Draw();

       //tex = new TLatex(4.0,6.15e-39,"DAMA/LIBRA phase-1, Na-recoil:  5 #sigma & 90%");
         tex = new TLatex(3.0,6.15e-39,"DAMA/LIBRA phase-1");
       tex->SetTextFont(42);
       tex->SetTextSize(0.04);
       tex->SetLineWidth(2);
       tex->SetTextColor(kBlue);
       //tex->Draw();

       tex = new TLatex(11.2,6.15e-42,"5-#sigma,");
       tex->SetTextFont(42);
       tex->SetTextSize(0.04);
       tex->SetLineWidth(2);
       tex->SetTextColor(kBlue);
       //tex->Draw();


       tex = new TLatex(13.2,6.15e-42,"90%");
       tex->SetTextFont(42);
       tex->SetTextSize(0.04);
       tex->SetLineWidth(2);
       tex->SetTextColor(kBlue);
       //tex->Draw();
     //
       //DAMA2009 allowed region
       TGraph *g_dama2009 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/dama2009.txt", "%lg  %lg");
       g_dama2009->SetName("g_dama2009");
       g_dama2009->SetFillColor(kBlue);
       g_dama2009->SetFillStyle(1001);
       g_dama2009->Draw("f");

         
       //DAMA2009 allowed region
       TGraph *g_dama1 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/dama_up_band.txt", "%lg  %lg");
       g_dama1->SetName("g_dama1");
       g_dama1->SetFillColor(kBlue);
       g_dama1->SetFillStyle(1001);
       g_dama1->Draw("f");

       TGraph *g_damal_90 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/DAMA_UP90_RL.txt", "%lg  %lg");
       g_damal_90->SetName("g_damal_90");
       g_damal_90->SetFillColor(kBlue);
       g_damal_90->SetFillStyle(1001);
       g_damal_90->Draw("f");

       //DAMA2009 allowed region
       TGraph *g_dama2 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/dama_down_band.txt", "%lg  %lg");
       g_dama2->SetName("g_dama2");
       g_dama2->SetFillColor(kBlue);
       g_dama2->SetFillStyle(1001);
       g_dama2->Draw("f");

       TGraph *g_damah_90 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/DAMA_Down90_RL.txt", "%lg  %lg");
       g_damah_90->SetName("g_damah_90");
       g_damah_90->SetFillColor(kBlue);
       g_damah_90->SetFillStyle(1001);
       g_damah_90->Draw("f");


        
       //CDMS-II Si legend
       pave = new TPave(3.2,1.2e-42,3.8,1.8e-42,4,"br");
       pave->SetFillColor(kAzure+1);
       pave->SetLineColor(kAzure+1);
       pave->SetFillStyle(3001);
       pave->SetShadowColor(0);
       //pave->Draw();

       tex = new TLatex(4.0,1.3e-42,"CDMS-II Si");
       tex->SetTextFont(42);
       tex->SetTextSize(0.03440367);
       tex->SetLineWidth(2);
       tex->SetTextColor(kAzure+1);
       //tex->Draw();
     //
       //CDMS-II Si 68 allowed
       double cdms_si_allow68_x[cdms_si_allow68_bin], cdms_si_allow68_y[cdms_si_allow68_bin];
       for(int i=0;i<cdms_si_allow68_bin;i++)
         { cdms_si_allow68_x[i] = cdms_si_allow68[i][1];
           cdms_si_allow68_y[i] = cdms_si_allow68[i][2];
         }
       TGraph *g_cdms_si_allow68 = new TGraph(cdms_si_allow68_bin,cdms_si_allow68_x,cdms_si_allow68_y);
       g_cdms_si_allow68->SetName("g_cdms_si_allow68");
       g_cdms_si_allow68->SetFillColor(kAzure+1);
       g_cdms_si_allow68->SetFillStyle(3001);
       //g_cdms_si_allow68->Draw("f");
       
       
     //===============================================================
       //CoGent AM legend

       pave = new TPave(3.2,0.103e-40,3.8,0.15e-40,4,"br");
       pave->SetFillColor(kGreen+2);
       pave->SetLineColor(kGreen+2);
       pave->SetShadowColor(0);
       //pave->Draw();

         tex = new TLatex(3.6,3.5e-42,"CoGeNT(2013)");
       tex->SetTextFont(42);
       tex->SetTextSize(0.04);
       tex->SetTextColor(kMagenta-4);
       tex->SetLineWidth(2);
       //tex->Draw();

       //CoGent AM allowed region
      double cogent2013_x[cogent2013_bin], cogent2013_y[cogent2013_bin];
     for(int i=0;i<cogent2013_bin;i++)
       { cogent2013_x[i] = cogent2013[i][1];
         cogent2013_y[i] = cogent2013[i][2];
         
       }
      
      TGraph *g_cogent2013 = new TGraph(cogent2013_bin,cogent2013_x,cogent2013_y);
      g_cogent2013->SetName("g1_cogent2013");
      g_cogent2013->SetFillColor(kMagenta-4);
      g_cogent2013->Draw("f");

      //===============================================================

      //texono2013 result
      TGraph *g_texono2013 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/texono2013.txt","%lg  %lg");
      g_texono2013->SetName("g_texono2013");
      g_texono2013->SetLineWidth(5);
      g_texono2013->SetLineColor(1);
      //g_texono2013->Draw("c");
      
      tex = new TLatex(10.0,4.7e-41,"/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO (2013)");
      tex->SetTextFont(42);
      tex->SetTextSize(0.03422619);
      tex->SetTextAngle(3.2);
      tex->SetLineWidth(2);
      //tex->Draw();

      
      //DAMIC2014
      TGraph *g_DAMIC = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/DAMIC.txt","%lg  %lg");
      g_DAMIC->SetName("g_DAMIC");
      g_DAMIC->SetLineColor(kRed-2);
      g_DAMIC->SetLineWidth(5);
      //g_DAMIC->Draw("l");

      tex = new TLatex(1.75,9.0e-39,"DAMIC");
      tex->SetTextColor(kRed-2);
      tex->SetTextFont(42);
      tex->SetTextSize(0.03422619);
      tex->SetTextAngle(338);
      tex->SetLineWidth(2);
      //tex->Draw();

      //superCDMS2014
      TGraph *g_superCDMS2014 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/superCDMS2014.txt","%lg  %lg");
      g_superCDMS2014->SetName("g_superCDMS2014");
      g_superCDMS2014->SetLineColor(kCyan-2);
      g_superCDMS2014->SetLineWidth(3);
      //g_superCDMS2014->Draw("c");

      tex = new TLatex(7.4,2.5e-42,"SuperCDMS (2014)");
      tex->SetTextColor(kCyan-2);
      tex->SetTextFont(42);
      tex->SetTextSize(0.03);
      tex->SetTextAngle(308);
      tex->SetLineWidth(2);
      //tex->Draw();


      //CDMS2015
      TGraph *g_CDMS2015 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDMS_2015.txt","%lg  %lg");
      g_CDMS2015->SetName("g_CDMS2015");
      g_CDMS2015->SetLineColor(kYellow+4);
      g_CDMS2015->SetLineWidth(3);
      //g_CDMS2015->Draw("c");
       
      tex = new TLatex(3.0,2e-41,"CDMSlite (2016)");
      tex->SetTextColor(kYellow+4);
      tex->SetTextFont(42);
      tex->SetTextSize(0.03);
      tex->SetTextAngle(338);
      tex->SetLineWidth(2);
      //tex->Draw();


      //Cresst2015
      TGraph *g_Cresst2015 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/Cresst_SI_2015.txt","%lg  %lg");
      g_Cresst2015->SetName("g_Cresst2015");
      g_Cresst2015->SetLineColor(kMagenta);
      g_Cresst2015->SetLineWidth(5);
      //g_Cresst2015->Draw("l");
      
      tex = new TLatex(2.15,7.5e-40,"CRESST-II (2016)");
      tex->SetTextColor(6);
      tex->SetTextFont(42);
      tex->SetTextSize(0.03422619);
      tex->SetTextAngle(305);
      tex->SetLineWidth(2);
      // tex->Draw();


      
      //PandaX2016
      TGraph *g_PandaX2016 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/PandaX_20170217_v2.txt","%lg  %lg");
      g_PandaX2016->SetName("g_PandaX2016");
      g_PandaX2016->SetLineColor(kGreen+3);
      g_PandaX2016->SetLineWidth(5);
      //g_PandaX2016->Draw("l");
      
      tex = new TLatex(4.3,7.0e-43,"PandaX (2016)");
      tex->SetTextColor(kGreen+3);
      tex->SetTextFont(42);
      tex->SetTextSize(0.03422619);
      tex->SetTextAngle(290);
      tex->SetLineWidth(2);
      //tex->Draw();



      //Lux2016
      TGraph *g_Lux2016 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/Lux_20170217_cm.txt","%lg  %lg");
      g_Lux2016->SetName("g_Lux2016");
      g_Lux2016->SetLineColor(kOrange+7);
      g_Lux2016->SetLineWidth(5);
      //g_Lux2016->Draw("l");
      
      tex = new TLatex(5.5,7.0e-43,"LUX (2016)");
      tex->SetTextColor(kOrange+7);
      tex->SetTextFont(42);
      tex->SetTextSize(0.03422619);
      tex->SetTextAngle(290);
      tex->SetLineWidth(2);
      //tex->Draw();

      //cdex2016
         /*
      TGraph *g_cdex2016 = new TGraph("CDEX_2016.txt","%lg  %lg");
      g_cdex2016->SetName("g_cdex2016");
      g_cdex2016->SetLineColor(kBlack);
      g_cdex2016->SetLineWidth(8);
      g_cdex2016->Draw("l");
          */
         /*
      tex = new TLatex(7.5,2.8e-43,"CDEX-1a");
      tex->SetTextColor(kBlack);
      tex->SetTextFont(62);
      tex->SetTextSize(0.03422619);
      tex->SetTextAngle(15);
      tex->SetLineWidth(2);
      tex->Draw();
          */
      //cdex2016 new BS
      TGraph *g_cdex2016new = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/c1_lihb_20170217.dat","%lg  %lg");
      g_cdex2016new->SetName("g_cdex2016new");
      g_cdex2016new->SetLineColor(kRed+1);
      g_cdex2016new->SetLineWidth(8);
      //g_cdex2016new->Draw("l");

      tex = new TLatex(3.5,5.0e-41,"CDEX-1 Modulation (this work)");
      tex->SetTextColor(kRed+1);
      tex->SetTextFont(40);
      tex->SetTextSize(0.04);
      tex->SetTextAngle(310);
      tex->SetLineWidth(2);
      //tex->Draw();


      TGraph *g_npc = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/Bounds_mass_cx_fit_npc_50eV.txt", "%lg %lg %*lg");
      g_npc->SetName("g_npc");
      g_npc->SetLineWidth(5);
      g_npc->SetLineColor(kRed);
      //g_npc->Draw("l");
      
      tex = new TLatex(4.5,5e-40,"NPC-TEXONO");
      tex->SetTextColor(2);
      tex->SetTextFont(42);
      tex->SetTextSize(0.03422619);
      tex->SetLineWidth(2);
      //tex->Draw();

     /////////////////////////
     // XMASS 2018
     /////////////////////////
      TGraph *xmass = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/xmass2018.txt", "%lg %lg");
      xmass->SetName("xmass");
      xmass->SetLineWidth(5);
      xmass->SetLineColor(kGray);
     // xmass->Draw("c");

      tex = new TLatex(5.6,7e-40,"XMASS 2018");
      tex->SetTextColor(kGray);
      tex->SetTextFont(42);
      tex->SetTextSize(0.05);
      tex->SetTextAngle(280);
      tex->SetLineWidth(2);
     // tex->Draw();

     /////////////////////////
     // XMASS 2018 v2
     /////////////////////////
      TGraph *xmass2 = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/xmass_2018_v2.txt", "%lg %lg");
      xmass2->SetName("xmass2");
      xmass2->SetLineWidth(5);
      xmass2->SetLineColor(kGray+2);
      //xmass2->Draw("c");

      tex = new TLatex(4.2,7e-40,"XMASS 2018");
      tex->SetTextColor(kGray+2);
      tex->SetTextFont(42);
      tex->SetTextSize(0.05);
      tex->SetTextAngle(283);
      tex->SetLineWidth(2);
      //tex->Draw();


     /////////////////////////
     // C1B spectrum
     /////////////////////////
         /*
      //TGraph *c1b_2017_cpc = new TGraph("c1b_2017_cpc.txt","%lg  %lg");
      TGraph *c1b_2017_cpc = new TGraph("CDEX1B-limit-SI.txt","%lg  %lg");
      c1b_2017_cpc->SetName("c1b_2017_cpc");
      c1b_2017_cpc->SetLineWidth(5);
      c1b_2017_cpc->SetLineColor(kBlack);
      c1b_2017_cpc->SetLineStyle(9);
      c1b_2017_cpc->Draw("l");

      tex = new TLatex(3.0,5.2e-41,"CDEX-1B unmodulated : Run 1");
      tex->SetTextColor(kGreen);
      tex->SetTextFont(62);
      tex->SetTextSize(0.03);
      tex->SetTextAngle(319);
      tex->SetLineWidth(2);
      //tex->Draw();

      TGraph *c1b_2017_cpc0 = new TGraph("CDEX1B-limit-SI_SIDM_all.txt","%lg  %lg");
      c1b_2017_cpc0->SetName("c1b_2017_cpc0");
      c1b_2017_cpc0->SetLineWidth(5);
      c1b_2017_cpc0->SetLineColor(kRed);
      c1b_2017_cpc0->SetLineStyle(1);
      c1b_2017_cpc0->Draw("l");

      tex = new TLatex(3.0,5.2e-42,"CDEX-1B cortrected");
      tex->SetTextColor(kGreen);
      tex->SetTextFont(62);
      tex->SetTextSize(0.03);
      tex->SetTextAngle(319);
      tex->SetLineWidth(2);
      //tex->Draw();
     */

     /////////////////////////
     // C1B annual modulation
     /////////////////////////
     /*
      tex = new TLatex(10,1.4e-41,"C1B-AM : old");
      tex->SetTextColor(kRed);
      tex->SetTextFont(62);
      tex->SetTextSize(0.035);
      tex->SetTextAngle(15);
      tex->SetLineWidth(2);
      //tex->Draw();

      tex = new TLatex(8.2,2e-41,"CDEX-1B Modulation (this work)");
      tex->SetTextColor(kRed);
      tex->SetTextFont(60);
      tex->SetTextSize(0.046);
      tex->SetTextAngle(13);
      tex->SetLineWidth(2);
      //tex->Draw();

     //
      TGraph *g_v60_53E_p1525 = new TGraph("simu_v60_53E_p1525.txt","%lg  %lg");
      g_v60_53E_p1525->SetName("g_v60_53E_p1525");
      g_v60_53E_p1525->SetLineWidth(6);
      g_v60_53E_p1525->SetLineColor(kGray);
      g_v60_53E_p1525->SetLineStyle(1);
      //g_v60_53E_p1525->Draw("c");

      TGraph *g_v60_sys_kl_53E_p1525 = new TGraph("simu_v60_sys_kl_53E_p1525.txt","%lg  %lg");
      g_v60_sys_kl_53E_p1525->SetName("g_v60_sys_kl_53E_p1525");
      g_v60_sys_kl_53E_p1525->SetLineWidth(6);
      g_v60_sys_kl_53E_p1525->SetLineColor(kBlue-1);
      g_v60_sys_kl_53E_p1525->SetLineStyle(2);
      //g_v60_sys_kl_53E_p1525->Draw("c");

     //
      TGraph *g_v60_sys_trimm_53E_p1525 = new TGraph("simu_v60_sys_trimm_53E_p1525.txt","%lg  %lg");
      g_v60_sys_trimm_53E_p1525->SetName("g_v60_sys_trimm_53E_p1525");
      g_v60_sys_trimm_53E_p1525->SetLineWidth(6);
      g_v60_sys_trimm_53E_p1525->SetLineColor(kRed-1);
      g_v60_sys_trimm_53E_p1525->SetLineStyle(2);
      //g_v60_sys_trimm_53E_p1525->Draw("c");

     //
      TGraph *g_v60_sys_trimp_53E_p1525 = new TGraph("simu_v60_sys_trimp_53E_p1525.txt","%lg  %lg");
      g_v60_sys_trimp_53E_p1525->SetName("g_v60_sys_trimp_53E_p1525");
      g_v60_sys_trimp_53E_p1525->SetLineWidth(6);
      g_v60_sys_trimp_53E_p1525->SetLineColor(kGreen-1);
      g_v60_sys_trimp_53E_p1525->SetLineStyle(2);
      //g_v60_sys_trimp_53E_p1525->Draw("c");

     //
      TGraph *g_v60_53E_p1525_nonflat = new TGraph("simu_v60_53E_p1525_nonflat.txt","%lg  %lg");
      g_v60_53E_p1525_nonflat->SetName("g_v60_53E_p1525_nonflat");
      g_v60_53E_p1525_nonflat->SetLineWidth(6);
      g_v60_53E_p1525_nonflat->SetLineColor(kGray);
      g_v60_53E_p1525_nonflat->SetLineStyle(1);
     //g_v60_53E_p1525_nonflat->Draw("c");

     //
      TGraph *g_v60_nosys_53E_p1525 = new TGraph("simu_v60_nosys_53E_p1525.txt","%lg  %lg");
      g_v60_nosys_53E_p1525->SetName("g_v60_nosys_53E_p1525");
      g_v60_nosys_53E_p1525->SetLineWidth(6);
      g_v60_nosys_53E_p1525->SetLineColor(kGray);
      g_v60_nosys_53E_p1525->SetLineStyle(1);
     //g_v60_nosys_53E_p1525->Draw("c");

     //
     const int reso_mx = 500;
     double mass[reso_mx];
     double data_v60_53E_p1525[reso_mx];
     double data_v60_sys_kl_53E_p1525[reso_mx];
     double data_v60_sys_trimm_53E_p1525[reso_mx];
     double data_v60_sys_trimp_53E_p1525[reso_mx];
     double data_v60_sys_nonflat_53E_p1525[reso_mx];
     double data_v60_sys_nosys_53E_p1525[reso_mx];

     double data_sys[reso_mx];

     double x, y0, y1, y2, y3, y4, y5;

     for(int i=0;i<reso_mx;i++)
     {
       g_v60_53E_p1525->GetPoint(i,x,y0);
       mass[i] = x;
       data_v60_53E_p1525[i] = y0;
     //
       g_v60_sys_kl_53E_p1525->GetPoint(i,x,y1);
       data_v60_sys_kl_53E_p1525[i] = y1;
     //
       g_v60_sys_trimm_53E_p1525->GetPoint(i,x,y2);
       data_v60_sys_trimm_53E_p1525[i] = y2;
     //
       g_v60_sys_trimp_53E_p1525->GetPoint(i,x,y3);
       data_v60_sys_trimp_53E_p1525[i] = y3;

       g_v60_53E_p1525_nonflat->GetPoint(i,x,y4);
       data_v60_sys_nonflat_53E_p1525[i] = y4;

       g_v60_nosys_53E_p1525->GetPoint(i,x,y5);
       data_v60_sys_nosys_53E_p1525[i] = y5;

     //printf("%d %f %e %e %e %e\n",i,mass[i],data_v60_53E_p1525[i],data_v60_sys_kl_53E_p1525[i],data_v60_sys_trimm_53E_p1525[i],data_v60_sys_trimp_53E_p1525[i]);
     }

     for(int i=0;i<reso_mx;i++)
     {
       data_sys[i] = 0.0;
       if(data_sys[i]<data_v60_53E_p1525[i]) { data_sys[i] = data_v60_53E_p1525[i]; }
       if(data_sys[i]<data_v60_sys_kl_53E_p1525[i]) { data_sys[i] = data_v60_sys_kl_53E_p1525[i]; }
       if(data_sys[i]<data_v60_sys_trimp_53E_p1525[i]) { data_sys[i] = data_v60_sys_trimp_53E_p1525[i]; }
       if(data_sys[i]<data_v60_sys_trimm_53E_p1525[i]) { data_sys[i] = data_v60_sys_trimm_53E_p1525[i]; }
       if(data_sys[i]<data_v60_sys_nonflat_53E_p1525[i]) { data_sys[i] = data_v60_sys_nonflat_53E_p1525[i]; }
       if(data_sys[i]<data_v60_sys_nosys_53E_p1525[i]) { data_sys[i] = data_v60_sys_nosys_53E_p1525[i]; }

       printf("%f %e\n",mass[i],data_sys[i]);
     }

     TGraph *g_sys = new TGraph(reso_mx,mass,data_sys);
     g_sys->SetName("g_sys");
     g_sys->SetLineWidth(5);
     g_sys->SetLineColor(kRed);
     g_sys->SetLineStyle(1);
     //g_sys->Draw("c");

     //
      TGraph *g_v51_sys_53E_p1525 = new TGraph("simu_v51_sys_53E_p1525.txt","%lg  %lg");
      g_v51_sys_53E_p1525->SetName("g_v51_sys_53E_p1525");
      g_v51_sys_53E_p1525->SetLineWidth(4);
      g_v51_sys_53E_p1525->SetLineColor(kRed-7);
      g_v51_sys_53E_p1525->SetLineStyle(9);
     // g_v51_sys_53E_p1525->Draw("c");
     //

     //g_sys->Draw("c");

     //
      tex = new TLatex(12.5,1.3e-41,"CDEX-1B time-integrated");
      tex->SetTextColor(kBlack);
      tex->SetTextFont(60);
      tex->SetTextSize(0.04);
      tex->SetTextAngle(9);
      tex->SetLineWidth(2);
      //tex->Draw();

      TGraph *nlimit_cjpl_linear = new TGraph("natural_limit_linear.txt","%lg  %lg");
      nlimit_cjpl_linear->SetName("nlimit_cjpl_linear");
     // nlimit_cjpl_linear->Draw("c");

      TGraph *nlimit_cjpl_real = new TGraph("natural_limit_real.txt","%lg  %lg");
      nlimit_cjpl_real->SetName("nlimit_cjpl_real");
     // nlimit_cjpl_real->Draw("c");


      TGraph *nlimit_cjpl_linear_np2 = new TGraph("natural_limit_linear_np2.txt","%lg  %lg");
      nlimit_cjpl_linear_np2->SetName("nlimit_cjpl_linear_np2");
     // nlimit_cjpl_linear_np2->Draw("c");

      TGraph *nlimit_cjpl_real_np2 = new TGraph("natural_limit_real_np2.txt","%lg  %lg");
      nlimit_cjpl_real_np2->SetName("nlimit_cjpl_real_np2");
     // nlimit_cjpl_real_np2->Draw("c");

         */
         /*
      TGraph *c1b_migdal_raw = new TGraph("c1b_migdal_raw.txt","%lg  %lg");
      c1b_migdal_raw->SetName("c1b_migdal_raw");
      c1b_migdal_raw->Draw("c");
         */
         
      TGraph *c1b_migdal_adjust = new TGraph("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/c1b_migdal_adjust.txt","%lg  %lg");
      c1b_migdal_adjust->SetName("c1b_migdal_adjust");
     c1b_migdal_adjust->SetLineWidth(5);
         c1b_migdal_adjust->SetLineColor(2);

      //c1b_migdal_adjust->Draw("l");
          
     /*
      TGraph *migdal_max = new TGraph("migdal_max.txt","%lg  %lg");
      migdal_max->SetName("migdal_max");
     // migdal_max->Draw("c");

      TGraph *np2_project = new TGraph("np2_limit.txt","%lg  %lg");
      np2_project->SetName("np2_project");
      np2_project->Draw("l");

      TGraph *atm_project = new TGraph();
      atm_project->SetName("atm_project");
      for(int i=0;i<9;i++)
      {
        np2_project->GetPoint(i,x,y0);
        atm_project->SetPoint(i,x,7.0*y0);
      }
      atm_project->Draw("l");
     */
           leg->AddEntry(gCMB,"CMB","f");
       leg->AddEntry(gXQC,"XQC","f");
           leg->AddEntry(gcresst_surf,"CRESST(2017) Surface","f");
       leg->AddEntry(gcresstII,"CRESST II","f");
     leg->AddEntry(g_damah_90,"DAMA2009","f");
    leg->AddEntry(g_cogent2013,"CoGeNT(2013)","f");
       leg->Draw();
    
      plot->Print("CDEX_Earth.pdf");
}
