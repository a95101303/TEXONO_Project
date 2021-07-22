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
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_real_Migdal_Lower_Original.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_real_Migdal_Upper.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/CDEX_real_Brem_Lower_Original.h"

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Lower.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Upper.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Lower_Original.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Brem_Lower.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Brem_Upper.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Migdal_Lower_All_Bins.h"

void SetLimits_CDEX_Comparison()
{
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
    TH2F *frame = new TH2F("frame","",100,0.05,2.07,100,1e-38,1e-31);
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
  frame->GetXaxis()->SetTitle("M_{#chi}[GeV]");
  frame->GetYaxis()->SetTitle("#sigma_{SI}[cm^{2}]");

    

     
    
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


//TEXONO
    TLegend *leg = new TLegend(0.15,0.14,0.4,0.35);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);

    int Work_Line_Width=5;

    TGraph *cdex_real_Migdal = new TGraph(); cdex_real_Migdal->SetName("cdex_real_MD");

    for(int i=0;i<16;i++)
    {
        cout << "CDEX_real_Migdal_Lower[i][0]: " << CDEX_real_Migdal_Lower[i][0] << endl;
        cout << "CDEX_real_Migdal_Lower[i][1]: " << CDEX_real_Migdal_Lower[i][1] << endl;

        cdex_real_Migdal->SetPoint((i), CDEX_real_Migdal_Lower[i][0],CDEX_real_Migdal_Lower[i][1]);
    }


    cdex_real_Migdal->SetLineStyle(1);
    cdex_real_Migdal->SetLineColor(2);
    cdex_real_Migdal->SetLineWidth(Work_Line_Width);
    cdex_real_Migdal->Draw("lsame");
    cdex_real_Migdal->SetFillStyle(3144);

    leg->AddEntry(cdex_real_Migdal,"Earth-constrained MD","l");

    TLatex *tex5;
    tex5 = new TLatex(0.2,1e-36,"CDEX-1b(Migdal Effect)");
    tex5->SetTextFont(42);
    tex5->SetTextSize(0.04);
    tex5->SetLineWidth(2);
    tex5->SetTextColor(kRed+3);
    //tex5->Draw();
    
     TGraph *cdex_vacuum_Migdal = new TGraph(); cdex_vacuum_Migdal->SetName("cdex_vacuum_MD");

     for(int i=0;i<16;i++)
     {
         cout << "CDEX_real_Migdal_Lower_Original[i][0]: " << CDEX_real_Migdal_Lower_Original[i][0] << endl;
         cout << "CDEX_real_Migdal_Lower_Original[i][1]: " << CDEX_real_Migdal_Lower_Original[i][1] << endl;

         cdex_vacuum_Migdal->SetPoint((i), CDEX_real_Migdal_Lower_Original[i][0],CDEX_real_Migdal_Lower_Original[i][1]);
     }

     cdex_vacuum_Migdal->SetLineStyle(5);
     cdex_vacuum_Migdal->SetLineWidth(5);
     cdex_vacuum_Migdal->SetLineColor(2);
     cdex_vacuum_Migdal->SetLineWidth(Work_Line_Width);
     cdex_vacuum_Migdal->Draw("Lsame");
     cdex_vacuum_Migdal->SetFillStyle(3144);

     leg->AddEntry(cdex_vacuum_Migdal,"Vacuum-constrained MD","l");

     TLatex *tex_CDEX_MD;
     tex_CDEX_MD = new TLatex(0.05,1e-33,"CDEX-1b(Migdal Effect)");
     tex_CDEX_MD->SetTextFont(42);
     tex_CDEX_MD->SetTextSize(0.04);
     tex_CDEX_MD->SetLineWidth(2);
     tex_CDEX_MD->SetTextColor(kOrange+10);
     tex_CDEX_MD->SetTextAngle(-15);
     //tex_CDEX_MD->Draw();

        

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
    
    
    /*
    TGraph *texono_real_Migdal = new TGraph(); texono_real_Migdal->SetName("texono_real");

    for(int i=0;i<16;i++)
    {
        cout << "TEXONO_real_Migdal_Lower[i][0]: " << TEXONO_real_Migdal_Lower[i][0] << endl;
        cout << "TEXONO_real_Migdal_Lower[i][1]: " << TEXONO_real_Migdal_Lower[i][1] << endl;

        texono_real_Migdal->SetPoint((i+1), TEXONO_real_Migdal_Lower[i][0],TEXONO_real_Migdal_Lower[i][1]);
    }
     */

    leg->Draw();
      plot->Print("CDEX_Comparison.pdf");
}
