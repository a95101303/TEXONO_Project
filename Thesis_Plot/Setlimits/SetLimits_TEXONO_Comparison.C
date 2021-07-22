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

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Brem_Lower_1st_Original.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Comparison/TEXONO_real_Brem_Lower_1st.h"

void SetLimits_TEXONO_Comparison()
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
    TH2F *frame = new TH2F("frame","",100,0.05,2.34,100,1e-37,1e-26);
  //frame->GetXaxis()->SetTitle("WIMP Mass (GeV/c^{2})");
  //frame->GetYaxis()->SetTitle("SI Corss section (cm^{2})");
  //frame->GetXaxis()->SetTitle("m_{#chi} (GeV/c^{2})");
//frame->GetYaxis()->SetTitle("#sigma_{SI} (cm^{2})");

  //frame->GetXaxis()->CenterTitle();
  //frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitle("M_{#chi}[GeV]");
    frame->GetYaxis()->SetTitle("#sigma_{SI}[cm^{2}]");
  frame->GetXaxis()->SetMoreLogLabels();
  frame->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

     

    

     
    


//TEXONO
    int Work_Line_Width=5;


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
    TLegend *leg = new TLegend(0.15,0.14,0.4,0.35);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);


    
    
     TGraph *texono_vacuum_brem = new TGraph(); texono_vacuum_brem->SetName("texono_vacuum_brem");

     for(int i=0;i<16;i++)
     {
         cout << "TEXONO_real_Brem_Lower_1st[i][0]: " << TEXONO_real_Brem_Lower_1st_Original[i][0] << endl;
         cout << "TEXONO_real_Brem_Lower_1st[i][1]: " << TEXONO_real_Brem_Lower_1st_Original[i][1] << endl;

         texono_vacuum_brem->SetPoint((i), TEXONO_real_Brem_Lower_1st_Original[i][0],TEXONO_real_Brem_Lower_1st_Original[i][1]);
     }
     
      
     texono_vacuum_brem->SetLineColor(kOrange+10);
     texono_vacuum_brem->SetLineStyle(5);
     texono_vacuum_brem->SetLineWidth(Work_Line_Width);
     texono_vacuum_brem->SetFillStyle(3144);
     texono_vacuum_brem->Draw("lsame");

     leg->AddEntry(texono_vacuum_brem,"Vacuum-constrained Brem","l");

     TLatex *tex_BREM;
     tex_BREM = new TLatex(0.05,5e-28,"TEXONO(Brem)");
     tex_BREM->SetTextFont(42);
     tex_BREM->SetTextSize(0.04);
     tex_BREM->SetLineWidth(2);
     tex_BREM->SetTextColor(kOrange+10);
     tex_BREM->SetTextAngle(-20);
     //tex_BREM->Draw();
     
     
     TGraph *texono_real_brem = new TGraph(); texono_real_brem->SetName("texono_real_BR");

     for(int i=0;i<19;i++)
     {
         cout << "TEXONO_real_Brem_Lower[i][0]: " << TEXONO_real_Brem_Lower_1st[i][0] << endl;
         cout << "TEXONO_real_Brem_Lower[i][1]: " << TEXONO_real_Brem_Lower_1st[i][1] << endl;

         texono_real_brem->SetPoint((i), TEXONO_real_Brem_Lower_1st[i][0],TEXONO_real_Brem_Lower_1st[i][1]);
     }
     
      
     texono_real_brem->SetLineColor(kOrange+10);
     texono_real_brem->SetLineStyle(1);
     texono_real_brem->SetLineWidth(Work_Line_Width);
     texono_real_brem->Draw("lsame");
     texono_real_brem->SetFillStyle(3144);
    // texono_real_brem->SetMarkerStyle(8);
    // texono_real_brem->SetMarkerSize(0.5);
    // texono_real_brem->SetMarkerColor(2);

    leg->AddEntry(texono_real_brem,"Earth-constrained Brem","l");

     TLatex *tex7;
     tex7 = new TLatex(0.1,1e-29,"TEXONO(Brem)");
     tex7->SetTextFont(42);
     tex7->SetTextSize(0.04);
     tex7->SetLineWidth(2);
     tex7->SetTextColor(kGreen+2);
     tex7->SetTextAngle(-15);
     //tex7->Draw();

     TGraph *texono_vacuum_Migdal = new TGraph(); texono_vacuum_Migdal->SetName("texono_vacuum_Migdal");

     for(int i=0;i<16;i++)
     {
         cout << "TEXONO_real_Migdal_Lower[i][0]: " << TEXONO_real_Migdal_Lower_Original[i][0] << endl;
         cout << "TEXONO_real_Migdal_Lower[i][1]: " << TEXONO_real_Migdal_Lower_Original[i][1] << endl;

         texono_vacuum_Migdal->SetPoint((i), TEXONO_real_Migdal_Lower_Original[i][0],TEXONO_real_Migdal_Lower_Original[i][1]);
     }
     texono_vacuum_Migdal->SetLineStyle(5);
     texono_vacuum_Migdal->SetLineColor(kBlue+1);
     texono_vacuum_Migdal->SetLineWidth(Work_Line_Width);
     texono_vacuum_Migdal->SetFillStyle(3144);
     texono_vacuum_Migdal->Draw("lsame");

    leg->AddEntry(texono_vacuum_Migdal,"Vacuum-constrained MD","l");

     TLatex *tex_CDEX_MD;
     tex_CDEX_MD = new TLatex(0.05,3e-32,"TEXONO(Migdal)");
     tex_CDEX_MD->SetTextFont(42);
     tex_CDEX_MD->SetTextSize(0.04);
     tex_CDEX_MD->SetLineWidth(2);
     tex_CDEX_MD->SetTextColor(kOrange+10);
     tex_CDEX_MD->SetTextAngle(-20);
     //tex_CDEX_MD->Draw();


    TGraph *texono_real_Migdal = new TGraph(); texono_real_Migdal->SetName("texono_real_MD");

    for(int i=0;i<16;i++)
    {
        cout << "TEXONO_real_Migdal_Lower[i][0]: " << TEXONO_real_Migdal_Lower[i][0] << endl;
        cout << "TEXONO_real_Migdal_Lower[i][1]: " << TEXONO_real_Migdal_Lower[i][1] << endl;

        texono_real_Migdal->SetPoint((i), TEXONO_real_Migdal_Lower[i][0],TEXONO_real_Migdal_Lower[i][1]);
    }
    
    texono_real_Migdal->SetLineStyle(1);
    texono_real_Migdal->SetLineColor(kBlue+1);
    texono_real_Migdal->SetLineWidth(Work_Line_Width);
    texono_real_Migdal->Draw("lfsame");
    texono_real_Migdal->SetFillStyle(3144);

    leg->AddEntry(texono_real_Migdal,"Earth-constrained MD","l");

    TLatex *tex6;
    tex6 = new TLatex(0.05,5e-33,"TEXONO(Migdal Effect)");
    tex6->SetTextFont(42);
    tex6->SetTextSize(0.04);
    tex6->SetLineWidth(2);
    tex6->SetTextColor(kCyan-4);
    tex6->SetTextAngle(-15);
    //tex6->Draw();


    leg->Draw();

    // DAMA legend
    //CoGent AM legend
      plot->Print("TEXONO_Comparison.pdf");
}
