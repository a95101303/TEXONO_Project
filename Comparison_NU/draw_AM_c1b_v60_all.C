#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "TCutG.h"
#include "TCut.h"
#include "cdms_si_allow68.h"
#include "superCDMS2014_band.h"
#include "cogent2013.h"

void draw_AM_c1b_v60_all()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gStyle->SetFrameLineWidth(2);	
  //gStyle->SetLabelOffset(0.007,"Y");
  //gStyle->SetLabelSize(0.04, "Y");
  //gStyle->SetNdivisions(510,"Y");
  //gStyle->SetLabelSize(0.04, "X");
  //gStyle->SetNdivisions(510,"X");
  //gStyle->SetTitleSize(0.045, "X" );
  //gStyle->SetTitleSize(0.05, "Y" );
  //gStyle->SetTitleFont(42,"X");
  //gStyle->SetTitleFont(42,"Y");
  //gStyle->SetLabelFont(42,"X");
  //gStyle->SetLabelFont(42,"Y");
  //gStyle->SetTitleOffset(1.0,"X");
  //gStyle->SetTitleOffset(1.5,"Y");
  TFile *fin=new TFile("/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/Codes/1.root");
    
  TGraph *tr=(TGraph*)fin->Get("Graph");
    
    TFile *fin1=new TFile("/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/Codes/2.root");
    
    TGraph *tr1=(TGraph*)fin1->Get("Graph");

    TFile *fin2=new TFile("/Users/ms08962476/Report/AS_GeIA/Analysis/CDEX_Analysis_method/Codes/Neutrino_floor.root");
    
    TGraph *tr2=(TGraph*)fin2->Get("Graph");

    
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

  TH2F *frame = new TH2F("frame","",100,1,1e04,100,1e-50,1e-36);
//  TH2F *frame = new TH2F("frame","",100,3.0,28,100,4e-42,9e-40);
  //frame->GetXaxis()->SetTitle("WIMP Mass (GeV/c^{2})");
  //frame->GetYaxis()->SetTitle("SI Corss section (cm^{2})");
  frame->GetXaxis()->SetTitle("m_{#chi} (GeV/c^{2})");
  frame->GetYaxis()->SetTitle("#sigma_{SI} (cm^{2})");
  //frame->GetXaxis()->CenterTitle();
  //frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetMoreLogLabels();
  frame->Draw();
  gPad->SetLogx();
  gPad->SetLogy();

  // DAMA legend
  //CoGent AM legend

  TLatex *tex;
  TPave *pave = new TPave(3.2,5.6e-38,3.8,8e-38,4,"br");
  pave->SetFillColor(kBlue-4);
  pave->SetLineColor(kBlue-4);
  pave->SetLineWidth(0);
  pave->SetFillStyle(3001);
  pave->SetShadowColor(0);
  pave->Draw();

  pave = new TPave(3.35,6.15e-38,3.635,7.44e-38,4,"br");
  pave->SetFillColor(kBlue);
  pave->SetLineColor(kBlue);
  pave->SetFillStyle(1001);
  pave->SetBorderSize(1);
  pave->SetShadowColor(0);
  pave->Draw();

  tex = new TLatex(4.0,6.15e-38,"DAMA/LIBRA ");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->SetTextColor(kBlue-4);
  tex->Draw();

  tex = new TLatex(10,6.15e-38,"90%");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->SetTextColor(kBlue);
  tex->Draw();
//
  //DAMA2009 allowed region
  TGraph *g_dama2009 = new TGraph("dama2009.txt", "%lg  %lg");
  g_dama2009->SetName("g_dama2009");
  g_dama2009->SetFillColor(kBlue-4);
  g_dama2009->SetFillStyle(3001);
  g_dama2009->Draw("f");

  //DAMA2009 allowed region
  TGraph *g_dama1 = new TGraph("dama_up_band.txt", "%lg  %lg");
  g_dama1->SetName("g_dama1");
  g_dama1->SetFillColor(kBlue-4);
  g_dama1->SetFillStyle(3001);
  g_dama1->Draw("f");

  TGraph *g_damal_90 = new TGraph("DAMA_UP90_RL.txt", "%lg  %lg");
  g_damal_90->SetName("g_damal_90");
  g_damal_90->SetFillColor(kBlue);
  g_damal_90->SetFillStyle(1001);
  g_damal_90->Draw("f");

  //DAMA2009 allowed region
  TGraph *g_dama2 = new TGraph("dama_down_band.txt", "%lg  %lg");
  g_dama2->SetName("g_dama2");
  g_dama2->SetFillColor(kBlue-4);
  g_dama2->SetFillStyle(3001);
  g_dama2->Draw("f");

  TGraph *g_damah_90 = new TGraph("DAMA_Down90_RL.txt", "%lg  %lg");
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
  
  
  //CoGent AM legend

  pave = new TPave(3.2,1.8e-38,3.8,2.2e-38,4,"br");
  pave->SetFillColor(kGreen+2);
  pave->SetLineColor(kGreen+2);
  pave->SetShadowColor(0);
  pave->Draw();

  tex = new TLatex(4.0,2e-38,"CoGeNT 90%");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetTextColor(kGreen+2);
  tex->SetLineWidth(2);
  tex->Draw();

  //CoGent AM allowed region
 double cogent2013_x[cogent2013_bin], cogent2013_y[cogent2013_bin];
for(int i=0;i<cogent2013_bin;i++) 
  { cogent2013_x[i] = cogent2013[i][1]; 
    cogent2013_y[i] = cogent2013[i][2]; 
    
  }
 
 TGraph *g_cogent2013 = new TGraph(cogent2013_bin,cogent2013_x,cogent2013_y);
 g_cogent2013->SetName("g1_cogent2013");
 g_cogent2013->SetFillColor(kGreen+2);
 g_cogent2013->SetLineColor(1);
 g_cogent2013->SetLineWidth(5);
 g_cogent2013->Draw("f");

 //texono2013 result
 TGraph *g_texono2013 = new TGraph("texono2013.txt","%lg  %lg");
 g_texono2013->SetName("g_texono2013");
 g_texono2013->SetLineWidth(5);
 g_texono2013->SetLineColor(1);
 //g_texono2013->Draw("c");
 
 tex = new TLatex(10.0,4.7e-41,"TEXONO (2013)");
 tex->SetTextFont(42);
 tex->SetTextSize(0.03422619);
 tex->SetTextAngle(3.2);
 tex->SetLineWidth(2);
 //tex->Draw();

 
 //DAMIC2014
 TGraph *g_DAMIC = new TGraph("DAMIC.txt","%lg  %lg");
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
 TGraph *g_superCDMS2014 = new TGraph("superCDMS2014.txt","%lg  %lg");
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
 TGraph *g_CDMS2015 = new TGraph("CDMS_2015.txt","%lg  %lg");
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
 TGraph *g_Cresst2015 = new TGraph("Cresst_SI_2015.txt","%lg  %lg");
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
 //tex->Draw();


 
 //PandaX2016
 TGraph *g_PandaX2016 = new TGraph("PandaX_20170217_v2.txt","%lg  %lg");
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
 TGraph *g_Lux2016 = new TGraph("Lux_20170217_cm.txt","%lg  %lg");
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
 TGraph *g_cdex2016 = new TGraph("CDEX_2016.txt","%lg  %lg");
 g_cdex2016->SetName("g_cdex2016");
 g_cdex2016->SetLineColor(kBlack);
 g_cdex2016->SetLineWidth(8);
 //g_cdex2016->Draw("l");
 
 tex = new TLatex(15.5,1e-42,"CDEX-1a");
 tex->SetTextColor(kBlack);
 tex->SetTextFont(62);
 tex->SetTextSize(0.03422619); 
 tex->SetTextAngle(15);
 tex->SetLineWidth(2);
 tex->Draw();

 //cdex2016 new BS
 TGraph *g_cdex2016new = new TGraph("c1_lihb_20170217.dat","%lg  %lg");
 g_cdex2016new->SetName("g_cdex2016new");
 g_cdex2016new->SetLineColor(kRed+1);
 g_cdex2016new->SetLineWidth(8);
 //g_cdex2016new->Draw("l");

 tex = new TLatex(3.5,5.0e-41,"CDEX-1 (this work)");
 tex->SetTextColor(kRed+1);
 tex->SetTextFont(42);
 tex->SetTextSize(0.04);
 tex->SetTextAngle(310);
 tex->SetLineWidth(2);
 //tex->Draw();


 TGraph *g_npc = new TGraph("Bounds_mass_cx_fit_npc_50eV.txt", "%lg %lg %*lg");
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
 TGraph *xmass = new TGraph("xmass2018.txt", "%lg %lg");
 xmass->SetName("xmass");
 xmass->SetLineWidth(5);
 xmass->SetLineColor(kGray);
 //xmass->Draw("c");

 tex = new TLatex(5.6,7e-40,"XMASS 2018");
 tex->SetTextColor(kGray);
 tex->SetTextFont(42);
 tex->SetTextSize(0.05);
 tex->SetTextAngle(280);
 tex->SetLineWidth(2);
 //tex->Draw();

/////////////////////////
// XMASS 2018 v2
/////////////////////////
 TGraph *xmass2 = new TGraph("xmass_2018_v2.txt", "%lg %lg");
 xmass2->SetName("xmass2");
 xmass2->SetLineWidth(5);
 xmass2->SetLineColor(kGray+2);
 //xmass2->Draw("c");

 tex = new TLatex(4.2,7e-40,"XMASS 2018 v2");
 tex->SetTextColor(kGray+2);
 tex->SetTextFont(42);
 tex->SetTextSize(0.05);
 tex->SetTextAngle(283);
 tex->SetLineWidth(2);
 //tex->Draw();


/////////////////////////
// C1B spectrum
/////////////////////////
    double x11;double x12; double y21;double y22;
 //TGraph *c1b_2017_cpc = new TGraph("c1b_2017_cpc.txt","%lg  %lg");
 TGraph *c1b_2017_cpc = new TGraph("CDEX1B-limit-SI.txt","%lg  %lg");
 c1b_2017_cpc->SetName("c1b_2017_cpc");
 c1b_2017_cpc->SetLineWidth(5);
 c1b_2017_cpc->SetLineColor(kBlack);
 c1b_2017_cpc->SetLineStyle(9);
 c1b_2017_cpc->Draw("l");
 c1b_2017_cpc->GetPoint(8,x11,y21);
    
    cout << "x1: " << x11 << endl;
    cout << "y21: " << y21 << endl;

 tex = new TLatex(3.0,5.2e-41,"CDEX-1B unmodulated : Run 1");
 tex->SetTextColor(kGreen);
 tex->SetTextFont(62);
 tex->SetTextSize(0.03);
 tex->SetTextAngle(319);
 tex->SetLineWidth(2);
 //tex->Draw();

/////////////////////////
// C1B annual modulation
/////////////////////////

 tex = new TLatex(10,1.4e-41,"C1B-AM : old");
 tex->SetTextColor(kRed);
 tex->SetTextFont(62);
 tex->SetTextSize(0.035);
 tex->SetTextAngle(15);
 tex->SetLineWidth(2);
 //tex->Draw();

 tex = new TLatex(1.2,1.9e-41,"CDEX-1B Modulation");
 tex->SetTextColor(kRed);
 tex->SetTextFont(62);
 tex->SetTextSize(0.046);
 tex->SetTextAngle(8);
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
const int reso_mx = 500;
double mass[reso_mx];
double data_v60_53E_p1525[reso_mx];
double data_v60_sys_kl_53E_p1525[reso_mx];
double data_v60_sys_trimm_53E_p1525[reso_mx];
double data_v60_sys_trimp_53E_p1525[reso_mx];
double data_sys[reso_mx];

double x, y0, y1, y2, y3;

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
//printf("%d %f %e %e %e %e\n",i,mass[i],data_v60_53E_p1525[i],data_v60_sys_kl_53E_p1525[i],data_v60_sys_trimm_53E_p1525[i],data_v60_sys_trimp_53E_p1525[i]);
}

for(int i=0;i<reso_mx;i++)
{
  data_sys[i] = 0.0;
  if(data_sys[i]<data_v60_sys_trimm_53E_p1525[i]) { data_sys[i] = data_v60_sys_trimm_53E_p1525[i]; }
  if(data_sys[i]<data_v60_53E_p1525[i]) { data_sys[i] = data_v60_53E_p1525[i]; }
  if(data_sys[i]<data_v60_sys_kl_53E_p1525[i]) { data_sys[i] = data_v60_sys_kl_53E_p1525[i]; }
  //if(data_sys[i]<data_v60_sys_trimm_53E_p1525[i]) { data_sys[i] = data_v60_sys_trimm_53E_p1525[i]; }
  if(data_sys[i]<data_v60_sys_trimp_53E_p1525[i]) { data_sys[i] = data_v60_sys_trimp_53E_p1525[i]; }
  //printf("%f %e\n",mass[i],data_sys[i]);
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
 //g_v51_sys_53E_p1525->Draw("c");
//

//g_sys->Draw("c");
    TLatex *tex1 = new TLatex(15,0.8e-40,"TEXONO");
    tex1->SetTextColor(kOrange);
    tex1->SetTextFont(62);
    tex1->SetTextSize(0.04);
    tex1->SetTextAngle(10);
    tex1->Draw();

    TLatex *tex2 = new TLatex(15,0.8e-49,"Neutrino Floor");
    tex2->SetTextColor(kBlue);
    tex2->SetTextFont(62);
    tex2->SetTextSize(0.04);
    tex2->Draw();

    tr->SetLineColor(kOrange);
    tr->SetLineWidth(3);
    tr->Draw("PLsame");
    
    tr1->SetLineColor(kRed);
    tr1->SetLineWidth(3);
    tr1->Draw("PLsame");

    tr2->SetLineColor(kBlue);
    tr2->SetLineWidth(3);
    tr2->Draw("PLsame");

    tr->GetPoint(1,x12,y22);
    cout << "x12: " << x12 << endl;
    cout << "y22: " << y22 << endl;

    TLine *line1 = new TLine(0,y21,x11,y21);
    TLine *line2 = new TLine(0,y22,x12,y22);
    line1->Draw("Lsame");
    line2->Draw("Lsame");
  
    plot->Print("112.png");

//
 tex = new TLatex(10.5,0.5e-41,"CDEX-1B unmodulated");
 tex->SetTextColor(kBlack);
 tex->SetTextFont(62);
 tex->SetTextSize(0.04);
 tex->SetTextAngle(10);
 tex->SetLineWidth(2);
 //tex->Draw();

}
