#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraph2D.h"

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"

void TGraph2D_Test()//CDEX:Threshold=160eV, TEXONO:Threshold=200eV
{
    TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,600,400);
    TGraph2D *dt = new TGraph2D();
    
    dt->SetPoint(0,1,1,2);
    dt->SetPoint(1,2,2,10);

    gStyle->SetPalette(1);
    dt->SetMarkerStyle(20);
    dt->Draw("pcol");
}

