#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"

void Sensitivity_Line_for_CRESST_Overlap()
{
    string Exp_Flux[2]={"3_CRESST_Flux","4_CRESST_Flux"};
    string Mass_Point[4]={"20","10","2","0P2"};
    double GetPoint_X_Bef[2];double GetPoint_Y_Bef[2];double Slope_Earth_Bef;
    double GetPoint_X_Aft[2];double GetPoint_Y_Aft[2];double Slope_Earth_Aft;

    for(int Exp=1; Exp<2; Exp++)
    {
        for(int Mass_INT=0; Mass_INT<4; Mass_INT++)
        {
            TGraph *NU_STS; TGraph *NU_STS_Bent;TGraph *NU_STS_Earth;TGraph *NU_STS_Earth_Bent;
            
        TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/NU_STS_.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS=(TGraph*)fin2->Get("Threshold_Plot");
                NU_STS->SetLineColor(2);
                NU_STS->SetLineStyle(1);

        TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/NU_STS_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS_Bent=(TGraph*)fin3->Get("Threshold_Plot");
                NU_STS_Bent->SetLineColor(2);
                NU_STS_Bent->SetLineStyle(5);
        TFile *fin4 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/NU_STS_Earth_.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS_Earth=(TGraph*)fin4->Get("Threshold_Plot");
                NU_STS_Earth->SetLineColor(4);
                NU_STS_Earth->SetLineStyle(1);

            for(Int_t i=0; i<2; i++) {

            NU_STS_Earth->GetPoint(i,GetPoint_X_Bef[i],GetPoint_Y_Bef[i]);
            cout<<i<<"th element of X array: "<<GetPoint_X_Bef[i]<<endl;
            cout<<i<<"th element of Y array: "<<GetPoint_Y_Bef[i]<<endl;
            }
            Slope_Earth_Bef = (GetPoint_Y_Bef[1]-GetPoint_Y_Bef[0])/(GetPoint_X_Bef[1]-GetPoint_X_Bef[0]);
            cout << "Slope_Earth_Bef: " << Slope_Earth_Bef << endl;
            cout << "X=0, Y= " << GetPoint_Y_Bef[0]-GetPoint_X_Bef[0]*Slope_Earth_Bef << endl;
            
        TFile *fin5 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/NU_STS_Earth_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS_Earth_Bent=(TGraph*)fin5->Get("Threshold_Plot");
                NU_STS_Earth_Bent->SetLineColor(4);
                NU_STS_Earth_Bent->SetLineStyle(5);

            for(Int_t i=0; i<2; i++) {

            NU_STS_Earth_Bent->GetPoint(i,GetPoint_X_Aft[i],GetPoint_Y_Aft[i]);
            cout<<i<<"th element of X array: "<<GetPoint_X_Aft[i]<<endl;
            cout<<i<<"th element of Y array: "<<GetPoint_Y_Aft[i]<<endl;
            }
            Slope_Earth_Aft = (GetPoint_Y_Aft[1]-GetPoint_Y_Aft[0])/(GetPoint_X_Aft[1]-GetPoint_X_Aft[0]);
            cout << "Slope_Earth_Aft: " << Slope_Earth_Aft << endl;
            cout << "X=0, Y= " << GetPoint_Y_Aft[0]-GetPoint_X_Aft[0]*Slope_Earth_Aft << endl;

             

            TCanvas *c3 = new TCanvas("c3");
            gStyle->SetOptFit(0);
            gStyle->SetOptStat(0);
            cout << "YES4" << endl;

            TLegend *leg = new TLegend(0.1,0.1,0.4,0.3);
            leg->SetFillColor(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            leg->SetBorderSize(0);
            leg->SetTextFont(22);

            /*
            MD->Draw("ALP");
             */
            NU_STS->Draw("ALP");
            NU_STS_Bent->Draw("LPsame");
            NU_STS_Earth->Draw("LPsame");
            NU_STS_Earth_Bent->Draw("LPsame");
            
            leg->Draw();
            leg->AddEntry("",Form("%sGeV",Mass_Point[Mass_INT].c_str()),"l");
            leg->AddEntry("","Solid Line: Without bending","l");
            leg->AddEntry("","Dash  Line: With    bending","l");

            TF1 *Linear_Line     = new TF1("Linear_Line","x",1e-40,1e-28);
            
            TF1 *Linear_Line_Bef= new TF1("Linear_Line_Bef","[1]*x+[2]",GetPoint_X_Bef[0],GetPoint_X_Bef[1]);
            Linear_Line_Bef->SetParameter(1,Slope_Earth_Bef);
            Linear_Line_Bef->SetParameter(2,GetPoint_Y_Bef[0]-GetPoint_X_Bef[0]*Slope_Earth_Bef );
            Linear_Line_Bef->SetLineColor(8);
            //Linear_Line_Bef->Draw("Lsame");

            
            TF1 *Linear_Line_Aft= new TF1("Linear_Line_Aft","[1]*x+[2]",GetPoint_X_Aft[0],GetPoint_X_Aft[1]);
            Linear_Line_Aft->SetParameter(1,Slope_Earth_Aft);
            Linear_Line_Aft->SetParameter(2,GetPoint_Y_Aft[0]-GetPoint_X_Aft[0]*Slope_Earth_Aft );
            Linear_Line_Aft->SetLineColor(9);
            //Linear_Line_Aft->Draw("Lsame");

            
            Linear_Line->SetLineColor(3);
            Linear_Line->Draw("Lsame");

            cout << "YES7" << endl;

            cout << "YES8" << endl;
            c3->SetLogy();
            cout << "YES9" << endl;
            c3->SetLogx();
            cout << "YES10" << endl;
            
            c3->Print(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/All_%sGeV_STS.pdf",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
        }
    }
}
    

