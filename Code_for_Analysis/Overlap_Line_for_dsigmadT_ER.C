#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLine.h"
#include "TLine.h"
#include "dsigma_dT2.h"

void Overlap_Line_for_dsigmadT_ER()
{
        /*
        cout << scientific;
        cout << "c_1(9e-42,1): " << c_1(9e-42,0.5) << endl;
        cout << "d_1(4e-34,0.5): " << d_1(4e-34,0.5) << endl;
    
        cout << "c_1(4e-34,0.5): " << c_1(4e-34,0.5) << endl;//c_1(4e-34,0.5) = 3.519178e+00
        cout << "CS_Try(3.519178e+00,0.5): " << CS_Try(3.519178e+00,0.5) << endl;//

        cout << "CS_Try(5.276072e-04,0.5): " << CS_Try(5.276072e-04,0.5) << endl;
    
        cout << "DS_Try(1e-9,0.5): " << DS_Try(1e-9,1) << endl;
        cout << "d_1(DS_Try(1e-9,1),0.5): " << d_1(DS_Try(1e-9,1),1) << endl;
        cout << "c_1(DS_Try(1e-9,1),0.5): " << c_1(DS_Try(1e-9,1),1) << endl;
        cout << "CS_Try(c_1(DS_Try(1e-9,1),1),0.5): " << CS_Try(c_1(DS_Try(1e-9,1),1),1) << endl;//
         */
    
        
        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        cout << "YES4" << endl;
        
        double Mx            = 1;//GeV/c^2
        double Atomic_Mass   = AXe;//
        double velocity_c    = 100.*1e3/3e8;//Velocity(c)
        double Cro_Sec       = DS_Try(1e-9,Mx);//(Sigma_e)
        cout << "Cro_Sec: " << Cro_Sec << endl;
        TF1 *fdsigma_dT_keV_FUNC = new TF1 ("fdsigma_dT_keV_FUNC", "fdsigma_dT_keV([0],[1],[2],[3],x)", 0., max_recoil_A_keV(Mx,velocity_c,Atomic_Mass));
        fdsigma_dT_keV_FUNC->SetParameter( 0, Mx);
        fdsigma_dT_keV_FUNC->SetParameter( 1, Cro_Sec);
        fdsigma_dT_keV_FUNC->SetParameter( 2, velocity_c);
        fdsigma_dT_keV_FUNC->SetParameter( 3, Atomic_Mass);
        fdsigma_dT_keV_FUNC->SetLineColor(4);
        cout << "N_Nucleus: " << 1e5*(1.8/((unified_atomic_mass_g*(ASi))))*total_Sigma(1,100,Cro_Sec,Mx,AXe) << endl;
    //cout << "fdsigma_dT_keV_FUNC->GetMean(): " << fdsigma_dT_keV_FUNC->GetMean() << endl;
    
        cout << "fdsigma_dT_keV_FUNC->Eval(): " << fdsigma_dT_keV_FUNC->Eval(0.12) << endl;
        cout << "max_recoil_A_keV(Mx,velocity_c,Atomic_Mass): " << max_recoil_A_keV(Mx,velocity_c,Atomic_Mass) << endl;
        TFile *fin2 = TFile::Open("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/SI_c1_XeData_Vel/c1_1_0GeV/DCS.root");
        TH1F  *c1_V1   = (TH1F*)fin2->Get("03333");
        c1_V1->Scale(7*1e-36*1e-15);
        c1_V1->GetYaxis()->SetRangeUser(1e-40,1e-10);
        c1_V1->SetLineColor(2);
    
        TFile *fin3 = TFile::Open("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/SI_d1_XeData_Vel/d1_1_0GeV/DCS.root");
        TH1F  *d1_V1   = (TH1F*)fin3->Get("03333");
        d1_V1->Scale(1e-18*1e-15);
        d1_V1->SetLineColor(3);

        TLegend *leg = new TLegend(0.1,0.1,0.4,0.3);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);

        leg->AddEntry(c1_V1,"c1","l");
        leg->AddEntry(d1_V1,"d1","l");
        leg->AddEntry("","100(km/s)","");

        c3->Draw();
        c3->SetLogy();
        c1_V1->Draw("HISTsame");
        d1_V1->Draw("HISTsame");
        fdsigma_dT_keV_FUNC->Draw("Lsame");

        leg->Draw("same");
         
}

        /*
            TGraph *NU_STS; TGraph *NU_STS_Bent;TGraph *NU_STS_Earth;TGraph *NU_STS_Earth_Bent;
            
        TFile *fin2 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/MD_STS.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS=(TGraph*)fin2->Get("Threshold_Plot");
                NU_STS->SetLineColor(2);
                NU_STS->SetLineStyle(1);

        TFile *fin3 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/MD_STS_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
                NU_STS_Bent=(TGraph*)fin3->Get("Threshold_Plot");
                NU_STS_Bent->SetLineColor(2);
                NU_STS_Bent->SetLineStyle(5);
        TFile *fin4 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/MD_STS_Earth.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
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
            
        TFile *fin5 = TFile::Open(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/%s/%sGeV/MD_STS_Earth_Bent.root",Exp_Flux[Exp].c_str(),Mass_Point[Mass_INT].c_str()));
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

             
            cout << "Difference: " << (Slope_Earth_Aft-Slope_Earth_Bef)/Slope_Earth_Bef << endl;
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

            
            MD->Draw("ALP");
             
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
         */
    

