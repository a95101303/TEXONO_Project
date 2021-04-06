#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"

void PRL_Check_For_Heavy_DM_Event_Number()//CDEX:Threshold=160eV, TEXONO:Threshold=200eV
{
    //FLux for the different experiment!
    int Exp_INDEX=0;//0 is for the MAJORANA experiment, 1 is for DAMA
    double Area[2]={2.5*2.5*3.14,10.2*10.2};
    double  WIMP_Mass_Array[14]={1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18};//WIMP
    double Flux_1[14];//Count/(day*m^2) for XENON1T and Count/(day*cm^2) for CRESST
    for(int kkk=0; kkk<14; kkk++)
    {
        //Flux_1[kkk]=365*Total_Flux(WIMP_Mass_Array[kkk])*1e4*Area[0];//Count/(day*m^2)
        Flux_1[kkk]=365*Total_Flux(WIMP_Mass_Array[kkk])*Area[Exp_INDEX];//Count/(day*cm^2)
    }
    TCanvas *c3 = new TCanvas("c3");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TGraph *Flux = new TGraph(14,WIMP_Mass_Array,Flux_1);
    Flux->GetXaxis()->SetTitle("m_{#chi}(GeV)");
    Flux->GetYaxis()->SetTitle("Evt/year*m^{2}");
    Flux->GetXaxis()->SetLimits(1e4,1e19);
    Flux->GetYaxis()->SetRangeUser(1,1e+15);
    Flux->SetLineColor(2);
    Flux->SetLineWidth(5);

    Flux->Draw("ALP");
    c3->SetLogy();
    c3->SetLogx();

    c3->Print("PRL_HEAVY_DM_Check/Flux_MAJORANA.pdf");
     
    //Puzzle of the lower bound of the cross-section for MIMP(Solved)
    
    /*
    int Exp_INDEX=0;//0 for XENON1T and 1 for TEXONO
    string Exp_Type[2]={"Xe","Ge"};
    double  WIMP_Mass_Array[9]={1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18};//WIMP
    int WIMP_Mass_Array_Int[9]={10,11,12,13,14,15,16,17,18};
    //double WIMP_Mass_Array[11]={1e10,1e11,1e12,1e13,1e14,1e15,1e16,1e17,1e18,1e19,1e20};//MIMP
    string Type_of_Model[4]={"NU","MD","BR","MDMPA"};
    double Sigma_SI_Array[9];
    
    //double Sigma_SI[11];//MIMP
for(int kkk=0; kkk<9;kkk++)
    {
    //double CPKKD_EXCLUSION[Number];
        //=======
        double Sigma_SI=7e-38 * pow(10,kkk);
    Sigma_SI_Array[kkk]=Sigma_SI;
        //=======
    double Event_Number=0;
    double Mass = WIMP_Mass_Array[kkk];//GeV
    cout << "Mass: " << Mass << endl;
    int Type_of_Model_INT=0;
    //=======================Recoil Spectrum set for three processes==============================
    double T_QF_Original_Bef_Array[reso_T];double Factor1_Original_Bef_Array[reso_T];
        
    double *T_QF_Original_Bef=RecoilX_Event(0,0,Mass,Sigma_SI,Type_of_Model_INT,2);
    for(int i=0;i<reso_T;i++){T_QF_Original_Bef_Array[i]=T_QF_Original_Bef[i];}
    double *Factor1_Original_Bef=RecoilX_Event(1,0,Mass,Sigma_SI,Type_of_Model_INT,2);
    for(int i=0;i<reso_T;i++){Factor1_Original_Bef_Array[i]=(Factor1_Original_Bef[i]);}
    
    for(int i=0;i<reso_T;i++)
    {
        if(Exp_INDEX==0 and T_QF_Original_Bef_Array[i]>1)
        {
            Event_Number = Event_Number + (906./2000.)*3.2*Factor1_Original_Bef_Array[i];
        }
        if(Exp_INDEX==1 and T_QF_Original_Bef_Array[i]>0.2)
        {
            Event_Number = Event_Number + (906./2000.)*1.06*Factor1_Original_Bef_Array[i];
        }
    }
    TCanvas *c3 = new TCanvas("c3");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TGraph *ER_Spectrum_Bef = new TGraph(reso_T,T_QF_Original_Bef_Array,Factor1_Original_Bef_Array);
    ER_Spectrum_Bef->GetXaxis()->SetTitle("Energy[keV]");
    ER_Spectrum_Bef->GetYaxis()->SetTitle("Evt/year*kg*keV");
    ER_Spectrum_Bef->GetXaxis()->SetLimits(1e-2,1e3);
    ER_Spectrum_Bef->GetYaxis()->SetRangeUser(1e-8,1e+15);
    ER_Spectrum_Bef->SetLineColor(2);
    ER_Spectrum_Bef->SetLineWidth(5);
        
    TLegend *leg = new TLegend(0.1,0.6,0.3,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    if(Exp_INDEX==0)leg->AddEntry("","Weight: 3.2kg",""); if(Exp_INDEX==1)leg->AddEntry("","Weight: 1.06kg","");
    leg->AddEntry("",Form("m_{#chi}=10^{%i}GeV",WIMP_Mass_Array_Int[kkk]),"");
    leg->AddEntry(ER_Spectrum_Bef,Form("%.3f #times 10^{-%.f}",CNFNV(0,Sigma_SI),CNFNV(1,Sigma_SI)),"l");
    if(Exp_INDEX==0)leg->AddEntry("",Form("Total Count/year*3.2 kg > 1keV for XENON1T: %.3f",Event_Number),"");
    if(Exp_INDEX==1)leg->AddEntry("",Form("Total Count/year*1.06kg > 0.2keV for TEXONO: %.3f",Event_Number),"");
    ER_Spectrum_Bef->Draw("ALP");

    if(Exp_INDEX==0){
    TLine *Threshold_of_XENON1T= new TLine(1,0,1,1e20);
    Threshold_of_XENON1T->SetLineColor(45);
    Threshold_of_XENON1T->Draw("Lsame");}

    if(Exp_INDEX==1){
    TLine *Threshold_of_XENON1T= new TLine(0.2,0,0.2,1e20);
    Threshold_of_XENON1T->SetLineColor(45);
    Threshold_of_XENON1T->Draw("Lsame");}

    leg->Draw();

    c3->SetLogy();
    c3->SetLogx();

    cout << "Event_Number: " << Event_Number << endl;
        c3->Print(Form("PRL_HEAVY_DM_Check/%s_%i.pdf",Exp_Type[Exp_INDEX].c_str(),kkk));

    }
     
    for(int kkk=0; kkk<9;kkk++)
    {
        cout << WIMP_Mass_Array[kkk] << "," << Sigma_SI_Array[kkk] << "," << endl;
    }
    */
}

