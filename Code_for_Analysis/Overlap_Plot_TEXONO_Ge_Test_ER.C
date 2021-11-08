#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
//#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "velocity_distribution_2000_Ave_ER.h"

#include "cpkkd_calculation_New.h"
void Overlap_Plot_TEXONO_Ge_Test_ER()
{

    const int Number=29;
    double Sigma_SI_Array[Number];
    //Method1 Threshold: 200eV
    double Sigma_SI_With_Threshold_M1[Number];double Sigma_SI_With_Threshold_earth_M1[Number];
    double Ratio_for_Average_M1[Number];double Ratio_for_PREM_M1[Number];

    //Method2 Threshold: 200km/s
    double Sigma_SI_With_Threshold_M2[Number];double Sigma_SI_With_Threshold_earth_M2[Number];
    double Ratio_for_Average_M2[Number];double Ratio_for_PREM_M2[Number];

    //Method1 Threshold: 200eV-300eV
    double Sigma_SI_With_Threshold_M3[Number];double Sigma_SI_With_Threshold_earth_M3[Number];
    double Ratio_for_Average_M3[Number];double Ratio_for_PREM_M3[Number];
    double Sigma_SI_With_Threshold_M3_Error[Number];double Error_X[Number];

    double CPKKD_EXCLUSION[Number];
    double Mass=1;//2.34
    int Take_Plot=0;
    string Type_of_Model="ER"; int Type_of_Model_INT=4;
    cout << "max_recoil_A_EM_keV(): " << max_recoil_A_EM_keV(2.34, 779.135*1000.0/2.99792458e8, AGe) << endl;
    
    /*
    Sigma_SI_Array[0]=1*(3e-30);
    Sigma_SI_With_Threshold_M1[0]=1e-42;
    Sigma_SI_With_Threshold_M3[0]=1e-42;
     */
    
    int Point_Number=0;
     
for(int kkk=36;kkk<37;kkk++)//for(int kkk=31;kkk<41;kkk++)
{
    for(int lll=9; lll<10; lll++)//for(int lll=1; lll<10; lll++)
    {
    float Multiply_Number=0;
        
    char fname[100];string Sfname;
    int YYY=0; int Addition_Option=0;

    //Sigma_SI_Array[Point_Number]=Multiply_Number*TMath::Power(10,-YYY);
    Sigma_SI_Array[Point_Number] = 1;
    TH1F *Flux_HIST_Random; TH1F *Flux_HIST_Aft_Collision_Earth; TH1F *Flux_HIST_Aft_Collision_EARTH;
    cout << "=======Right======== " << endl;
    cout <<"Multiply_Number:" << Multiply_Number << endl;
    cout << "=======Right======== " << endl;

    TFile *fin = TFile::Open("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/1_CDEX_Flux/1GeV/2.root");
    Flux_HIST_Random=(TH1F*)fin->Get("Flux_HIST_Random");
    Flux_HIST_Aft_Collision_EARTH=(TH1F*)fin->Get("Flux_HIST_Aft_Collision_EARTH");

    /*
    double *RE_DATA_Aft    =Hist_SetLimit_Plot_v2_Extract_Peak(0);
    double *RE_Rate_Aft    =Hist_SetLimit_Plot_v2_Extract_Peak(1);
    double *RE_DATA_Err_Aft=Hist_SetLimit_Plot_v2_Extract_Peak(2);
    double *RE_Rate_Err_Aft=Hist_SetLimit_Plot_v2_Extract_Peak(3);
    double FACTOR = cpkkd_calculation_New(Flux_HIST_Aft_Collision_EARTH,Mass,Sigma_SI_Array[Point_Number],RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft);
    
    cout << "YES!: " << FACTOR << endl;
    cout << "SIGMASI*FACTOR: " << Sigma_SI_Array[Point_Number]*FACTOR << endl;
        
    CPKKD_EXCLUSION[Point_Number]=Sigma_SI_Array[Point_Number]*FACTOR;
     */
    double T_QF_Original_Bef_Array[reso_T]; double T_QF_Original_Aft_Array[reso_T];
    double Factor1_Original_Bef_Array[reso_T]; double Factor1_Original_Aft_Array[reso_T];

    double *T_QF_Original_Bef=RecoilX_Event(0,Flux_HIST_Random,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        T_QF_Original_Bef_Array[i]=T_QF_Original_Bef[i];}
        cout << "======================================" << endl;
    double *T_QF_Original_Aft=RecoilX_Event(0,Flux_HIST_Aft_Collision_EARTH,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        T_QF_Original_Aft_Array[i]=T_QF_Original_Aft[i];}
        cout << "======================================" << endl;

    double *Factor1_Original_Bef=RecoilX_Event(1,Flux_HIST_Random,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        Factor1_Original_Bef_Array[i]=(Factor1_Original_Bef[i]);}
        cout << "======================================" << endl;

    double *Factor1_Original_Aft=RecoilX_Event(1,Flux_HIST_Aft_Collision_EARTH,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,1);
    for(int i=0;i<reso_T;i++){
        Factor1_Original_Aft_Array[i]=(Factor1_Original_Aft[i]);}
        cout << "Sigma_SI_Array[Point_Number]: " << Sigma_SI_Array[Point_Number] << endl;
        cout << "======================================" << endl;

        /*
        for(int i=0;i<reso_T;i++){
            cout << "T_QF_Original_Bef_Array[i]: " << T_QF_Original_Bef_Array[i] << endl;
            cout << "Factor1_Original_Bef_Array[i]: " << Factor1_Original_Bef_Array[i] << endl;
        }
        for(int i=0;i<reso_T;i++){
            cout << "T_QF_Original_Aft_Array[i]: " << T_QF_Original_Aft_Array[i] << endl;
            cout << "Factor1_Original_Aft_Array[i]: " << Factor1_Original_Aft_Array[i] << endl;
        }
         */
    double RecoilX_Event_Original_M1=0; double RecoilX_Event_Aft_EARTH_M1=0;
    double RecoilX_Event_Original_M3=0; double RecoilX_Event_Aft_EARTH_M3=0;
        int Number_of_Bin=0;
    for(int i=0;i<reso_T;i++)
    {
        //cout << "Factor1_Original_Bef_Array[i]: " << Factor1_Original_Bef_Array[i] << endl;
        //cout << "Factor1_Original_Aft_Array[i]: " << Factor1_Original_Aft_Array[i] << endl;
        if(T_QF_Original_Bef_Array[i]>0.7 and T_QF_Original_Bef_Array[i]<1.7)
        {
            Number_of_Bin = Number_of_Bin + 1;
            RecoilX_Event_Original_M1 = RecoilX_Event_Original_M1 + Factor1_Original_Bef_Array[i];
        }
        //cout << "===========================================================================";
        if(T_QF_Original_Aft_Array[i]>0.7 and T_QF_Original_Bef_Array[i]<0.75)
        {
            RecoilX_Event_Aft_EARTH_M1 = RecoilX_Event_Aft_EARTH_M1 + Factor1_Original_Aft_Array[i];
        }
        if(T_QF_Original_Bef_Array[i]>0.7 and T_QF_Original_Bef_Array[i]<0.75) RecoilX_Event_Original_M3 = RecoilX_Event_Original_M3 + Factor1_Original_Bef_Array[i];
        if(T_QF_Original_Aft_Array[i]>0.7 and T_QF_Original_Aft_Array[i]<0.75) RecoilX_Event_Aft_EARTH_M3 = RecoilX_Event_Aft_EARTH_M3 + Factor1_Original_Aft_Array[i];
    }
        cout << "RecoilX_Event_Original_M1: " << RecoilX_Event_Original_M1 << endl;
        cout << "Number_of_Bin: " << Number_of_Bin << endl;
        cout << "RecoilX_Event_Original_M1/Number_of_Bin: " << RecoilX_Event_Original_M1/(double)Number_of_Bin << endl;
        cout << "60/(RecoilX_Event_Original_M1/(double)Number_of_Bin): " << 0.5/(RecoilX_Event_Original_M1/(double)Number_of_Bin) << endl;
        double Scaling = 0.5/(RecoilX_Event_Original_M1/(double)Number_of_Bin);
        cout << "DS_Try(1e-9,0.5): " << CS_Try(1,0.5) << endl;
        cout << "Final: " << CS_Try(1*sqrt(Scaling),0.5) << endl;
        cout << "CS_Try: " << CS_Try(5.28*1e-4,0.5) << endl;
    double EARTH_Original=0;
    double EARTH_Bigger_Than_Threshold=0;
    double Earth_Bigger_Than_Threshold=0;
    double Event_Number = Flux_HIST_Aft_Collision_EARTH->GetEntries();


    //Method2 Threhold: 200eV->1.009keV->201km/s
    for(int jjj=0; jjj<2000; jjj++)
    {//Double_t xcenter = h3->GetZaxis()->GetBinCenter(27);
        double EARTH_Original_Bin=0; double BINX_EARTH=0;double BINX_Earth=0;
        EARTH_Original_Bin=Flux_HIST_Random->GetXaxis()->GetBinCenter(jjj);
        BINX_EARTH=Flux_HIST_Aft_Collision_EARTH->GetXaxis()->GetBinCenter(jjj);
        //cout << "EARTH_Original_Bin: " << EARTH_Original_Bin << "Flux_HIST_Random->GetBinContent(jjj): " << Flux_HIST_Random->GetBinContent(jjj) << endl;
        //cout << "BINX_EARTH: " << BINX_EARTH << "(Flux_HIST_Aft_Collision_EARTH->GetBinContent(jjj)): " << (Flux_HIST_Aft_Collision_EARTH->GetBinContent(jjj)) << endl;

        if(EARTH_Original_Bin>278)
        {
            EARTH_Original = EARTH_Original + (Flux_HIST_Random->GetBinContent(jjj));
        }
        if(BINX_EARTH>278)
        {
            EARTH_Bigger_Than_Threshold  = EARTH_Bigger_Than_Threshold  + (Flux_HIST_Aft_Collision_EARTH->GetBinContent(jjj));
        }
    }
    Ratio_for_PREM_M2[Point_Number]   =(EARTH_Bigger_Than_Threshold/EARTH_Original);

    Sigma_SI_With_Threshold_M2[Point_Number]      =Sigma_SI_Array[Point_Number]*(EARTH_Bigger_Than_Threshold/EARTH_Original);

        cout << "EARTH_Bigger_Than_Threshold: " << EARTH_Bigger_Than_Threshold << endl;
        cout << "EARTH_Original: " << EARTH_Original << endl;

    cout << "EARTH_Bigger_Than_Threshold/EARTH_Original: " << EARTH_Bigger_Than_Threshold/EARTH_Original << endl;

    Sigma_SI_With_Threshold_M1[Point_Number]      =Sigma_SI_Array[Point_Number]*(RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1);
        
    cout << "RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1: " << RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1 << endl;
    cout << "Sigma_SI_With_Threshold_M1[Point_Number]: " << Sigma_SI_With_Threshold_M1[Point_Number] << endl;
    
    Sigma_SI_With_Threshold_M3[Point_Number]      =Sigma_SI_Array[Point_Number]*(RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3);
    
    cout << "RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3: " << RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3<< endl;
    cout << "Sigma_SI_With_Threshold_M3[Point_Number]: " << Sigma_SI_With_Threshold_M3[Point_Number]  << endl;
 
   Sigma_SI_With_Threshold_M3_Error[Point_Number] =Sigma_SI_Array[Point_Number]*Binomial_Error(Event_Number, (RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3) );
    cout << " Sigma_SI_With_Threshold_M3_Error[Point_Number]: " <<  Sigma_SI_With_Threshold_M3_Error[Point_Number]  << endl;
    Error_X[Point_Number]=0;
        
    TCanvas *c3 = new TCanvas("c3");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    double SetOutline_X[1]={0};double SetOutline_Y[1]={0};
    TGraph *SetOutLine = new TGraph(1,SetOutline_X,SetOutline_Y);
    SetOutLine->GetXaxis()->SetRangeUser(0,3);
    SetOutLine->GetYaxis()->SetRangeUser(0,1e+15);

    const int data_Bin_ER = reso_T;
        //const int data_Bin_ER = data_bin;

    //Energy recoil Spectrum
    TGraph *ER_Spectrum_Bef = new TGraph(data_Bin_ER,T_QF_Original_Bef_Array,Factor1_Original_Bef_Array);
    //TGraph *ER_Spectrum_Bef = new TGraph(data_bin,T_QF_Original_Bef_Array,Factor1_Original_Bef_Array);
    ER_Spectrum_Bef->GetXaxis()->SetTitle("Energy[keV]");
    //ER_Spectrum_Bef->GetYaxis()->SetTitle("Count");
    ER_Spectrum_Bef->GetYaxis()->SetTitle("d<#sigma>/dT");

    ER_Spectrum_Bef->GetXaxis()->SetLimits(1e-2,3e+0);
    //ER_Spectrum_Bef->GetXaxis()->SetLimits(1e-2,1e1);
    ER_Spectrum_Bef->GetYaxis()->SetRangeUser(1e-21,1e-13);
    //ER_Spectrum_Bef->GetYaxis()->SetRangeUser(1e-7,1e+7);
     
    //ER_Spectrum_Bef->GetXaxis()->SetLimits(0,1);
    //ER_Spectrum_Bef->GetYaxis()->SetRangeUser(1,12);

    ER_Spectrum_Bef->SetLineColor(2);
    ER_Spectrum_Bef->SetLineWidth(5);
    TGraph *ER_Spectrum_Aft = new TGraph(data_bin,T_QF_Original_Aft_Array,Factor1_Original_Aft_Array);
    ER_Spectrum_Aft->GetXaxis()->SetTitle("Energy[keV]");
    ER_Spectrum_Aft->GetYaxis()->SetTitle("Count");
    ER_Spectrum_Aft->GetXaxis()->SetLimits(1e-2,1e+2);
    ER_Spectrum_Aft->GetYaxis()->SetRangeUser(1e-18,1e-13);
    ER_Spectrum_Aft->SetLineColor(3);


 /*
    double *Data_X       =Hist_SetLimit_Plot_v2_Extract_Peak(0);
    double *Data_Y       =Hist_SetLimit_Plot_v2_Extract_Peak(1);
    double *Data_X_Error =Hist_SetLimit_Plot_v2_Extract_Peak(2);
    double *Data_Y_Error =Hist_SetLimit_Plot_v2_Extract_Peak(3);
    double Data_X_Point[257]; double Data_Y_Point[257]; double Data_X_Error_Point[257]; double Data_Y_Error_Point[257];
        for(int i=0; i<257; i++)
        {
        Data_X_Point[i]=Data_X[i]; Data_Y_Point[i]=(Data_Y[i]); Data_X_Error_Point[i]=Data_X_Error[i]; Data_Y_Error_Point[i]=(Data_Y_Error[i]);
            
        }
    TGraphErrors *Data_Point = new TGraphErrors(257,Data_X_Point,Data_Y_Point,Data_X_Error_Point,Data_Y_Error_Point);
 */

        double Number[2]={0,0};int Scale[2]={0,0};
        for(int kkk=0; kkk<15; kkk++)
        {
            if(Number[0]!=0)break;
            else if(TMath::Power(10,kkk)*(RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1)<1)continue;
            else
            {
                Number[0]=TMath::Power(10,kkk)*(RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1);
                Scale[0]=kkk;
            }
        }
        for(int kkk=0; kkk<15; kkk++)
        {
            if(Number[1]!=0)break;
            else if(TMath::Power(10,kkk)*(RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3)<1)continue;
            else
            {
                Number[1]=TMath::Power(10,kkk)*(RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3);
                Scale[1]=kkk;
            }
        }
    TLegend *leg = new TLegend(0.3,0.6,0.5,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
        /*
    leg->AddEntry("",Form("%s:#sigma_{SI}:%.2fe-%i, M_{#chi}=%.3fGeV",Type_of_Model.c_str(),Multiply_Number,YYY,Mass),"");
    leg->AddEntry(ER_Spectrum_Bef,"ER_Spectrum_Bef","l");
    leg->AddEntry(ER_Spectrum_Aft,"ER_Spectrum_Aft","l");
    leg->AddEntry("",Form("Ratio of >200eV: %.3f #times 10^{-%i}",Number[0],Scale[0]),"l");
    leg->AddEntry("",Form("Ratio of >200eV and <250eV: %.3f #times 10^{-%i}",Number[1],Scale[1]),"l");
        */
        
    ER_Spectrum_Bef->Draw("ALP");
    //ER_Spectrum_Aft->Draw("ALP");
    //Data_Point->Draw("same");
    ER_Spectrum_Aft->Draw("LPsame");
    leg->Draw();
    c3->SetLogy();
    c3->SetLogx();
    if(Take_Plot==1)c3->Print(Form("%s_%s_Aft.pdf",Type_of_Model.c_str(),Sfname.c_str()));
 

     //Velocity Distributions
     /*
        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

    int ReBin_Number=20;
    //Flux_HIST_Random->Rebin(ReBin_Number);
    //Flux_HIST_Aft_Collision_Earth->Rebin(ReBin_Number);
    //Flux_HIST_Aft_Collision_EARTH->Rebin(ReBin_Number);

    Flux_HIST_Random->Draw("hist");
    //Flux_HIST_Aft_Collision_Earth->Draw("histsame");
    Flux_HIST_Aft_Collision_EARTH->Draw("histsame");
    
    Flux_HIST_Random->GetXaxis()->SetTitle("Velocity[km/s]");
    Flux_HIST_Random->GetYaxis()->SetTitle("Possibility");
    Flux_HIST_Random->GetYaxis()->SetRangeUser(0,1);
    Flux_HIST_Random->GetXaxis()->SetRangeUser(0,14);

    //clbr_max_he->GetXaxis()->SetRangeUser(0,0.5*1e7);
    
    TLegend *leg = new TLegend(0.3,0.6,0.5,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->AddEntry("",Form("#sigma_{SI}:%.2fe-%i, M_{#chi}=%iGeV",Multiply_Number,YYY,Mass),"");
    leg->AddEntry(Flux_HIST_Random,"WIMP Flux","l");
    //leg->AddEntry(Flux_HIST_Aft_Collision_Earth,"Through Average Density","l");
    leg->AddEntry(Flux_HIST_Aft_Collision_EARTH,"Through PREM Model","l");
    //leg->AddEntry("",Form("Ratio of [Velocity>201km/s]: %.5f",EARTH_Bigger_Than_Threshold/EARTH_Original),"l");
    leg->Draw();

     //   cout <<"Multiply_Number:" << lll << endl;
     //   cout <<"Sigma_SI:" << kkk << endl;

    c3->Print(Form("Sigma_SI_Flux/%s.pdf",Sfname.c_str()));
    */
        
    Point_Number = Point_Number+1;
    cout << "Point_Number: " << Point_Number << endl;
    }
}
/*
    TGraph *Threshold_Plot_2 = new TGraph(Point_Number,Sigma_SI_Array,CPKKD_EXCLUSION);
    Threshold_Plot_2->SetLineColor(5);
    Threshold_Plot_2->SetMarkerColor(5);
    Threshold_Plot_2->SetName("Threshold_Plot");
    Threshold_Plot_2->SetLineColor(2);
    Threshold_Plot_2->SetMarkerColor(2);
    Threshold_Plot_2->GetXaxis()->SetLimits(1e-40,1e-28);
    Threshold_Plot_2->GetYaxis()->SetRangeUser(1e-40,1e-29);
    Threshold_Plot_2->SetTitle("");
    Threshold_Plot_2->GetXaxis()->SetTitle("Sigma_SI");
    Threshold_Plot_2->GetYaxis()->SetTitle("(Sigma_SI) * FACTOR");
    Threshold_Plot_2->Draw("AL");

    TF1 *Linear_Line = new TF1("Linear_Line","x",1e-40,1e-28);
    Linear_Line->SetLineColor(4);
    Linear_Line->Draw("Lsame");
    
    gPad->SetLogx();
    gPad->SetLogy();
    cout << "Yes?!3" << endl;
    c4->Print("EXCLUSION_Plot.pdf");*/



/*
    Threshold_Plot_0->SetName("Threshold_Plot");
    Threshold_Plot_0->SetLineColor(2);
    Threshold_Plot_0->SetMarkerColor(2);
    Threshold_Plot_0->GetXaxis()->SetLimits(1e-40,1e-28);
    Threshold_Plot_0->GetYaxis()->SetRangeUser(1e-40,1e-25);
    Threshold_Plot_0->SetTitle("");
    Threshold_Plot_0->GetXaxis()->SetTitle("Sigma_SI");
    Threshold_Plot_0->GetYaxis()->SetTitle("(Sigma_SI) #times (Event_Fraction)");
*/
    
    TCanvas *c4 = new TCanvas("c4");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    cout << "Yes?!1" << endl;
    
    for(int kkk=0; kkk<Point_Number;kkk++)
    {
        cout << "Sigma_SI_Array: " << Sigma_SI_Array[kkk] << endl;
        cout << "Sigma_SI_With_Threshold_M1[Point_Number] : " << Sigma_SI_With_Threshold_M1[kkk]  << endl;
        cout << "Sigma_SI_With_Threshold_M3[Point_Number] : " << Sigma_SI_With_Threshold_M3[kkk]  << endl;
    }
    cout << "Yes?!2" << endl;

    TGraph *Threshold_Plot_2 = new TGraph(Point_Number,Sigma_SI_Array,Sigma_SI_With_Threshold_M1);
    Threshold_Plot_2->SetLineColor(5);
    Threshold_Plot_2->SetMarkerColor(5);
    Threshold_Plot_2->SetName("Threshold_Plot");
    Threshold_Plot_2->SetLineColor(2);
    Threshold_Plot_2->SetMarkerColor(2);
    Threshold_Plot_2->GetXaxis()->SetLimits(1e-42,1e-27);
    Threshold_Plot_2->GetYaxis()->SetRangeUser(1e-42,1e-29);
    Threshold_Plot_2->SetTitle("");
    Threshold_Plot_2->GetXaxis()->SetTitle("Sigma_SI");
    Threshold_Plot_2->GetYaxis()->SetTitle("(Sigma_SI) #times (Event_Fraction)");

    
    TGraph *Threshold_Plot_4 = new TGraph(Point_Number,Sigma_SI_Array,Sigma_SI_With_Threshold_M3);
    Threshold_Plot_2->SetLineColor(7);
    Threshold_Plot_2->SetMarkerColor(7);

    /*
    TGraphErrors *Threshold_Plot_4 = new TGraphErrors(Point_Number,Sigma_SI_Array,Sigma_SI_With_Threshold_M3,Error_X,Sigma_SI_With_Threshold_M3_Error);
    Threshold_Plot_2->SetLineColor(7);
    Threshold_Plot_2->SetMarkerColor(7);
    */

    TF1 *Linear_Line = new TF1("Linear_Line","x",1e-42,1e-27);
    Linear_Line->SetLineColor(4);
    
    //Threshold_Plot_0->Draw("A");
    //Threshold_Plot_1->Draw("Lsame");
    Threshold_Plot_2->Draw("AL");
    //Threshold_Plot_3->Draw("Lsame");
    //Threshold_Plot_4->Draw("Lsame");
    //Threshold_Plot_4->Draw("Lsame");
    //Threshold_Plot_5->Draw("Lsame");
    Linear_Line->Draw("Lsame");
    
    TLegend *leg = new TLegend(0.5,0.2,0.8,0.5);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->AddEntry("",Form("M_{#chi}=%.2fGeV",Mass),"");
    //leg->AddEntry(Threshold_Plot_0,"Threshold:201km/s","lP");
   //leg->AddEntry(Threshold_Plot_1,"Through Average Density","lP");
    leg->AddEntry(Threshold_Plot_2,"Count> 200eV","lP");
   //leg->AddEntry(Threshold_Plot_3,"Through Average Density","lP");
    leg->AddEntry(Threshold_Plot_4,"Count: 200-300eV","lP");
   // leg->AddEntry(Threshold_Plot_5,"Through Average Density","lP");
    leg->AddEntry(Linear_Line ,"X=Y","l");
    leg->Draw();

    
    gPad->SetLogx();
    gPad->SetLogy();
    cout << "Yes?!3" << endl;
    cout << "Point_Number: " << Point_Number << endl;

    //c4->Print("Recoil_Spectrum_CDEX/2P5GeV_MD_EXCLUSION_Plot_1_CDEX.pdf");
    if(Take_Plot==1)
    {
        char fout_name[100];
        sprintf(fout_name,Form("Recoil_Spectrum_CDEX_0P9GeV/%s_EXCLUSION_Plot_1_CDEX.root",Type_of_Model.c_str()));
        TFile *fout=new TFile(fout_name,"recreate");
        Threshold_Plot_2->Write();
    }
}
    /*
    TGraph *Threshold_Plot_Ratio = new TGraph(11,Sigma_SI_Array,Ratio_for_Average);
    Threshold_Plot_Ratio->SetName("Threshold_Plot_Ratio");
    Threshold_Plot_Ratio->SetLineColor(2);
    Threshold_Plot_Ratio->GetXaxis()->SetRangeUser(1e-42,1e-28);
    Threshold_Plot_Ratio->GetYaxis()->SetRangeUser(0,1);
    Threshold_Plot_Ratio->SetTitle("Threshold_Plot_Ratio");
    Threshold_Plot_Ratio->GetXaxis()->SetTitle("Sigma_SI");
    Threshold_Plot_Ratio->GetYaxis()->SetTitle("Ratio");
    Threshold_Plot_Ratio->SetMarkerStyle(20);
    Threshold_Plot_Ratio->SetMarkerColor(2);
    Threshold_Plot_Ratio->SetMarkerSize(0.5);

    TGraph *Threshold_Plot_Ratio_1 = new TGraph(11,Sigma_SI_Array,Ratio_for_PREM);
    Threshold_Plot_Ratio_1->SetMarkerStyle(20);
    Threshold_Plot_Ratio_1->SetMarkerColor(3);
    Threshold_Plot_Ratio_1->SetMarkerSize(0.5);

    Threshold_Plot_Ratio->Draw("AP");
    Threshold_Plot_Ratio_1->Draw("Psame");
    */
    

