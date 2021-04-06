#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"

void Overlap_Plot_CDEX_Ge()//CDEX:Threshold=160eV, TEXONO:Threshold=200eV
{
    //double CPKKD_EXCLUSION[Number];
    const int Number=29;int Take_Plot=1;//Plot or not
    //==============Exp and Model=============//
    int Experiment_Type=0;double Threshold[2]={0.175,0.2};//0 for CDEX and 1 for TEXONO, (keV)
    int Type_of_Model_INT=2;string Type_of_Model[4]={"NU","MD","BR","MDMPA"};
    //==============Mass==============//
    string Mass_Point[17]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10"};
    //================================//
    double Sigma_SI_Array[Number];
    double Sigma_SI_With_Threshold_M1[Number];//>200eV
    double Sigma_SI_With_Threshold_M3[Number];//200-300eV
    
    double Mass=0;
for(int Mass_INT=10; Mass_INT<16; Mass_INT++){//Open1
    int Check_ZERO=0;
    int Point_Number=0;//Initialize
    cout << "Mass_INT: " << Mass_INT << endl;
    for(int FILE=1; FILE<40; FILE++){//Open1
        string path = Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/%sGeV/%i.root",Mass_Point[Mass_INT].c_str(),FILE);
        ifstream fin(path);
        if(fin.is_open() and Check_ZERO==0){//Open
        //===============Input the ROOTFILE for all==============
            cout << "path: " << path << endl;
        TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/%sGeV/%i.root",Mass_Point[Mass_INT].c_str(),FILE));
        TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
        Double_t mx,sigma_si;
        T1_TREE->SetBranchAddress("mx",&mx);T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
        T1_TREE->GetEntry(0);
        Mass=mx;Sigma_SI_Array[Point_Number]=sigma_si;
            cout << "Sigma_SI_Array[Point_Number]: " << Sigma_SI_Array[Point_Number] << endl;
        //==============Input STandard Flux as well as the "attenuated" flux==============
        TH1F *Flux_HIST_Random;TH1F *Flux_HIST_Aft_Collision_EARTH;
        Flux_HIST_Random=(TH1F*)ROOT_FILE->Get("Flux_HIST_Random");
        Flux_HIST_Aft_Collision_EARTH=(TH1F*)ROOT_FILE->Get("Flux_HIST_Aft_Collision_EARTH");
        //=======================Recoil Spectrum set for three processes==============================
        double T_QF_Original_Bef_Array[reso_T];    double T_QF_Original_Aft_Array[reso_T];
        double Factor1_Original_Bef_Array[reso_T]; double Factor1_Original_Aft_Array[reso_T];

          
        double *T_QF_Original_Bef=RecoilX_Event(0,Flux_HIST_Random,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,0);
        for(int i=0;i<reso_T;i++){T_QF_Original_Bef_Array[i]=T_QF_Original_Bef[i];}
           
        double *T_QF_Original_Aft=RecoilX_Event(0,Flux_HIST_Aft_Collision_EARTH,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,0);
        for(int i=0;i<reso_T;i++){T_QF_Original_Aft_Array[i]=T_QF_Original_Aft[i];}
            
        double *Factor1_Original_Bef=RecoilX_Event(1,Flux_HIST_Random,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,0);
        for(int i=0;i<reso_T;i++){Factor1_Original_Bef_Array[i]=(Factor1_Original_Bef[i]);}
             
        double *Factor1_Original_Aft=RecoilX_Event(1,Flux_HIST_Aft_Collision_EARTH,Mass,Sigma_SI_Array[Point_Number],Type_of_Model_INT,0);
        for(int i=0;i<reso_T;i++){Factor1_Original_Aft_Array[i]=(Factor1_Original_Aft[i]);}
        //===========================Calculate the Ratio passing through====================//
        double RecoilX_Event_Original_M1=0; double RecoilX_Event_Aft_EARTH_M1=0;
        double RecoilX_Event_Original_M3=0; double RecoilX_Event_Aft_EARTH_M3=0;
        
        for(int i=0;i<reso_T;i++)
        {
            if(T_QF_Original_Bef_Array[i]>Threshold[Experiment_Type])//>Threshold
            {RecoilX_Event_Original_M1  = RecoilX_Event_Original_M1 + Factor1_Original_Bef_Array[i];}
            if(T_QF_Original_Aft_Array[i]>Threshold[Experiment_Type])//>Threshold
            {RecoilX_Event_Aft_EARTH_M1 = RecoilX_Event_Aft_EARTH_M1 + Factor1_Original_Aft_Array[i];}
            if(T_QF_Original_Bef_Array[i]>Threshold[Experiment_Type] and T_QF_Original_Bef_Array[i]<Threshold[Experiment_Type]+0.05)//>Threshold and <Threshold+0.1
            {RecoilX_Event_Original_M3  = RecoilX_Event_Original_M3 + Factor1_Original_Bef_Array[i];}
            if(T_QF_Original_Aft_Array[i]>Threshold[Experiment_Type] and T_QF_Original_Aft_Array[i]<Threshold[Experiment_Type]+0.05)//>Threshold and <Threshold+0.1
            {RecoilX_Event_Aft_EARTH_M3 = RecoilX_Event_Aft_EARTH_M3 + Factor1_Original_Aft_Array[i];}
        }

        //Possion for >threshold
        Sigma_SI_With_Threshold_M1[Point_Number]      =Sigma_SI_Array[Point_Number]*(RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1);
        cout << "RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1: " << RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1 << endl;
        cout << "Sigma_SI_With_Threshold_M1[Point_Number]: " << Sigma_SI_With_Threshold_M1[Point_Number] << endl;
        //Possion for >threshold and <threshold+0.1
        Sigma_SI_With_Threshold_M3[Point_Number]      =Sigma_SI_Array[Point_Number]*(RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3);
        cout << "RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3: " << RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3<< endl;
        cout << "Sigma_SI_With_Threshold_M3[Point_Number]: " << Sigma_SI_With_Threshold_M3[Point_Number]  << endl;
     
        if(Sigma_SI_With_Threshold_M1[Point_Number]==0 or Sigma_SI_With_Threshold_M3[Point_Number]==0)
        {
            Check_ZERO=1;
        }
        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);

            double Data_RE_CDEX[4]= {0.208,0.30,0.35,0.4};//CDEX: 0.208
            double Data_RATE_CDEX[4]= {7.3,7.6,7.4,7.9};
            double Data_RE_Err[4]= {0,0,0,0};
            double Data_RE_Rate_Err[4]= {2.7,3.5,1.5,2};

            TGraphErrors *CDEXData = new TGraphErrors(4,Data_RE_CDEX,Data_RATE_CDEX,Data_RE_Err,Data_RE_Rate_Err);

        //Energy recoil Spectrum
        TGraph *ER_Spectrum_Bef = new TGraph(reso_T,T_QF_Original_Bef_Array,Factor1_Original_Bef_Array);
        ER_Spectrum_Bef->GetXaxis()->SetTitle("Energy[keVee]");
        ER_Spectrum_Bef->GetYaxis()->SetTitle("Count[Evts/keV/kg/day]");
        ER_Spectrum_Bef->GetXaxis()->SetLimits(1e-2,1e2);
        ER_Spectrum_Bef->GetYaxis()->SetRangeUser(1e-8,1e+15);

        ER_Spectrum_Bef->SetLineColor(2);
        ER_Spectrum_Bef->SetLineWidth(5);
        TGraph *ER_Spectrum_Aft = new TGraph(reso_T,T_QF_Original_Aft_Array,Factor1_Original_Aft_Array);
        ER_Spectrum_Aft->GetXaxis()->SetTitle("Energy[keVee]");
        ER_Spectrum_Aft->GetYaxis()->SetTitle("Count[Evts/keV/kg/day]");
        ER_Spectrum_Aft->SetLineColor(3);

        TLegend *leg = new TLegend(0.3,0.6,0.5,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
                    
            cout << "Mass: " << Mass << endl;
        leg->AddEntry("",Form("%s:#sigma_{SI}:%.2f #times 10^{-%.f}, M_{#chi}=%.2fGeV",Type_of_Model[Type_of_Model_INT].c_str(),CNFNV(0,Sigma_SI_Array[Point_Number]),CNFNV(1,Sigma_SI_Array[Point_Number]),mx),"");
        leg->AddEntry(ER_Spectrum_Bef,"ER_Spectrum_Bef","l");
        leg->AddEntry(ER_Spectrum_Aft,"ER_Spectrum_Aft","l");
        //leg->AddEntry("",Form("Ratio of >200eV: %.3f #times 10^{-%.f}",CNFNV(0,RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1),CNFNV(1,RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1)),"l");
        leg->AddEntry("",Form("Ratio of >0.16keV and <0.21eV: %.3f #times 10^{-%.f}",CNFNV(0,RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3),CNFNV(1,RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3)),"l");
            
        ER_Spectrum_Bef->Draw("ALP");
        ER_Spectrum_Aft->Draw("LPsame");
        CDEXData->Draw("LPsame");
        leg->Draw();
        c3->SetLogy();
        c3->SetLogx();
        if(Take_Plot==1)c3->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/%sGeV/Recoil_Spectrum/%s_%i.pdf",Mass_Point[Mass_INT].c_str(),Type_of_Model[Type_of_Model_INT].c_str(),FILE));
     
            
        Point_Number = Point_Number+1;
        cout << "Point_Number: " << Point_Number << endl;
        }//Close2
    }
        TCanvas *c4 = new TCanvas("c4");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        
        for(int kkk=0; kkk<Point_Number;kkk++)
        {
            cout << "Sigma_SI_Array: " << Sigma_SI_Array[kkk] << endl;
            cout << "Sigma_SI_With_Threshold_M1[Point_Number] : " << Sigma_SI_With_Threshold_M1[kkk]  << endl;
            cout << "Sigma_SI_With_Threshold_M3[Point_Number] : " << Sigma_SI_With_Threshold_M3[kkk]  << endl;
        }

        TGraph *Threshold_Plot_2 = new TGraph(Point_Number,Sigma_SI_Array,Sigma_SI_With_Threshold_M3);
        Threshold_Plot_2->SetLineColor(5);
        Threshold_Plot_2->SetMarkerColor(5);
        Threshold_Plot_2->SetName("Threshold_Plot_Bin_Possion");
        Threshold_Plot_2->SetLineColor(2);
        Threshold_Plot_2->SetMarkerColor(2);
        Threshold_Plot_2->GetXaxis()->SetLimits(1e-42,1e-27);
        Threshold_Plot_2->GetYaxis()->SetRangeUser(1e-42,1e-29);
        Threshold_Plot_2->SetTitle("");
        Threshold_Plot_2->GetXaxis()->SetTitle("Sigma_SI");
        Threshold_Plot_2->GetYaxis()->SetTitle("(Sigma_SI) #times (Event_Fraction)");

        TGraph *Threshold_Plot_4 = new TGraph(Point_Number,Sigma_SI_Array,Sigma_SI_With_Threshold_M1);
        Threshold_Plot_4->SetLineColor(7);
        Threshold_Plot_4->SetMarkerColor(7);

        TF1 *Linear_Line = new TF1("Linear_Line","x",1e-42,1e-27);
        Linear_Line->SetLineColor(4);
        
        Threshold_Plot_2->Draw("AL");
        //Threshold_Plot_4->Draw("Lsame");
    
        Linear_Line->Draw("Lsame");
        
        TLegend *leg = new TLegend(0.1,0.5,0.4,0.8);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry("",Form("M_{#chi}=%.2fGeV",Mass),"");
        leg->AddEntry(Threshold_Plot_2,"Count >0.16keV and <0.21keV","lP");
        //leg->AddEntry(Threshold_Plot_4,Form("Count: -300eV"),"lP");
        leg->AddEntry(Linear_Line ,"X=Y","l");
        leg->Draw();

        
        gPad->SetLogx();
        gPad->SetLogy();
        if(Take_Plot==1)
        {
            char fout_name[300];
            sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/%sGeV/%s_Line.root",Mass_Point[Mass_INT].c_str(),Type_of_Model[Type_of_Model_INT].c_str()));
            TFile *fout=new TFile(fout_name,"recreate");
            Threshold_Plot_2->Write();
            c4->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/1_CDEX_Flux/%sGeV/%s_Line.pdf",Mass_Point[Mass_INT].c_str(),Type_of_Model[Type_of_Model_INT].c_str()));
        }
    }
    
}
    

