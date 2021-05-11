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

void Overlap_Plot_TEXONO_Ge3_CRESST()//CDEX:Threshold=160eV, TEXONO:Threshold=200eV
{
    //double CPKKD_EXCLUSION[Number];
    const int Number=40;int Take_Plot=1;//Plot or not
    //int NU_Bins[12]={7,7,7,7,7,2,2,2,2,2,0,0};int MD_Bins[12]={35,35,35,35,35,7,7,7,2,2,2,2};int BR_Bins[12]={35,35,35,35,35,35,35,7,7,7,2,2};
    int NU_Bins[12]={0,0,0,0,0,0,0,0,0,0,0,0};int MD_Bins[12]={0,0,0,0,0,0,0,0,0,0,0,0};int BR_Bins[12]={0,0,0,0,0,0,0,0,0,0,0,0};
    //==============Exp and Model=============//
    int Experiment_Type=0;double Threshold[2]={0.2,0.3};//0 for CDEX and 1 for TEXONO, (keV)
    int Type_of_Model_INT=0;string Type_of_Model[4]={"NU","MD","BR","MDMPA"};
    //==============Mass==============//
for(int kkk=1;kkk<6;kkk++)
    {
        int Mass_INT=kkk;
        //string Mass_Point[19]={"2","1","0P9","0P8","0P7","0P6","0P5","0P4","0P3","0P2","0P1","0P09","0P08","0P07","0P06","0P05","10","5","7"};
        //string Mass_Point[10]={"0P2","0P19","0P18","0P17","0P16","0P15","0P14","0P13","0P12","0P11"};
        //string Mass_Point[12]={"20","19","17","15","13","11","9","7","5","4","3","2P35"};
        string Mass_Point[3]={"20","10","2","0P2"};
        //================================//
        double Sigma_SI_Array[Number];
        double Sigma_SI_With_Threshold_M1[Number];//>200eV
        double Sigma_SI_With_Threshold_M3[Number];//200-300eV
        int Point_Number=0;//Initialize
        int Check_ZERO=0;

        double Mass=0;
    for(int FILE=1; FILE<5; FILE++){//Open1
        string path = Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/%i_STS_Bent.root",Mass_Point[Mass_INT].c_str(),FILE);
        //cout << "path: " << path << endl;
        ifstream fin(path);
        if(fin.is_open() and Check_ZERO==0){//Open
        //===============Input the ROOTFILE for all==============
        TFile *ROOT_FILE = TFile::Open(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/%i_STS_Bent.root",Mass_Point[Mass_INT].c_str(),FILE));
        TTree *T1_TREE = (TTree*)ROOT_FILE->Get("t1");
        Double_t mx,sigma_si;
        T1_TREE->SetBranchAddress("mx",&mx);T1_TREE->SetBranchAddress("sigma_si",&sigma_si);
        T1_TREE->GetEntry(0);
        Mass=mx;Sigma_SI_Array[Point_Number]=sigma_si;
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
            
            if(T_QF_Original_Bef_Array[i]>0.020)//>Threshold
            {RecoilX_Event_Original_M1  = RecoilX_Event_Original_M1 + Factor1_Original_Bef_Array[i];}
            if(T_QF_Original_Aft_Array[i]>0.020)//>Threshold
            {RecoilX_Event_Aft_EARTH_M1 = RecoilX_Event_Aft_EARTH_M1 + Factor1_Original_Aft_Array[i];}
            
            if(Type_of_Model_INT==0)
            {
                if(T_QF_Original_Bef_Array[i]>0.020)
                {RecoilX_Event_Original_M3  = RecoilX_Event_Original_M3 + Factor1_Original_Bef_Array[i];}
                if(T_QF_Original_Aft_Array[i]>0.020)
                {RecoilX_Event_Aft_EARTH_M3 = RecoilX_Event_Aft_EARTH_M3 + Factor1_Original_Aft_Array[i];}
            }

            /*
            if(Type_of_Model_INT==0)
            {
                if(T_QF_Original_Bef_Array[i]>p103_le_VrV_ON_NaI1_50eV[NU_Bins[kkk]][0]-0.025 and T_QF_Original_Bef_Array[i]<p103_le_VrV_ON_NaI1_50eV[NU_Bins[kkk]][0]+0.025)
                {RecoilX_Event_Original_M3  = RecoilX_Event_Original_M3 + Factor1_Original_Bef_Array[i];}
                if(T_QF_Original_Aft_Array[i]>p103_le_VrV_ON_NaI1_50eV[NU_Bins[kkk]][0]-0.025 and T_QF_Original_Aft_Array[i]<p103_le_VrV_ON_NaI1_50eV[NU_Bins[kkk]][0]+0.025)
                {RecoilX_Event_Aft_EARTH_M3 = RecoilX_Event_Aft_EARTH_M3 + Factor1_Original_Aft_Array[i];}
            }
            if(Type_of_Model_INT==1)
            {
                if(T_QF_Original_Bef_Array[i]>p103_le_VrV_ON_NaI1_50eV[MD_Bins[kkk]][0]-0.025 and T_QF_Original_Bef_Array[i]<p103_le_VrV_ON_NaI1_50eV[MD_Bins[kkk]][0]+0.025)
                {RecoilX_Event_Original_M3  = RecoilX_Event_Original_M3 + Factor1_Original_Bef_Array[i];}
                if(T_QF_Original_Aft_Array[i]>p103_le_VrV_ON_NaI1_50eV[MD_Bins[kkk]][0]-0.025 and T_QF_Original_Aft_Array[i]<p103_le_VrV_ON_NaI1_50eV[MD_Bins[kkk]][0]+0.025)
                {RecoilX_Event_Aft_EARTH_M3 = RecoilX_Event_Aft_EARTH_M3 + Factor1_Original_Aft_Array[i];}
            }
            if(Type_of_Model_INT==2)
            {
                if(T_QF_Original_Bef_Array[i]>p103_le_VrV_ON_NaI1_50eV[BR_Bins[kkk]][0]-0.025 and T_QF_Original_Bef_Array[i]<p103_le_VrV_ON_NaI1_50eV[BR_Bins[kkk]][0]+0.025)
                {RecoilX_Event_Original_M3  = RecoilX_Event_Original_M3 + Factor1_Original_Bef_Array[i];}
                if(T_QF_Original_Aft_Array[i]>p103_le_VrV_ON_NaI1_50eV[BR_Bins[kkk]][0]-0.025 and T_QF_Original_Aft_Array[i]<p103_le_VrV_ON_NaI1_50eV[BR_Bins[kkk]][0]+0.025)
                {RecoilX_Event_Aft_EARTH_M3 = RecoilX_Event_Aft_EARTH_M3 + Factor1_Original_Aft_Array[i];}
            }
             */
        }

            cout << "RecoilX_Event_Aft_EARTH_M1: " << RecoilX_Event_Aft_EARTH_M1 << endl;
        //Possion for >threshold
        Sigma_SI_With_Threshold_M1[Point_Number]      =Sigma_SI_Array[Point_Number]*(RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1);
        cout << "RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1: " << RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1 << endl;
        cout << "Sigma_SI_With_Threshold_M1[Point_Number]: " << Sigma_SI_With_Threshold_M1[Point_Number] << endl;
        //Possion for >threshold and <threshold+0.1
            cout << "RecoilX_Event_Original_M3: " <<RecoilX_Event_Original_M3 << endl;
            cout << "RecoilX_Event_Aft_EARTH_M3: " << RecoilX_Event_Aft_EARTH_M3 << endl;

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

        //Energy recoil Spectrum
        TGraph *ER_Spectrum_Bef = new TGraph(reso_T,T_QF_Original_Bef_Array,Factor1_Original_Bef_Array);
        ER_Spectrum_Bef->GetXaxis()->SetTitle("Energy[keVee]");
        ER_Spectrum_Bef->GetYaxis()->SetTitle("Count[Evts/kg/day]");
        ER_Spectrum_Bef->GetXaxis()->SetLimits(1e-2,1e2);
        ER_Spectrum_Bef->GetYaxis()->SetRangeUser(1e-8,1e+15);
        ER_Spectrum_Bef->SetName("ER_Spectrum_Bef");
        ER_Spectrum_Bef->SetLineColor(2);
        ER_Spectrum_Bef->SetMarkerColor(2);
        ER_Spectrum_Bef->SetLineWidth(5);
            
        TGraph *ER_Spectrum_Aft = new TGraph(reso_T,T_QF_Original_Aft_Array,Factor1_Original_Aft_Array);
        ER_Spectrum_Aft->GetXaxis()->SetTitle("Edet[keVee]");
        ER_Spectrum_Aft->GetYaxis()->SetTitle("Count[Evts/kg/day]");
        ER_Spectrum_Aft->SetName("ER_Spectrum_Aft");
        ER_Spectrum_Aft->SetLineColor(3);
        ER_Spectrum_Aft->SetMarkerColor(3);
            /*
        double Data_RE_TEXONO[4]= {0.225,0.275,0.325,0.375};//CDEX Brem: 0.265(for threshold 0.25), CDEX Migdal: 0.2, CDEX Brem: 0.206
        double Data_RATE_TEXONO[4]= {49.432,33.209,8.993,13.404};
        double Data_RE_Err[4]= {0.025,0.025,0.025,0.025};
        double Data_RE_Rate_Err[4]= {12.451,7.239,4.753,3.445};
            */
            
        double *RE_DATA_1    =Hist_SetLimit_Plot_v2_Extract_Peak(0);
        double *RE_Rate_1    =Hist_SetLimit_Plot_v2_Extract_Peak(1);
        double *RE_DATA_Err_1=Hist_SetLimit_Plot_v2_Extract_Peak(2);
        double *RE_Rate_Err_1=Hist_SetLimit_Plot_v2_Extract_Peak(3);
             
            /*
            for(int kkk=0; kkk<257; kkk++)
            {
                cout << "RE_DATA_1: " << RE_DATA_1[kkk] << endl;
                cout << "RE_Rate_1: " << RE_Rate_1[kkk] << endl;
                cout << "RE_DATA_Err_1: " << RE_DATA_Err_1[kkk] << endl;
                cout << "RE_Rate_Err_1: " << RE_Rate_Err_1[kkk] << endl;
            }
             */
            
        TGraphErrors *TEXONOData = new TGraphErrors(257,RE_DATA_1,RE_Rate_1,RE_DATA_Err_1,RE_Rate_Err_1);
        TEXONOData->SetName("TEXONOData");
             
        TLegend *leg = new TLegend(0.3,0.6,0.5,0.9);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
                    
            cout << "Mass: " << Mass << endl;
        leg->AddEntry("",Form("%s:#sigma_{SI}:%.2f #times 10^{-%.f}, M_{#chi}=%.2fGeV",Type_of_Model[Type_of_Model_INT].c_str(),CNFNV(0,Sigma_SI_Array[Point_Number]),CNFNV(1,Sigma_SI_Array[Point_Number]),mx),"");
        leg->AddEntry(ER_Spectrum_Bef,"Recoil Spectrum Before the attenuation","l");
        leg->AddEntry(ER_Spectrum_Aft,"Recoil Spectrum After the attenuation","l");
            /*
        leg->AddEntry("",Form("Ratio of >200eV: %.3f #times 10^{-%.f}",CNFNV(0,RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1),CNFNV(1,RecoilX_Event_Aft_EARTH_M1/RecoilX_Event_Original_M1)),"");
        leg->AddEntry("",Form("Ratio of >200eV and <250eV: %.3f #times 10^{-%.f}",CNFNV(0,RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3),CNFNV(1,RecoilX_Event_Aft_EARTH_M3/RecoilX_Event_Original_M3)),"");
            */
        ER_Spectrum_Bef->Draw("ALP");
        ER_Spectrum_Aft->Draw("LPsame");
        TEXONOData->Draw("Psame");
        leg->Draw();
        c3->SetLogy();
        c3->SetLogx();
                cout << "=================================================================================" << endl;
            char fout_name1[300];
            sprintf(fout_name1,Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/Recoil_Spectrum/%s_%i_STS_Bent.root",Mass_Point[Mass_INT].c_str(),Type_of_Model[Type_of_Model_INT].c_str(),FILE));
            TFile *fout1=new TFile(fout_name1,"recreate");
            ER_Spectrum_Bef->Write();
            ER_Spectrum_Aft->Write();
            TEXONOData->Write();
                cout << "=================================================================================" << endl;

        if(Take_Plot==1){
        c3->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/Recoil_Spectrum/%s_%i_STS_Bent.pdf",Mass_Point[Mass_INT].c_str(),Type_of_Model[Type_of_Model_INT].c_str(),FILE));
            cout << "=================================================================================" << endl;
        }
     
            
        Point_Number = Point_Number+1;
        cout << "Point_Number: " << Point_Number << endl;
        }//Close2
    }//Close1

        Sigma_SI_Array[4] = 3.5e-29;
        Sigma_SI_With_Threshold_M3[4] =0;
        
        TCanvas *c4 = new TCanvas("c4");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        
        for(int kkk=0; kkk<Point_Number;kkk++)
        {
            cout << "Sigma_SI_Array: " << Sigma_SI_Array[kkk] << endl;
            cout << "Sigma_SI_With_Threshold_M1[Point_Number] : " << Sigma_SI_With_Threshold_M1[kkk]  << endl;
            cout << "Sigma_SI_With_Threshold_M3[Point_Number] : " << Sigma_SI_With_Threshold_M3[kkk]  << endl;
        }

        TGraph *Threshold_Plot_2 = new TGraph(Point_Number+1,Sigma_SI_Array,Sigma_SI_With_Threshold_M3);
        Threshold_Plot_2->SetLineColor(5);
        Threshold_Plot_2->SetMarkerColor(5);
        Threshold_Plot_2->SetName("Threshold_Plot");
        Threshold_Plot_2->GetXaxis()->SetLimits(1e-42,1e-26);
        Threshold_Plot_2->GetYaxis()->SetRangeUser(1e-42,1e-28);
        Threshold_Plot_2->SetTitle("");
        Threshold_Plot_2->GetXaxis()->SetTitle("Sigma_SI");
        Threshold_Plot_2->GetYaxis()->SetTitle("(Sigma_SI) #times (Event_Fraction)");

        TGraph *Threshold_Plot_4 = new TGraph(Point_Number+1,Sigma_SI_Array,Sigma_SI_With_Threshold_M1);
        Threshold_Plot_2->SetLineColor(7);
        Threshold_Plot_2->SetMarkerColor(7);

        TF1 *Linear_Line = new TF1("Linear_Line","x",1e-42,1e-27);
        Linear_Line->SetLineColor(4);
        
        Threshold_Plot_2->Draw("AL");
        Linear_Line->Draw("Lsame");
        
        TLegend *leg = new TLegend(0.1,0.5,0.4,0.8);
        leg->SetFillColor(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetTextFont(22);
        leg->AddEntry("",Form("M_{#chi}=%.2fGeV",Mass),"");
        if(Type_of_Model_INT==0)leg->AddEntry(Threshold_Plot_2,Form("Bin: %i",NU_Bins[kkk]+1),"lP");
        if(Type_of_Model_INT==1)leg->AddEntry(Threshold_Plot_2,Form("Bin: %i",MD_Bins[kkk]+1),"lP");
        if(Type_of_Model_INT==2)leg->AddEntry(Threshold_Plot_2,Form("Bin: %i",BR_Bins[kkk]+1),"lP");

        //leg->AddEntry(Threshold_Plot_4,Form("Count: -300eV"),"lP");
        leg->AddEntry(Linear_Line ,"X=Y","l");
        leg->Draw();


        gPad->SetLogx();
        gPad->SetLogy();
        if(Take_Plot==1)
        {
            char fout_name[300];
            sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/%s_STS_Bent.root",Mass_Point[Mass_INT].c_str(),Type_of_Model[Type_of_Model_INT].c_str()));
            TFile *fout=new TFile(fout_name,"recreate");
            Threshold_Plot_2->Write();
            c4->Print(Form("/Users/yehchihhsiang/Desktop/Analysis/CDEX_Analysis_method/Codes/3_CRESST_Flux/%sGeV/%s_STS_Bent.pdf",Mass_Point[Mass_INT].c_str(),Type_of_Model[Type_of_Model_INT].c_str()));
        }
    }
}
    

