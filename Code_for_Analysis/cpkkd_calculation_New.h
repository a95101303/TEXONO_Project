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

#include "TCutG.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TAxis.h"
//#include "B_L_Henke_data_PE_f1_f2.h"

//Point_Number


const int Number_for_Brem = 251;
const int n_point_Brem = Number_for_Brem ;

/*
const int reso_T=Number_for_Brem ;
const int dm_spec_resolution=Number_for_Brem;
const int n_point = Number_for_Brem ;
*/


const int reso_T=2510;
//const int reso_T=251;

const int dm_spec_resolution=2510;
const int n_point = 251;


const int interp = 50;
const int n_point1 = 251*interp;
const int dotnumber = 1e5;


//Recoil_Energy_array_extension
//Fixed_parameter
//For BREM
/*
const int reso_T=25100;
const int dm_spec_resolution=25100;
const int n_point_Brem = 25100;
*/
//For MIGDAL,NUCLEAR
/*
const int reso_T=2000;
const int dm_spec_resolution=2000;
const int n_point1 = 251;
*/
const int Data_element=257;

/*
const double a0_dE = 30.0/1000.0;
const double a1_dE = 50.0/(1000.0*sqrt(10.36));
 */
const double a0_dE = 33.4992/1000.0;
const double a1_dE = 13.2145/1000.0;

//for Lindhard k = 0.22
const double p0 = 0.888572;
const double p1 = 1.35418;

double antiquenching(double x)  //input f*ER, output ER
{
  double f;
    double lnf = p0*TMath::Log(x)+p1;  // x[0] = f*ER
    f = TMath::Exp(lnf);
  return f;
}
double counttrans(double x)  //translate counts (Qnr)
{
  double f;
  double lnf = p0*TMath::Log(x)+p1;  // x[0] = f*ER
  double trans = TMath::Exp(lnf)*p0/x;
  f = trans;

  return f;
}

//

float Find_minimum_WIMP_mass(float Minimum_energy_recoil, double A = AGe)//keV
{
    float Find_WIMP_mass;
    cout << "1: " << endl;
    for(int k=0 ; k<100 ; k++)
        {
            double WIMP_max_T = 1000.0*max_recoil_A(1+0.01*k, 779.135*1000.0/2.99792458e8, A);
            if(Find_WIMP_mass==0 and WIMP_max_T>=Minimum_energy_recoil)
            {
                Find_WIMP_mass = 1+0.01*k;
            }
        }
    return Find_WIMP_mass;
}

double test(double *array_test)
{
    return array_test[0];
}
//=============================
double *T_QF_Array(double WIMP_mx=10)//X
{
    static double T_QF[reso_T];double T[reso_T];double A = AGe;
    double WIMP_max_T = 1000.0*max_recoil_A(WIMP_mx, 779.135*1000.0/2.99792458e8, A)+1.8; //keV
    
    
    for(int i=0;i<reso_T;i++)
    {
        T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
        T_QF[i] = TQF(T[i]);
    }
    return T_QF;
}
double Total_Flux(double WIMP_mx)
{
    double sum=0;
    for(int j=0;j<2000;j++)
    {
        sum = sum + velo_dist_Ave[j][3];
    }
    double Total_Flux=0;
    for(int j=0;j<2000;j++)
    {
        float v_cm_day = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_cmday1;
        Total_Flux = Total_Flux + (rohx/WIMP_mx)*v_cm_day*(1/(sum))*velo_dist_Ave[j][3];
        /*
        cout << "(1/(sum))*velo_dist_Ave[j][3]: " << (1/(sum))*velo_dist_Ave[j][3] << endl;
        cout << "v: " << 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])<< endl;
        cout << "v_cm_day: " << v_cm_day<< endl;
        cout << "N_atom_Xe_1kg*(rohx/WIMP_mx)*v_cm_day*(1/(sum))*velo_dist_Ave[j][3]: " << N_atom_Xe_1kg*(rohx/WIMP_mx)*v_cm_day*(1/(sum))*velo_dist_Ave[j][3] << endl;
         */
    }
    return Total_Flux;
}
static vector<double> Possibilities_array;
static int check=0;//For running for once

double Possibility()//v(km/s)
{
    if(check==0)
    {
        double Vecolity[2000];double Pre_Possiblity[2000];
        double Pre_Sum=0;
        for(int JJJ=0; JJJ<2000; JJJ++){Pre_Sum = Pre_Sum + velo_dist_Ave[JJJ][3];}
        cout << "Pre_Sum: " << Pre_Sum  << endl;
        for(int JJJ=0; JJJ<2000; JJJ++){Pre_Possiblity[JJJ] = velo_dist_Ave[JJJ][3]/Pre_Sum;}

        double sum_for_Possibility=0;
        double Possiblity[2000];
        //cout << "sum_for_Possibility: " << sum_for_Possibility << endl;
        for(int j=0;j<2000;j++)
        {
            float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2]);//(km/s)
            Vecolity[j] = v;
            sum_for_Possibility = sum_for_Possibility + Pre_Possiblity[j]*0.395*1e5;
        }
        for(int j=0;j<2000;j++)
        {
            Possiblity[j] = Pre_Possiblity[j]/sum_for_Possibility;

            //cout << "Possiblity[j]: " << Possiblity[j] << endl;
        }
        //cout << "sum_for_Possibility: " << sum_for_Possibility << endl;
        //cout << "1/sum_for_Possibility: " << 1/sum_for_Possibility << endl;

        TH1F   *Flux_HIST = new TH1F("Flux_HIST","Flux_HIST",2000,0,791);
        for(int kkk=0;kkk<2000;kkk++){Flux_HIST->SetBinContent(kkk+1,Possiblity[kkk]);}

        double Total_Prob_Temp=0;
        for(int j=0; j<DM_Beta_Right.size()-1; j++)
        {
            double Pob_temp= 0;
            double LLL=0;

            int binx_start    = Flux_HIST->GetXaxis()->FindBin(DM_Beta_Right[j]*(3e8/1e3)*(1e-3));
            int binx_end      = Flux_HIST->GetXaxis()->FindBin(DM_Beta_Right[j+1]*(3e8/1e3)*(1e-3));
            //cout << "binx_start: " << binx_start << endl;
            //cout << "binx_end: "   << binx_end << endl;

            for(int kkk=binx_start; kkk<binx_end+1; kkk++)
            {
                Pob_temp = Pob_temp + Flux_HIST->GetBinContent(binx_end);
                LLL = LLL + 1;
            }
            Total_Prob_Temp = Total_Prob_Temp + Pob_temp;
            Possibilities_array.push_back(Pob_temp/(double)LLL);
        }

        /*
        for(int j=0; j<bin_for_ER-1; j++)
        {
            double Pob_temp= 0;
            int binx_start    = Flux_HIST->GetXaxis()->FindBin(vel_dist_ER_kms[j]);
            int binx_end      = Flux_HIST->GetXaxis()->FindBin(vel_dist_ER_kms[j+1]);
            for(int kkk=binx_start; kkk<binx_end+1; kkk++)
            {
                //cout << "Flux_HIST->GetBinContent(kkk): " << Flux_HIST->GetBinContent(kkk) << endl;
                Pob_temp = Pob_temp + Flux_HIST->GetBinContent(binx_end);
            }
            cout << "1/sum: " << 1/sum << endl;
            cout << "vel_dist_ER_kms[j]: " << vel_dist_ER_kms[j] << endl;
            cout << "vel_dist_ER_kms[j+1]: " << vel_dist_ER_kms[j+1] << endl;
            cout << "Pob_temp: " << Pob_temp << endl;
            Total_Prob_Temp = Total_Prob_Temp + Pob_temp;
            Possibilities_array.push_back(Pob_temp);
        }
         */
        
            check=1;
    }
    /*
    for(int kkk=0; kkk<Possibilities_array.size(); kkk++)
    {
        cout << "Possibilities_array_2: " << Possibilities_array[kkk] << endl;
        cout << "kkk: " << kkk << endl;
    }
     */
    return 0;
}


double *RecoilX_Event(int Option, TH1F *Flux,double WIMP_mx,double Sigma_SI,int Model_of_Interaction, int Conventional_or_not, int File_index)
{ //Model_of_Interaction=0 for Nuclear, Model_of_Interaction=1 for Migdal(N=2), Model_of_Interaction=3 for Brem, Model_of_Interaction=4 for Migdal MPA
  //Conventional_or_not=0 for real distributions, Conventional_or_not=1 for conventional distributions
    static double Mass_TEMP=0.09;
    
    double sum=0;
    for(int j=0;j<2000;j++)
    {
        sum = sum + velo_dist_Ave[j][3];
    }
    std::cout<<"-----------start-------------"<<std::endl;
    //======================For the migdal Effect======================//
    //const int interp = 50;
    //const int interp = 20;
    const int interp = 2;
    const int n_point1 = 251*interp;
    const int dotnumber = 1e5;

    cout << " n_point1: " <<  n_point1 << endl;
    static int check=0;
    static double dpdE[n_point1][lnumber];
    static double Eem_keV[n_point1];  // in keV
    double transvalue = 3.e08*0.85/(2.*TMath::Pi());//(2.*TMath::Pi())*0.8;/////////////////////
    const int lnumber = 8;
    int nl[5] = {1,2,3,3,2};
    double Enl_Series[lnumber] = {11,1.4,1.2,0.17,0.12,0.035,0.015,0.0065};  //energy nl in keV


    //=====================For the Brem Effect======================//
    
    const int Ew_point = 520;
    //const int Ew_point_interp = 520*100;
    //const int Ew_point_interp = 520*50;
    //const int Ew_point_interp = 520*50;
    //const int Ew_point_interp = 520*20;
    const int Ew_point_interp = 520*5;
    
    const double Ewstart = 1e-2;  // keV
    const double Ewfinish = 30; // keV
    const double Edetstart = 1e-2;  // keV
    const double Edetfinish = 30; // keV
    static double Ew[Ew_point];static double fw1[Ew_point];static double fw2[Ew_point];static double fw[Ew_point];
    static double Ew_interp[Ew_point_interp];static double fw_interp[Ew_point_interp];
    static double Edett[n_point];


    if(check==0){//Open1
    //======================For the migdal Effect======================//
    std::ifstream ifs;
    ifs.open("Ge.dat",std::ifstream::in);
    double Ee[n_point+1][lnumber],prob[n_point+1][lnumber];
    double Ee_keV[n_point+1][lnumber],prob_trans[n_point+1][lnumber];
    
  for (int i=0;i<lnumber;i++) {
    for (int j=0;j<n_point+1;j++) {
      ifs >> Ee[j][i] >> prob[j][i]; // energy in eV
    }}

  ifs.close();

  //interp dpdE
  for (int nli=0;nli<lnumber;nli++) {
    for (int i=0;i<n_point1;i++) {
      if (i<100) {
        Eem_keV[i] = Ee[1][0]/1000./100.*(double)(i+1);
        dpdE[i][nli] = prob[1][nli]*transvalue;
      }
      else {
        int in = i/interp;
        int inn = i-in*interp;
        Eem_keV[i] = (Ee[in+1][0]-Ee[in][0])/1000./((double)interp)*(double)(inn+1)+Ee[in][0]/1000.;
        dpdE[i][nli] = prob[in][nli]+(prob[in][nli]-prob[in+1][nli])/(Ee[in][0]/1000.-Ee[in+1][0]/1000.)*(Eem_keV[i]-Ee[in][0]/1000.);
        dpdE[i][nli] = dpdE[i][nli]*transvalue;
          if(nli==0)
          {
              //cout << "Eem_keV[i]: " << Eem_keV[i] << endl;
              //cout << "dpdE[i][nli]: " << dpdE[i][nli] << endl;
          }
      }
    }
  }
        //=====================For the Brem Effect======================//
              //read f(w)
        std::ifstream ifss;
        ifss.open("ge_nist.nff",std::ifstream::in);
          for (int i=0;i<Ew_point;i++) {
            ifss >> Ew[i] >> fw1[i] >> fw2[i]; // energy in eV
            Ew[i] = Ew[i];  // energy in keV
        //    std::cout << Ew[i] <<std::endl;
            if (fw1[i] <-100.) fw1[i] = 0.;
            fw[i] = TMath::Sqrt(fw1[i]*fw1[i]+fw2[i]*fw2[i]);
          }
          ifss.close();

          for (int i=0;i<Ew_point_interp;i++) {
            double mean = log(Ewstart)+(log(Ewfinish)-log(Ewstart))/((double)Ew_point_interp)*((double)i);
            Ew_interp[i] = TMath::Exp(mean);
          }

          //interp fw
          int index = 0;
          for (int i=0;i<Ew_point_interp;i++) {
            if (Ew_interp[i]<=Ew[0]) { fw_interp[i] = fw[0]; }
            else if (Ew_interp[i]>=Ew[Ew_point-1]) { fw_interp[i] = fw[Ew_point-1]; }
            else {
              while (Ew_interp[i]>=Ew[index]) {index++;}
              fw_interp[i] = (fw[index]-fw[index-1])/(Ew[index]-Ew[index-1])*(Ew_interp[i]-Ew[index])+fw[index];
            }
          }
        
          for (int i=0;i<n_point;i++) {
              double mean = log(Edetstart)+(log(Edetfinish)-log(Edetstart))/((double)n_point)*((double)i);
              Edett[i] = TMath::Exp(mean);
  }


    }  //Close1
    check = check + 1;

    for (int i=0;i<Ew_point_interp;i++) {
        //cout << "Ew_interp[i]: " << Ew_interp[i] << "fw_interp[i]: " << fw_interp[i] << endl;
    }
/*
    for (int i=0;i<lnumber;i++)
    {
        for (int j=0;j<n_point+1;j++)
        {
        ifs >> Ee[j][i] >> prob[j][i]; // energy in eV
        Ee_keV[j][i] = Ee[j][i]*1e-3; //keV
        if(j==0){prob_trans[j][i]=0;}
        else{prob_trans[j][i] = prob[j][i]*transvalue;} //Prob
        if(i==1 or i==2)
            {
                cout << "Ee_keV[j][i]: " << Ee_keV[j][i] << endl;
                cout << "prob_trans[j][i]: " << prob_trans[j][i] << endl;
            }
        }
    }
     */
    //==================================================================//
    cout << "WIMPMass: " << WIMP_mx << endl;
    cout << "Option: " << Option << endl;
     


    //==============Constant_for_calculation==================
    double A = AGe;//Only for CRESST
    //double A = 20;
    
    static double T[reso_T];static double T_QF[reso_T];double recoil[reso_T];
    double v;double  MaxV=0;
    if(Conventional_or_not==0){
        MaxV = Flux->GetXaxis()->GetBinCenter(Flux->FindLastBinAbove());
        cout << "MaxV: " << MaxV << endl;}
    if(Conventional_or_not==1){
        MaxV = 779.135;
        cout << "MaxV: " << MaxV << endl;}


    
    double WIMP_max_T =0;
    if(Model_of_Interaction==0 and Conventional_or_not==0){//Nuclear Recoil
        WIMP_max_T = 1000.0*max_recoil_A(WIMP_mx, MaxV*1000.0/2.99792458e8, A)+0.2;} //keV
    if(Model_of_Interaction==1 and Conventional_or_not==0){//**Migdal Effect**
        WIMP_max_T = max_recoil_A_EM_keV(WIMP_mx, MaxV*1000.0/2.99792458e8, A)+0.2;} //keV
    if(Model_of_Interaction==2 and Conventional_or_not==0){//Brem
        WIMP_max_T = max_recoil_A_EM_keV(WIMP_mx, MaxV*1000.0/2.99792458e8, A)+0.2;} //keV
    if(Model_of_Interaction==3 and Conventional_or_not==0){//Migdal Effect: MPA
        WIMP_max_T = max_recoil_A_EM_keV(WIMP_mx, MaxV*1000.0/2.99792458e8, A)+0.2;} //keV

    if(Conventional_or_not==1){//Conventional Distribution
        WIMP_max_T = max_recoil_A_EM_keV(WIMP_mx, MaxV*1000.0/2.99792458e8, A)+0.2;} //keV
    if(Conventional_or_not==2){//Conventional Distribution
        WIMP_max_T = 1000.0*max_recoil_A(WIMP_mx, 779.135*1000.0/2.99792458e8, A)+0.2;} //keV
    if(Model_of_Interaction==4){//electron Recoil
        WIMP_max_T = Energy_DM(WIMP_mx,779*kms1_to_c)+0.02;
    } //keV
    double WIMP_max_T_NU = 1000.0*max_recoil_A(WIMP_mx, MaxV*1000.0/2.99792458e8, A);
    double WIMP_max_T_EM = max_recoil_A_EM_keV(WIMP_mx, MaxV*1000.0/2.99792458e8, A);

    //cout << "(WIMP_max_T_NU): " << (1000.0*max_recoil_A(WIMP_mx, MaxV*1000.0/2.99792458e8, A)) << endl;
    //cout << "WIMP_max_T_EM: " << max_recoil_A_EM_keV(WIMP_mx, MaxV*1000.0/2.99792458e8, A) << endl;

    /*
    cout << "(max_recoil_A_keV(WIMP_mx, v, A)): " << (max_recoil_A_keV(WIMP_mx, MaxV*kms1_to_c, A)) << endl;
    cout << "TQF(T[i]): " << TQF(1) << endl;
     */
    //==========Maximum_energy_recoil_based_on_Dark_Matter==========
    // Data_implemented_no_rebin
    double Event_Number_accumulated=0;
    //Min_Point_found
    //===============================================================
    //Theoretical line of the rate of WIMP Method 2
    static double recoilX[reso_T];
    
    //Fuck!!!!!Please remember to initiate
    for(int i=0;i<reso_T;i++)
    {
        recoilX[i] = 0.0;
    }

    if(Model_of_Interaction==0)//Nuclear-recoil Only
    {
        cout << "Nuclear Recoil" << endl;
        for(int i=0;i<reso_T;i++)
        {
            T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
            T_QF[i] = TQF(T[i]); //Normal Case for Ge
            //T_QF[i] = (T[i]);//CRESST CASE QF=1
            recoilX[i] = 0.0;
            //cout << "T[i]: " << T[i] << endl;
            for(int j=0;j<2000;j++)
            {
                float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_c;
                float v_cm_day = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_cmday1;

                    if((max_recoil_A_keV(WIMP_mx, v, A))>T[i] and Model_of_Interaction==0)
                    {
                        cout << "fdsigma_dT_keV(WIMP_mx, Sigma_SI, v, A, T[i]): " << fdsigma_dT_keV(WIMP_mx, Sigma_SI, v, A, T[i]) << endl;
                        if(Conventional_or_not==0)recoilX[i] = recoilX[i] + rate_scale_QF(T[i])*fdsigma_dT_keV(WIMP_mx, Sigma_SI, v, A, T[i])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(Flux->GetBinContent(j))/(Flux->GetEntries());
                        if(Conventional_or_not==1)recoilX[i] = recoilX[i] + rate_scale_QF(T[i])*fdsigma_dT_keV(WIMP_mx, Sigma_SI, v, A, T[i])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1/(sum))*velo_dist_Ave[j][3];
                        //if(Conventional_or_not==2)recoilX[i] = recoilX[i] + 365*rate_scale_QF(T[i])*AAAA_keV(WIMP_mx, Sigma_SI, v, 131, T[i])*N_atom_Xe_1kg*(rohx/WIMP_mx)*v_cm_day*(1/(sum))*velo_dist_Ave[j][3];//For Heavy DM
                        if(Conventional_or_not==2)recoilX[i] = recoilX[i] + 365*rate_scale_QF(T[i])*AAAA_keV(WIMP_mx, Sigma_SI, v, 72.64, T[i])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1/(sum))*velo_dist_Ave[j][3];//For Heavy DM
                        //cout << "T_QF[i]: " << T_QF[i] << endl;
                        //if(Conventional_or_not==0)recoilX[i] = recoilX[i] + (T[i])*fdsigma_dT_keV(WIMP_mx, Sigma_SI, v, A, T[i])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(Flux->GetBinContent(j))/5000;//Only for CRESST
                    }
                
            }
        }
    }
    static double Velocity_Min=799;

    if(Mass_TEMP!=WIMP_mx)
    {
        Velocity_Min=799;
        Mass_TEMP=WIMP_mx;
    }
    
    double Minimum=0.16;

    if(Model_of_Interaction==1)//Migdal Effect
    {//Open9
        cout << "Migdal Effect" << endl;
        for(int Number_of_Level=1;Number_of_Level<6;Number_of_Level++)
        {//Open8
            cout << "Number_of_Level: " << Number_of_Level << endl;
            for(int i=0;i<reso_T;i++)
            {//Open7
                T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
                
                T_QF[i] = TQF(T[i]);
                //recoilX[i] = 0.0;
                double Edet=T[i];
                if(Edet<WIMP_max_T_EM and Option>0)
                {//Open6
                    //cout << "WIMP_max_T_EM: " << WIMP_max_T_EM << endl;
                    for(int kkk=0; kkk<n_point1; kkk++)
                    {//Open5
                        //For Electronic Recoil
                        double Eem = Eem_keV[kkk];//Energy after the ionization (keV)
                        double Enl = Enl_Series[Number_of_Level];//Binding Energy (keV)
                        double Eel = Eem + Enl;//Total Energy
                        //For Nuclear Recoil
                        double Eee = Edet - Eel;
                        double ER  = antiquenching(Eee);
                        double NU_Min_V = min_v_keV(WIMP_mx,ER,A);//beta
                        
                        //cout << "counttrans(Eee): " << counttrans(Eee) << endl;
                        if(ER > 0 and ER < max_recoil_A_keV(WIMP_mx,MaxV*kms1_to_c, A) and Option==1)
                        {//Open1
                            //cout << "NU_Min_V: " << NU_Min_V*3e8/1e3 << endl;
                            //cout << "min_v_EM(WIMP_mx,A,ER,Eel): " << min_v_EM(WIMP_mx,A,ER,Eel) << endl;
                            for(int j=0;j<2000;j++)
                            {//Open2
                                    float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_c;
                                    float v_cm_day = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_cmday1;
                                

                        if((max_recoil_A_EM_keV(WIMP_mx, v, A))>T[i] and v>min_v_EM(WIMP_mx,A,ER,Eel)*kms1_to_c)//T is for Edet
                            { //Open3
                            if( Eel > Edet-TQF(max_recoil_A_keV(WIMP_mx,v, A)) and Eel < max_recoil_A_EM_keV(WIMP_mx, v, A))
                                { //Open4
                                    
                                        if(Conventional_or_not==0)recoilX[i] = recoilX[i] +  counttrans(Eee)*(Eem_keV[kkk+1]-Eem_keV[kkk])*fdsigma_dERdEEM_keV(dpdE[kkk][Number_of_Level]*(2*ER/(72.64*0.95)),WIMP_mx, Sigma_SI, v, A, ER)*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(Flux->GetBinContent(j))/(Flux->GetEntries());//For the real case
                                        if(Conventional_or_not==1)recoilX[i] = recoilX[i] + counttrans(Eee)*(Eem_keV[kkk+1]-Eem_keV[kkk])*fdsigma_dERdEEM_keV(dpdE[kkk][Number_of_Level]*(2*ER/(72.64*0.95)),WIMP_mx, Sigma_SI, v, A, ER)*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1/sum)*velo_dist_Ave[j][3];//For the normal distribution
                                    if(v*3e8/1e3<Velocity_Min)
                                    {
                                        cout << "WIMP_mx: " << WIMP_mx << endl;
                                        cout << "Velocity_Min: " << Velocity_Min << endl;
                                        Velocity_Min=v*3e8/1e3;
                                    }
                                }//Close4
                            }//Close3
                            }//Close2
                        } //Close1
                        
                    }//Close5
                }//Close6
            }//Close7
        }//close8
    }//Close9
    
    if(Model_of_Interaction==2)//Brem
    {//Open9
        cout << "BREM" << endl;
        cout << "Sigma_SI: " << Sigma_SI << endl;
        cout << "WIMPMASS: " << WIMP_mx << endl;
        for(int i=0;i<reso_T;i++)
        {//Open8
          //  T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
        }//Close8
        //for(int i=0;i<reso_T;i++)
        for(int i=0;i<n_point_Brem;i++)
        {//Open7
            //cout << "T[i]: " << T[i] << endl;
            //double Edet=T[i];
            double Edet=Edett[i];//Edet = ER*Q+EB
            T[i] = Edett[i];
            if(Edet<WIMP_max_T_EM and Option>0)
            {//Open6
                //cout << "WIMP_max_T_EM: " << WIMP_max_T_EM << endl;
                //for(int kkk=0; kkk<BIN_henke_f; kkk++)
                for(int kkk=0; kkk<Ew_point_interp; kkk++)
                {//Open5
                    //For Brem Energy
                    //double EB = f_henke[kkk][0]/1e3;//Energy of Brem (keV)
                    double EB = Ew_interp[kkk];//Energy of Brem (keV)
                    //For Nuclear Recoil
                    double Eee = Edet - EB;
                    double ER  = antiquenching(Eee);
                    double NU_Min_V = min_v_keV(WIMP_mx,ER,A);//beta
                    
                    if(ER > 0  and ER < max_recoil_A_keV(WIMP_mx,MaxV*kms1_to_c, A) and (WIMP_max_T_EM - EB)>0 and Option==1)
                    {//Open1
                        for(int j=0;j<2000;j++)
                        {//Open2
                            float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_c;//Beta
                            float v_cm_day = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_cmday1;
                            
                            //if((max_recoil_A_EM_keV(WIMP_mx, v, A))>T[i] and v>min_v_BM(WIMP_mx,A,EB)*kms1_to_c and v>NU_Min_V)//T is for Edet
                                if((max_recoil_A_EM_keV(WIMP_mx, v, A))>T[i])//T is for Edet
                                {//Open3
                                    //double ER_Max_NU=max_recoil_A_keV(WIMP_mx,v, A);
                                    //double ER_Max_BE=max_min_recoil_A_BM_keV(0,WIMP_mx,v,A,EB);
                                    //double ER_Min_BE=max_min_recoil_A_BM_keV(1,WIMP_mx,v,A,EB);
                                    double EB_Min_BE=Edet-TQF(max_min_recoil_A_BM_keV(0,WIMP_mx,v,A,EB));
                                    double EB_Max_BE=Edet-TQF(max_min_recoil_A_BM_keV(1,WIMP_mx,v,A,EB));

                                    //if(ER<ER_Max_BE and ER>ER_Min_BE and ER<ER_Max_NU)
                                    if(EB>EB_Min_BE and EB<EB_Max_BE)
                                    {//Open4
                                        
                                        int Option=0;
                                        if(Conventional_or_not==0 and kkk==0)recoilX[i] = recoilX[i] +  counttrans(Eee)*(Ew_interp[kkk])*fdsigma_dERdw_keV(Option,WIMP_mx,v,Sigma_SI,A,ER,EB,fw_interp[kkk])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(Flux->GetBinContent(j))/(Flux->GetEntries());//For the real case
                                        if(Conventional_or_not==0 and kkk!=0)recoilX[i] = recoilX[i] +  counttrans(Eee)*(Ew_interp[kkk+1]-Ew_interp[kkk])*fdsigma_dERdw_keV(Option,WIMP_mx,v,Sigma_SI,A,ER,EB,fw_interp[kkk])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(Flux->GetBinContent(j))/(Flux->GetEntries());//For the real case
                                        //=====================
                                        if(Conventional_or_not==1 and kkk==0)recoilX[i] = recoilX[i] +  counttrans(Eee)*(Ew_interp[kkk])*fdsigma_dERdw_keV(Option,WIMP_mx,v,Sigma_SI,A,ER,EB,fw_interp[kkk])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(1/sum)*velo_dist_Ave[j][3];//For the real case
                                        if(Conventional_or_not==1 and kkk!=0)recoilX[i] = recoilX[i] +  counttrans(Eee)*(Ew_interp[kkk+1]-Ew_interp[kkk])*fdsigma_dERdw_keV(Option,WIMP_mx,v,Sigma_SI,A,ER,EB,fw_interp[kkk])*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(1/sum)*velo_dist_Ave[j][3];//For the real case
                                         
                                        if(v*3e8/1e3<Velocity_Min)
                                        {
                                            //cout << "Velocity_Min: " << Velocity_Min << endl;
                                            Velocity_Min=v*3e8/1e3;
                                        }
                                }//Close4
                            }//Close3
                        }//Close2
                    } //Close1
                }//Close5
            }//Close6
        }//Close7
    }//Close9
    
    if(Model_of_Interaction==3)//Migdal Effect: MPA
    {//Open9
        cout << "Migdal Effect: MPA" << endl;
        for(int Number_of_Level=1;Number_of_Level<6;Number_of_Level++)
        {//Open8
            cout << "Number_of_Level: " << Number_of_Level << endl;
            for(int i=0;i<reso_T;i++)
            {//Open7
                T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
                //cout << "T[i]: " << T[i] << endl;
                T_QF[i] = TQF(T[i]);
                double Edet=T[i];
                if(Edet<WIMP_max_T_EM and Option>0)
                {//Open6
                    for(int kkk=0; kkk<BIN_henke_f; kkk++)
                    //for(int kkk=0; kkk<n_point; kkk++)
                    {//Open5
                        //For Electronic Recoil
                        double Eel = (f_henke[kkk][0])*(1e-3);//Total Energy (keV)
                        double Enl = Enl_Series[Number_of_Level];//Binding Energy (keV)
                        double Eem = Eel - Enl;//
                        //For Nuclear Recoil
                        double Eee = Edet - Eel;
                        double ER  = antiquenching(Eee);
                        double NU_Min_V = min_v_keV(WIMP_mx,ER,A);//beta
                                        
                        if(Eem > 0 and ER > 0 and ER < max_recoil_A_keV(WIMP_mx,MaxV*kms1_to_c, A) and Option==1)
                        {//Open1
                            //cout << "Eel: " << Eel << endl;
                            for(int j=0;j<2000;j++)
                            {//Open2
                                    float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_c;
                                    float v_cm_day = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_cmday1;

                        if((max_recoil_A_EM_keV(WIMP_mx, v, A))>T[i] and v>min_v_EM(WIMP_mx,A,ER,Eel)*kms1_to_c)//T is for Edet
                            { //Open3
                            if( Eel > Edet-TQF(max_recoil_A_keV(WIMP_mx,v, A)) and Eel < max_recoil_A_EM_keV(WIMP_mx, v, A) and v>NU_Min_V)
                                { //Open4
                                    recoilX[i] = recoilX[i] +  ((1e-3)*(f_henke[kkk+1][0]-f_henke[kkk][0]))*(9e16)*fdsigma_dERdEr_MPA(WIMP_mx,v,Sigma_SI,A,ER,Eel,Cross_Section_MPA(Eel*1e3,f_henke[kkk][2]))*N_atom_Ge_1kg*(rohx/WIMP_mx)*v_cm_day*(1)*(Flux->GetBinContent(j))/10000;//For the real case
                                }//Close4
                            }//Close3
                            }//Close2
                        } //Close1
                        
                    }//Close5
                }//Close6
            }//Close7
        }//close8
    }//Close9

    double dm_mass= WIMP_mx;
      const double percm_GeV = 1.97326971780039025e-14;
  float  me = 0.510998928E-3;
  float red_mass = ((me*dm_mass)/(me+dm_mass));
  float alpha = (1E-37*TMath::Pi())/(pow(red_mass,2)*pow(percm_GeV,2));

    /*
    if(Model_of_Interaction==4)//Electronic-recoil Only
    {
        cout << "Electronic Recoil//(From my Code)" << endl;
        for(int i=0;i<reso_T;i++)
        {
            T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
            recoilX[i] = 0.0;//No need to attenuate
            //cout << "T[i]: " << T[i] << endl;

            for(int j=0;j<2000;j++)
            {
                float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2]);
                float v_cm_day = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2])*kms1_to_cmday1;
                
                //if((max_recoil_A_for_ER_keV(v,WIMP_mx))>T[i] and v<544 and v>sqrt(T[i]/(WIMP_mx*1e6)))
                if( Energy_DM(WIMP_mx,v*kms1_to_c)>T[i] and (v*kms1_to_c)>(sqrt(T[i]/(WIMP_mx*1e6))) and v<776)
                    {
                        if(v>50 and v<650 and T[i]<1.8)
                        {
                        //if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 1e3*1e-1*1e-18*N_atom_Xe_1kg*(rohx/WIMP_mx)*1e-15*fdsigma_dT_ER_New(filename,v,T[i])*v_cm_day*(1/(sum))*velo_dist_Ave[j][3]*0.395;
                        //if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 1e-15*fdsigma_dT_ER_New(v,T[i],WIMP_mx)*(1/(sum))*velo_dist_Ave[j][3]*v*kms1_to_c*(0.395*kms1_to_c);
                        //if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 1e-18*1e-15*fdsigma_dT_ER_New(v,T[i],WIMP_mx)*(1/(sum))*velo_dist_Ave[j][3]*v*kms1_to_c*(0.395*kms1_to_c);
                        if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 1e-36*1e-15*fdsigma_dT_ER_New(v,T[i],WIMP_mx)*(1/(sum))*velo_dist_Ave[j][3]*v*1e5*(0.395*1e5);

                        //if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 1e-15*fdsigma_dT_ER_New(v,T[i],WIMP_mx)*(1/(sum))*velo_dist_Ave[j][3]*v*kms1_to_c*(0.395*kms1_to_c);
                        
                        //if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 1e-18*1e-15*fdsigma_dT_ER_New(v,T[i],WIMP_mx)*(1/(sum))*velo_dist_Ave[j][3]*v*1e5*(0.395*1e5);
                        //if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 1e-18*alpha*1e-15*fdsigma_dT_ER_New(v,T[i],WIMP_mx)*(1/(sum))*velo_dist_Ave[j][3]*v*1e5*(0.395*1e5);

                        //if(T[i]>2.477480e-01)cout << "fdsigma_dT_ER(filename,v*kms1_to_c,T[i]): " << fdsigma_dT_ER(filename,v*kms1_to_c,T[i]) << endl;
                        }
                    }
            }
            //cout << "recoilX[i]: " << recoilX[i] << endl;
        }
    }
     */

    //const int data_Bin_ER = data_bin;
    const int data_Bin_ER = reso_T;
    if(Model_of_Interaction==4)//Electronic-recoil Only for Lahkwinder
    {
        cout << "Electronic Recoil//(From Lahkwinder's Code)" << endl;
        
        fdsigma_dT_ER_New(File_index,0.1,0.1,WIMP_mx);//For filling the
        Possibility();

        for(int kkk=0; kkk<Possibilities_array.size(); kkk++)
        {
            cout << "Possibilities_array_2: " << Possibilities_array[kkk] << endl;
            cout << "kkk: " << kkk << endl;
        }

        for(int i=1;i<data_Bin_ER;i++)
        //for(int i=0;i<data_bin;i++)
        //for(int i=1000;i<1001;i++)
        {
            T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
            //T[i] = reference[i][0]*1e-3; // keV ==>Used to do the check with Lakh's codes
            recoilX[i] = 0.0;//No need to attenuate

            double A_Check=0;double B_Check=0;
            //for(int j=0;j<bin_for_ER-1;j++)
            //for(int j=0;j<1;j++)
            for(int j=0;j<1999;j++)
            {
                //Definetly true==>Used to do the check with Lakh's codes
                /*
                double v_cms_now = vel_dist_ER_cms[j];//cm/s
                double v_cms_next = vel_dist_ER_cms[j+1];//cm/s
                double v_kms_now = vel_dist_ER_cms[j]/(1e5);// km/s
                double dv=(v_cms_next)-(v_cms_now);
                if(Conventional_or_not==1)recoilX[i] = recoilX[i] + 86400*1e-36*1e-15*fdsigma_dT_ER_New(v_kms_now,T[i],WIMP_mx)*Possibilities_array[j]*v_cms_next*(dv);
                cout << "v_kms_now: " << v_kms_now << endl;
                cout << "fdsigma_dT_ER_New(v_kms_now,T[i],WIMP_mx): " << fdsigma_dT_ER_New(v_kms_now,T[i],WIMP_mx) << endl;
                 */
                
                double v_cms_now  = (velo_dist_Ave[j][1])*1e5;//cm/s
                double v_cms_next = (velo_dist_Ave[j][2])*1e5;//cm/s
                double v_kms_now  = (velo_dist_Ave[j][1]);// km/s
                double dv=(v_cms_next)-(v_cms_now);
                double dR_Factor  = N_atom_1kg_Ge_Electron*DM_density_ER/(WIMP_mx);
                //double dR_Factor  = 1;
                
                //if(Conventional_or_not==1)recoilX[i] = recoilX[i] + dR_Factor*86400*1e-36*1e-15*fdsigma_dT_ER_New(v_kms_now,T[i],WIMP_mx)*(1/sum)*velo_dist_Ave[j][3]*v_cms_next*(dv);
                //if(v_kms_now>vel_dist_ER_kms[0] and v_kms_now<vel_dist_ER_kms[1]) A_Check = A_Check + fdsigma_dT_ER_New(v_kms_now,T[i],WIMP_mx)*(1/sum)*velo_dist_Ave[j][3];
                
                /*
                cout << "=======================================================" << endl;
                cout << "T[i]: " << T[i] << endl;
                cout << "(1/sum)*velo_dist_Ave[j][3]:  " << (1/sum)*velo_dist_Ave[j][3] << endl;
                 */

                double Vel_Temp=0;int Vel_block=0;
                for(int LMN=0;LMN<DM_Beta_Right.size()-1;LMN++)
                {
                    if(v_kms_now>DM_Beta_Now[LMN]*(3e8/1e3)*(1e-3)  and v_kms_now<DM_Beta_Now[LMN+1]*(3e8/1e3)*(1e-3) )
                    {
                        Vel_Temp=DM_Beta_Now[LMN+1]*(3e8/1e3)*(1e-3);//km/s
                        Vel_block=LMN;
                    }
                }
                double dsigma_dT_ER_temp = 0;
                if(File_index==0 or File_index==1)
                {
                    if(T[i]<0.012)//Xe doesn't have the data of the DCS for recoil energy below 0.012keV
                    {
                        dsigma_dT_ER_temp=0;
                    }
                    else
                    {
                        dsigma_dT_ER_temp= fdsigma_dT_ER_New(File_index,Vel_Temp,T[i],WIMP_mx);
                    }

                }
                if(File_index==2 or File_index==3)
                {
                    if(T[i]<0.08)//Xe doesn't have the data of the DCS for recoil energy below 0.08keV
                    {
                        dsigma_dT_ER_temp=0;
                    }
                    else
                    {
                        dsigma_dT_ER_temp= fdsigma_dT_ER_New(File_index,Vel_Temp,T[i],WIMP_mx);
                    }

                }
                if(File_index==0 or File_index==2)recoilX[i] = recoilX[i] + dR_Factor*86400*1e-36*1e-15*dsigma_dT_ER_temp*Possibilities_array[Vel_block]*v_cms_next*(dv);//c1
                if(File_index==1 or File_index==3)recoilX[i] = recoilX[i] + dR_Factor*86400*1e-18*1e-15*dsigma_dT_ER_temp*Possibilities_array[Vel_block]*v_cms_next*(dv);//d1
                //  cm3/(keV*day)

            }
            /*
            cout << "T[i]: " << T[i] << endl;
            cout << "recoilX[i]: " << recoilX[i] << endl;
            cout << "check_A: " << A_Check << endl;
            cout << "check_B: " << B_Check << endl;
             */
        }
    }
    for(int kkk=0; kkk<Possibilities_array.size(); kkk++)
    {
        //cout << "vel_dist_ER_cms[kkk]: " << vel_dist_ER_cms[kkk] << endl;
        //cout << "vel_dist_ER_cms[kkk+1]: " << vel_dist_ER_cms[kkk+1] << endl;
        //cout << "Possibilities_array: " << Possibilities_array[kkk] << endl;
    }
    //===============================================================
    double sig_E; double dEx;
    static double Factor1[dm_spec_resolution];
    
    //===============================================================
    if(Model_of_Interaction==0 and Option==1)
    {
        double pEx = T_QF[0];
        //Find out the Full energy spectrum with the resolution of the detector based on the theory
        for(int jjj=0 ; jjj<dm_spec_resolution ; jjj++)
        {
            Factor1[jjj]=0;
            for(int j=0 ; j<dm_spec_resolution ; j++)
            {
                //Find the theoretical line of the rate of WIMP
                if (recoilX[j]>0)
                {
                    sig_E = a0_dE + a1_dE*sqrt(T_QF[j]);
                    
                    if(j==0) { dEx = T_QF[0]; } else { dEx = T_QF[j] - pEx; }
                    
                    if( (recoilX[j]>0) )
                    {
                        Factor1[jjj] = Factor1[jjj] + recoilX[j]*exp(-pow((T_QF[j]-T_QF[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                    }
                    pEx = T_QF[j];
                }
                
            }
        }
    }
    if(Model_of_Interaction>=1 and Option==1)
    {
        double pEx = T[0];
        //Find out the Full energy spectrum with the resolution of the detector based on the theory
        for(int jjj=0 ; jjj<dm_spec_resolution ; jjj++)
        {
            Factor1[jjj]=0;
            for(int j=0 ; j<dm_spec_resolution ; j++)
            {
                //Find the theoretical line of the rate of WIMP
                if (recoilX[j]>0)
                {
                    sig_E = a0_dE + a1_dE*sqrt(T[j]);
                    
                    if(j==0) { dEx = T[0]; } else { dEx = T[j] - pEx; }
                    
                    if( (recoilX[j]>0) )
                    {
                        double Accumulated = recoilX[j]*exp(-pow((T[j]-T[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                        if(Accumulated>0)Factor1[jjj] = Factor1[jjj] + Accumulated;
                    }
                    pEx = T[j];
                }
                
            }
        }
    }
    
    cout << "2: " << endl;

    for(int i=0; i<reso_T; i++)
    {
        
            //cout << "i: " << i << endl;
        //cout << "WIMP_mx: " << WIMP_mx << endl;
        //cout << "Energy_DM(WIMP_mx,779*kms1_to_c): " << Energy_DM(WIMP_mx,779*kms1_to_c) << endl;

            cout << "T[i]: " << T[i] << endl;
            //cout << "T_QF[i]: " << T_QF[i] << endl;
            cout << "RecoilX: " << recoilX[i] << endl;
            //cout << "Factor1: " << Factor1[i] << endl;
         
    }
    
    cout << "WIMP_max_T: " << WIMP_max_T << endl;
    cout << "WIMP_mx_Check: " << WIMP_mx << endl;
    cout << "Velocity_Min_Check: " << Velocity_Min << endl;
    cout << "================================================" << endl;
    if(Option==0 and Model_of_Interaction==0)return (T_QF);
    if(Option==0 and Model_of_Interaction>0)return (T);
    if(Option==1 and Conventional_or_not==0)return Factor1;
    if(Option==1 and Conventional_or_not==1)return recoilX;
    if(Option==1 and Conventional_or_not==2)return Factor1;
}
//=======Scaling_Factor_Get
double *cpkkd_calculation_Scaling_Factor(double Sigma_SI,int Option, double *T, double *recoilX, double *RE_DATA, double *RE_Rate, double *RE_DATA_Err, double *RE_Rate_Err, double WIMP_mx=10)//Y
{
    cout << "=====================cpkkd_calculation_Scaling_Factor===================" << endl;
    cout << "WIMP_mx: " <<  WIMP_mx << endl;
    const int Data_element_N =257;
    static double VALUE_1[2+Data_element_N];
    //===============================================================
    double sig_E; double dEx;
    static double Factor[Data_element_N];
    static double Factor1[dm_spec_resolution];

    //===============================================================
    double pEx = T[0];
    //Find out the Full energy spectrum with the resolution of the detector based on the detector
    for(int jjj=0 ; jjj<Data_element_N ; jjj++)
    {
            Factor[jjj]=0;
            for(int j=0 ; j<dm_spec_resolution ; j++)
            {
                //Find the theoretical line of the rate of WIMP
                if (recoilX[j]>0)
                {
                    sig_E = a0_dE + a1_dE*sqrt(T[j]);
                    
                    if(j==0) { dEx = T[0]; } else { dEx = T[j] - pEx; }
                    
                    if( (recoilX[j]>0)&&(T[j]>=0.7e-3) )
                    {
                        //Factor[jjj] = Factor[jjj] + recoilX[j]*exp(-pow((T[j]-RE_DATA[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                        double Accumulated = recoilX[j]*exp(-pow((T[j]-RE_DATA[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                        if(Accumulated>0)Factor[jjj] = Factor[jjj] + Accumulated;

                    }
                    pEx = T[j];
                }
            }
    }
    /*
    for(int jjj=0 ; jjj<Data_element_N ; jjj++)
    {
        cout << "RE_DATA[jjj]: " << RE_DATA[jjj] << endl;
        cout << "RE_Rate[jjj]: " << RE_Rate[jjj] << endl;
        cout << "RE_Rate_Err[iii]: " << RE_Rate_Err[jjj] << endl;
        cout << "Factor[jjj]: " << Factor[jjj] << endl;
    }
     */
    
    //Find out the scaling factor
    //===============================================================
    double Initial_Factor; double Ratio; double bin_Lowest=0;
    Initial_Factor= Factor[0]/(RE_Rate[0]+RE_Rate_Err[0]);
    //Initial_Factor= Factor[2]/(RE_Rate[2]+RE_Rate_Err[2]);
    //bin_Lowest = 2;
    
    /*
    for(int iii=0 ; iii<50 ; iii++)
    {
        double Ratio1 = (Factor[iii]/(RE_Rate[iii]+RE_Rate_Err[iii]));
        cout << "Ratio1: " << Ratio1 << endl;
        if(iii>0 && abs(Ratio1)>Initial_Factor && abs(1/Ratio1)!=0)
        {
            
            Initial_Factor= Ratio1;
            bin_Lowest = iii;
        }
    }
    */
    
    cout << "1/Ratio: " << 1/Initial_Factor << endl;
    
    
    //===============================================================
    /*
    for(int jjj=0 ; jjj<dm_spec_resolution ; jjj++)
    {
        Factor2[jjj]=0;
        Factor2[jjj] = Factor1[jjj]*(1/Initial_Factor);
    }
     */
    //cout << "Factor1[jjj]: " << Factor1[int(bin_Lowest)] << endl;

    pEx = T[0];
    //Find out the Full energy spectrum with the resolution of the detector based on the theory
    for(int jjj=0 ; jjj<dm_spec_resolution ; jjj++)
    {
        Factor1[jjj]=0;
        for(int j=0 ; j<dm_spec_resolution ; j++)
        {
            //Find the theoretical line of the rate of WIMP
            if (recoilX[j]>0)
            {
                sig_E = a0_dE + a1_dE*sqrt(T[j]);
                
                if(j==0) { dEx = T[0]; } else { dEx = T[j] - pEx; }
                
                if( (recoilX[j]>0) )
                {
                    double Accumulated = recoilX[j]*exp(-pow((T[j]-T[jjj]),2)/(2*pow(sig_E,2)))*dEx/(sig_E*sqrt(2*PI));
                    if(Accumulated>0)Factor1[jjj] = Factor1[jjj] + Accumulated;
                }
                pEx = T[j];
            }
            
        }
    }

    /*
    for(int i=0; i<reso_T; i++)
    {
            cout << "RecoilX: " << recoilX[i] << endl;
            cout << "Factor: " << Factor[i] << endl;
    }*/

    
    for(int jjj=0 ; jjj<dm_spec_resolution ; jjj++)
    {
        /*
        cout << "T[jjj]" << T[jjj] << endl;
        cout << "Factor1: " << Factor1[jjj]*(1/Initial_Factor) << endl;
         */
    }
     

    
    for(int jjj=0 ; jjj<2+dm_spec_resolution; jjj++)
    {
        VALUE_1[jjj+2] = Factor1[jjj]*(1/Initial_Factor);
    }
    VALUE_1[0] = (1/Initial_Factor);
    VALUE_1[1] = (bin_Lowest);
    cout << " (1/Initial_Factor): " <<  VALUE_1[0] << endl;
    cout << " (bin_Lowest): " <<  VALUE_1[1] << endl;

    cout << "1/Initial_Factor * Sigma_SI: " << (1/Initial_Factor) * Sigma_SI << endl;
    return( VALUE_1 );
    
    //===============================================================
}

