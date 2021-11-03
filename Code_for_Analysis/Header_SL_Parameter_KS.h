const int DataBin = 255;

#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_06GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_07GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_08GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_09GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_10GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_12GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_20GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_25GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_30GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_40GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_50GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_60GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_70GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_80GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_0_90GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_1_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_1_3GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_1_5GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_2_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_2_5GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_3_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_3_5GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_4_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_5_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_10_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/DataJune_ShortRange/DataJune_20_0GeV_c1_dcs.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_LR_DCS_June.h"
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Chi-e/AM_DATA/Header_LR_DCS_Dec.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <string>
#include <filesystem>
namespace fs = std::filesystem;

const double GeV = 1.0;
const double eV     = 1.0e-9 * GeV;

//These parameters are for calculating the sigma_e(Unit transformation!)
const double hbar_eV_s = 6.582119569E-16;//(eV*s)
const double Light_Speed = 2.99792458E8;//(m/s)
const double GeV_for_transform = 1E9;//(eV)
const double hbarC_divide_GeV_m  = hbar_eV_s*Light_Speed/(GeV_for_transform);//(m)
const double hbarC_divide_GeV_cm = 1e2*hbarC_divide_GeV_m;//(cm)
const double hbarC_divide_GeV_cm_Square = pow(hbarC_divide_GeV_cm,2);//(cm^2)

//Other parameters that should be used
const double avogadro_number = 6.02214129e+23;
const double density = 0.3; //GeV cm^{-3} for (DM)
const double Averaged_Velocity=232.;//(km/s)
const double Electron_Mass_MeV=0.511;//MeV
const double Me_times_alpha=3.7;//keV/c^2

const double V_to_C = 1e3/3e8;//km/s to beta
const double MeV_to_GeV = 1E-3;//
const double GeV_to_keV = 1E+6;//

double Final_Energy(int Times, double Initial_Energy)
{
    double A=Initial_Energy; double B;
    for(int kkk=0; kkk<Times ; kkk++)
    {
        B = A - (A*0.5);
        A = B;
    }
    return A;
}

double RMe(double mx)//mx(GeV/c^2),Reduce_Mass_e(RMe)
{
    //cout << "RMe: " << 1e-3*Electron_Mass_MeV*1e-3*1000*mx/((1e-3*Electron_Mass_MeV)+1e-3*1000*mx) << endl;
    return 1e-3*Electron_Mass_MeV*mx/((1e-3*Electron_Mass_MeV)+mx);//GeV/c^2
}
double max_recoil_A_for_ER_keV(double velocity, double mx)//mx(GeV/c^2),Velocity(km/s)
{
  double max_recoil_A_0 = 0.5*mx*1e6*(velocity*V_to_C)*(velocity*V_to_C); //(keV)
    cout << "max_recoil_A_0: " << max_recoil_A_0 << endl;
  return max_recoil_A_0;
}
double c_1(double CS, double mx)//Cross-Section(CS) to c1
{
    return sqrt( ( CS*PI ) / ( RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square )  );
}
double d_1(double CS, double mx)//Cross-Section(CS) to d1
{
    return sqrt( (CS*PI*pow(Me_times_alpha,4)) /( RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square*1e24 ) );
}
double electron_number(double A = AGe)
{
    return (1000.0*avogadro_number)/(A);//kg^{-1}
}

/*
double dsigma_dT_keV_ER(double CS, double Month, double T, double mx, double A)//dsigma_dT_keV for Electron Recoil, Month1: June, Month2: Dec
{
    double rate_factor = pow(d_1(CS,mx),2)*1E-18*1E+3*86400*3E+10*(density*electron_number(A))/(mx);
    if(Month==1)return rate_factor*LongRangeDcs_June(T*1E+3 , mx )*1e-29;
    if(Month==2)return rate_factor*LongRangeDcs_Dec(T*1E+3 , mx )*1e-29;
}
*/
double Max_Recoil_A_keV_ER(double V, double mx)//Electron Recoil Max recoil energy, mx(GeV/c^2), velocity(V)
{
    double Beta = (V*1e3)/(3E+8);//Beta(c)
    return 0.5*mx*(1e6)*(Beta)*(Beta);//Energy of DM basically
    //return 0.5*RMe(mx)*(1e6)*(Beta)*(Beta);//Energy of DM basically

}
double Max_Recoil_A_keV_ER_Free(double V, double mx)//Electron Recoil Max recoil energy, mx(GeV/c^2), velocity(V)
{
    double Beta = (V*1e3)/(3E+8);//Beta(c)
    return 0.5*0.511*(1e3)*(Beta)*(Beta);//Energy of DM basically
    //return 0.5*RMe(mx)*(1e6)*(Beta)*(Beta);//Energy of DM basically

}
/*
double
{
    TF1 *fa2 = new TF1("fa2","dsigma_dT_keV_ER([0],x,[1])",95.0*1E-3,95.0*1E-3+5);
    fa2->SetParameter(0,2);
    fa2->SetParameter(1,mx);
}
*/

double CS_Try(double c_1, double mx)
{
    return c_1*c_1*RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square/(PI);
}
double DS_Try(double d_1, double mx)
{
    return d_1*d_1*RMe(mx)*RMe(mx)*hbarC_divide_GeV_cm_Square*1e24/(PI*pow(Me_times_alpha,4));
}
//================From paper 1905.06348v2===============//
/*
double Max_Recoil_A_keV_ER_Free(double mx, double V_Initial)//Electron Recoil Free electron, Max recoil energy, mx(GeV/c^2), velocity(V)
{
    //cout << "V_Initial: " << V_Initial << endl;
    double V_Final = (V_Initial)*(2*mx/(mx+Electron_Mass_MeV*MeV_to_GeV));
    //cout << "V_Final" << V_Final << endl;
    double Beta    = (V_Final*1e3)/(3E+8);//Beta(c)
    //cout << "Beta: " << Beta << endl;
    //cout << "0.5*Electron_Mass_MeV*1e3*(Beta)*(Beta): " << 0.5*Electron_Mass_MeV*1e3*(Beta)*(Beta) << endl;
    return 0.5*Electron_Mass_MeV*1e3*(Beta)*(Beta);//Energy of DM basically
}
 */


const double MeV       = 1.0e-3 * 1;//MeV to GeV
const double mElectron = 0.5109989461 * MeV;//GeV
const double aEM       = 1.0 / 137.035999139;
const double qRef      = aEM*mElectron;//GeV/c^2

double F_DM(int Option, double q, double mMediator)//
{
    //Contact interactions
        if(Option==0) return 1.0;
    //General dark photon
        else if(Option==1)    return (qRef*qRef+mMediator*mMediator)/(q*q+mMediator*mMediator);
    //Long range interaction
        else if (Option==2)    return qRef*qRef/q/q;
    //Electric dipole interaction
        else if (Option==3)        return qRef/q;
    //Error
        else
        {
            cout << "Something wrong from F_DM" << endl;
            return 0;
        }
}
//=====================================================//
double v_min_DM(double Ee, double q, double Mx, double v_int)//All in GeV (Not totally right)
{
    double Delta_Ee = Ee;// eV (Free-electron Energy + Binding Energy)
    double v_min_beta = (Delta_Ee/(q)) + (q/(2*Mx));//All in GeV
    double v_min      = v_min_beta*(3e8/1e3);
    //cout << "v_min: " << v_min << endl;
    double Bool_for_truncate = 0;
    if(v_int>v_min) Bool_for_truncate=1;
    else Bool_for_truncate=0;
    return Bool_for_truncate;
}
//=====================================================//


double dE_dX_Crystal(double Cross_Section, double mx, double velocity)//velocity(km/s)
{
    const int Ei_Number=500;const int qi_Number=900;
    double v_beta  = velocity*1e3/3e8;
    //static double form_factor_table[qi_Number][Ei_Number];
    vector<double> aux_list;
    vector<vector<double>> form_factor_table(900, vector<double>(500, 0.0));
    vector<double> Ee_List(500,0);vector<double> q_List(900,0);
    //Get the form factor
    std::ifstream input("Test.txt");//Input the auxiliary file
    int Element_Number=-1;double data;
    while(Element_Number<Ei_Number*qi_Number)//
    {
        Element_Number = Element_Number + 1;
        input >> data;
        aux_list.push_back(data);
    }
    //All the necessary parameters
    double prefactor = 2.0 * eV;
    double wk        = 2.0 / 137.0;
    double dE        = 0.1 * eV;
    unsigned int i = 0;

    const double MeV       = 1.0e-3 * 1;//MeV to GeV
    const double mElectron = 0.5109989461 * MeV;//GeV
    const double aEM       = 1.0 / 137.035999139;
    const double dq        = aEM * mElectron;

    TH2F   *HIST_q_E = new TH2F("HIST_q_E","HIST_q_E",500,0,50,900,0,9);
    HIST_q_E->GetZaxis()->SetRangeUser(1e-3,100);
    HIST_q_E->GetYaxis()->SetRangeUser(0,8);

    for(int Ei = 0; Ei < 500; Ei++)
        for(int qi = 0; qi < 900; qi++)
        {
            form_factor_table[qi][Ei] = prefactor * (qi + 1) / dE * wk / 4.0 * aux_list[i++];
            double Ee = (Ei + 1) * 1e-1;//The energy of electron recoil(eV)
            double Unit_Y = 0.02 * (qi + 1) ;//?(Alpha*Me)
            double q  = (Unit_Y+0.5) * dq;// momentum transfer(GeV/c)
            //cout << "q: " << q << endl;
            Ee_List[Ei] = (Ee+0.5)*1e-9;//GeV
            q_List[qi]  = q ;
            /*
            cout << "Ei*1e-1: " << Ei*1e-1 << endl;
            cout << "qi*1e-2: " << 2*qi*1e-2 << endl;
            cout << "form_factor_table: " << form_factor_table[qi][Ei]*form_factor_table[qi][Ei] << endl;
             */
            HIST_q_E->Fill(Ee,Unit_Y,form_factor_table[qi][Ei]*form_factor_table[qi][Ei]);//
        }
    
    
    double reduce_mass_Si_N = (mx)*(1e-3*unified_atomic_mass_MeV*28)/((mx)+(1e-3*unified_atomic_mass_MeV*28));// (GeV/c^2)
    double reduce_mass_Si_e = (mx)*(0.511*1e-3)/((mx)+(0.511*1e-3)); // (GeV/c^2)
    
    double v_rel                 = 1.;
    double q_max                 = 2*reduce_mass_Si_N*(velocity*1e3/3e8);// (GeV/c)
    double Max_Energy_e          = 0.5*(reduce_mass_Si_N)*(v_beta)*(v_beta)*1e9;//eV
    double Max_Recoil_Energy     = 0.5*(mx)*1e9*(v_beta)*(v_beta);//GeV
    //cout << "q_max: " << q_max << endl;
    //double m_cell = 4*72./ 6.02e23;    //each unit cell has 4 Fe and 4 O
    //double n_cell = 2.7/m_cell;
    //double n_cell = 4.9e+22;

    double Max_beta = 784.*1e3/3e8;
    //========================Prefactor============================
    //Parameter: Si_Number_density
    //The number density of Silicon: 5e22(atoms/cm^3)
    //cout << "Prefactor: " << Prefactor << endl;
    //=======================The rest of the integration=======================
    const double Ethr = 1.1*eV;
    double sum_L=0.0;
    double DM_Max_Energy = 0.5*mx*Max_beta*Max_beta;//GeV
    int   Edm_Bin_Number = 500;
    double dE_DM         = DM_Max_Energy/Edm_Bin_Number;
    double dE_dX         = 0;
    
    //for(int i=Ethr/dE;i<Ei_Number;i++)
    /*   // Energy as the integrator
    for(int kkk=0; kkk<Edm_Bin_Number; kkk++)//Energy of DM
    {
        //cout << "kkk: " << kkk << endl;
        double E_DM_now = dE_DM*kkk;//GeV
        double V_DM_now = sqrt((E_DM_now*2)/(mx))*(3e8/1e3);//km/s
        double v_beta   = sqrt((E_DM_now*2)/(mx));
        //double v_beta   = 220.*1e3/3e8;

        for(int i=0;i<Ei_Number;i++)//Energy of electrons
        {
            double V_min_T = sqrt((2*Ee_List[i])/(mx))*(3e8/1e3);
            for(int qi=0;qi<qi_Number;qi++)//Momentum transfer
            {
                if( V_DM_now > v_min_DM(Ee_List[i],q_List[qi],mx,0) and V_DM_now >V_min_T and E_DM_now>1.1*1e-9 and q_List[qi]<q_max)
                {
                    if( Ee_List[i]<DM_Max_Energy )
                    {
                    double Prefactor = (n_cell*Cross_Section*aEM*mElectron*mElectron)/(reduce_mass_Si_e*reduce_mass_Si_e*v_beta*v_rel);
                    double sum_L = ( (dE*Ee_List[i])* ( 0.02*dq*(1.0/(q_List[qi]*q_List[qi])) )*1*form_factor_table[qi][i] );
                    dE_dX = dE_dX + Prefactor*sum_L;
                    }
                }
            }
        }
    }
     
     
    cout << "dE_dX " << (dE_dX) << endl;
    cout << "1/dE_dX " << (1./dE_dX) << endl;
    cout << "(1/dE_dX)*dE_DM " << (1./dE_dX)*dE_DM*1e4 << endl;
     */
    
    const int   dv_Bin_Number = 50;
    double dv_DM         = 784.*(V_to_C)/dv_Bin_Number;
    
    const double N_Avo = 6.02e23;            //N_avo
    double rho =2.7;        //g cm-3
    double m_cell = 4*72./N_Avo;    //each unit cell has 4 Fe and 4 O
    double n_cell = (rho/m_cell);
     
    //double n_cell = 6e21/2;
    cout << "mx: " << mx << endl;
    cout << "n_cell:" << n_cell << endl;
    double sum_array[dv_Bin_Number];
    
    for(int lll=0; lll<dv_Bin_Number; lll++)
    {
        double sum           = 0;
        double V_DM_now = dv_DM*(lll+0.5);
        
        for(int i=0;i<Ei_Number;i++)//Energy of electrons
        {
            for(int qi=0;qi<qi_Number;qi++)//Momentum transfer
            {
                if( V_DM_now > v_min_DM(Ee_List[i],q_List[qi],mx,0)*(V_to_C) and V_DM_now>sqrt(2*Ee_List[i]/mx) and V_DM_now>sqrt(2*(1.1*1e-9)/mx) and q_List[qi]<2*reduce_mass_Si_N*(V_DM_now) )
                {
                    double Prefactor = (n_cell*aEM*mElectron*mElectron)/(reduce_mass_Si_e*reduce_mass_Si_e*v_rel);
                    double sum_L = ( (dE*Ee_List[i])* ( 0.02*dq*(1.0/(q_List[qi]*q_List[qi])) )*1*form_factor_table[qi][i] );
                    if(sum_L!=0 and Prefactor!=0)
                    {
                        sum = sum + (Prefactor*sum_L);
                    }
                }
            }
        }
        sum_array[lll] = sum;
    }
    for(int lll=0; lll<dv_Bin_Number; lll++)
    {
        double V_DM_now = dv_DM*(lll+0.5);
        if(sum_array[lll]>0)dE_dX = dE_dX + mx*pow(V_DM_now,2)*dv_DM*1/sum_array[lll];
    }
    
    cout << "sigma " << (dE_dX)/2e5 << endl;

    //cout << "Energy_Loss(keV): " << (dEdX_M)*(dEdX_M)/(2*mx)*1e5*1e6 << endl;//keV
    //cout << "Energy_Loss(keV): " << (dEdX_M)*(dEdX_M)/(2*mx)*2e5*1e6 << endl;//keV
    //cout << "Energy_Loss(keV): " << (dEdX_M)*2e5*1e6 << endl;//keV

    cout << "Energy_DM(keV): " << 0.5*(mx*1e6)*(v_beta)*(v_beta);//keV
    //Check the form facotr
    
    TCanvas *c3 = new TCanvas("c3");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    HIST_q_E->Draw("colz");
    c3->SetLogz();
    c3->Print("Check_qE.pdf");
    
    return 0;
}
//y=A*x+B
double Slope_A(double x1, double y1, double x2, double y2){return (y2-y1)/(x2-x1);}
double Constant_B(double slope, double x1, double y1){return y1-x1*slope;}
double Expected_Y(double x, double A, double B){return A*x+B;}
//
double fdsigma_dT_ER(string mx, double v_int, double T)//mx(GeV/c^2),v(c), dsigma_dT, T(KeV)
{//Vecloty-dependent dsigma_dT
    
       static string Pre_mx="";
    double v = v_int*(1e3/3e8)*1e3;

    //if(T>0.3 and v_int>300)cout << "T: " << T << endl;
    //if(T>0.3 and v_int>300)cout << "v_int: " << v_int << endl;
    //if(T>0.3 and v_int>300)cout << "v: " << v << endl;
    
       const int DM_Beta_N=10;
       string DM_Beta[DM_Beta_N]={"1P000","1P167","1P333","1P500","1P667","1P833","2P000","2P167","2P333","2P500"};
        //=================Check the cross sections used===========================
       int Now_File=0;
       double DM_Beta_Exact_double[DM_Beta_N+1]={0.000,1.000,1.167,1.333,1.500,1.667,1.833,2.000,2.167,2.333,2.500};
       for(int N=0; N<DM_Beta_N; N++){if(v>DM_Beta_Exact_double[N] and v<DM_Beta_Exact_double[N+1]) Now_File=N;}
       if(v>DM_Beta_Exact_double[9])Now_File=9;
        //====================================================================
       static vector<vector<double>> Energy_Transfer_array(DM_Beta_N);
       static vector<vector<double>> dsigma_dT_array(DM_Beta_N);

       static vector<vector<double>> Energy_Transfer_array_Hist(DM_Beta_N);
       static vector<vector<double>> dsigma_dT_array_Hist(DM_Beta_N);

       static vector<vector<double>> Slope(DM_Beta_N);static vector<vector<double>> Constant(DM_Beta_N);
       static vector<TH1F*>      Hist_DCS_array;
       static vector<TGraph*> TG_Hist_DCS_array;

       int Number_Exe=10;
    
       if(Pre_mx!=mx)
       {
           Pre_mx = mx;
           //fill(Energy_Transfer_array.begin(), Energy_Transfer_array.end(), 0);
           //fill(dsigma_dT_array.begin(), dsigma_dT_array.end(), 0);
           //Make sure that there is no element in.
           //Make sure that there is no element in.
           for(int N=0; N<Number_Exe; N++)
           {
                   string filename(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%sGeV/%s.txt",mx.c_str(),DM_Beta[N].c_str()));
                   cout << "filename: " << filename << endl;
                   vector<char> bytes;
                   vector<char> Energy_Transfer;
                   vector<char> dsigma_dT;
                   string Power="E";
                   int counter_E=0;int counter_sigma=0;int Counter=0;
                   string Energy_Transfer_String="";
                   string dsigma_dT_String="";
                   char byte = 0;
                   ifstream input_file(filename);
                   if (!input_file.is_open())
                   {
                       cerr << "Could not open the file - '"
                            << filename << "'" << endl;
                       return EXIT_FAILURE;
                   }
                   while (input_file.get(byte))
                   {
                       bytes.push_back(byte);
                       if(byte==' ')continue;
                       //^<^Start for energy transfer^<^
                       if(byte==(char)123){counter_E=1; continue;}//Start counting for energy transfer^<^, (char)123 = "{"
                       if(counter_E!=0 and byte!=(char)123){Energy_Transfer_String = Energy_Transfer_String + byte;}//Counting for energy transfer^0^
                       if(counter_E!=0 and byte==','){//End counting for energy transfer^_^
                           Counter = Counter + 1;
                           Energy_Transfer_array[N].push_back(stod(Energy_Transfer_String));
                           Energy_Transfer_String="",counter_E=0;counter_sigma=1;continue;
                       }
                       
                       //^<^Start for differential cross sections^<^
                       if(counter_sigma!=0 and byte=='D'){dsigma_dT_String = dsigma_dT_String + Power;}//Replace E with D^0^
                       else if(counter_sigma!=0 and byte==(char)125){//Replace E with D^0^, (char)125 = "}"
                           dsigma_dT_array[N].push_back(stod(dsigma_dT_String));
                           dsigma_dT_String="";counter_sigma=0;
                       }
                       else if(counter_sigma!=0){dsigma_dT_String = dsigma_dT_String + byte;}//Counting for dsigma_dT
                        
                   }
               
                    input_file.close();
                    
                   //^^TGraph Plot
                   Float_t TG_Lower_eV[DM_Beta_N][Counter];
                   Float_t TG_dsigma_dT[DM_Beta_N][Counter];
                   Float_t TG_Slope[DM_Beta_N][Counter];
                   Float_t TG_Constant[DM_Beta_N][Counter];


                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       TG_Lower_eV[N][kkk]  = Energy_Transfer_array[N][kkk];//TGraph
                       TG_dsigma_dT[N][kkk] = TMath::Log10(dsigma_dT_array[N][kkk]);//TGraaph
                   }//For extrapolating the cross sections between 0 to 80eV
                    
                   for(int kkk=0; kkk<Counter-1; kkk++)//Add the slope from the second interval ( from 80eV ) into the bin from the second element
                   {
                       TG_Slope[N][kkk+1]    = Slope_A(TG_Lower_eV[N][kkk],TG_dsigma_dT[N][kkk],TG_Lower_eV[N][kkk+1],TG_dsigma_dT[N][kkk+1]);
                       TG_Constant[N][kkk+1] = Constant_B(TG_Slope[N][kkk+1],TG_Lower_eV[N][kkk],TG_dsigma_dT[N][kkk]);
                       //cout << "Check: " << TG_Constant[N][kkk+1]+TG_Slope[N][kkk+1]*TG_Lower_eV[N][kkk+1] << "==?" << TG_dsigma_dT[N][kkk+1] << endl;
                       //cout << "hist_DFCS->GetBinContent(kkk+1): " << hist_DFCS->GetBinContent(kkk+1) << endl;
                   }
                   //Aussume the first line between 0 to 80eV is the same as the second one
                   TG_Slope[N][0]    = Slope_A(TG_Lower_eV[N][0],TG_dsigma_dT[N][0],TG_Lower_eV[N][1],TG_dsigma_dT[N][1]);
                   TG_Constant[N][0] = Constant_B(TG_Slope[N][0],TG_Lower_eV[N][0],TG_dsigma_dT[N][0]);
                   //============================Find the extrapolated line============================
                   Float_t TG_Lower_eV_F[DM_Beta_N][Counter+1];//Add one point (0) to the array
                   Float_t TG_dsigma_dT_F[DM_Beta_N][Counter+1];
                   TG_Lower_eV_F[N][0]=0;TG_dsigma_dT_F[N][0]=TG_Constant[N][0];
                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       TG_Lower_eV_F[N][kkk+1]  = Energy_Transfer_array[N][kkk];//TGraph
                       TG_dsigma_dT_F[N][kkk+1] = TMath::Log10(dsigma_dT_array[N][kkk]);//TGraaph
                   }//For extrapolating the cross sections between 0 to 80eV

                   TGraph* gr = new TGraph(Counter+1,TG_Lower_eV_F[N],TG_dsigma_dT_F[N]);
                   TG_Hist_DCS_array.push_back(gr);
                   //============================Final values==========================================
                   for(int kkk=0; kkk<Counter+1; kkk++)//Put all of them into another array for histogram
                   {
                       Energy_Transfer_array_Hist[N].push_back(TG_Lower_eV_F[N][kkk]);
                       dsigma_dT_array_Hist[N].push_back(TG_dsigma_dT_F[N][kkk]);
                   }
                   //==================================================================================
                   Float_t DIff_Lower_eV_Hist[DM_Beta_N][Counter];
                   Double_t Diff_Total[DM_Beta_N];
                   const double Bin_Number=1e+4;
                   Double_t Bin_Number_interval[DM_Beta_N][Counter];
                   Double_t Width_Per_Bin[DM_Beta_N][Counter];

                   double Bin_Number_total;
                   double checkkk;
                   for(int kkk=0; kkk <Counter; kkk++)
                   {
                       DIff_Lower_eV_Hist[N][kkk]=Energy_Transfer_array_Hist[N][kkk+1]-Energy_Transfer_array_Hist[N][kkk];//Hist
                       Diff_Total[N] = Diff_Total[N] + DIff_Lower_eV_Hist[N][kkk];
                   }
                   for(int kkk=0; kkk <Counter; kkk++)
                   {
                       checkkk = checkkk + DIff_Lower_eV_Hist[N][kkk]/Diff_Total[N];
                       Bin_Number_interval[N][kkk]=Bin_Number*(DIff_Lower_eV_Hist[N][kkk]/Diff_Total[N]);
                       Width_Per_Bin[N][kkk]=(DIff_Lower_eV_Hist[N][kkk])/(Bin_Number_interval[N][kkk]);
                       Bin_Number_total = Bin_Number_total + Bin_Number_interval[N][kkk];
                   }
                   cout.precision(10);
                   for(int kkk=0; kkk <Counter; kkk++)
                   {
                      // cout << "Max: "      << Energy_Transfer_array_Hist[N][kkk]+(Bin_Number_interval[N][kkk])*Width_Per_Bin[N][kkk] << endl;

                   }
                   //cout << "Diff_Total: " << Diff_Total[N] << endl;
                   vector<vector<double>> Lower_eV_array_for_Hist(DM_Beta_N);
                   vector<vector<double>> Lower_eV_center_array_for_Hist(DM_Beta_N);
                   vector<vector<double>> dsigma_dT_for_Hist(DM_Beta_N);

                   int test_bin_Number=0;
                  //==================================================================================
                   for(int kkk=0; kkk< Counter; kkk++)//Tell about the interval designated
                   {
                       /*
                       cout << "=========================================================" << endl;
                       cout << "Bin_Number_interval[N][kkk]: " << Bin_Number_interval[N][kkk] << endl;
                       cout << "kkk: " << kkk << endl;
                       cout << "Energy_Transfer_array_Hist[N][kkk]: " << Energy_Transfer_array_Hist[N][kkk] << endl;
                        
                       cout << "TG_Slope[N][kkk]: " << TG_Slope[N][kkk] << endl;cout << "TG_Constant[N][kkk]: " << TG_Constant[N][kkk] << endl;
                        */
                       //cout << "=========================================================" << endl;
                            for(int jjj=0; jjj<Bin_Number_interval[N][kkk]; jjj++)
                            {
                                Lower_eV_array_for_Hist[N].push_back(Energy_Transfer_array_Hist[N][kkk]+(jjj)*Width_Per_Bin[N][kkk]);
                                Lower_eV_center_array_for_Hist[N].push_back(Energy_Transfer_array_Hist[N][kkk]+(jjj+0.5)*Width_Per_Bin[N][kkk]);
                                dsigma_dT_for_Hist[N].push_back(TG_Slope[N][kkk]*Lower_eV_center_array_for_Hist[N][test_bin_Number]+TG_Constant[N][kkk]);
                                /*
                                cout << "Lower_eV_center_array_for_Hist[N]: " << Lower_eV_center_array_for_Hist[N][test_bin_Number] << endl;
                                cout << "dsigma_dT_for_Hist[N]: " << dsigma_dT_for_Hist[N][test_bin_Number] << endl;
                                cout << "TMath::Power(10,dsigma_dT_for_Hist[N][test_bin_Number]): " << TMath::Power(10,dsigma_dT_for_Hist[N][test_bin_Number]) << endl;
                                 */
                                test_bin_Number = test_bin_Number + 1;
                            }
                       //cout << "=========================================================" << endl;
                   }
               
                    
                   Float_t    Lower_eV_Hist_Final[DM_Beta_N][test_bin_Number];
                   Float_t    dsigma_dT_Hist_Final[DM_Beta_N][test_bin_Number];

               for(int kkk=0; kkk<test_bin_Number; kkk++)
               //for(int kkk=0; kkk<3; kkk++)
               {Lower_eV_Hist_Final[N][kkk]=Lower_eV_array_for_Hist[N][kkk];}//Hist
                TH1F* hist_DFCS = new TH1F("hist","",test_bin_Number-1, Lower_eV_Hist_Final[N]);//Differential Cross-section
               for(int kkk=0; kkk<test_bin_Number; kkk++)
               //for(int kkk=0; kkk<3; kkk++)
               {hist_DFCS->SetBinContent(kkk+1, TMath::Power(10,dsigma_dT_for_Hist[N][kkk]));}//Hist
               
                   Hist_DCS_array.push_back(hist_DFCS);

                   /*
                    //^^Produce the Hist
                   Float_t    Lower_eV[DM_Beta_N][];
                   Lower_eV[N][0]=0;
                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       Lower_eV[N][kkk+1]=Energy_Transfer_array[N][kkk];//Hist
                       //cout << "Energy_transfer: " << Lower_eV[N][kkk+1] << endl;
                   }
                   TH1F* hist_DFCS = new TH1F("hist","",Counter, Lower_eV[N]);//Differential Cross-section
                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       hist_DFCS->SetBinContent(kkk+1, TMath::Log10(dsigma_dT_array[N][kkk]));//Hist
                       //cout << "dsigma_dT_array[kkk+1]: " << dsigma_dT_array[N][kkk] << endl;
                   }
                   Hist_DCS_array.push_back(hist_DFCS);
                    */

           }
       }
        Int_t binx = Hist_DCS_array[Now_File]->GetXaxis()->FindBin(T*1e3);
        //cout << "binx: " << binx << endl;
        //if(T>0.3 and v_int>300)cout << "binx: " << binx << endl;
        double DCS = Hist_DCS_array[Now_File]->GetBinContent(binx);
        //cout << "DCS: " << DCS << endl;
        //if(T>0.3 and v_int>300)cout << "DCS: " << DCS << endl;
        /*
         
         //gr->Draw("ALP");
         //c3->Print("Check_ECS.pdf");

        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
    
        c3->SetLogy();
        hist_DFCS->GetXaxis()->SetRangeUser(0,435);
        hist_DFCS->GetYaxis()->SetRangeUser(1e-8,1);
         
        hist_DFCS->Draw("HIST");
         */
    
        //Give out the root file to check the content of the DCS
        /*
        char fout_name[100];
        sprintf(fout_name,Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%sGeV/DCS.root",mx.c_str()));
        TFile *fout=new TFile(fout_name,"recreate");
        for(int N=0; N<Number_Exe; N++){   Hist_DCS_array[N]->Write(DM_Beta[N].c_str());}
        //for(int N=0; N<DM_Beta_N; N++){TG_Hist_DCS_array[N]->Write(DM_Beta[N].c_str());}

        fout->Close();
        */
    
         return DCS;
}

static vector<string> FileName;//
static vector<int> Mass_Int;//Real Mass*100 for storage and comparison

int Find_File()
{
    int kkk=0;
    std::string path = "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/SI_c1_XeData_Vel";
    for (const auto & entry : fs::directory_iterator(path))
    {
        if ( std::find(std::begin(FileName), std::end(FileName), entry.path()) != std::end(FileName) ){cout << "OK" ;}
        else{FileName.push_back(entry.path());}
        
        int Start_file=0;
        char byte = 0;
        ifstream input_file(entry.path());
        while(input_file.get(byte))
            {
                if(byte='l')Start_file=1;
                
        }
        kkk = kkk + 1;
    }
    for(int kkk=0; kkk< FileName.size(); kkk++){cout << "FileName: " << FileName[kkk] << endl;}
    return 0;
}

double fdsigma_dT_ER_New(string mx, double v_int, double T)//mx(GeV/c^2),v(c), dsigma_dT, T(KeV)
{//Vecloty-dependent dsigma_dT
    
       static string Pre_mx="";
       double v = v_int*(1e3/3e8)*1e3;
       Find_File();
       const int DM_Beta_N=12;
       // string DM_Beta[DM_Beta_N]={"06667","08333","10000","11670","13330","15000","16670","18330","20000","21670","23330","25000"};//Ge 0.5GeV
       string DM_Beta[DM_Beta_N]={"03333","05000","06667","08333","10000","11670","13330","15000","16670","18330","20000","21670"};//Xe 0.5GeV
       //const int DM_Beta_N=7;
       //string DM_Beta[DM_Beta_N]={"11670","13330","15000","16670","18330","20000","21670"};//Xe 0.2GeV
        //=================Check the cross sections used===========================
       int Now_File=0;
       //double DM_Beta_Exact_double[DM_Beta_N+1]={0.000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670,2.3330,2.5000};//Ge 0.5GeV
       double DM_Beta_Exact_double[DM_Beta_N+1]={0.000,0.3333,0.5000,0.6667,0.8333,1.0000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};//Xe 0.5GeV
      // double DM_Beta_Exact_double[DM_Beta_N+1]={0.000,1.1670,1.3330,1.5000,1.6670,1.8330,2.0000,2.1670};//Xe 0.2GeV
       for(int N=0; N<DM_Beta_N; N++){if(v>DM_Beta_Exact_double[N] and v<DM_Beta_Exact_double[N+1]) Now_File=N;}
       if(v>DM_Beta_Exact_double[DM_Beta_N-1])Now_File=DM_Beta_N-1;
        //====================================================================
       static vector<vector<double>> Energy_Transfer_array(DM_Beta_N);
       static vector<vector<double>> dsigma_dT_array(DM_Beta_N);

       static vector<vector<double>> Energy_Transfer_array_Hist(DM_Beta_N);
       static vector<vector<double>> dsigma_dT_array_Hist(DM_Beta_N);

       static vector<vector<double>> Slope(DM_Beta_N);static vector<vector<double>> Constant(DM_Beta_N);
       static vector<TH1F*>      Hist_DCS_array;
       static vector<TGraph*> TG_Hist_DCS_array;

       int Number_Exe=12;
    
       if(Pre_mx!=mx)
       {
           Pre_mx = mx;
           //fill(Energy_Transfer_array.begin(), Energy_Transfer_array.end(), 0);
           //fill(dsigma_dT_array.begin(), dsigma_dT_array.end(), 0);
           //Make sure that there is no element in.
           //Make sure that there is no element in.
           for(int N=0; N<Number_Exe; N++)
           {
                   
                   //string filename(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/%sGeV/%s.txt",mx.c_str(),DM_Beta[N].c_str()));
                   //string filename(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/SI_d1_XeData_Vel/d1_%sGeV/DM_%sGeV_%sV.h",mx.c_str(),mx.c_str(),DM_Beta[N].c_str()));
                   string filename(Form("/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/SI_d1_XeData_Vel/Xe_d1_0_500GeV/DM_0_50GeV_%sV.h",DM_Beta[N].c_str()));
                  //cout << "filename: " << filename << endl;
                   vector<char> bytes;
                   vector<char> Energy_Transfer;
                   vector<char> dsigma_dT;
                   string Power="E";
                   int counter_E_start=0;int counter_E=0;
                   int counter_sigma=0;int Counter=0;
                   string Energy_Transfer_String="";
                   string dsigma_dT_String="";
                   char byte = 0;
                   ifstream input_file(filename);
                   if (!input_file.is_open())
                   {
                       cerr << "Could not open the file - '"
                            << filename << "'" << endl;
                       return EXIT_FAILURE;
                   }
                   while (input_file.get(byte))
                   {
                       bytes.push_back(byte);
                       //cout << byte;
                       if(byte==' ')continue;
                       //^<^Start for energy transfer^<^
                       if(byte==(char)123 and counter_E_start==0){counter_E_start=counter_E_start+1; continue;}//Start counting for energy transfer^<^, (char)123 = "{"
                       if(counter_E_start!=0 and byte==(char)123){counter_E=1;continue;}//Start counting in this array
                       if(counter_E==1 and byte!=','){Energy_Transfer_String = Energy_Transfer_String + byte;}//Counting for energy transfer^0^
                       if(counter_E==1 and byte==','){//End counting for energy transfer^_^
                           Counter = Counter + 1;
                           Energy_Transfer_array[N].push_back(stod(Energy_Transfer_String));
                           Energy_Transfer_String="",counter_E=0;counter_sigma=1;continue;
                       }
                       if(counter_sigma==1 and byte!=(char)125){dsigma_dT_String = dsigma_dT_String + byte;}//Counting for energy transfer^0^
                       if(counter_sigma==1 and byte==(char)125){//
                           if(dsigma_dT_String=="0.000000e+00"){Counter = Counter - 1;break;}
                           dsigma_dT_array[N].push_back(stod(dsigma_dT_String));
                           dsigma_dT_String="";counter_sigma=0;
                       }
                       

                   }
               
                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       //cout << "Energy_Transfer_array[N][kkk]: " << Energy_Transfer_array[N][kkk] << endl;
                       //cout << "dsigma_dT_array[N][kkk]: " << dsigma_dT_array[N][kkk] << endl;
                   }
                    input_file.close();
                    
                   //^^TGraph Plot
                   Float_t TG_Lower_eV[DM_Beta_N][Counter];
                   Float_t TG_dsigma_dT[DM_Beta_N][Counter];
                   Float_t TG_Slope[DM_Beta_N][Counter];
                   Float_t TG_Constant[DM_Beta_N][Counter];


                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       TG_Lower_eV[N][kkk]  = Energy_Transfer_array[N][kkk];//TGraph
                       TG_dsigma_dT[N][kkk] = TMath::Log10(dsigma_dT_array[N][kkk]);//TGraaph
                       //cout << "TG_Lower_eV[N][kkk]: " << TG_Lower_eV[N][kkk] << endl;
                       //cout << "TG_dsigma_dT[N][kkk]: " << TG_dsigma_dT[N][kkk] << endl;
                   }//For extrapolating the cross sections between 0 to 80eV
                    
                   for(int kkk=0; kkk<Counter-1; kkk++)//Add the slope from the second interval ( from 80eV ) into the bin from the second element
                   {
                       TG_Slope[N][kkk+1]    = Slope_A(TG_Lower_eV[N][kkk],TG_dsigma_dT[N][kkk],TG_Lower_eV[N][kkk+1],TG_dsigma_dT[N][kkk+1]);
                       TG_Constant[N][kkk+1] = Constant_B(TG_Slope[N][kkk+1],TG_Lower_eV[N][kkk],TG_dsigma_dT[N][kkk]);
                       //cout << "Check: " << TG_Constant[N][kkk+1]+TG_Slope[N][kkk+1]*TG_Lower_eV[N][kkk+1] << "==?" << TG_dsigma_dT[N][kkk+1] << endl;
                       //cout << "hist_DFCS->GetBinContent(kkk+1): " << hist_DFCS->GetBinContent(kkk+1) << endl;
                   }
                   //Aussume the first line between 0 to 80eV is the same as the second one
                   TG_Slope[N][0]    = Slope_A(TG_Lower_eV[N][0],TG_dsigma_dT[N][0],TG_Lower_eV[N][1],TG_dsigma_dT[N][1]);
                   TG_Constant[N][0] = Constant_B(TG_Slope[N][0],TG_Lower_eV[N][0],TG_dsigma_dT[N][0]);
                   //============================Find the extrapolated line============================
                   Float_t TG_Lower_eV_F[DM_Beta_N][Counter+1];//Add one point (0) to the array
                   Float_t TG_dsigma_dT_F[DM_Beta_N][Counter+1];
                   TG_Lower_eV_F[N][0]=0;TG_dsigma_dT_F[N][0]=TG_Constant[N][0];
                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       TG_Lower_eV_F[N][kkk+1]  = Energy_Transfer_array[N][kkk];//TGraph
                       TG_dsigma_dT_F[N][kkk+1] = TMath::Log10(dsigma_dT_array[N][kkk]);//TGraaph
                   }//For extrapolating the cross sections between 0 to 80eV

                   TGraph* gr = new TGraph(Counter+1,TG_Lower_eV_F[N],TG_dsigma_dT_F[N]);
                   TG_Hist_DCS_array.push_back(gr);
                   //============================Final values==========================================
                   for(int kkk=0; kkk<Counter+1; kkk++)//Put all of them into another array for histogram
                   {
                       Energy_Transfer_array_Hist[N].push_back(TG_Lower_eV_F[N][kkk]);
                       dsigma_dT_array_Hist[N].push_back(TG_dsigma_dT_F[N][kkk]);
                   }
                   //==================================================================================
                   Float_t DIff_Lower_eV_Hist[DM_Beta_N][Counter];
                   Double_t Diff_Total[DM_Beta_N];
                   const double Bin_Number=1e+4;
                   Double_t Bin_Number_interval[DM_Beta_N][Counter];
                   Double_t Width_Per_Bin[DM_Beta_N][Counter];

                   double Bin_Number_total;
                   double checkkk;
                   for(int kkk=0; kkk <Counter; kkk++)
                   {
                       DIff_Lower_eV_Hist[N][kkk]=Energy_Transfer_array_Hist[N][kkk+1]-Energy_Transfer_array_Hist[N][kkk];//Hist
                       Diff_Total[N] = Diff_Total[N] + DIff_Lower_eV_Hist[N][kkk];
                   }
                   for(int kkk=0; kkk <Counter; kkk++)
                   {
                       checkkk = checkkk + DIff_Lower_eV_Hist[N][kkk]/Diff_Total[N];
                       Bin_Number_interval[N][kkk]=Bin_Number*(DIff_Lower_eV_Hist[N][kkk]/Diff_Total[N]);
                       Width_Per_Bin[N][kkk]=(DIff_Lower_eV_Hist[N][kkk])/(Bin_Number_interval[N][kkk]);
                       Bin_Number_total = Bin_Number_total + Bin_Number_interval[N][kkk];
                   }
                   cout.precision(10);
                   for(int kkk=0; kkk <Counter; kkk++)
                   {
                      // cout << "Max: "      << Energy_Transfer_array_Hist[N][kkk]+(Bin_Number_interval[N][kkk])*Width_Per_Bin[N][kkk] << endl;

                   }
                   //cout << "Diff_Total: " << Diff_Total[N] << endl;
                   vector<vector<double>> Lower_eV_array_for_Hist(DM_Beta_N);
                   vector<vector<double>> Lower_eV_center_array_for_Hist(DM_Beta_N);
                   vector<vector<double>> dsigma_dT_for_Hist(DM_Beta_N);

                   int test_bin_Number=0;
                  //==================================================================================
                   for(int kkk=0; kkk< Counter; kkk++)//Tell about the interval designated
                   {
                       /*
                       cout << "=========================================================" << endl;
                       cout << "Bin_Number_interval[N][kkk]: " << Bin_Number_interval[N][kkk] << endl;
                       cout << "kkk: " << kkk << endl;
                       cout << "Energy_Transfer_array_Hist[N][kkk]: " << Energy_Transfer_array_Hist[N][kkk] << endl;
                        
                       cout << "TG_Slope[N][kkk]: " << TG_Slope[N][kkk] << endl;cout << "TG_Constant[N][kkk]: " << TG_Constant[N][kkk] << endl;
                        */
                       //cout << "=========================================================" << endl;
                            for(int jjj=0; jjj<Bin_Number_interval[N][kkk]; jjj++)
                            {
                                Lower_eV_array_for_Hist[N].push_back(Energy_Transfer_array_Hist[N][kkk]+(jjj)*Width_Per_Bin[N][kkk]);
                                Lower_eV_center_array_for_Hist[N].push_back(Energy_Transfer_array_Hist[N][kkk]+(jjj+0.5)*Width_Per_Bin[N][kkk]);
                                dsigma_dT_for_Hist[N].push_back(TG_Slope[N][kkk]*Lower_eV_center_array_for_Hist[N][test_bin_Number]+TG_Constant[N][kkk]);
                                /*
                                cout << "Lower_eV_center_array_for_Hist[N]: " << Lower_eV_center_array_for_Hist[N][test_bin_Number] << endl;
                                cout << "dsigma_dT_for_Hist[N]: " << dsigma_dT_for_Hist[N][test_bin_Number] << endl;
                                cout << "TMath::Power(10,dsigma_dT_for_Hist[N][test_bin_Number]): " << TMath::Power(10,dsigma_dT_for_Hist[N][test_bin_Number]) << endl;
                                 */
                                test_bin_Number = test_bin_Number + 1;
                            }
                       //cout << "=========================================================" << endl;
                   }
               
                    
                   Float_t    Lower_eV_Hist_Final[DM_Beta_N][test_bin_Number];
                   Float_t    dsigma_dT_Hist_Final[DM_Beta_N][test_bin_Number];

               for(int kkk=0; kkk<test_bin_Number; kkk++)
               //for(int kkk=0; kkk<3; kkk++)
               {Lower_eV_Hist_Final[N][kkk]=Lower_eV_array_for_Hist[N][kkk];}//Hist
                TH1F* hist_DFCS = new TH1F("hist","",test_bin_Number-1, Lower_eV_Hist_Final[N]);//Differential Cross-section
               for(int kkk=0; kkk<test_bin_Number; kkk++)
               //for(int kkk=0; kkk<3; kkk++)
               {hist_DFCS->SetBinContent(kkk+1, TMath::Power(10,dsigma_dT_for_Hist[N][kkk]));}//Hist
               
                   Hist_DCS_array.push_back(hist_DFCS);

                   /*
                    //^^Produce the Hist
                   Float_t    Lower_eV[DM_Beta_N][];
                   Lower_eV[N][0]=0;
                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       Lower_eV[N][kkk+1]=Energy_Transfer_array[N][kkk];//Hist
                       //cout << "Energy_transfer: " << Lower_eV[N][kkk+1] << endl;
                   }
                   TH1F* hist_DFCS = new TH1F("hist","",Counter, Lower_eV[N]);//Differential Cross-section
                   for(int kkk=0; kkk<Counter; kkk++)
                   {
                       hist_DFCS->SetBinContent(kkk+1, TMath::Log10(dsigma_dT_array[N][kkk]));//Hist
                       //cout << "dsigma_dT_array[kkk+1]: " << dsigma_dT_array[N][kkk] << endl;
                   }
                   Hist_DCS_array.push_back(hist_DFCS);
                    */

           }
       }
        Int_t binx = Hist_DCS_array[Now_File]->GetXaxis()->FindBin(T*1e3);
        //cout << "binx: " << binx << endl;
        //if(T>0.3 and v_int>300)cout << "binx: " << binx << endl;
        double DCS = Hist_DCS_array[Now_File]->GetBinContent(binx);
        //cout << "DCS: " << DCS << endl;
        //if(T>0.3 and v_int>300)cout << "DCS: " << DCS << endl;
    
        /*
         
         //gr->Draw("ALP");
         //c3->Print("Check_ECS.pdf");

        TCanvas *c3 = new TCanvas("c3");
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
    
        c3->SetLogy();
        hist_DFCS->GetXaxis()->SetRangeUser(0,435);
        hist_DFCS->GetYaxis()->SetRangeUser(1e-8,1);
         
        hist_DFCS->Draw("HIST");
         */
    
        //Give out the root file to check the content of the DCS
        /*
        char fout_name[100];
        sprintf(fout_name,"/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/ER_cross_section/SI_d1_XeData_Vel/Xe_d1_0_500GeV/DCS.root");
        TFile *fout=new TFile(fout_name,"recreate");
        for(int N=0; N<Number_Exe; N++){   Hist_DCS_array[N]->Write(DM_Beta[N].c_str());}
        //for(int N=0; N<DM_Beta_N; N++){TG_Hist_DCS_array[N]->Write(DM_Beta[N].c_str());}

        fout->Close();
        */
    
         return DCS;
}

/*
double dE_dX_crystal_Alex()
{
    
}
 */
/*
double From_Factor_Si_e()
{
    const int reso_mx = 150;
    double Ee[reso_mx];
    double P_f1[reso_mx];double P_f2[reso_mx];
    
    double Ee_f1_V, P_f1_V;
    double Ee_f2_V, P_f2_V;

    //=======================
    TGraph *Si_f1 = new TGraph("Si_f1.txt","%lg  %lg");
    TGraph *Si_f2 = new TGraph("Si_f2.txt","%lg  %lg");
    //=======================F1 and F2=======================//
    for(int i=0;i<reso_mx;i++)
    {
      Si_f1->GetPoint(i,Ee_f1_V,P_f1_V);
      Ee[i]   = Ee_f1_V;
      P_f1[i]    = P_f1_V;
      cout << "P_f1[i]: " << P_f1[i] << endl;
    }
    for(int i=0;i<reso_mx;i++)
    {
      Si_f2->GetPoint(i,Ee_f2_V,P_f2_V);
      P_f2[i]    = P_f2_V;;
      cout << "P_f2[i]: " << P_f1[i] << endl;
    }
    //=========================F_total=========================//
    for(int i=0; i<reso_mx;i++)
    {
        cout << "Ee: " << Ee[i]*1e3 << endl;
        cout << "F: "  << F_total(P_f1[i],P_f2[i]) << endl;
    }
    return 0;
}
 */
//Find the dsigma_dT
/*
double dsigma_dT_keV_ERSS(double Sigma_SI, double mx, double V)//ER_Square
{
    //return Sigma_SI*(1./(4*RMe(mx)*RMe(mx)*V*V_to_C*V*V_to_C))*F_DM(mx,V)*F_DM(mx,V);
    return Sigma_SI*(1./(4*RMe(mx)*1e6*RMe(mx)*1e6*V*V_to_C*V*V_to_C));//(cm^2*c^2/keV^2)
}
double dsigma_dT_keV_ER(double Sigma_SI, double mx, double V, double T, double M_DP)//Sigma_SI(cm^2), mx(GeV/c^2), V(km/s), q(momentum from above), M_DP(keV/c^2)[Dark Photon Mass], ER
{
    //return Sigma_SI*(1./(4*RMe(mx)*RMe(mx)*V*V_to_C*V*V_to_C))*F_DM(mx,V)*F_DM(mx,V);
    double Result = 2*Electron_Mass_MeV*1e3*(dsigma_dT_keV_ERSS(Sigma_SI,mx,V)*F_DM(T, M_DP)*F_DM(T, M_DP));
    if( Max_Recoil_A_keV_ER_Free(V,mx) < T)
    {
        Result = 0;
    }
    return Result;//
}


double Total_Sigma_ER(double CS, double V, double mx)//Total Sigma for Electron Recoil, V(km/s)
{
    int reso_T=1000;double T[reso_T];double total_Sigma=0;
     //double WIMP_max_T   = max_recoil_A_for_ER_keV(V,mx); //keV
    double WIMP_max_T   = Max_Recoil_A_keV_ER_Free(V,mx); //keV
    cout << "WIMP_max_T: " << WIMP_max_T << endl;
    //======
    for(int i=0;i<reso_T;i++)
    {
        T[i] = ((double)i+0.5)*((WIMP_max_T)/(double)reso_T); // keV
    }
    //======
    double dEx=0;
    double pEx = T[0];
    for(int i=0;i<reso_T;i++)
    {
        if(i==0) { dEx = T[0]; }
        else { dEx = T[i] - pEx; }
        total_Sigma = total_Sigma + (dsigma_dT_keV_ER(CS, mx, V, T[i], 1e9)*dEx);
        pEx = T[i];
    }
    cout << "total_Sigma: " << total_Sigma << endl;
    return total_Sigma;
}
double F_total(double F1, double F2)
{
    return sqrt(F1*F1+F2*F2);
}
*/
