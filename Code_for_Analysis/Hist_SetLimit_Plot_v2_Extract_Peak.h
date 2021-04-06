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
#include "VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"

/*
#include "VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"
#include "dsigma_dT2.h"
#include "cpkkd_chi_N_v3.h"
#include "cpkkd_calculation_New.h"
 */

//=============================New_Error============================//
float New_Error(float Original_error, float Platform_Error, float fitting_error)
{
    return sqrt( (Original_error*Original_error) + (Platform_Error*Platform_Error) + (fitting_error*fitting_error)  );
}
//=================FInd the resolution of the energy peak=================//
const Double_t* Get_Error(TGraphErrors *cdexdata, double Lower_bound, double Upper_bound, double Parameter0, double Parameter1, double Parameter2)
{
    const Double_t *par_Error_latter;
    TF1 *fitfun=new TF1("fitfun","[0]*TMath::Exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))",Lower_bound,Upper_bound);
    fitfun->SetParameters(0,Parameter0); // Changable constant ==> Find the error band
    fitfun->SetParLimits(1,Parameter1+0.001,Parameter1+0.004);  // fixed mean
    fitfun->SetParLimits(2,Parameter2-0.000001,Parameter2+0.000001);  // fixed sigma
    cdexdata->Fit(fitfun,"","",Lower_bound,Upper_bound);
    par_Error_latter = fitfun->GetParErrors();
    return par_Error_latter;
}
TF1* Draw_Uncertainty(double Lower_bound, double Upper_bound, int Color_Used, double Parameter0, double Parameter1, double Parameter2)
{
    TF1 *Drawfun=new TF1("Drawfun","[0]*TMath::Exp(-((x-[1])*(x-[1]))/(2*[2]*[2]))",Lower_bound,Upper_bound);
    Drawfun->SetLineColor(Color_Used);
    Drawfun->SetParameter(0,Parameter0);  // fixed constant
    Drawfun->SetParameter(1,Parameter1);  // fixed mean
    Drawfun->SetParameter(2,Parameter2);  // fixed sigma
    return Drawfun;
}
float Resolution_of_Signal(float Energy_keV)
{
    double a0_dE = 30.0/1000.0;
    double a1_dE = 50.0/(1000.0*sqrt(10.36));
    return a0_dE + a1_dE*sqrt(Energy_keV);
}
//=================For finding the parameter of Gaussian ( Area/Y/constant )=================//
float Gaussian_Constant(float Area_of_Gaussian, float Uncertainty)
{
    return Area_of_Gaussian/(Uncertainty*sqrt(2*3.14159));
}
float Gaussian_Area(float Constant, float Uncertainty)
{
    return Constant*Uncertainty*sqrt(2*3.14159);
}
float Extract_K_L_shell(float Target_X, vector<float> Constant, vector<float> Mean, vector<float> Error)
{
    float YYY=0;
    for(int kkk=0 ; kkk<10 ; kkk++)
    {
        YYY = YYY + Constant[kkk]*TMath::Exp(-((Target_X-Mean[kkk])*(Target_X-Mean[kkk]))/(2*Error[kkk]*Error[kkk]));
    }
    return YYY;
}
//===========================================================================================//
float Find_Gaussian_Y(float Target_X, float constant, float mean, float error)
{
    return constant*TMath::Exp(-((Target_X-mean)*(Target_X-mean))/(2*error*error));
}
float Find_Gaussian_Y_over_all(float Target_X, vector<float> par)
{
    return par[1]*TMath::Exp(-((Target_X-par[2])*(Target_X-par[2]))/(2*par[3]*par[3]))+par[4]*TMath::Exp(-((Target_X-par[5])*(Target_X-par[5]))/(2*par[6]*par[6]))+par[7]*TMath::Exp(-((Target_X-par[8])*(Target_X-par[8]))/(2*par[9]*par[9]))+par[10]*TMath::Exp(-((Target_X-par[11])*(Target_X-par[11]))/(2*par[12]*par[12]))+par[13]*TMath::Exp(-((Target_X-par[14])*(Target_X-par[14]))/(2*par[15]*par[15]))+par[16]*TMath::Exp(-((Target_X-par[17])*(Target_X-par[17]))/(2*par[18]*par[18]));
}
//=============================================For Fitting===========================================//
Double_t peak_ran(Double_t *x,Double_t *par)
{
    return par[0]*TMath::Exp(-((x[0]-par[1])*(x[0]-par[1]))/(2*par[2]*par[2]));
}
Double_t xray_peak_over_all(Double_t *x,Double_t *par)
{
    return par[0]+par[1]*TMath::Exp(-((x[0]-par[2])*(x[0]-par[2]))/(2*par[3]*par[3]))+par[4]*TMath::Exp(-((x[0]-par[5])*(x[0]-par[5]))/(2*par[6]*par[6]))+par[7]*TMath::Exp(-((x[0]-par[8])*(x[0]-par[8]))/(2*par[9]*par[9]))+par[10]*TMath::Exp(-((x[0]-par[11])*(x[0]-par[11]))/(2*par[12]*par[12]))+par[13]*TMath::Exp(-((x[0]-par[14])*(x[0]-par[14]))/(2*par[15]*par[15]))+par[16]*TMath::Exp(-((x[0]-par[17])*(x[0]-par[17]))/(2*par[18]*par[18]));
}
//================================================//

double *Hist_SetLimit_Plot_v2_Extract_Peak(int Option=0)
{
TLegend *leg= new TLegend(0.4,0.5,0.8,0.9);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.03);
leg->SetBorderSize(0);
leg->SetTextFont(22);
//68%->90%
const float SE_to_NT_P = 1.64458/0.994;
const int Data_element = 257;
double RE_DATA_Original[Data_element]; double RE_Rate_Original[Data_element]; double RE_DATA_Err_Original[Data_element]; double RE_Rate_Err_Original[Data_element];
double RE_DATA[Data_element]; double RE_Rate[Data_element]; double RE_DATA_Err[Data_element]; double RE_Rate_Err[Data_element];
static double RE_DATA_Aft[Data_element]; static double RE_Rate_Aft[Data_element]; static double RE_DATA_Err_Aft[Data_element]; static double RE_Rate_Err_Aft[Data_element];
for(int jjj=0 ; jjj<257 ; jjj++)
{
        RE_DATA_Original[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
        RE_Rate_Original[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err_Original[jjj]= 0;
        RE_Rate_Err_Original[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*SE_to_NT_P;
    //if(jjj==0) cout << "RE_Rate_Aft[jjj]+RE_Rate_Err_Aft[jjj]=2: " << RE_Rate_Original[jjj]+RE_Rate_Err_Original[jjj] << endl;

}
    
for(int jjj=0 ; jjj<257 ; jjj++) //Difference: 0.0501216 => 50.1216(eV/bin)  Minimim => 0.200486keV
{//0.997983
    if(jjj==0){
        RE_DATA[jjj]= (p103_le_VrV_ON_NaI1_50eV[jjj][0]/0.997617)/0.999957;
        //cout <<"RE_DATA[0]: " << (p103_le_VrV_ON_NaI1_50eV[jjj][0]/0.997617)/0.999957 << endl;
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1]*0.997617*0.999957;
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*SE_to_NT_P;}
    else if(jjj==256){
        RE_DATA[jjj]= (p103_le_VrV_ON_NaI1_50eV[jjj][0]/0.997617)/0.999957;
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1]*0.997617*0.999957;
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*SE_to_NT_P;}
    else{
        RE_DATA[jjj]= (p103_le_VrV_ON_NaI1_50eV[jjj][0]/0.997617)/0.999957;
        RE_Rate[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1]*0.997617*0.999957;
        RE_DATA_Err[jjj]= 0;
        RE_Rate_Err[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*SE_to_NT_P;}
}
vector<float> Two_Region_Fit={0.2,3.97,13.025};
vector<float> K_Shell={10.369,9.66,8.98,6.54,4.97};
vector<float> L_Shell={1.2977,1.1936,1.0961,0.769,0.564};
vector<float> K_shell_Exp; vector<float> K_shell_Exp_Upper; vector<float> K_shell_Exp_Lower;
vector<float> Transformation={0.133,0.111,0.119,0.106,0.106};
vector<float> L_shell_Area_upper; vector<float> L_shell_Area_middle; vector<float> L_shell_Area_lower;
Double_t Expe_Mean[5]; Double_t Expe_Error[5]; Double_t Expe_Constant[5];

TGraphErrors *cdexdata = new TGraphErrors(257,RE_DATA,RE_Rate,RE_DATA_Err,RE_Rate_Err);
cdexdata->SetName("cdexdata");
cdexdata->SetLineColor(1);
cdexdata->GetXaxis()->SetRangeUser(0,2.5);
//cdexdata->GetXaxis()->SetRangeUser(0,2.025);
cdexdata->GetYaxis()->SetRangeUser(0,180);
cdexdata->SetTitle("The rate of WIMP");
cdexdata->GetXaxis()->SetTitle("Recoil Energy[keV]");
cdexdata->GetYaxis()->SetTitle("Rate[Count/day*kg*keV]");
cdexdata->SetMarkerStyle(21);

Double_t par_latter[19];
Double_t par_latter_Full[19];
const Double_t *par_Error_latter;
    
Double_t par_previous[19];
Double_t par_previous_Full[19];

vector<float> Total_Constant;
vector<float> Total_Constant_Upper;
vector<float> Total_Constant_Lower;
vector<float> Total_Mean;
vector<float> Total_Error;
//==============================================================================================================================//
// The fitting in the latter region
//==> Find means and errors of the peaks of the Gaussians
//==>*Important because those parameters should be fixed in the next round
//Find the means and errors of data-driven gaussians
//==============================================================================================================================//
//Use the already-known parameters to fit for the best fitting as getting the platform of the whole-range peaks
//====High_energy_fitting
//10.369keV peak
TF1 *Gaus1    =new TF1("gaus1","gaus",K_Shell[0]-0.3,K_Shell[0]+0.3);
cdexdata->Fit(Gaus1,"R");
Gaus1->GetParameters(&par_latter[0]);
    Expe_Mean[0]=(par_latter[1]); Expe_Error[0]=(par_latter[2]);
par_Error_latter = Get_Error(cdexdata,K_Shell[0]-0.3,K_Shell[0]+0.3,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun01=Draw_Uncertainty(K_Shell[0]-0.3,K_Shell[0]+0.3,1,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun02=Draw_Uncertainty(K_Shell[0]-0.3,K_Shell[0]+0.3,2,par_latter[0]+par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
TF1 *fitfun03=Draw_Uncertainty(K_Shell[0]-0.3,K_Shell[0]+0.3,3,par_latter[0]-par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
K_shell_Exp.push_back(par_latter[0]);K_shell_Exp.push_back(par_latter[1]);K_shell_Exp.push_back(par_latter[2]);
K_shell_Exp_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
K_shell_Exp_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
    Total_Constant.push_back(par_latter[0]);Total_Mean.push_back(par_latter[1]);Total_Error.push_back(par_latter[2]);
    Total_Constant_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
    Total_Constant_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
    
//9.66keV peak
TF1 *Gaus2    =new TF1("gaus2","gaus",K_Shell[1]-0.3,K_Shell[1]+0.3);
cdexdata->Fit(Gaus2,"R+");
Gaus2->GetParameters(&par_latter[0]);
    Expe_Mean[1]=(par_latter[1]); Expe_Error[1]=(par_latter[2]);
    par_Error_latter = Get_Error(cdexdata,K_Shell[1]-0.3,K_Shell[1]+0.3,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun05=Draw_Uncertainty(K_Shell[1]-0.3,K_Shell[1]+0.3,1,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun06=Draw_Uncertainty(K_Shell[1]-0.3,K_Shell[1]+0.3,2,par_latter[0]+par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
TF1 *fitfun07=Draw_Uncertainty(K_Shell[1]-0.3,K_Shell[1]+0.3,3,par_latter[0]-par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
K_shell_Exp.push_back(par_latter[0]);K_shell_Exp.push_back(par_latter[1]);K_shell_Exp.push_back(par_latter[2]);
K_shell_Exp_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
K_shell_Exp_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
    Total_Constant.push_back(par_latter[0]);Total_Mean.push_back(par_latter[1]);Total_Error.push_back(par_latter[2]);
    Total_Constant_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
    Total_Constant_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
//8.98 keV peak
TF1 *Gaus3    =new TF1("gaus3","gaus",K_Shell[2]-0.3,K_Shell[2]+0.3);
cdexdata->Fit(Gaus3,"R+");
Gaus3->GetParameters(&par_latter[0]);
    Expe_Mean[2]=(par_latter[1]); Expe_Error[2]=(par_latter[2]);
    par_Error_latter = Get_Error(cdexdata,K_Shell[2]-0.3,K_Shell[2]+0.3,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun09 =Draw_Uncertainty(K_Shell[2]-0.3,K_Shell[2]+0.3,1,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun010=Draw_Uncertainty(K_Shell[2]-0.3,K_Shell[2]+0.3,2,par_latter[0]+par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
TF1 *fitfun011=Draw_Uncertainty(K_Shell[2]-0.3,K_Shell[2]+0.3,3,par_latter[0]-par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
K_shell_Exp.push_back(par_latter[0]);K_shell_Exp.push_back(par_latter[1]);K_shell_Exp.push_back(par_latter[2]);
K_shell_Exp_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
K_shell_Exp_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
    Total_Constant.push_back(par_latter[0]);Total_Mean.push_back(par_latter[1]);Total_Error.push_back(par_latter[2]);
    Total_Constant_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
    Total_Constant_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
//6.54 keV peak
TF1 *Gaus4    =new TF1("gaus4","gaus",K_Shell[3]-0.3,K_Shell[3]+0.3);
cdexdata->Fit(Gaus4,"R+");
Gaus4->GetParameters(&par_latter[0]);
    Expe_Mean[3]=(par_latter[1]); Expe_Error[3]=(par_latter[2]);
    par_Error_latter = Get_Error(cdexdata,K_Shell[3]-0.3,K_Shell[3]+0.3,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun012=Draw_Uncertainty(K_Shell[3]-0.3,K_Shell[3]+0.3,1,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun013=Draw_Uncertainty(K_Shell[3]-0.3,K_Shell[3]+0.3,2,par_latter[0]+par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
TF1 *fitfun014=Draw_Uncertainty(K_Shell[3]-0.3,K_Shell[3]+0.3,3,par_latter[0]-par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
K_shell_Exp.push_back(par_latter[0]);K_shell_Exp.push_back(par_latter[1]);K_shell_Exp.push_back(par_latter[2]);
K_shell_Exp_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
K_shell_Exp_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
    Total_Constant.push_back(par_latter[0]);Total_Mean.push_back(par_latter[1]);Total_Error.push_back(par_latter[2]);
    Total_Constant_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
    Total_Constant_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
    //4.97 keV peak
TF1 *Gaus5    =new TF1("gaus5","gaus",K_Shell[4]-0.3,K_Shell[4]+0.3);
cdexdata->Fit(Gaus5,"R+");
Gaus5->GetParameters(&par_latter[0]);
    Expe_Mean[4]=(par_latter[1]); Expe_Error[4]=(par_latter[2]);
    par_Error_latter = Get_Error(cdexdata,K_Shell[4]-0.3,K_Shell[4]+0.3,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun015=Draw_Uncertainty(K_Shell[4]-0.3,K_Shell[4]+0.3,1,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun016=Draw_Uncertainty(K_Shell[4]-0.3,K_Shell[4]+0.3,2,par_latter[0]+par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
TF1 *fitfun017=Draw_Uncertainty(K_Shell[4]-0.3,K_Shell[4]+0.3,3,par_latter[0]-par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
K_shell_Exp.push_back(par_latter[0]);K_shell_Exp.push_back(par_latter[1]);K_shell_Exp.push_back(par_latter[2]);
K_shell_Exp_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
K_shell_Exp_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
    Total_Constant.push_back(par_latter[0]);Total_Mean.push_back(par_latter[1]);Total_Error.push_back(par_latter[2]);
    Total_Constant_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
    Total_Constant_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);

for(int kkk=0 ; kkk<5 ; kkk++)
{
    Expe_Constant[kkk]  =(K_shell_Exp[0+(3*kkk)]);
    Expe_Error[kkk]     =(K_shell_Exp[2+(3*kkk)]);
    L_shell_Area_middle.push_back( abs(Gaussian_Area(Expe_Constant[kkk],Expe_Error[kkk]))*Transformation[kkk] );
    L_shell_Area_upper.push_back( abs(Gaussian_Area(K_shell_Exp_Upper[kkk],Expe_Error[kkk]))*Transformation[kkk] );
    L_shell_Area_lower.push_back( abs(Gaussian_Area(K_shell_Exp_Lower[kkk],Expe_Error[kkk]))*Transformation[kkk]  );
    cout << "Constant: " << Expe_Constant[kkk] << "Sigma: " << Expe_Error[kkk] << endl;
    cout << "K_shell_Area: " << abs(Gaussian_Area(K_shell_Exp_Upper[kkk],Expe_Error[kkk]))*Transformation[kkk] << endl;
    cout << "L_shell_Area_theoretical: " << L_shell_Area_middle[kkk] << endl;
}
    
float Calculated_constant_upper;
float Calculated_constant_middle;
float Calculated_constant_lower;
//====Low_energy_fitting
Calculated_constant_upper  = Gaussian_Constant(L_shell_Area_upper[0],Resolution_of_Signal(L_Shell[0]));
Calculated_constant_middle = Gaussian_Constant(L_shell_Area_middle[0],Resolution_of_Signal(L_Shell[0]));
Calculated_constant_lower  = Gaussian_Constant(L_shell_Area_lower[0],Resolution_of_Signal(L_Shell[0]));

TF1 *fitfun018=Draw_Uncertainty(L_Shell[0]-0.2,L_Shell[0]+0.2,4,Calculated_constant_middle,L_Shell[0],Resolution_of_Signal(L_Shell[0]));
TF1 *fitfun019=Draw_Uncertainty(L_Shell[0]-0.2,L_Shell[0]+0.2,2,Calculated_constant_upper,L_Shell[0],Resolution_of_Signal(L_Shell[0]));
TF1 *fitfun020=Draw_Uncertainty(L_Shell[0]-0.2,L_Shell[0]+0.2,3,Calculated_constant_lower,L_Shell[0],Resolution_of_Signal(L_Shell[0]));
Total_Constant.push_back(Calculated_constant_middle);Total_Mean.push_back(L_Shell[0]-0.006);Total_Error.push_back(Resolution_of_Signal(L_Shell[0]));
Total_Constant_Upper.push_back(Calculated_constant_upper);
Total_Constant_Lower.push_back(Calculated_constant_lower);

//===============================
Calculated_constant_upper  = Gaussian_Constant(L_shell_Area_upper[1],Resolution_of_Signal(L_Shell[1]));
Calculated_constant_middle = Gaussian_Constant(L_shell_Area_middle[1],Resolution_of_Signal(L_Shell[1]));
Calculated_constant_lower  = Gaussian_Constant(L_shell_Area_lower[1],Resolution_of_Signal(L_Shell[1]));
    
TF1 *fitfun021=Draw_Uncertainty(L_Shell[1]-0.2,L_Shell[1]+0.2,30,Calculated_constant_middle,L_Shell[1],Resolution_of_Signal(L_Shell[1]));
TF1 *fitfun022=Draw_Uncertainty(L_Shell[1]-0.2,L_Shell[1]+0.2,2,Calculated_constant_upper,L_Shell[1],Resolution_of_Signal(L_Shell[1]));
TF1 *fitfun023=Draw_Uncertainty(L_Shell[1]-0.2,L_Shell[1]+0.2,3,Calculated_constant_lower,L_Shell[1],Resolution_of_Signal(L_Shell[1]));
Total_Constant.push_back(Calculated_constant_middle);Total_Mean.push_back(L_Shell[1]);Total_Error.push_back(Resolution_of_Signal(L_Shell[1]));
Total_Constant_Upper.push_back(Calculated_constant_upper);
Total_Constant_Lower.push_back(Calculated_constant_lower);

//===============================
Calculated_constant_upper  = Gaussian_Constant(L_shell_Area_upper[2],Resolution_of_Signal(L_Shell[2]));
Calculated_constant_middle = Gaussian_Constant(L_shell_Area_middle[2],Resolution_of_Signal(L_Shell[2]));
Calculated_constant_lower  = Gaussian_Constant(L_shell_Area_lower[2],Resolution_of_Signal(L_Shell[2]));

TF1 *fitfun024=Draw_Uncertainty(L_Shell[2]-0.2,L_Shell[2]+0.2,45,Calculated_constant_middle,L_Shell[2],Resolution_of_Signal(L_Shell[2]));
TF1 *fitfun025=Draw_Uncertainty(L_Shell[2]-0.2,L_Shell[2]+0.2,2,Calculated_constant_upper,L_Shell[2],Resolution_of_Signal(L_Shell[2]));
TF1 *fitfun026=Draw_Uncertainty(L_Shell[2]-0.2,L_Shell[2]+0.2,3,Calculated_constant_lower,L_Shell[2],Resolution_of_Signal(L_Shell[2]));
Total_Constant.push_back(Calculated_constant_middle);Total_Mean.push_back(L_Shell[2]);Total_Error.push_back(Resolution_of_Signal(L_Shell[2]));
Total_Constant_Upper.push_back(Calculated_constant_upper);
Total_Constant_Lower.push_back(Calculated_constant_lower);

//===============================
//===============================
Calculated_constant_upper  = Gaussian_Constant(L_shell_Area_upper[3],Resolution_of_Signal(L_Shell[3]));
Calculated_constant_middle = Gaussian_Constant(L_shell_Area_middle[3],Resolution_of_Signal(L_Shell[3]));
Calculated_constant_lower  = Gaussian_Constant(L_shell_Area_lower[3],Resolution_of_Signal(L_Shell[3]));

TF1 *fitfun030=Draw_Uncertainty(L_Shell[3]-0.2,L_Shell[3]+0.2,9,Calculated_constant_middle,L_Shell[3],Resolution_of_Signal(L_Shell[3]));
TF1 *fitfun031=Draw_Uncertainty(L_Shell[3]-0.2,L_Shell[3]+0.2,2,Calculated_constant_upper,L_Shell[3],Resolution_of_Signal(L_Shell[3]));
TF1 *fitfun032=Draw_Uncertainty(L_Shell[3]-0.2,L_Shell[3]+0.2,3,Calculated_constant_lower,L_Shell[3],Resolution_of_Signal(L_Shell[3]));
Total_Constant.push_back(Calculated_constant_middle);Total_Mean.push_back(L_Shell[3]);Total_Error.push_back(Resolution_of_Signal(L_Shell[3]));
Total_Constant_Upper.push_back(Calculated_constant_upper);
Total_Constant_Lower.push_back(Calculated_constant_lower);

//===============================
    
Calculated_constant_upper  = Gaussian_Constant(L_shell_Area_upper[4],Resolution_of_Signal(L_Shell[4]));
Calculated_constant_middle = Gaussian_Constant(L_shell_Area_middle[4],Resolution_of_Signal(L_Shell[4]));
Calculated_constant_lower  = Gaussian_Constant(L_shell_Area_lower[4],Resolution_of_Signal(L_Shell[4]));

TF1 *fitfun033=Draw_Uncertainty(L_Shell[4]-0.2,L_Shell[4]+0.2,12,Calculated_constant_middle,L_Shell[4],Resolution_of_Signal(L_Shell[4]));
TF1 *fitfun034=Draw_Uncertainty(L_Shell[4]-0.2,L_Shell[4]+0.2,2,Calculated_constant_upper,L_Shell[4],Resolution_of_Signal(L_Shell[4]));
TF1 *fitfun035=Draw_Uncertainty(L_Shell[4]-0.2,L_Shell[4]+0.2,3,Calculated_constant_lower,L_Shell[4],Resolution_of_Signal(L_Shell[4]));
Total_Constant.push_back(Calculated_constant_middle);Total_Mean.push_back(L_Shell[4]);Total_Error.push_back(Resolution_of_Signal(L_Shell[4]));
Total_Constant_Upper.push_back(Calculated_constant_upper);
Total_Constant_Lower.push_back(Calculated_constant_lower);

//===============================
//====Platform_fitting
float Platform_constant;
float Platform_Error;
TF1 *Constant0    =new TF1("Constant0","[0]",11.2,13);
cdexdata->Fit(Constant0,"R+");
Platform_constant = Constant0->GetParameter(0);
Platform_Error    = Constant0->GetParError(0);
TF1 *Constant_Draw    =new TF1("Constant1","[0]",0,2.5);
Constant_Draw->SetParameter(0,Platform_constant);
Constant_Draw->SetLineColor(14);
cout << "Platform_Error" << Platform_Error << endl;

// The fitting in the previous region
//==> Find means and errors of the peaks of the Gaussians
//==>*Important because those parameters should be fixed in the next round
    
    
for(int jjj=0 ; jjj<257 ; jjj++)
{
    float Error_New_fitting = 0;
    RE_DATA_Aft[jjj]= RE_DATA[jjj];
    //cout << "RE_DATA_Aft[jjj]: " << RE_DATA_Aft[jjj] << endl;
    //Fitting_Error
    Error_New_fitting = Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant_Upper,Total_Mean,Total_Error)-Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant_Lower,Total_Mean,Total_Error);
    RE_Rate_Aft[jjj]= (RE_Rate[jjj]-Platform_constant)-Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant,Total_Mean,Total_Error);
    if(jjj==7)
    {
        cout << "=============================================================" << endl;
        cout << "RE_Rate[jjj]: " << RE_Rate[jjj] << endl;
        cout << "Platform_constant: " << Platform_constant << endl;
        cout << "Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant,Total_Mean,Total_Error): " << Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant,Total_Mean,Total_Error) << endl;
        cout << "RE_Rate_Aft[jjj]: " << RE_Rate_Aft[jjj] << endl;
        cout << "=============================================================" << endl;
     }
    RE_DATA_Err_Aft[jjj]= 0.025;
    RE_Rate_Err_Aft[jjj] = New_Error(p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994,Platform_Error*1.64458/0.994,Error_New_fitting);
    //if(jjj==0) cout << "RE_Rate_Aft[jjj]+RE_Rate_Err_Aft[jjj]=2: " << RE_Rate_Aft[jjj]+RE_Rate_Err_Aft[jjj] << endl;
    
}
/*
TGraphErrors *cdexdata1 = new TGraphErrors(257,RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft);
cdexdata1->SetName("cdexdata1");
cdexdata1->SetLineColor(6);
cdexdata1->GetXaxis()->SetRangeUser(0,12);
cdexdata1->GetYaxis()->SetRangeUser(-10,100);
cdexdata1->SetTitle("The rate of WIMP");
cdexdata1->GetXaxis()->SetTitle("Recoil Energy[keV]");
cdexdata1->GetYaxis()->SetTitle("Rate[Count/day*kg*keV]");
cdexdata1->SetMarkerStyle(21);
cdexdata1->SetMarkerColor(6);
cdexdata1->Draw();
    c1->Print("See.pdf");
    */
//double *CPKKD_Aft = cpkkd_calculation_New(RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft,3);
if(Option==0) return RE_DATA_Aft;
if(Option==1) return RE_Rate_Aft;
if(Option==2) return RE_DATA_Err_Aft;
if(Option==3) return RE_Rate_Err_Aft;
    
    /*
     //Normally distributed Velocity without the sturctral world
     double Scaling_old[29]; double Scaling_new[29]; double WIMP_mass_array[29];
     for(int k=2; k<31 ; k++)
     {
     WIMP_mass_array[k-2] = (double)k;
     Scaling_new[k-2] = 1e-40*(cpkkd_calculation_Scaling_Factor(RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft,(double)k))[0];
     Scaling_old[k-2] = 1e-40*(cpkkd_calculation_Scaling_Factor(RE_DATA_Original,RE_Rate_Original,RE_DATA_Err_Original,RE_Rate_Err_Original,(double)k))[0];
     {
     cout << "(cpkkd_calculation_Scaling_Factor(RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft,(double)k))[1]: " << (cpkkd_calculation_Scaling_Factor(RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft,(double)k))[1] << endl;
     cout << "Scaling_new[k-2]: " << Scaling_new[k-2] << endl;
     cout << "(cpkkd_calculation_Scaling_Factor(RE_DATA_Original,RE_Rate_Original,RE_DATA_Err_Original,RE_Rate_Err_Original,(double)k))[1]: " << (cpkkd_calculation_Scaling_Factor(RE_DATA_Original,RE_Rate_Original,RE_DATA_Err_Original,RE_Rate_Err_Original,(double)k))[1] << endl;
     cout << "Scaling_old[k-2]: " << Scaling_old[k-2] << endl;
     }
     }
     
     
     double *Dark_10GeV_cpkkd = cpkkd_calculation_New(RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft,3);
     double *Dark_10GeV_T_QF  = T_QF_Array(3);
     TGraph *gcpkkdX2 = new TGraph(reso_T,Dark_10GeV_T_QF,Dark_10GeV_cpkkd);
     gcpkkdX2->SetName("gcpkkdX2");
     gcpkkdX2->SetLineColor(7);
     
     TGraph *Set_Limits_plot = new TGraph(29,WIMP_mass_array,Scaling_old);
     Set_Limits_plot->SetLineColor(2);
     Set_Limits_plot->SetMarkerStyle(21);
     Set_Limits_plot->SetMarkerColor(2);
     
     Set_Limits_plot->SetMarkerSize(1.2);
     Set_Limits_plot->SetTitle("Set_Limits_plot");
     Set_Limits_plot->GetXaxis()->SetTitle("WIMP mass[GeV]");
     Set_Limits_plot->GetYaxis()->SetTitle("");
     Set_Limits_plot->GetXaxis()->SetRangeUser(0,40);
     Set_Limits_plot->GetYaxis()->SetRangeUser(1e-44,1e-33);
     
     TGraph *Set_Limits_plot_1 = new TGraph(29,WIMP_mass_array,Scaling_new);
     Set_Limits_plot_1->SetLineColor(3);
     Set_Limits_plot_1->SetMarkerColor(3);
     Set_Limits_plot_1->SetMarkerStyle(21);
     Set_Limits_plot_1->SetMarkerSize(1.2);
     
     
     leg->AddEntry(Set_Limits_plot,"Original case","pl");
     leg->AddEntry(Set_Limits_plot_1 ,"Subtracted case","pl");
     
     Set_Limits_plot->Draw("apl");
     Set_Limits_plot_1->Draw("plsame");
     
     //cdexdata1->Draw("ap");
     //gcpkkdX2->Draw("plsame");
     leg->Draw();
     
     c1->SetLogy();
     c1->Print("TEXONO_setlimit_3GeV.png");
     //68% (1sigma) => 90% (2sigma)
     */
    
    
    
    // y = 0.997983 x
/*
//=======================Scaling_factor_for_energy=================//
 //Check the Energy fitting parameters with the theoretical values
 Double_t K_Shell_sort[5]={10.369,9.66,8.98,6.54,4.97};
 Double_t K_Shell_sort_Error[5]={0,0,0,0,0};
 
 Double_t par_Scaling[2];
TGraphErrors *Peak_Position           = new TGraphErrors(5,K_Shell_sort,Expe_Mean,K_Shell_sort_Error,Expe_Error);
TF1 *Energy_transformed = new TF1("f1","[0]*x",0,12);
Peak_Position->SetName("Peak_Position");
Peak_Position->SetLineColor(1);
Peak_Position->GetXaxis()->SetLimits(0,12);
Peak_Position->GetYaxis()->SetRangeUser(0,12);
Peak_Position->SetTitle("Energy_transformed");
Peak_Position->GetXaxis()->SetTitle("Theoretical K-shell Energies[keV]");
Peak_Position->GetYaxis()->SetTitle("Experimental K-shell Energies[keV]");
Peak_Position->SetMarkerStyle(21);
Peak_Position->Fit(Energy_transformed,"R");
Energy_transformed->GetParameters(&par_Scaling[0]);
leg->AddEntry(Energy_transformed,Form("y=%fx",par_Scaling[0]),"lp");
Peak_Position->Draw("apE");
*/
    
/*
Double_t Theoretical_Resolution[5];Double_t Experimental_Resolution[5]; Double_t Experimental_Resolution_Error[5];
for(int kkk=0 ; kkk<5 ; kkk++)
{
    Theoretical_Resolution[kkk] = Resolution_of_Signal(K_Shell_sort[kkk]);
    Experimental_Resolution[kkk]= abs(Par_Use_latter[3+(3*kkk)]);
    Experimental_Resolution_Error[kkk] = abs(fitfun04->GetParError(3+(3*kkk)));
    cout << "Theoretical_Resolution[kkk]: " << Theoretical_Resolution[kkk] << endl;
    cout << "Experimental_Resolution[kkk]: " << Experimental_Resolution[kkk] << endl;
    cout << "Experimental_Resolution_Error[kkk]: " << Experimental_Resolution_Error[kkk] << endl;
}

TGraphErrors *Peak_Width_Theoretical  = new TGraphErrors(5,K_Shell_sort,Theoretical_Resolution,K_Shell_sort_Error,K_Shell_sort_Error);
Peak_Width_Theoretical->SetName("Peak_Width_Theoretical");
Peak_Width_Theoretical->SetLineColor(1);
Peak_Width_Theoretical->GetXaxis()->SetRangeUser(4,12);
Peak_Width_Theoretical->GetYaxis()->SetRangeUser(0,0.2);
Peak_Width_Theoretical->SetTitle("Energy_transformed");
Peak_Width_Theoretical->GetXaxis()->SetTitle("Theoretical K-shell Energies[keV]");
Peak_Width_Theoretical->GetYaxis()->SetTitle("Resolution");
Peak_Width_Theoretical->SetMarkerStyle(21);

TGraphErrors *Peak_Width_Experiment   = new TGraphErrors(5,K_Shell_sort,Experimental_Resolution,K_Shell_sort_Error,Experimental_Resolution_Error);
Peak_Width_Experiment->SetLineColor(2);
Peak_Width_Experiment->SetMarkerStyle(21);
Peak_Width_Experiment->SetMarkerColor(2);
TF1 *Peak_Width_Experiment_Fit = new TF1("f1","[0]+[1]*sqrt(x)",6,11);
Peak_Width_Experiment->Fit(Peak_Width_Experiment_Fit,"R");

Peak_Width_Theoretical->Draw("apE");
Peak_Width_Experiment->Draw("pEsame");
leg->AddEntry(Peak_Width_Theoretical,"Theoretical Resolution","p");
leg->AddEntry(Peak_Width_Experiment,"Experimental Resolution","p");

TLatex *tex = new TLatex(Total_Mean[5],160,"^{68}Ge");
tex->SetTextColor(4);
tex->SetTextFont(42);
tex->SetTextSize(0.03422619);
tex->SetLineWidth(2);

TLatex *tex1 = new TLatex(Total_Mean[6],10,"^{68}Ga");
tex1->SetTextColor(30);
tex1->SetTextFont(42);
tex1->SetTextSize(0.03422619);
tex1->SetLineWidth(2);

TLatex *tex2 = new TLatex(Total_Mean[7]-0.1,40,"^{65}Zn");
tex2->SetTextColor(45);
tex2->SetTextFont(42);
tex2->SetTextSize(0.03422619);
tex2->SetLineWidth(2);

TLatex *tex3 = new TLatex(Total_Mean[8],30,"^{55}Fe");
tex3->SetTextColor(9);
tex3->SetTextFont(42);
tex3->SetTextSize(0.03422619);
tex3->SetLineWidth(2);

TLatex *tex4 = new TLatex(Total_Mean[9],30,"^{49}V");
tex4->SetTextColor(12);
tex4->SetTextFont(42);
tex4->SetTextSize(0.03422619);
tex4->SetLineWidth(2);

TLatex *tex5 = new TLatex(0.05,10,"Platform");
tex5->SetTextColor(14);
tex5->SetTextFont(42);
tex5->SetTextSize(0.03422619);
tex5->SetLineWidth(2);

cdexdata->Draw("apE");
cdexdata1->Draw("pEsame");
//Gaus1->Draw("Lsame");
    
fitfun01->Draw("same");
fitfun02->Draw("Lsame");
fitfun03->Draw("Lsame");
fitfun05->Draw("Lsame");
fitfun06->Draw("Lsame");
fitfun07->Draw("Lsame");
fitfun09->Draw("Lsame");
fitfun010->Draw("Lsame");
fitfun011->Draw("Lsame");
fitfun012->Draw("Lsame");
fitfun013->Draw("Lsame");
fitfun014->Draw("Lsame");
fitfun015->Draw("Lsame");
fitfun016->Draw("Lsame");
fitfun017->Draw("Lsame");
fitfun018->Draw("Lsame");
fitfun019->Draw("Lsame");
fitfun020->Draw("Lsame");
fitfun021->Draw("Lsame");
fitfun022->Draw("Lsame");
fitfun023->Draw("Lsame");
fitfun024->Draw("Lsame");
fitfun025->Draw("Lsame");
fitfun026->Draw("Lsame");
fitfun030->Draw("Lsame");
fitfun031->Draw("Lsame");
fitfun032->Draw("Lsame");
fitfun033->Draw("Lsame");
fitfun034->Draw("Lsame");
fitfun035->Draw("Lsame");
Constant_Draw->Draw("Lsame");
tex->Draw("same");
tex1->Draw("same");
tex2->Draw("same");
tex3->Draw("same");
tex4->Draw("same");
tex5->Draw("same");
*/
   
}
