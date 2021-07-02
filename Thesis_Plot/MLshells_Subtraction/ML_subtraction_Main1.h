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
#include "/Users/yehchihhsiang/Desktop/GITHUB_TEXONO/Code_for_Analysis/VrV_le_PPC103_Xenon_Subtracted_ON_NaI1_180418_190830_50eV_date20200420.h"

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
{//
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
float Resolution_of_Signal(float Energy_keV)//Detector Resolution==>Real TEXONO detector
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
float Gaussian_Area(float Constant, float Uncertainty)//A*C*sqrt(2*pi)=Gaussian Area
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

double *ML_subtraction_Main(int Option=0)
{
TCanvas * c1 = new TCanvas("c", "c", 0,0,1000,1000);
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.04,"XY");
gStyle->SetTitleFont(62,"XY");
gStyle->SetLegendFont(62);

TLegend *leg= new TLegend(0.1,0.7,0.3,0.9);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.04);
leg->SetBorderSize(0);
leg->SetTextFont(22);

const float SE_to_NT_P = 1.64458/0.994;//68%->90%
const int Data_element = 257; // I have 257 Points for my dataset
//*Coding Issue here: You should remember that the type of the number of the element in an array, it should be the "const".
    
double RE_DATA_Original[Data_element];//Original Energy from data
double RE_Rate_Original[Data_element];//Original Rate(count/kg*keV*day) from data
double RE_DATA_Err_Original[Data_element];//Oringal energy range for each bin
double RE_Rate_Err_Original[Data_element];//Oringal Rate range for each bin
    
double RE_DATA[Data_element];//After
double RE_Rate[Data_element];//After
double RE_DATA_Err[Data_element];//After
double RE_Rate_Err[Data_element];//After
    
static double RE_DATA_Aft[Data_element];
static double RE_Rate_Aft[Data_element];
static double RE_DATA_Err_Aft[Data_element];
static double RE_Rate_Err_Aft[Data_element];
//If you would like to return the array, then you should use static function.
    
/*
for(int jjj=0 ; jjj<257 ; jjj++)
{
        RE_DATA_Original[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][0];
        RE_Rate_Original[jjj]= p103_le_VrV_ON_NaI1_50eV[jjj][1];
        RE_DATA_Err_Original[jjj]= 0;
        RE_Rate_Err_Original[jjj] = p103_le_VrV_ON_NaI1_50eV[jjj][2]*SE_to_NT_P;

}
 */
    
for(int jjj=0 ; jjj<257 ; jjj++) //Difference: 0.0501216 => 50.1216(eV/bin)  Minimim => 0.200486keV
{//We can tune the bin width, and if the bin width is tuned to (Original_Value)*A, then the rate should time (1/A).
    if(jjj==0){
        RE_DATA[jjj]= (p103_le_VrV_ON_NaI1_50eV[jjj][0]/0.997617)/0.999957;
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
//*The middle energies of K_shell and L_shell
vector<float> K_Shell={10.369,9.66,8.98,6.54,4.97};
vector<float> L_Shell={1.2977,1.1936,1.0961,0.769,0.564};
//*First we will get the areas of K Shell
vector<float> K_shell_Exp;
vector<float> K_shell_Exp_Upper;
vector<float> K_shell_Exp_Lower;
//*Then, we have the transformation between two areas(the K shell areas * Transformation = L shell areas)
vector<float> Transformation={0.133,0.111,0.119,0.106,0.106};
//*Second, we will get the areas of L shell
vector<float> L_shell_Area_upper;
vector<float> L_shell_Area_middle;
vector<float> L_shell_Area_lower;
//Get the fitting parameters
Double_t Expe_Mean[5];
Double_t Expe_Error[5];
Double_t Expe_Constant[5];

//Drawing the plot! First, plotting the data with the modification version
TGraphErrors *cdexdata = new TGraphErrors(257,RE_DATA,RE_Rate,RE_DATA_Err,RE_Rate_Err);
cdexdata->SetName("cdexdata");
cdexdata->SetLineColor(1);
    cdexdata->GetXaxis()->SetRangeUser(4.5,11);//(4.5,11)(0,2)
cdexdata->GetYaxis()->SetRangeUser(0,800);//(0,800),(0,200)
cdexdata->SetTitle("");
cdexdata->GetXaxis()->SetTitle("Recoil Energy[keV]");
cdexdata->GetYaxis()->SetTitle("Rate[Count/day*kg*keV]");
cdexdata->SetMarkerStyle(8);
cdexdata->Draw("AP");
    
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
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%High_energy_fitting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//10.369keV peak
TF1 *Gaus1    =new TF1("gaus1","gaus",K_Shell[0]-0.3,K_Shell[0]+0.3);
cdexdata->Fit(Gaus1,"R");
Gaus1->GetParameters(&par_latter[0]);
Expe_Mean[0]=(par_latter[1]); Expe_Error[0]=(par_latter[2]);
//par_latter[0]=constant par_latter[1]=mean par_latter[2]=uncertainty
par_Error_latter = Get_Error(cdexdata,K_Shell[0]-0.3,K_Shell[0]+0.3,par_latter[0],par_latter[1],par_latter[2]);
//Draw the uncertainty for the Gaussians
TF1 *fitfun01=Draw_Uncertainty(K_Shell[0]-0.3,K_Shell[0]+0.3,1,par_latter[0],par_latter[1],par_latter[2]);
TF1 *fitfun02=Draw_Uncertainty(K_Shell[0]-0.3,K_Shell[0]+0.3,2,par_latter[0]+par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
TF1 *fitfun03=Draw_Uncertainty(K_Shell[0]-0.3,K_Shell[0]+0.3,3,par_latter[0]-par_Error_latter[0]*SE_to_NT_P,par_latter[1],par_latter[2]);
K_shell_Exp.push_back(par_latter[0]);K_shell_Exp.push_back(par_latter[1]);K_shell_Exp.push_back(par_latter[2]);
K_shell_Exp_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);//Constant error
K_shell_Exp_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);
Total_Constant.push_back(par_latter[0]);Total_Mean.push_back(par_latter[1]);Total_Error.push_back(par_latter[2]);
Total_Constant_Upper.push_back(par_latter[0]+par_Error_latter[0]*SE_to_NT_P);
Total_Constant_Lower.push_back(par_latter[0]-par_Error_latter[0]*SE_to_NT_P);

//9.66keV peak
TF1 *Gaus2    =new TF1("gaus2","gaus",K_Shell[1]-0.3,K_Shell[1]+0.3);
cdexdata->Fit(Gaus2,"R+");
Gaus2->GetParameters(&par_latter[0]);
Expe_Mean[1]=(par_latter[1]); Expe_Error[1]=(par_latter[2]);//Get the parameters including the mean and error
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Transformation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(int kkk=0 ; kkk<5 ; kkk++)
{//Calculate the areas of the K shell
    Expe_Constant[kkk]  =(K_shell_Exp[0+(3*kkk)]);
    Expe_Error[kkk]     =(K_shell_Exp[2+(3*kkk)]);
    L_shell_Area_middle.push_back( abs(Gaussian_Area(Expe_Constant[kkk],Expe_Error[kkk]))*Transformation[kkk] );
    L_shell_Area_upper.push_back( abs(Gaussian_Area(K_shell_Exp_Upper[kkk],Expe_Error[kkk]))*Transformation[kkk] );
    L_shell_Area_lower.push_back( abs(Gaussian_Area(K_shell_Exp_Lower[kkk],Expe_Error[kkk]))*Transformation[kkk]  );
    cout << "Constant: " << Expe_Constant[kkk] << "Sigma: " << Expe_Error[kkk] << endl;
    cout << "K_shell_Area: " << abs(Gaussian_Area(K_shell_Exp_Upper[kkk],Expe_Error[kkk]))*Transformation[kkk] << endl;
    cout << "L_shell_Area_theoretical: " << L_shell_Area_middle[kkk] << endl;
}
    
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Low_energy_fitting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
//Because of the high peak near Ge68, we use the latter part of the energy spectrum to find the platform constant and error
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
    //Fitting_Error
    Error_New_fitting = Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant_Upper,Total_Mean,Total_Error)-Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant_Lower,Total_Mean,Total_Error);
    RE_Rate_Aft[jjj]= (RE_Rate[jjj]-Platform_constant)-Extract_K_L_shell(RE_DATA_Aft[jjj],Total_Constant,Total_Mean,Total_Error);
    RE_DATA_Err_Aft[jjj]= 0.025;
    RE_Rate_Err_Aft[jjj] = New_Error(p103_le_VrV_ON_NaI1_50eV[jjj][2]*1.64458/0.994,Platform_Error*1.64458/0.994,Error_New_fitting);
}

TGraphErrors *cdexdata1 = new TGraphErrors(257,RE_DATA_Aft,RE_Rate_Aft,RE_DATA_Err_Aft,RE_Rate_Err_Aft);
cdexdata1->SetName("cdexdata1");
cdexdata1->SetLineColor(6);
//cdexdata1->GetXaxis()->SetRangeUser(0,2.5);
cdexdata1->GetXaxis()->SetRangeUser(0,11);
cdexdata1->GetYaxis()->SetRangeUser(-10,800);
cdexdata1->GetXaxis()->SetTitle("Recoil Energy[keV]");
cdexdata1->GetYaxis()->SetTitle("Rate[Count/day*kg*keV]");
cdexdata1->SetMarkerStyle(8);
cdexdata1->SetMarkerColor(6);
//cdexdata1->Draw("PEsame");

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    float Shift=0.2;
TLatex *tex11 = new TLatex(K_shell_Exp[1]-Shift,750,"^{68, 71}Ge");
tex11->SetTextColor(4);
tex11->SetTextFont(42);
tex11->SetTextSize(0.03422619);
tex11->SetLineWidth(2);

TLatex *tex22 = new TLatex(K_shell_Exp[4]-Shift,45,"^{68}Ga");
tex22->SetTextColor(30);
tex22->SetTextFont(42);
tex22->SetTextSize(0.03422619);
tex22->SetLineWidth(2);

TLatex *tex33 = new TLatex(K_shell_Exp[7]-Shift,130,"^{65}Zn");
tex33->SetTextColor(45);
tex33->SetTextFont(42);
tex33->SetTextSize(0.03422619);
tex33->SetLineWidth(2);

TLatex *tex44 = new TLatex(K_shell_Exp[10]-Shift,50,"^{55}Fe");
tex44->SetTextColor(7);
tex44->SetTextFont(42);
tex44->SetTextSize(0.03422619);
tex44->SetLineWidth(2);

TLatex *tex55 = new TLatex(K_shell_Exp[13]-Shift,50,"^{49}V");
tex55->SetTextColor(9);
tex55->SetTextFont(42);
tex55->SetTextSize(0.03422619);
tex55->SetLineWidth(2);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TLatex *tex = new TLatex(Total_Mean[5],160,"^{68,71}Ge");
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
//cdexdata1->Draw("pEsame");
//Gaus1->Draw("Lsame");
    
leg->AddEntry(cdexdata,"Raw data","P");
//leg->AddEntry(cdexdata1,"L/M shells X-ray subtracted data","P");
    
Gaus1->SetLineColor(4);
Gaus1->Draw("Lsame");
Gaus2->SetLineColor(30);
Gaus2->Draw("Lsame");
Gaus3->SetLineColor(45);
Gaus3->Draw("Lsame");
Gaus4->SetLineColor(7);
Gaus4->Draw("Lsame");
Gaus5->SetLineColor(9);
Gaus5->Draw("Lsame");
     

fitfun01->Draw("Lsame");
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
 
    
/*
fitfun018->SetLineColor(4);
fitfun018->Draw("Lsame");
fitfun021->SetLineColor(30);
fitfun021->Draw("Lsame");
fitfun024->SetLineColor(45);
fitfun024->Draw("Lsame");
fitfun030->SetLineColor(9);
fitfun030->Draw("Lsame");
fitfun033->SetLineColor(12);
fitfun033->Draw("Lsame");
  */
    /*
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
    */
Constant_Draw->Draw("Lsame");
tex->Draw("same");
tex1->Draw("same");
tex2->Draw("same");
tex3->Draw("same");
tex4->Draw("same");
tex5->Draw("same");
tex11->Draw("same");
tex22->Draw("same");
tex33->Draw("same");
tex44->Draw("same");
tex55->Draw("same");
leg->Draw("same");
    
c1->Print("Resolution_Plot11.pdf");
}
