#include <iostream>
#include <fstream>
#include "TH2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include <string>
#include <cstring>

std::string lineData = "";
using namespace std;
const double alpha = 1/137.;
const double m_e = 0.5e6; 		//eV
const double light_c=3.e8;
const double N_Avo = 6.02e23;			//N_avo

const int num_q_bins = 900;
const int num_er_bins = 500;
double er_binsize = 0.1; 			//eV
double q_binsize = 0.02;			//alpha m_e [eV]

double v_min_esc_E = 784*1.e3/light_c;
double v_min_typical = 300*1.e3/light_c;

char ExpName[25] = "SNOLAB_EF_Ge";
char FileName[3] = "LR";

double E_th = 1.1;			//eV
double depth = 2000*100;		//cm

const int MaxNo_massDM = 26;

double massDM[MaxNo_massDM] =  {0.06, 0.07,  0.08, 0.09,  0.10, 0.12, 0.20,
			0.25, 0.30,  0.40, 0.50,  0.60, 0.70, 0.80,
			0.90, 1.00,  1.30, 1.50,  2.00, 2.50, 3.00,
			3.50, 4.00,  5.00, 10.00, 20.00};

#include "crystal_f2.h"

double plotf2(int er_rescale,int q_rescale,double m_DM){

	char in_fileName[300];
	sprintf(in_fileName,"/home/alex/Desktop/Code/Earth_Effect/AnalyticalMethod/C_Ge.dat");
	std::ifstream file(in_fileName);
	getline(file, lineData); // Get first line.
	
	
	double f2[num_er_bins][num_q_bins];
	int er_bin=1, q_bin=1;
	double Ee, q; 				//eV
	double dEe = er_binsize;
	double dq = q_binsize*(m_e*alpha);	//eV
	double f2_data;
	double F2_DM=0 ;
	gStyle->SetPalette(kSunset);
	TColor::InvertPalette();
	TH2D* h2d_f2 = new TH2D("","",num_er_bins,0,num_er_bins*dEe, num_q_bins, 0, num_q_bins*dq);
		
	double S_e[num_er_bins][num_q_bins];
	double v_min[num_er_bins][num_q_bins];
	double rho =2.7;		//g cm-3
	double m_cell = 4*72./N_Avo;	//each unit cell has 4 Fe and 4 O
	double n_cell = rho/m_cell;
	double mu_e = m_e*m_DM/(m_e+m_DM);
	double Z_eff =1;

	double sigma;
	

	double v_d = sqrt(2*E_th/m_DM);
	double v_0 = 784.*1000/light_c;
	
	int v_bin =50;
	double dv = (double)abs(v_d-v_0)/v_bin;
	double v_DM;
	double v_rel =1;
	
	while(file >> f2_data){
		h2d_f2->SetBinContent(er_bin,q_bin,f2_data);
		q_bin++;		
		if(q_bin==num_q_bins+1) {er_bin++; q_bin = 1;}
	}
	
//	int er_rescale=7;
//	int q_rescale=11;
	
	dEe *= er_rescale;
	dq *= q_rescale;
	
	TH2D* h2d_f2_new = new TH2D("","",num_er_bins/er_rescale,0,num_er_bins/er_rescale*dEe, num_q_bins/q_rescale, 0, num_q_bins/q_rescale*dq/alpha/m_e);
	TH2D* h2d_S_e = new TH2D("","",num_er_bins/er_rescale,0,num_er_bins/er_rescale*dEe, num_q_bins/q_rescale, 0, num_q_bins/q_rescale*dq/alpha/m_e);
	
	double old_bin_value=0;
	
	int no_bin[num_er_bins/er_rescale][num_q_bins/q_rescale];
	
	for(int i=0;i<num_er_bins/er_rescale;i++){
		for(int j=0;j<num_q_bins/q_rescale;j++){
		no_bin[i][j]=0;
		v_min[i][j]=0;
		}
	}
	
	for(int i=0;i<(int)(num_er_bins/er_rescale)*(int)er_rescale;i++){
		for(int j=0;j<(int)(num_q_bins/q_rescale)*(int)q_rescale;j++){
			old_bin_value = h2d_f2_new->GetBinContent(i/er_rescale+1, j/q_rescale+1);
			h2d_f2_new->SetBinContent(i/er_rescale+1, j/q_rescale+1, old_bin_value + h2d_f2->GetBinContent(i,j));
			if(h2d_f2->GetBinContent(i,j)!=0) no_bin[i/er_rescale][j/q_rescale]++;
		}	
	}
	
	for(int i=0;i<num_er_bins/er_rescale;i++){
		for(int j=0;j<num_q_bins/q_rescale;j++){
				old_bin_value = h2d_f2_new->GetBinContent(i+1,j+1);
				if(no_bin[i][j]!=0) h2d_f2_new->SetBinContent(i+1,j+1,old_bin_value/no_bin[i][j]);
				else h2d_f2_new->SetBinContent(i+1,j+1,0);
				
		}
	}

	for(int i=0;i<num_er_bins/er_rescale;i++){
		for(int j=0;j<num_q_bins/q_rescale;j++){
		

		Ee = (i+0.5)*dEe;
		q = (j+0.5)*dq;

		if(strcmp(FileName, "SR")==0) F2_DM = 1;
		if(strcmp(FileName, "LR")==0) F2_DM = pow(alpha*m_e/q,4);		//F2= 1 / pow(alpha*m_e/q,4);	
		
		h2d_S_e->SetBinContent(i+1,j+1,h2d_f2_new->GetBinContent(i+1,j+1) * F2_DM * dq/pow(q,2) *dEe*Ee);
//		else h2d_S_e->SetBinContent(i,j,0);
		
//		cerr<<i<<' '<<j<<' '<<h2d_f2_new->GetBinContent(i+1,j+1)<<' '<<h2d_S_e->GetBinContent(i+1,j+1)<<endl;
		
		v_min[i][j] = Ee/q + q / (2*m_DM);
		
		}

	}
	
	
	double sum[v_bin];
	double total_sum=0;
	
	for(int i=0;i<v_bin;i++){
		v_DM = min(v_d,v_0)+(i+0.5)*dv;
		sum[i]=0;
		for(int j=0;j<num_er_bins/er_rescale;j++){
			for(int k=0;k<num_q_bins/q_rescale;k++){
				if(v_DM > v_min[j][k]){		
					sum[i] += (n_cell*pow(m_e,2)*alpha) / (pow(mu_e,2)*(v_rel)) * h2d_S_e->GetBinContent(j+1,k+1);
//					sum[i] += (n_cell*pow(m_e,2)*alpha) / (mu_e*Z_eff*alpha*m_e) * h2d_S_e->GetBinContent(j,k);
				}		
			}		
		}
//		cerr<<i<<' '<<v_DM<<' '<<sum[i]<<endl;
//		if(sum[i]!=0) total_sum += m_DM * pow(v_DM,2) * dv/sum[i];
	}
	
	for(int i=0;i<v_bin;i++){
		v_DM = min(v_d,v_0)+(i+0.5)*dv;
		
		if(sum[i]!=0) total_sum += m_DM * pow(v_DM,2) * dv/sum[i]; //pow(v_DM,2) because denominator has a v_DM
		cerr<<m_DM/1.e9<<' '<<v_DM*3e8/1000. <<' '<<sum[i]/v_DM<<' '<<endl;	
	}
	
	//cerr<<(n_cell*pow(m_e,2)*alpha) / (pow(mu_e,2)*v_rel)<<endl;
//	cerr<<total_sum<<endl; 
	sigma = total_sum/depth;
//	TCanvas* c1= new TCanvas("c1");
//	c1->SetLogz();
//	h2d_f2_new->Draw("COLZ");
//	h2d_f2_new->SetMinimum(1e-3);
//	h2d_f2_new->SetMaximum(1e2);
	cerr<<m_DM<<' '<<sigma<<endl;
	delete h2d_f2 ;
	delete h2d_f2_new;
	delete h2d_S_e;
return sigma;
}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
void get_sigma(){
char out_fileName[100];
TFile *f1=NULL;
sprintf(out_fileName,"/home/alex/Desktop/Code/Earth_Effect/AnalyticalMethod/%s%s.root",ExpName,FileName); //extract data from tree
f1 = new TFile( out_fileName, "RECREATE"); 

double mass, sigma;
TTree tr("tr","a Tree with upper limit by earth effect");
tr.Branch("mass",&mass,"mass/D"); 			
tr.Branch("sigma",&sigma,"sigma/D");

	for(int i=0;i<4;i++){
		mass = 1e6*pow(10,i);
//		mass = massDM[i] * 1e9;
		sigma = plotf2(10,1,mass);
		tr.Fill();
	}
	
f1->cd();
tr.SetMarkerStyle(4);
tr.SetLineWidth(3);
tr.Write();
f1->Write();
f1->Close();
}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//TH2D* f2_hist(){

//	char in_fileName[300];
//	sprintf(in_fileName,"/home/alex/Desktop/Code/Earth_Effect/AnalyticalMethod/C.dat");
//	std::ifstream file(in_fileName);
//	getline(file, lineData); // Get first line.
//	
//	int er_bin=1, q_bin=1;
//	double Ee, q; 				//eV
//	double dEe = er_binsize;
//	double dq = q_binsize*(m_e*alpha);	//eV
//	double f2_data;
//	double F2_DM=0 ;
//	gStyle->SetPalette(kSunset);
//	TColor::InvertPalette();
//	TH2D* h2d_f2 = new TH2D("","",num_er_bins,0,num_er_bins*dEe, num_q_bins, 0, num_q_bins*dq);
//	
//	while(file >> f2_data){
//		h2d_f2->SetBinContent(er_bin,q_bin,f2_data);
//		q_bin++;		
//		if(q_bin==num_q_bins+1) {er_bin++; q_bin = 1;}
//	}
//h2d_f2->Draw("COLZ");
//return h2d_f2;
//}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////

TGraph* f_diff(double m_DM){

double sigma = 1e-40;
double m_e = 0.5e6;
double mu_e = m_e*m_DM/(m_e+m_DM);
double v_rel = 1./137;
double v_DM , v_min;
int num_v_bins=100;
v_DM = 784.*1000./3e8;
double dv = (double)v_DM/num_v_bins;
double dEe = er_binsize*10.;
double q;
double dq = q_binsize*10.*(m_e*alpha);	//eV
TH2D* h2 = new TH2D();
h2 = rescale_f2_hist(10,10);

double E;
double sum_E=0;

double m_cell = 4*72./ 6.02e23;	//each unit cell has 4 Fe and 4 O
double n_cell = 2.7/m_cell;

TGraph* g1 = new TGraph(num_v_bins);
for(int k=0;k<num_v_bins;k++){
v_DM = k*dv;
cerr<<v_DM*3.e8/1000.<<' '<<0.5*m_DM*pow(v_DM,2)<<endl;
for(int i=1;i<num_er_bins/10.;i++){
		E = dEe*i;
		sum_E=0;
		for(int j=1;j<num_q_bins/10.;j++){
			q = j*dq;
			v_min = E/q;
			
			if(v_DM > v_min && q !=0 && E!=0) sum_E += E*dEe * dq/pow(q,2) * h2->GetBinContent(i,j) * n_cell * sigma * alpha *pow(m_e,2)/(pow(mu_e,2)*v_DM*v_rel);
		}
	}
g1->SetPoint(k,v_DM*3.e8/1000.,sum_E);
}


g1->Print();
g1->Draw("APL");
return g1;

}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
TF1* f_diff_free(double m_DM){

double sigma = 1e-40;
double m_e = 0.5e6;
double mu_e = m_e*m_DM/(m_e+m_DM);
double v_DM = 800*1000./3e8;

TF1* f1 = new TF1("f1","[0]/x",1, 10000);
f1->SetParameter(0,m_DM*sigma/(2*mu_e*mu_e*v_DM));

f1->Draw("same");

return f1;
}
