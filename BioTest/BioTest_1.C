#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"

/*
TF1 *fitting_Line = new TF1("fitting_Line","[0]*x+[1]",12,velocity_TH1F[Applied_Hist_1]->GetXaxis()->GetBinCenter( Maximum_Bin-1 ));//Fitting
g->Fit(fitting_Line,"R");
Double_t par[9];
fitting_Line->GetParameters(&par[0]);
Fitting_a.push_back(par[0]);
Fitting_b.push_back(par[1]);
cout << "par[0]: " << par[0] << endl;
cout << "par[1]: " << par[1] << endl;
*/
double Get_the_P_Value()
{
    
}


void BioTest_1()
{
    
    double alpha = 0.05;//Type 1 error( Like the dollar 0.05 cent you use to bet )
    
    
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);

    TF1 *Random_F = new TF1("Random_F","1",0,1);//Fitting

    double Energy_Loss_eV = Random_F->GetRandom();//Energy_Loss(eV)
    cout << "Energy_Loss_eV: " << Energy_Loss_eV << endl;

        

}//End_Main

     


