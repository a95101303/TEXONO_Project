#include "TMath.h"

const double e_coulomb=1.602e-19;
const double e_mass_kg=9.1e-31;
const double e_mass_eV=0.5e6;

double crossProduct(double vect_A[], double vect_B[], double cross_P[])
{
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    cout << "cross_P[0]: " << cross_P[0] << endl;
    cout << "cross_P[1]: " << cross_P[1] << endl;
    cout << "cross_P[2]: " << cross_P[2] << endl;

    return 0;
}

double Accelerate(double vect_B[], double cross_P[])
{
    double Force = e_coulomb*(vect_B,cross_P);
    double
}

double velocity_kms(double E_eV)//Energy(eV)
{
    return 3e8*sqrt((2*E_eV)/(e_mass_eV));
}

void Timing_Calculation()
{
    //Initiate the function
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);
    Double_t par[3];
    TRandom *eventGenerator = new TRandom(0);//
    eventGenerator->GetSeed();
    //Two scales of energies(keV,MeV) for electron
    double vect_B[3]={0.,0.,5.};
    double cross_P[3]={0,0,0};
    double e_V       = velocity_kms(1e2);
    eventGenerator->Sphere(par[0],par[1],par[2],e_V);
    double P_V
    
    crossProduct(vect_A,vect_B,cross_P);
    
}
    

