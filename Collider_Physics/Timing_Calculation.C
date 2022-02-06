#include "TMath.h"

const double e_coulomb=1.602e-19;
const double e_mass_kg=9.1e-31;
const double e_mass_eV=0.5e6;

double *crossProduct(double vect_A[], double vect_B[], double cross_P[])
{
    static double Array_1[3];
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    Array[0]=cross_P[0];Array[1]=cross_P[1];Array[2]=cross_P[2];
    return Array_1;
}

double *Accelerate(double vect_B[], double vect_V[])
{
    static double Array_2[3];
    double cross_P[3] = {0,0,0};
    double *VcrossB   = crossProduct(vect_B,vect_V,cross_P);
    Vcross[i++]=e_coulomb*Vcross[i++];
    return Array_2;
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
    double e_V       = velocity_kms(1e2);
    eventGenerator->Sphere(par[0],par[1],par[2],e_V);
    double vect_B[3]={0.,0.,5.};//Tesla
    double vect_V[3] = {par[0],par[1],par[2]};
    cout << "vect_V[0]: " << vect_V[0] << "vect_V[1]: " << vect_V[1] << "vect_V[2]: " << vect_V[2] << endl;
    double *vect_A=Accelerate(vect_B,vect_V);
    cout << "vect_A: " << vect_A[0] << endl;
    cout << "vect_A: " << vect_A[1] << endl;
    cout << "vect_A: " << vect_A[2] << endl;

   // crossProduct(vect_A,vect_B,cross_P);
}
    

