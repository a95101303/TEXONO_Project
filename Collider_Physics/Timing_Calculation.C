#include "TMath.h"

const double e_coulomb=1.602e-19;
const double e_mass_kg=9.1e-31;
const double e_mass_eV=0.5e6;


double Scalar(double vect_A[]){return sqrt(vect_A[0]*vect_A[0]+vect_A[1]*vect_A[1]+vect_A[2]*vect_A[2]);}

double *crossProduct(double vect_A[], double vect_B[], double cross_P[])
{
    static double Array_1[3];
    cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
    for(int axis=0; axis<3; axis++)
    {
        Array_1[axis]=cross_P[axis];
    }
    return Array_1;
}

double *Accelerate(double vect_B[], double vect_V[])
{
    static double Array_2[3];
    double cross_P[3] = {0,0,0};
    double *VcrossB   = crossProduct(vect_B,vect_V,cross_P);
    for(int axis=0; axis<3; axis++)
    {
        VcrossB[axis]=VcrossB[axis]*(e_coulomb/e_mass_kg);
        Array_2[axis]=VcrossB[axis];
    }
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
    double e_V       = velocity_kms(5e4);
    cout << "e_V(c): " << e_V/3e8 << endl;
    eventGenerator->Sphere(par[0],par[1],par[2],e_V);
    double vect_B[3]={0.,0.,5.};//Tesla
    double vect_V[3] = {par[0],par[1],par[2]};
    cout << "vect_V[0]: " << vect_V[0] << "vect_V[1]: " << vect_V[1] << "vect_V[2]: " << vect_V[2] << endl;
    double *vect_A=Accelerate(vect_B,vect_V);
    double  Scal_A=Scalar(vect_A);
    double  Radius=e_V*e_V/(Scal_A);
    cout << "Radius: " << Radius << endl;
    //======================Run the loop======================
    double Time_End=8e-9;//(s)
    int  Bin_Number=1000;
    double Tx=0;
    /*
    for(int Bin=0; Bin<Bin_Number; Bin++)
    {
        double Time=Time_End*(Bin+1.)/double(Bin_Number);
        
    }
     */
    
}
    

