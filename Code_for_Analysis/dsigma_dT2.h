#include "math.h"
#ifndef INCLUDE_CONSTANT_H
#define INCLUDE_CONSTANT_H
#include "constant.h"
#include "cjpl_dome.h"
//#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2_Bent.h"
#include <iomanip>      // std::setprecision
 #include "B_L_Henke_data_PE_f1_f2.h"

#endif

//const double Minimum_Velocity[17]={30.6725,42.9415,45.3161,48.9865,51.2527,55.2104,60.3555,67.4794,77.7696,95.1836,134.761,141.885,150.592,161.674,174.734,190.169,0};
//Energy threshold Brem >0.01

const double Weighted_Atomic_Number=(AFe*Fe_Percentage+AO*O_Percentage+ASi*Si_Percentage);
const double Weighted_Atomic_Number_Cement=(AH2O*H2O_Percentage+ACa*Others);
const double Weighted_Atomic_Number_Shielding = (APb*0.3+AFe*0.1+AB*0.5+ACu*0.1);
double Ground_Level=0.01;//0.01(km)==10(m)
double Ground_Level_m=10;//10m

//===================Possion_YES_OR_NO=====================//
double Possion_GetRandom(int Random_or_Possibility,double Lamda)
{
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);
    double Random_Value=0;
    TF1 *f3 = new TF1("f3","TMath::PoissonI(x,[0])",0,5);
    f3->SetParameter(0,Lamda);
    Random_Value = f3->GetRandom();
    //cout << "f3->GetRandom()" << Random_Value << endl;
    //cout << "Out of 1: " << 1-TMath::PoissonI(0,Lamda)- TMath::PoissonI(1,Lamda)<< endl;
    if(Random_or_Possibility==0)return Random_Value;
    if(Random_or_Possibility==1)return TMath::PoissonI(0,Lamda);
    if(Random_or_Possibility==2)return TMath::PoissonI(1,Lamda);
    if(Random_or_Possibility==3)return TMath::PoissonI(2,Lamda);
}
int Possion_GetRandom_Full(double Lamda)
{
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);
    double Random_Value=0;
    TF1 *f3 = new TF1("f3","TMath::PoissonI(x,[0])",0,5);
    f3->SetParameter(0,Lamda);
    int Smaller_than_1=0;
    
    while(f3->GetRandom()<1)
    {
        Smaller_than_1 = Smaller_than_1 + 1;
    }
    
    
    return Smaller_than_1;
}
//====For MPA calculation===//
//Photon Energy(PE):eV
double Lamda_of_Photon(double PE)
{//E=hv, h=4.1*1e15(eV.s),v=c/lamda, E=h*(c/lamda), lamda = (h*c/E)
    double Lamda = (h_eV_s*3e8)/(PE);//m
    return Lamda*1e2;//cm
}
double Cross_Section_MPA(double PE, double f2)
{
    return 2*r_e_cm*Lamda_of_Photon(PE)*f2;
}
/*
//====For Migdal Effect===//
  const int interp = 10;
  const int n_point1 = 251*interp;
  const int dotnumber = 1e5;
  const int n_point = 251;
*/
//Set up
const double transvalue = 3.e08*0.85/(2.*TMath::Pi());//(2.*TMath::Pi())*0.8;/////////////////////
const int lnumber = 8;
const double Enl_Series[lnumber] = {11,1.4,1.2,0.17,0.12,0.035,0.015,0.0065};  //energy nl in keV

//====For Brem===//

double Atomic_Scattering_Factor(int Index)
{//f1(PE)+if2(PE)=f(PE),[abs(f)]^2=[(f1)^2+(f2)^2]
    double f1=f_henke[Index][1];
    double f2=f_henke[Index][2];
    return 4*f1*f1+4*f2*f2;
}

//for Lindhard k = 0.22
const double p0_ = 0.888572;
const double p1_ = 1.35418;

/*
double quenching(double x)
{
  double f;
  double ER = x;
  f = TMath::Exp((TMath::Log(ER)-p1_)/p0_)/ER;

  return f;
}
*/
double antiquenching_Def(double x)  //input f*ER, output ER
{
    

  double f;
    double lnf = p0_*TMath::Log(x)+p1_;  // x[0] = f*ER
    f = TMath::Exp(lnf);
  return f;
}

//=====For Nuclear Recoil=====//
 const double alpha = 1.0;
 const double qf0 = alpha*0.19816;
 const double qf1 = alpha*0.05052;
 const double qf2 = alpha*0.00378;
 const double qf3 = alpha*0.00192;
 const double qf4 = alpha*0.0016;
//keV
//T:keV
//QF of Ge
double QF(double T)
{ return (qf0+qf1*log10(T)+qf2*pow(log10(T),2)+qf3*pow(log10(T),3)+qf4*pow(log10(T),4)); }
//keV
//T:keV
double TQF(double T)
{ return T*QF(T); }

bool double_equals(double a, double b, double epsilon = 0.001)
{
    return std::abs(a - b) < epsilon;
}
double CNFNV(int A_P, double N)
{
    double Number[2]={0,0};//The first one is for the value < 10, other one is for the scale
    for(int kkk=0; kkk<45; kkk++)
    {//1
        if(Number[0]!=0)break;
        else if( double_equals(((TMath::Power(10,kkk))*N),1.0,0.00001) )
        {
            Number[0]=1.0;
            Number[1]=1.0*(kkk);
        }

        else if( ((TMath::Power(10,kkk))*N) < 1.0)continue;
        else
        {//2
            Number[0]=((TMath::Power(10,kkk))*N);
            Number[1]=1.0*(kkk);
        }//2
    }//1
    return Number[A_P];
}

double Length_of_vector(double A=0, double B=0, double C=0)
{
    return sqrt( A*A + B*B + C*C )  ;
}

double Binomial_Error(double N, double P)
{
    return sqrt( P*(1-P)/N );
}

double Volume_for_ball(double Radius)//Radius(km)
{
    return (4/3)*TMath::Pi()*TMath::Power(Radius*1e5,3);//(cm^3)
}

double Two_points_difference(double *The_first_point, double *The_second_point)
{
    double X_Difference=(The_first_point[0]-The_second_point[0]);
    double Y_Difference=(The_first_point[1]-The_second_point[1]);
    double Z_Difference=(The_first_point[2]-The_second_point[2]);
    double Total_Difference = sqrt(TMath::Power(X_Difference,2)+TMath::Power(Y_Difference,2)+TMath::Power(Z_Difference,2));
    return Total_Difference;
}
double Solution_for_Linear_Equation(int Sign, double a, double b, double c) // a(X^2) + b(X) + c = 0
{//Check OK!
    double Denominator=0; //Down
    double Numerator=0; //Up
    double Intermedium=0;
    Denominator = 2*a;
    if(Sign==1)Numerator = -b + sqrt(b*b-4*a*c);
    if(Sign==2)Numerator = -b - sqrt(b*b-4*a*c);
    
     //cout << "2*a: " << 2*a << endl;
     //cout << "Numerator: " << Numerator << endl;
     //cout << "(Numerator/Denominator): " << (Numerator/Denominator) << endl;
    return (Numerator/Denominator);
}

double Choose_solution_for_Linear_Equation(double a, double b, double c) // a(X^2) + b(X) + c = 0
{//Check OK!
    double A = Solution_for_Linear_Equation(1, a, b, c);
    double B = Solution_for_Linear_Equation(2, a, b, c);
    double Final_Choose=0;
    if(b*b-4*a*c>=0 and b*b+4*a*c>=0)
    {
        TF1 *f1 = new TF1("f1","1",-1,1);
        double Random = f1->GetRandom();
        if(Random>0)  Final_Choose = A;
        if(Random<=0) Final_Choose = B;
    }
    if(b*b+4*a*c>=0 and b*b-4*a*c<0) {Final_Choose = A;}
    if(b*b+4*a*c< 0 and b*b-4*a*c>=0){Final_Choose = B;}
    else{Final_Choose = 0;}
    
    return Final_Choose;
}
double Scale_for_scattering_process(int air_or_earth, double Layer_Place, double *Track_Point, double *Direction)
{
    double A = 1;
    double B = 2*Track_Point[0]*Direction[0]+2*Track_Point[1]*Direction[1]+2*Track_Point[2]*Direction[2];
    double C = -Layer_Place*Layer_Place+Track_Point[0]*Track_Point[0]+Track_Point[1]*Track_Point[1]+Track_Point[2]*Track_Point[2];
    double S1 = Solution_for_Linear_Equation(1,A,B,C);
    double S2 = Solution_for_Linear_Equation(2,A,B,C);
    cout << "S1: " << S1 << endl;
    cout << "S2: " << S2 << endl;
    if(S1<(1e-5) and S1>-(1e-5)) S1=0;
    if(S2<(1e-5) and S2>-(1e-5)) S2=0;
    if(B*B-4*A*C<0){S1=0;S2=0;}
    
    if(air_or_earth==0){
    if(abs(S1)<=1000)return S1;
    if(abs(S2)<=1000)return S2;}
    
    if(air_or_earth==1){
        if(S1>0 and S2>0){return min(S1,S2);}
        else if(S1>0 and S2==0){return S1;}
        else if(S1==0 and S2>0){return S2;}
        else{return max(S1,S2);}
    }
    
}
double Selection_for_STS(double A, double B)
{
    if(A>0 and B>0){return min(A,B);}
    else if(A>0 and B==0){return A;}
    else if(A==0 and B>0){return B;}
    else{return max(A,B);}
}

double *Vertical_vector_Function(double *Vector_Norm)//
{
    double Length = sqrt(Square_Sum(Vector_Norm,Vector_Norm));
    double NC[3]={Vector_Norm[0]/Length,Vector_Norm[1]/Length,Vector_Norm[2]/Length};
    
    //cout << "NC[0]*NC[0]+NC[1]*NC[1]+NC[2]*NC[2]: " << NC[0]*NC[0]+NC[1]*NC[1]+NC[2]*NC[2] << endl;
    TF1 *f1 = new TF1("f1","1",-1,1);
    double Vy =0;double a=0;double b=0;double c=0;
    while(b*b-4*a*c<=0)
    {
        Vy = f1->GetRandom();
        a = (NC[2]/NC[0])*(NC[2]/NC[0])+1;
        b = 2*Vy*(NC[1]/NC[0])*(NC[2]/NC[0]);
        c = -( 1-Vy*Vy*( 1 + (NC[1]/NC[0])*(NC[1]/NC[0]) ) );
    }
    //cout << "Vy: " << Vy << endl;
    
    double Vz_1 = Solution_for_Linear_Equation(1,a,b,c);
    double Vz_2 = Solution_for_Linear_Equation(2,a,b,c);
    //cout << "Vz_1: " << Vz_1 << endl;cout << "Vz_2: " << Vz_2 << endl;
    //cout << "Vz_1*Vz_1*a+Vz_1*b+c=0 : " << Vz_1*Vz_1*a+Vz_1*b+c << endl;
    //cout << "Vz_2*Vz_2*a+Vz_2*b+c=0 : " << Vz_2*Vz_2*a+Vz_2*b+c << endl;

    double Vz = 0;
    double Random_Vz=f1->GetRandom();
    if(Random_Vz>0) Vz=Vz_1;
    if(Random_Vz<=0)Vz=Vz_2;

    //cout << "Vz: " << Vz << endl;

    double Vx_1 =  sqrt(1-Vy*Vy-Vz*Vz);
    //cout << "Vx_1: " << Vx_1 << endl;
    //cout << "Check: " << -Vy*(NC[1]/NC[0])-Vz*(NC[2]/NC[0]) << endl;
    double Vx_2 = -sqrt(1-Vy*Vy-Vz*Vz);
    double Vx = -Vy*(NC[1]/NC[0])-Vz*(NC[2]/NC[0]);
    /*
    double Random_Vx=f1->GetRandom();
    if(Random_Vx>0) Vx=Vx_1;
    if(Random_Vx<=0)Vx=Vx_2;
     */
    //cout << "Vx: " << Vx << endl;

    //cout << "Check1==1 " << Vx*Vx+Vy*Vy+Vz*Vz << endl;
    //cout << "Check2==0 " << NC[0]*Vx+NC[1]*Vy+NC[2]*Vz << endl;

    static double V[3];
    V[0]=Vx;V[1]=Vy;V[2]=Vz;
    
    return V;
    /*
    TF1 *f1 = new TF1("f1","1",-1,1);
    double Vx = f1->GetRandom();
    double Constraint = sqrt(1-Vx*Vx);
    TF1 *f2 = new TF1("f2","1",-Constraint,Constraint);
    double Vy = f2->GetRandom();
    double Vz = -(Normalized_vector[0]/Normalized_vector[0])*Vx-(Normalized_vector[1]/Normalized_vector[0])*Vy

    cout << "Vx*Vx+Vy*Vy+Vz*Vz" << Vx*Vx+Vy*Vy+Vz*Vz << endl;
    cout << "Vector_Norm" << Vx*Normalized_vector[0]+Vy*Normalized_vector[1]+Vz*Normalized_vector[2] << endl;
     */
}
double *Aft_scatterd_Direction(int Earth_or_Air, double DM_mx, double DM_Velocity_Aft_Colliding, double *Direction_Bef, double Ratio_of_Energy_Loss_to_Atom)
{
    static double Direction[3];double Direction_VT[3];
    
    TF1 *Angle_CM = new TF1("f1","1",0,360);
    double Angle_Lab     = Phi_Rest_Frame(Earth_or_Air,DM_mx,Ratio_of_Energy_Loss_to_Atom);
    cout << "Angle_Lab: " << Angle_Lab << endl;
    
    double MDV = DM_Velocity_Aft_Colliding*TMath::Cos(Angle_Lab*Degree_to_Radian);//  Mother_direction_Value
    double VDV = DM_Velocity_Aft_Colliding*TMath::Sin(Angle_Lab*Degree_to_Radian);//Vertical_direction_Value
    
    double *Vertical_vector=Vertical_vector_Function(Direction_Bef);
    
    Direction[0]    =MDV*Direction_Bef[0];     Direction[1]   =MDV*Direction_Bef[1];  Direction[2]    =MDV*Direction_Bef[2];
    Direction_VT[0] =VDV*Vertical_vector[0];Direction_VT[1] =VDV*Vertical_vector[1];Direction_VT[2] =VDV*Vertical_vector[2];

    Direction[0]=Direction[0]+Direction_VT[0];
    Direction[1]=Direction[1]+Direction_VT[1];
    Direction[2]=Direction[2]+Direction_VT[2];
    
    cout << sqrt(Direction[0]*Direction[0]+Direction[1]*Direction[1]+Direction[2]*Direction[2]) << "==" << DM_Velocity_Aft_Colliding << "Yes?" << endl;

    double Length_Scattering = sqrt(Square_Sum(Direction,Direction));
    Direction[0] = Direction[0]/Length_Scattering;
    Direction[1] = Direction[1]/Length_Scattering;
    Direction[2] = Direction[2]/Length_Scattering;

    cout <<  sqrt(Direction[0]*Direction[0]+Direction[1]*Direction[1]+Direction[2]*Direction[2]) << "~1 Yes?" << endl;
    
    return Direction;

}
double Check_PREM_density()//Past has already checked!
{
    double Volume_For_Different_layers[12];
    double Density_for_Different_Layer[12];
    double Mass=0;
    for(int kkk=0; kkk<13;kkk++)
    {
        Volume_For_Different_layers[kkk]= Volume_for_ball(earth_table[kkk+1][0])-Volume_for_ball(earth_table[kkk][0]);
        Density_for_Different_Layer[kkk]=(earth_table[kkk+1][1]+earth_table[kkk][1])/2;
        cout << "Volume_For_Different_layers[kkk]: " << Volume_For_Different_layers[kkk] << endl;
        cout << "Density_for_Different_Layer[kkk]: " << Density_for_Different_Layer[kkk] << endl;
        cout << "earth_table[kkk][0]: " << earth_table[kkk+1][0] << endl;
        Mass = Mass + Volume_For_Different_layers[kkk]*Density_for_Different_Layer[kkk];
    }
    return (Mass/Volume_for_ball(6371.0));
}

//============================================For the real case ( Our experiment )=========================================
double KS_Earth_Path_Length(double Path_Vector_Unit_Y, double Path_Vector_Unit_Z,double Deep_of_Surface_Level,double Previous_Length)//OK!Checked 2020/02/02
{//   abs(A(Xbar,Ybar,Zbar)+(0,-0.029,6370.96))=(a) a depends on the layer of the earth
    double PY=Path_Vector_Unit_Y;double PZ=Path_Vector_Unit_Z;
    double Ground_Level=0.01;//0.01(km)==10(m)
    double Lab_Cor[3]={0,0.028,0};double Central_of_Earth[3]={0,0,Earth_Radius-Ground_Level};//LC and CE
    double LCtoCE[3]={Central_of_Earth[0]-Lab_Cor[0],Central_of_Earth[1]-Lab_Cor[1],Central_of_Earth[2]-Lab_Cor[2]};//(0,-0.029,r-0.01)
    
    
    double Solution_for_it=0;double Solution_for_it_1=0;double Solution_for_it_2=0;double Solution_for_it_3=0;
    Solution_for_it_1=Solution_for_Linear_Equation(1, 1, 2*PY*LCtoCE[1]+2*PZ*LCtoCE[2], TMath::Power(LCtoCE[1],2)+TMath::Power(LCtoCE[2],2)-TMath::Power(Deep_of_Surface_Level,2) );
    Solution_for_it_2=Solution_for_Linear_Equation(2, 1, 2*PY*LCtoCE[1]+2*PZ*LCtoCE[2], TMath::Power(LCtoCE[1],2)+TMath::Power(LCtoCE[2],2)-TMath::Power(Deep_of_Surface_Level,2) );
    Solution_for_it_3=Solution_for_Linear_Equation(1, 1, 2*(-PY)*LCtoCE[1]+2*(-PZ)*LCtoCE[2], TMath::Power(LCtoCE[1],2)+TMath::Power(LCtoCE[2],2)-TMath::Power(Deep_of_Surface_Level,2) );

    
    //Two parts of the length
    //Previous Length: From the inner to the outer length
    if( PZ < 0 )
    {
        if(Deep_of_Surface_Level<=6368){
            Solution_for_it = Solution_for_it_1 - Solution_for_it_2 - Previous_Length;}
        if(Deep_of_Surface_Level>6368){
            Solution_for_it = Solution_for_it_1 - Solution_for_it_2 - Solution_for_it_3 - Previous_Length;}
    }
    if( PZ >= 0 )
    {
        //cout << "Solution_for_it: " << Solution_for_it << endl;
        Solution_for_it = Solution_for_it_1;
    }
    return Solution_for_it;//km
}

//============================================For the real case ( Our experiment )=========================================
double KS_Earth_Path_Length_for_STS(int Option, double Path_Vector_Unit_Y, double Path_Vector_Unit_Z,double Deep_of_Surface_Level,double Previous_Length)
{//   abs(A(Xbar,Ybar,Zbar)+(0,-0.029,6370.96))=(a) a depends on the layer of the earth
    double PY=Path_Vector_Unit_Y;double PZ=Path_Vector_Unit_Z;
    double Ground_Level=0.01;//0.01(km)==10(m)
    double Lab_Cor[3]={0,0.028,0};double Central_of_Earth[3]={0,0,Earth_Radius-Ground_Level};//LC and CE
    double LCtoCE[3]={Central_of_Earth[0]-Lab_Cor[0],Central_of_Earth[1]-Lab_Cor[1],Central_of_Earth[2]-Lab_Cor[2]};//(0,-0.029,r-0.01)
    //cout << "LCtoCE[0]: "<< LCtoCE[0] << endl;cout << "LCtoCE[1]: "<< LCtoCE[1] << endl;cout << "LCtoCE[2]: "<< LCtoCE[2] << endl;
    
    
    double Solution_for_it=0;double Solution_for_it_1=0;double Solution_for_it_2=0;double Solution_for_it_3=0;
    Solution_for_it_1=Solution_for_Linear_Equation(1, 1, 2*PY*LCtoCE[1]+2*PZ*LCtoCE[2], TMath::Power(LCtoCE[1],2)+TMath::Power(LCtoCE[2],2)-TMath::Power(Deep_of_Surface_Level,2) );
    Solution_for_it_2=Solution_for_Linear_Equation(2, 1, 2*PY*LCtoCE[1]+2*PZ*LCtoCE[2], TMath::Power(LCtoCE[1],2)+TMath::Power(LCtoCE[2],2)-TMath::Power(Deep_of_Surface_Level,2) );
    Solution_for_it_3=Solution_for_Linear_Equation(1, 1, 2*(-PY)*LCtoCE[1]+2*(-PZ)*LCtoCE[2], TMath::Power(LCtoCE[1],2)+TMath::Power(LCtoCE[2],2)-TMath::Power(Deep_of_Surface_Level,2) );
    //cout << "Solution_for_it_1: " << Solution_for_it_1 << endl;
    //cout << "Solution_for_it_2: " << Solution_for_it_2 << endl;
    //cout << "Solution_for_it_3: " << Solution_for_it_3 << endl;

    //Two parts of the length
    if( PZ < 0 )
    {
        if(Deep_of_Surface_Level<=6368){
            Solution_for_it = Solution_for_it_1 - Solution_for_it_2 - Previous_Length;}
        if(Deep_of_Surface_Level>6368){
            //cout << "Solution_for_it_3: " << Solution_for_it_3 << endl;
            Solution_for_it = Solution_for_it_1 - Solution_for_it_2 - Previous_Length;}
    }
    if( PZ >= 0 )
    {
        //cout << "Solution_for_it: " << Solution_for_it << endl;
        Solution_for_it = Solution_for_it_1;
    }
    if(Option==0)return Solution_for_it;//km
    if(Option==1)return Solution_for_it_3;//km

}

double KS_Air_Path_Length(double Path_Vector_Unit_Y, double Path_Vector_Unit_Z,double High_of_Sea_Level,double Previous_Length)//OK!Checked 2020/02/02
{//   abs(A(Xbar,Ybar,Zbar)+(0,-0.029,6370.96))=(r+a) a depends on the layer of the air
    double PY=Path_Vector_Unit_Y;double PZ=Path_Vector_Unit_Z;
    double Lab_Cor[3]={0,0.028,0};double Central_of_Earth[3]={0,0,Earth_Radius-Ground_Level};//LC and CE
    double LCtoCE[3]={Central_of_Earth[0]-Lab_Cor[0],Central_of_Earth[1]-Lab_Cor[1],Central_of_Earth[2]-Lab_Cor[2]};//(0,-0.028,r-0.01)
    
    double Solution_for_it=0;Solution_for_it=Solution_for_Linear_Equation(1, 1, 2*PY*LCtoCE[1]+2*PZ*LCtoCE[2], TMath::Power(LCtoCE[1],2)+TMath::Power(LCtoCE[2],2)-TMath::Power(Earth_Radius+High_of_Sea_Level,2) );
    
    Solution_for_it = Solution_for_it - Previous_Length;
    return Solution_for_it;
}


double KS_Reactor_Path_Length(double Rea_R,double Path_Vector_Unit_X, double Path_Vector_Unit_Y, double Path_Vector_Unit_Z)//OK
{//Rea_R=26.5m (Water) Rea_R=27.5m (Total)
    if(Path_Vector_Unit_Z>=0)
    {
    double PX=Path_Vector_Unit_X;double PY=Path_Vector_Unit_Y;double PZ=Path_Vector_Unit_Z;
    double The_high_of_the_upper2=50;double Lab_Cor[3]={0,28,0};double Particle_Outward[3]={-PX,-PY,-PZ};
    double The_First_Point_on_Cylinder[3];double The_Second_Point_on_Cylinder[3];

    double CT_Point_on_the_Cylinder[2];
    for(int kkk=0; kkk<2; kkk++)
    {
        CT_Point_on_the_Cylinder[kkk]=Solution_for_Linear_Equation(kkk+1,(-PX)*(-PX)+(-PY)*(-PY),2*(-PY)*Lab_Cor[1],Lab_Cor[1]*Lab_Cor[1]-Rea_R*Rea_R);//+
        if(CT_Point_on_the_Cylinder[kkk]*Particle_Outward[2]<=-(The_high_of_the_upper2))
        {
            CT_Point_on_the_Cylinder[kkk] = (The_high_of_the_upper2)/(Path_Vector_Unit_Z);
        }
    }
    if(CT_Point_on_the_Cylinder[0]>=0 and CT_Point_on_the_Cylinder[1]>=0)
    {
    The_First_Point_on_Cylinder[0] =CT_Point_on_the_Cylinder[0]*Particle_Outward[0];
    The_First_Point_on_Cylinder[1] =CT_Point_on_the_Cylinder[0]*Particle_Outward[1]+Lab_Cor[1];
    The_First_Point_on_Cylinder[2] =CT_Point_on_the_Cylinder[0]*Particle_Outward[2];
    The_Second_Point_on_Cylinder[0]=CT_Point_on_the_Cylinder[1]*Particle_Outward[0];
    The_Second_Point_on_Cylinder[1]=CT_Point_on_the_Cylinder[1]*Particle_Outward[1]+Lab_Cor[1];
    The_Second_Point_on_Cylinder[2]=CT_Point_on_the_Cylinder[1]*Particle_Outward[2];

    double Normalize=sqrt( (The_Second_Point_on_Cylinder[0]-The_First_Point_on_Cylinder[0])*(The_Second_Point_on_Cylinder[0]-The_First_Point_on_Cylinder[0]) +(The_Second_Point_on_Cylinder[1]-The_First_Point_on_Cylinder[1])*(The_Second_Point_on_Cylinder[1]-The_First_Point_on_Cylinder[1])+(The_Second_Point_on_Cylinder[2]-The_First_Point_on_Cylinder[2])*(The_Second_Point_on_Cylinder[2]-The_First_Point_on_Cylinder[2]) );
        //double Normalize=1;
        /*
        cout << "The_Second_Point_on_Cylinder[0]-The_First_Point_on_Cylinder[0]: " << (The_Second_Point_on_Cylinder[0]-The_First_Point_on_Cylinder[0])/Normalize << endl;
        cout << "The_Second_Point_on_Cylinder[1]-The_First_Point_on_Cylinder[1]: " << (The_Second_Point_on_Cylinder[1]-The_First_Point_on_Cylinder[1])/Normalize << endl;
        cout << "The_Second_Point_on_Cylinder[2]-The_First_Point_on_Cylinder[2]: " << (The_Second_Point_on_Cylinder[2]-The_First_Point_on_Cylinder[2])/Normalize << endl;
         */
    return Two_points_difference(The_First_Point_on_Cylinder,The_Second_Point_on_Cylinder)/1e3;//km
    }
    else
    {
    return 0;
    }
        
    }
    if(Path_Vector_Unit_Z<0)
    {
        return 0;
    }
}

double KS_Cement_Path_Length(double Path_Vector_Unit_X, double Path_Vector_Unit_Y, double Path_Vector_Unit_Z)//OK!Checked 2020/02/02
{//
    double PX=Path_Vector_Unit_X;double PY=Path_Vector_Unit_Y;double PZ=Path_Vector_Unit_Z;
    if(PZ >= 0){

    double The_high_of_the_upper1=70;      double The_radius_of_the_upper1=57;//(m)
    double The_high_of_the_upper2=60;     double The_radius_of_the_upper2=56;//(m)

    double CT=0;double CT3=0;//Constant Term
    double Radius_of_WIMP_upper1=0;double Radius_of_WIMP_upper2=0;//Criteria for "Upper or Side"
    double Radius_on_XYPlane_Unit=sqrt(TMath::Power(Path_Vector_Unit_X,2)+TMath::Power(Path_Vector_Unit_Y,2));
    //
    double Lab_Cor[3]={0,28,0};
    CT  = (The_high_of_the_upper1/Path_Vector_Unit_Z);
    CT3 = (The_high_of_the_upper2/Path_Vector_Unit_Z);

    Radius_of_WIMP_upper1=sqrt(TMath::Power(Path_Vector_Unit_X*CT,2)+TMath::Power(Path_Vector_Unit_Y*CT,2));
    Radius_of_WIMP_upper2=sqrt(TMath::Power(Path_Vector_Unit_X*CT3,2)+TMath::Power(Path_Vector_Unit_Y*CT3,2));

    //
        //cout << "YES!1: " << endl;
    double CT1=0;double CT2=0;double Upper1_Point[3];double Upper2_Point[3];
    if(Radius_of_WIMP_upper1<The_radius_of_the_upper1)
    {//Upper case
        //cout << "Upper_Case1" << endl;
        CT1 = (The_high_of_the_upper1/Path_Vector_Unit_Z);
        Upper1_Point[0]=-Path_Vector_Unit_X*CT1+Lab_Cor[0];
        Upper1_Point[1]=-Path_Vector_Unit_Y*CT1+Lab_Cor[1];
        Upper1_Point[2]=-Path_Vector_Unit_Z*CT1+Lab_Cor[2];
    }
    else if(Path_Vector_Unit_Z==0 or Radius_of_WIMP_upper1>=The_radius_of_the_upper1)//C^2-2*C*(Y_Unit)*28+28^2=Radius^2
    {//Side case
        //cout << "Side_Case1" << endl;
        CT1 = Solution_for_Linear_Equation(1,(-PX)*(-PX)+(-PY)*(-PY),-2*Path_Vector_Unit_Y*Lab_Cor[1],Lab_Cor[1]*Lab_Cor[1]-The_radius_of_the_upper1*The_radius_of_the_upper1);
        //cout << "CT1: " << CT1 << endl;
        Upper1_Point[0]=-Path_Vector_Unit_X*CT1+Lab_Cor[0];
        Upper1_Point[1]=-Path_Vector_Unit_Y*CT1+Lab_Cor[1];
        Upper1_Point[2]=-Path_Vector_Unit_Z*CT1+Lab_Cor[2];
        //cout << "sqrt(CT1*CT1-2CPY*28+28^2)" << sqrt(CT1*CT1-2*CT1*PY*28+28*28) << endl;
    }
   // cout << "YES!2: " << endl;

    if(Radius_of_WIMP_upper2<=The_radius_of_the_upper2)
    {//Upper case
        //cout << "Upper_Case2" << endl;
        CT2 = (The_high_of_the_upper2/Path_Vector_Unit_Z);
        Upper2_Point[0]=-Path_Vector_Unit_X*CT2+Lab_Cor[0];
        Upper2_Point[1]=-Path_Vector_Unit_Y*CT2+Lab_Cor[1];
        Upper2_Point[2]=-Path_Vector_Unit_Z*CT2+Lab_Cor[2];
    }
    else if(Path_Vector_Unit_Z==0 or Radius_of_WIMP_upper2>=The_radius_of_the_upper2)
    {//Side case
        //cout << "Side_Case2" << endl;
        CT2 = Solution_for_Linear_Equation(1,(-PX)*(-PX)+(-PY)*(-PY),-2*Path_Vector_Unit_Y*Lab_Cor[1],Lab_Cor[1]*Lab_Cor[1]-The_radius_of_the_upper2*The_radius_of_the_upper2);
        Upper2_Point[0]=-Path_Vector_Unit_X*CT2+Lab_Cor[0];
        Upper2_Point[1]=-Path_Vector_Unit_Y*CT2+Lab_Cor[1];
        Upper2_Point[2]=-Path_Vector_Unit_Z*CT2+Lab_Cor[2];
    }
    /*
    cout << "Upper1_Point[0]: " << Upper1_Point[0] << endl;
    cout << "Upper1_Point[1]: " << Upper1_Point[1] << endl;
    cout << "Upper1_Point[2]: " << Upper1_Point[2] << endl;
    cout << "Upper2_Point[0]: " << Upper2_Point[0] << endl;
    cout << "Upper2_Point[1]: " << Upper2_Point[1] << endl;
    cout << "Upper2_Point[2]: " << Upper2_Point[2] << endl;

    cout << "sqrt(Upper1_Point[0]*Upper1_Point[0]+Upper1_Point[1]*Upper1_Point[1])<=57: " << sqrt(Upper1_Point[0]*Upper1_Point[0]+Upper1_Point[1]*Upper1_Point[1]) << endl;
    cout << "sqrt(Upper2_Point[0]*Upper2_Point[0]+Upper2_Point[1]*Upper2_Point[1])<=56: " << sqrt(Upper2_Point[0]*Upper2_Point[0]+Upper2_Point[1]*Upper2_Point[1]) << endl;
     */
    return (Two_points_difference(Upper1_Point,Upper2_Point))/1e3;//(km)
    }
    else{return 0;}
}

double Bool_If_Earth(double PX, double PY, double PZ)//OK!Checked 2020/02/02
{
    //cout << "PX: " << PX << "PY: " << PY << "PZ: " << PZ << endl;
    double The_radius_of_the_upper1=57;//m
    double Constant_HIGH=0;    double YES_or_NO=0;//Pass Earth or not ( 0 means no and non-0 means yes, and the value means the length (m)
    double Lab_Cor[3]={0,28,0};
   if(PZ > 0)
   {
    Constant_HIGH = (Ground_Level_m/PZ);//Constant will be positive
    double X_Vector = Constant_HIGH*(-PX);double Y_Vector = Constant_HIGH*(-PY);
    double Criteria_of_Bool = Length_of_vector(X_Vector+Lab_Cor[0],Y_Vector+Lab_Cor[1],0);//See if the radius is out of the cement
    //cout << "Criteria_of_Bool: " << Criteria_of_Bool << " < or >= ?" << The_radius_of_the_upper1 << endl;
    if(Criteria_of_Bool<The_radius_of_the_upper1)
       {
           YES_or_NO=0;
       }
    else if(Criteria_of_Bool>=The_radius_of_the_upper1)
      {
    Constant_HIGH
    =Solution_for_Linear_Equation(1,(-PX)*(-PX)+(-PY)*(-PY),-2*PY*Lab_Cor[1],Lab_Cor[1]*Lab_Cor[1]-The_radius_of_the_upper1*The_radius_of_the_upper1);
        YES_or_NO = Length_of_vector(PX*Constant_HIGH,PY*Constant_HIGH,PZ*Constant_HIGH);
        //cout << "Length_of_vector(PX*Constant_HIGH,PY*Constant_HIGH,0): " << Length_of_vector(-PX*Constant_HIGH,-PY*Constant_HIGH+28,0) << endl;
         // cout << "YES_or_NO: " << YES_or_NO << endl;

      }
       return YES_or_NO/1e3;//km
   }
   else if(PZ == 0)
   {
       Constant_HIGH
       =Solution_for_Linear_Equation(1,(-PX)*(-PX)+(-PY)*(-PY),-2*PY*Lab_Cor[1],Lab_Cor[1]*Lab_Cor[1]-The_radius_of_the_upper1*The_radius_of_the_upper1);
        YES_or_NO = Length_of_vector(PX*Constant_HIGH,PY*Constant_HIGH,PZ*Constant_HIGH);
        //cout << "Length_of_vector(PX*Constant_HIGH,PY*Constant_HIGH,0): " << Length_of_vector(-PX*Constant_HIGH,-PY*Constant_HIGH+28,0) << endl;
        //cout << "YES_or_NO: " << YES_or_NO << endl;
       return YES_or_NO/1e3;//km
   }
   
   else{return 0;}
}



//============================================================================================================================//

double kg_perm3_to_g_percm3(double kg_perm3)
{
    return (kg_perm3)*(1e3/1e6);
}
//=====================================Lab on the surface ( Ignore the real case )=====================================
double Path_Length_Air(double Path_Vector_Unit_Z, double High_of_Sea_Level,double Path_Length_Earth_Calculated,double Previous_Length)
{// From the surface of earth
    //cout << "Previous: " << Previous_Length << endl;
    double Solution_for_it=0;
    double A=0; double B=0; double C=0;
    A = 1;
    B = 2*(Path_Vector_Unit_Z)*Earth_Radius;
    C= -(2*Earth_Radius*High_of_Sea_Level+High_of_Sea_Level*High_of_Sea_Level);
    Solution_for_it = Solution_for_Linear_Equation(1,A,B,C);
    Solution_for_it = Solution_for_it - Path_Length_Earth_Calculated - Previous_Length;
    return Solution_for_it;
}
double Path_Length_Earth(double Path_Vector_Unit_Z, double Deep_of_Earth_Level,double Previous_Length)
{// From the surface of earth
    //cout << "Previous: " << Previous_Length << endl;
    double Solution_for_it=0; double Solution_for_it_1=0;double Solution_for_it_2=0;
    double A=0; double B=0; double C=0;
    A = 1;
    B = 2*(Path_Vector_Unit_Z)*Earth_Radius;
    C=  Earth_Radius*Earth_Radius-Deep_of_Earth_Level*Deep_of_Earth_Level;
    Solution_for_it_1 = Solution_for_Linear_Equation(1,A,B,C);
    Solution_for_it_2 = Solution_for_Linear_Equation(2,A,B,C);
    //cout << "Solution_for_it_1: " << Solution_for_it_1 << endl;
    //cout << "Solution_for_it_2: " << Solution_for_it_2 << endl;
    Solution_for_it = Solution_for_it_1-Solution_for_it_2-Previous_Length;
    //cout << "Solution_for_it: " << Solution_for_it << endl;
    return Solution_for_it;
}

double Path_Length_Earth_Rough(double Path_Vector_Unit_Z)
{//Check OK!
    double Path_Length=0;
    if(Path_Vector_Unit_Z<0)
    {
        Path_Length = 2*abs(Path_Vector_Unit_Z)*Earth_Radius;
    }
    else
    {
        Path_Length=0;
    }
    return Path_Length;
}
//============================================================================================================================//
double *Array_from_SPHC_to_CARC(double r, double theta, double Phi) //r > 0, 0 <= theta < pi(), 0 <= Phi < 2*pi()
{
    static double Array_With_CARC[3];
    if( (r<=0) or (theta<0) or (theta>TMath::Pi()) or (Phi<0) or (Phi>2*TMath::Pi()) ) // Debug on the angle
    {
        cout << "Something wrong! DO CHECK!" << endl;
        cout << "r: " << r << endl;
        cout << "theta: " << theta << endl;
        cout << "Phi: " << Phi << endl;
    }
    
    Array_With_CARC[0]= r * TMath::Sin(theta) * TMath::Cos(Phi); //X
    Array_With_CARC[1]= r * TMath::Sin(theta) * TMath::Sin(Phi);//Y
    Array_With_CARC[2]= r * TMath::Cos(theta);//Z = r*cos(theta)
    return Array_With_CARC;
}


double N_atom_1kg(double Atomic_mass_g)
{
    return 1000.0/(Atomic_mass_g);
}

double Mass_Density_of_Earth_New(double Average_Density, double Percentage, double Atomic_mass)
{
    return (Average_Density*Percentage)/(Atomic_mass);
}


double Collision_Time_Earth(double Pass_Length, double Mean_Free_Path)
{
    return (Pass_Length)/(Mean_Free_Path);
}

//0.5*M_Chi*(v/c)^2=Energy mx(GeV), Velocity(beta)
double Energy_DM(double mx, double beta)
{
    return 0.5*(mx*1e6)*(beta)*(beta);//keV
}
//Inverse the Energy of DM to the velocity of DM with the formula in the previous function
double Velocity_DM(double mx, double Energy)//Outcome:[km/s]
{
    //cout << "WHATWHAT???" << endl;
    //cout << "2*Energy:" << 2*Energy << endl;
    //cout << "mx*1e6: " << mx*1e6 << endl;
    //cout << "sqrt(2*Energy/(mx*1e6)): " << sqrt(2*Energy/(mx*1e6)) << endl;
    
    return sqrt(2*Energy/(mx*1e6))*(3e8/1e3); //[km/s]
}
double Function_of_material_possibility(int Earth_or_air) // 1 is earth 2 is air
{//比例依樣，所以密度高低沒有關係，設為1
    if(Earth_or_air==1)
    {  //cout << "EARTH: " << Weighted_Atomic_Number << endl ;
        return Weighted_Atomic_Number;
    }
    else if(Earth_or_air==2)
    {   //cout << "AIR: " << 14.99 << endl;;
        return 14.99;}
    else if(Earth_or_air==3)
    {  //cout << "EARTH: ";
        return Weighted_Atomic_Number_Cement;
    }
    else if(Earth_or_air==4)
    {  //cout << "EARTH: ";
        return AFe;
    }
    else if(Earth_or_air==5)
    {  //cout << "EARTH: ";
        return AH2O;
    }
    else if(Earth_or_air==6)
    {  //cout << "EARTH: ";
        return Weighted_Atomic_Number_Shielding;
    }

    
    /*
    cout <<"HIGH?" << endl;
    double Material_Possibility_Array[3];
    double Atom_Mass[3];
    if(Earth_or_air==1)//EARTH
    {
        Atom_Mass[0]        =AFe ; Atom_Mass[1]        =AO;Atom_Mass[2]        =ASi;
        double Fe_Fraction = Mass_Density_of_Earth_New(1,Fe_Percentage,atom_mass_Fe_g);
        double O_Fraction  = Mass_Density_of_Earth_New(1,O_Percentage,atom_mass_O_g);
        double Si_Fraction = Mass_Density_of_Earth_New(1,Si_Percentage,atom_mass_Si_g);
        double Total_Density = Fe_Fraction+O_Fraction+Si_Fraction;
        Material_Possibility_Array[0]= (Fe_Fraction/Total_Density);
        Material_Possibility_Array[1]= (O_Fraction/Total_Density);
        Material_Possibility_Array[2]= (Si_Fraction/Total_Density);
    }
    
    if(Earth_or_air==2)//AIR
    {
        Atom_Mass[0]        =AN ; Atom_Mass[1]        =AO;
        double N_Fraction = Mass_Density_of_Earth_New(1,N_Percentage_ATM,atom_mass_N_g);
        double O_Fraction  = Mass_Density_of_Earth_New(1,O_Percentage_ATM,atom_mass_O_g);
        double Total_Density = N_Fraction+O_Fraction;
        Material_Possibility_Array[0]= (N_Fraction/Total_Density);
        Material_Possibility_Array[1]= (O_Fraction/Total_Density);
        Material_Possibility_Array[2]= 0;
    }
    
    TH1F   *Material_Possibility = new TH1F("Material_Possibility","Material_Possibility",3,0,3);
    for(int kkk=0;kkk<3;kkk++){Material_Possibility->SetBinContent(kkk+1,Material_Possibility_Array[kkk]);}
    double Random_Energy= Material_Possibility->GetRandom();
    int Random_Energy_int = Random_Energy;
    
    cout << "Which atom you prefer? : " << Atom_Mass[Random_Energy_int] << endl;
    
    return Atom_Mass[Random_Energy_int];*/
    
    
}
//mx: mass of dark matter, unit: GeV
//velocity: natural unit
//sigma_SI: in unit of cm^2
//T: recoil nuclear energy, in unit of MeV
double max_recoil_A(double mx, double velocity, double atomic_mass)
{
    double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));

    double r_A = 4.0*pow(reduce_mass_A,2.0)/(mx*1000.0*unified_atomic_mass_MeV*atomic_mass);
    
    double max_recoil_A_0 = 0.5*mx*1000.0*velocity*velocity*r_A; //MeV
    
    return max_recoil_A_0;
}

double max_recoil_A_keV(double mx, double velocity, double atomic_mass)
{ return 1000.0*max_recoil_A(mx, velocity, atomic_mass); }
//

///////////////////////////////////////////////////////////////////////////////
// no unit
// T: MeV
double F2(double atomic_mass, double T)
{
    double s_fm = 0.9;//Thickness, fm
    double s = s_fm*fm_MeV1; //unit MeV^-1
    
    double q = sqrt(2.0*unified_atomic_mass_MeV*atomic_mass*T); //unit MeV
    double rn2_fm = pow(((1.23*pow(atomic_mass,(1.0/3.0)))-0.6),2.0)+(7.0/3.0)*pow((0.52*TMath::Pi()),2.0)-(5.0*s_fm*s_fm); //unit fm^2
    double rn = sqrt(rn2_fm)*fm_MeV1; //unit MeV^-1
    double form_factor = 3.0*((sin(q*rn)-(q*rn)*cos(q*rn))/pow((q*rn),3.0))*exp(-1.0*pow((q*s),2.0)/2.0);
    
    return form_factor*form_factor;
}

//cm^2/MeV
//mx: GeV
//sigma_SI: cm^2
//velocity: ratio of speed of light
//T: MeV
double fdsigma_dT(double mx, double sigma_SI, double velocity, double atomic_mass, double T)
{
    double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));
    double reduce_mass_n = mx*1000.0*unified_atomic_mass_MeV/((mx*1000.0)+unified_atomic_mass_MeV);
    
    double sigma_chiN0 = sigma_SI*pow(reduce_mass_A,2.0)*pow(atomic_mass,2.0)/pow(reduce_mass_n,2.0);
    double dsigma_dq2 = sigma_chiN0*F2(atomic_mass, T)/(4.0*pow(reduce_mass_A,2.0)*pow(velocity,2.0));
    double dsigma_dT_0;
    
    
    if(max_recoil_A(mx, velocity, atomic_mass)>T)
    {
        dsigma_dT_0 = 2.0*unified_atomic_mass_MeV*atomic_mass*dsigma_dq2;
    }
    else {
        dsigma_dT_0 = 0.0; }
     
    
    return dsigma_dT_0;
}
// cm^2/keV
//mx: GeV
//sigma_SI: cm^2
//velocity: ratio of speed of light
//T: keV
double fdsigma_dT_keV(double mx=10, double sigma_SI=1e-40, double velocity=0, double atomic_mass=0, double T=0)
{
    return (1.0/1000.0)*fdsigma_dT(mx, sigma_SI, velocity, atomic_mass, (T/1000.0));
}

double AAAA(double mx, double sigma_SI, double velocity, double atomic_mass, double T)//check for high-mass range above 10^10GeV
{
    double reduce_mass_A = (unified_atomic_mass_MeV*atomic_mass)*atomic_mass;//muA ~ mN*A
    double reduce_mass_n = mx*1000.0*unified_atomic_mass_MeV/((mx*1000.0)+unified_atomic_mass_MeV);
    double sigma_chiN0 = sigma_SI*pow(reduce_mass_A,2.0)*pow(atomic_mass,2.0)/pow(reduce_mass_n,2.0);
    double dsigma_dq2 = sigma_chiN0*F2(atomic_mass, T)/(4.0*pow(reduce_mass_A,2.0)*pow(velocity,2.0));
    double dsigma_dT_0;
    
    if(max_recoil_A(mx, velocity, atomic_mass)>T)
    {
        dsigma_dT_0 = 2.0*unified_atomic_mass_MeV*atomic_mass*dsigma_dq2;
    }
    else {
        dsigma_dT_0 = 0.0; }
     
    //dsigma_dT_0 = 2.0*unified_atomic_mass_MeV*atomic_mass*dsigma_dq2;
    return dsigma_dT_0;
}
// cm^2/keV
//mx: GeV
//sigma_SI: cm^2
//velocity: ratio of speed of light
//T: keV
double AAAA_keV(double mx=10, double sigma_SI=1e-40, double velocity=0, double atomic_mass=0, double T=0)
{
    return (1.0/1000.0)*AAAA(mx, sigma_SI, velocity, atomic_mass, (T/1000.0));
}


///////////////////////////////////////////////////////////////////////////////
// for Migdal effect
double max_recoil_A_EM(double mx, double velocity, double atomic_mass)
{
  double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));

  double max_recoil_A_0 = 0.5*reduce_mass_A*velocity*velocity; //MeV

  return max_recoil_A_0;
}

double max_recoil_A_EM_keV(double mx, double velocity, double atomic_mass)
{ return 1000.0*max_recoil_A_EM(mx, velocity, atomic_mass); }

///////////////////////////////////////////////////////////////////////////////
// for Bremsstrahlung
//PE(keV),mx(GeV)
double max_min_recoil_A_BM_keV(int MaxorMin, double mx, double velocity, double atomic_mass, double PE)
{//Photon Energy(PE),MaxorMin=0 Max, MaxorMin=1 Min,
  double Mn=unified_atomic_mass_GeV*atomic_mass;//GeV/c^2
  double un=mx*Mn/(mx+Mn);//GeV/c^2
    //cout << "PE: " << PE << endl;
  double Front_Factor = ( (un*un)*(velocity*velocity)/Mn );
  double Factor_1     = (  1 - (PE/(un*velocity*velocity*1e6))   );
  double Factor_2     = sqrt(  1 - (2*PE/(un*1e6*velocity*velocity))   );
    
  if(1 - (2*PE/(un*1e6*velocity*velocity))<0) Factor_2=0;
    
 double Return_Factor = 0;
 if(MaxorMin==0)Return_Factor = Front_Factor* ( Factor_1 + Factor_2 );//Max
 if(MaxorMin==1)Return_Factor = Front_Factor* ( Factor_1 - Factor_2 );//Min

 //cout << "Return_Factor: " << Return_Factor << endl;
  return Return_Factor*1e6;//keV
}



//For Migdal Effect
double fdsigma_dERdEEM_keV(double Probability, double mx=10, double sigma_SI=1e-40, double velocity=0, double atomic_mass=0, double T=0)
{
    double Cross_Section = fdsigma_dT_keV(mx, sigma_SI, velocity, atomic_mass, T)*Probability;
     return Cross_Section;
}

double RRR()
{
    static double check=0;
    check = check + 1;
    
    return check;
}


/*
double fdsigma_dT_keV_MD(double mx, double sigma_SI, double velocity, double atomic_mass, double ER, int Ionizaed_Energy)//T=ER
{
    //The fixed possibilities
    static double dpdE[n_point1][lnumber];
    static double Eem_keV[n_point1];  // in keV
    static double dpdE_Sum[n_point1];
    
    static int check=0;
    //======Scanning the files=====//
    if(check==0)
    {
        cout << "check: " << check << endl;
      std::ifstream ifs;
      ifs.open("Ge.dat",std::ifstream::in);
     double Ee[n_point+1][lnumber],prob[n_point+1][lnumber];

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
        }
      }
    }
        double WIMP_max_T = max_recoil_A_EM_keV(mx, velocity, atomic_mass); //keV
        //double dsigma_dERdEEM_keV=0;
    for(int jjj=0; jjj<n_point1;jjj++)
    {//Open0
        double Eem = Eem_keV[jjj];//Energy after the ionization (keV)
        for(int kkk=0; kkk<n_point1; kkk++)
        {//Open1
            for(int Number_of_Level=1;Number_of_Level<6;Number_of_Level++)
            {//Open3
                    //For Electronic Recoil
                    if(kkk!=jjj) continue;
                    double Enl = Enl_Series[Number_of_Level];//Binding Energy (keV)
                    double Eel = Eem + Enl;//Total Energy
                    if( (Eel) < WIMP_max_T)
                    {//Open4
                        //cout << "Number_of_Level: " << Number_of_Level << endl;
                        //cout << "dpdE[kkk][Number_of_Level]*(2*ER/(72.64*0.95)): " << dpdE[kkk][Number_of_Level]*(2*ER/(72.64*0.95)) << endl;
                        //dsigma_dERdEEM_keV = 0;
                        //cout << "dpdE_Sum[jjj] : " << dpdE_Sum[jjj]  << endl;
                        dpdE_Sum[jjj] = dpdE_Sum[jjj]+ dpdE[jjj][Number_of_Level];
//                    fdsigma_dERdEEM_keV(dpdE[kkk][Number_of_Level]*(2*ER/(72.64*0.95)),mx, sigma_SI, velocity, atomic_mass, ER);
                    }//Close4
                    else
                    {//Open5
                        continue;
                    }//Close5
            } //Close3
        }//Close1
    }//Close0

        
    }
    check = check + 1;
    
    for (int i=0;i<n_point1;i++) {
        for (int j=0;j<lnumber;j++) {
            //cout << "i: " << i << "j: " << j << endl;
            //cout << "Ee[i]" << Eem_keV[i] << endl;
              //cout << "Prob" << dpdE[i][j] << endl; // energy in eV
          }}

    //=======
    //Check Point
     //cout << "============================================" << endl;
    return dpdE_Sum[Ionizaed_Energy]*fdsigma_dT_keV(mx, sigma_SI, velocity, atomic_mass, ER)*(2*ER/(72.64*0.95));
}
*/

//mx:GeV/c^2, PE:keV
double min_v_BM(double mx, double atomic_mass, double PE)
{
  double Mn=unified_atomic_mass_GeV*atomic_mass;//GeV/c^2
  double un=mx*Mn/(mx+Mn);//GeV/c^2
    
  return sqrt( 2* (PE*1e-6) /un )*(3e8/1e3);//km/s
}

//ER,Eel:keV === Mn,un:GeV/c^2
double min_v_EM(double mx, double atomic_mass,double ER, double Eel)
{//(Mn*ER+un*Eel)/un/TMath::Sqrt(2.*Mn*ER)*300.
  double Mn=unified_atomic_mass_GeV*atomic_mass;//GeV/c^2
  double un=mx*Mn/(mx+Mn);//GeV/c^2
    
  return (Mn*ER+un*Eel)/un/TMath::Sqrt(2.*Mn*ER)*300.;//km/s
}


double fdsigma_dERdw_keV(int Option, double mx, double velocity, double sigma_SI, double atomic_mass, double ER, double EB, double ASF)
{
    double Mn=unified_atomic_mass_GeV*atomic_mass*1e6;//keV/c^2
    double First_Term  = (4*Fine_Structure/ (3*TMath::Pi()*EB) )*(ER/Mn);
    double Second_Term = ASF*ASF;
    double Third_Term  = fdsigma_dT_keV(mx, sigma_SI, velocity, atomic_mass, ER);
    return First_Term*Second_Term*Third_Term;
}

double fdsigma_dERdw_keV_approximation(int Option, double mx, double velocity, double sigma_SI, double atomic_mass, double ER, double PE)
{
    //
    double First_Term  = (4*Fine_Structure/ (3*TMath::Pi()*PE) )*(Atomic_Scattering_Factor(PE)*Atomic_Scattering_Factor(PE));
    //
    double Mn=unified_atomic_mass_GeV*atomic_mass*1e6;//keV/c^2
    double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));//MeV/c^2
    double reduce_mass_n = mx*1000.0*unified_atomic_mass_MeV/((mx*1000.0)+unified_atomic_mass_MeV);//MeV/c^2
    double sigma_chiN0 = sigma_SI*pow(reduce_mass_A,2.0)*pow(atomic_mass,2.0)/pow(reduce_mass_n,2.0);
    double Second_Term = (reduce_mass_A*reduce_mass_A*1e6*velocity*velocity*sigma_chiN0/(Mn*Mn));
    //
    double un=mx*Mn/(mx+Mn);//GeV/c^2
    double Factor_1     = (  1 - (PE/(un*velocity*velocity*1e6))   );
    double Factor_2     = sqrt(  1 - (2*PE/(un*1e6*velocity*velocity))   );
    if(1 - (2*PE/(un*1e6*velocity*velocity))<0) Factor_2=0;

    double Third_Term=Factor_1*Factor_2;
    
    return First_Term*Second_Term*Third_Term;
}
//mx:GeV,velocity:beta,ER:keV,Er:keV,Cross-section
double fdsigma_dERdEr_MPA(double mx, double velocity, double sigma_SI, double atomic_mass, double ER, double Er, double Cross_Section)//Er = E_{F}-E_{I}
{
    double Mn=unified_atomic_mass_GeV*atomic_mass;//GeV/c^2
    double un=mx*Mn/(mx+Mn);//GeV/c^2=1e3MeV/c^2
    
    double First_Term = (Me*Me)/((un*1e3)*(un*1e3)*(velocity)*(velocity));//No_Unit
    double Second_Term = (ER/Er);//No_Unit
    double Third_Term = ( (1)/(4*(TMath::Pi()*TMath::Pi())*Fine_Structure) );

    double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));//MeV/c^2
    double reduce_mass_n = mx*1000.0*unified_atomic_mass_MeV/((mx*1000.0)+unified_atomic_mass_MeV);//MeV/c^2
    double sigma_chiN0 = sigma_SI*pow(reduce_mass_A,2.0)*pow(atomic_mass,2.0)/pow(reduce_mass_n,2.0);

    return First_Term*sigma_chiN0*Second_Term*Cross_Section*Third_Term;
}
//T: MeV
double min_v(double mx, double T, double atomic_mass)
{
    double reduce_mass_A = mx*1000.0*unified_atomic_mass_MeV*atomic_mass/((mx*1000.0)+(unified_atomic_mass_MeV*atomic_mass));
    double r_A = 4.0*pow(reduce_mass_A,2.0)/(mx*1000.0*unified_atomic_mass_MeV*atomic_mass);
    
    double velocity0 = sqrt(2.0*T/(mx*1000.0*r_A));
    
    return velocity0;
}

//T: keV
double min_v_keV(double mx, double T, double atomic_mass)
{ return min_v(mx, T/1000.0, atomic_mass); }

/*
double *Velocity_Aft_collision_Coupled(int Collision_Time=0, double mx=10, double Sigma_SI_Default=1e-40, double Initial_Velocity=0, int Earth_or_air_S=0)
{
    static double RETURN_VALUE[3];//(1)Final_Velocity(2)Energy_Difference(T1+T2+...)(3)Energy_Difference(E_Final-E_initial)
    static double dpdE[n_point1][lnumber];
    static double Eem_keV[n_point1];  // in keV
    static int check=0;
    if(check==0)
    {
      std::ifstream ifs;
      ifs.open("Ge.dat",std::ifstream::in);
     double Ee[n_point+1][lnumber],prob[n_point+1][lnumber];

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
        }
      }
    }
    }
    check = check + 1;
    
    double Original_Energy=Energy_DM(mx,Initial_Velocity*1e3/3e8);
    if(Initial_Velocity!=0)
    {
        int Collision_Time_Temp=Collision_Time;
        double DM_Energy_Aft_Colliding=Energy_DM(mx,Initial_Velocity*1e3/3e8);
        //cout << "Initial_Energy: " << DM_Energy_Aft_Colliding << endl;
        double DM_Velocity_Aft_Colliding=Initial_Velocity;
        double Energy_Lost_Total=0;
        if(Collision_Time!=0 and DM_Velocity_Aft_Colliding>1)
        {
            for(int kkk=0 ; kkk<Collision_Time_Temp ; kkk++)
            {
                double ATOM_KIND = Function_of_material_possibility(Earth_or_air_S);
                //cout << "Earth_or_air_S: " << Earth_or_air_S << endl;
                cout << "Collision_TIME: kkk " << kkk << endl;
                //===================================Nuclear Recoil===================================
                TF1 *f3 = new TF1("f3","fdsigma_dT_keV([0],[1],[2],[3],x)",0,max_recoil_A_keV(mx,DM_Velocity_Aft_Colliding*1e3/3e8,ATOM_KIND));
                f3->SetParameter(0,mx);f3->SetParameter(1,Sigma_SI_Default);f3->SetParameter(2,(DM_Velocity_Aft_Colliding*1e3/3e8));f3->SetParameter(3,ATOM_KIND);
                double Random_Energy_NU= f3->GetRandom();
                cout << "=========================Nuclear Recoil=========================================" << endl;
                cout << "Random_Energy_NU: " << Random_Energy_NU << endl;
                //===================================Migdal_Effect===================================
                cout << "============================Migdal_Effect=======================================" << endl;
                TF1 *f4 = new TF1("f4","fdsigma_dT_keV_MD([0],[1],[2],[3],[4],x)",0,n_point1);
                f4->SetParameter(0,mx);f4->SetParameter(1,Sigma_SI_Default);f4->SetParameter(2,(DM_Velocity_Aft_Colliding*1e3/3e8));f4->SetParameter(3,ATOM_KIND);f4->SetParameter(4,Random_Energy_NU);
                cout << "===========================================================================================================" << endl;
                double Random_Energy_MD= Eem_keV[int(f4->GetRandom())];
                if(Random_Energy_MD<=1e-5) Random_Energy_MD=0;
                //====================================================================================
                cout << "DM_Energy_Aft_Colliding: " << DM_Energy_Aft_Colliding << endl;
                cout << "Max_Recoil_NU: " << max_recoil_A_keV(mx,DM_Velocity_Aft_Colliding*1e3/3e8,ATOM_KIND) << endl;
                cout << "Random_Energy_NU: " << Random_Energy_NU << endl;
                cout << "Max_Recoil_MD: " << max_recoil_A_EM_keV(mx,DM_Velocity_Aft_Colliding*1e3/3e8, ATOM_KIND) << endl;
                cout << "Random_Energy_MD: " << Random_Energy_MD << endl;
                Energy_Lost_Total = Energy_Lost_Total + (Random_Energy_NU + Random_Energy_MD);
                DM_Energy_Aft_Colliding = (DM_Energy_Aft_Colliding - Random_Energy_NU-Random_Energy_MD);
                if(DM_Velocity_Aft_Colliding<1 or DM_Energy_Aft_Colliding<0)
                {
                    cout << "BREAK: " << endl;
                    break;
                }

                cout << "DM_Energy_Aft_Colliding: " << DM_Energy_Aft_Colliding << endl;
                 
                DM_Velocity_Aft_Colliding = Velocity_DM(mx,DM_Energy_Aft_Colliding);
                cout << "DM_Velocity_Aft_Colliding: " << DM_Velocity_Aft_Colliding << endl;
                //====================================================================================
            }
        }
        else
        {
            //cout << "GHGH: " << endl;
            DM_Velocity_Aft_Colliding=Initial_Velocity;
        }
        RETURN_VALUE[0]=DM_Velocity_Aft_Colliding;
        RETURN_VALUE[1]=Energy_Lost_Total;
        RETURN_VALUE[2]=Energy_DM(mx,Initial_Velocity*1e3/3e8)-Energy_DM(mx,DM_Velocity_Aft_Colliding*1e3/3e8);
        return RETURN_VALUE;
    }
    else
    { RETURN_VALUE[0]=0; RETURN_VALUE[1]=0; RETURN_VALUE[2]=0;
        return RETURN_VALUE;}
}

*/
//Energy_Aft_multiple_collision(GOOD)
double *Velocity_Aft_collision(int Collision_Time=0, double mx=10, double Sigma_SI_Default=1e-40, double Initial_Velocity=0, int Earth_or_air_S=0)
{
    static double RETURN_VALUE[4];//(1)Final_Velocity(2)Energy_Difference(T1+T2+...)(3)Energy_Difference(E_Final-E_initial)
    
    //if(Earth_or_air_S==2){cout << "Earth_or_air_S: AIR: " << Collision_Time << endl; cout << "Velocity: " << Initial_Velocity << endl;}
    //if(Earth_or_air_S==1){cout << "Earth_or_air_S: EARTH: " << Collision_Time << endl;cout << "Velocity: " << Initial_Velocity << endl;}
    //cout <<"Velocity_Aft_collision_Collision_TIme: " << Collision_Time << endl;
    // Percentage of the material of the earth
    if(Initial_Velocity!=0)
    {
        int Collision_Time_Temp=Collision_Time;
        //cout << "Initial_Velocity: " << Initial_Velocity << endl;
        double DM_Energy_Aft_Colliding=Energy_DM(mx,Initial_Velocity*1e3/3e8);
        //cout << "DM_Energy_Bef: " << DM_Energy_Aft_Colliding << endl;
        double DM_Velocity_Aft_Colliding=Initial_Velocity;
        double Energy_Lost_Total=0;
        int Count=0;
        if(Collision_Time!=0 and DM_Velocity_Aft_Colliding>1)
        {
            for(int kkk=0 ; kkk<Collision_Time_Temp ; kkk++)
            {
                double ATOM_KIND = Function_of_material_possibility(Earth_or_air_S);// 1 is earth 2 is air
                //cout << "Earth_or_air_S: " << Earth_or_air_S << endl;
                //cout << "ATOM_KIND: " << ATOM_KIND << endl;
                TF1 *f3 = new TF1("f3","fdsigma_dT_keV([0],[1],[2],[3],x)",0,max_recoil_A_keV(mx,DM_Velocity_Aft_Colliding*1e3/3e8,ATOM_KIND));
                f3->SetParameter(0,mx);f3->SetParameter(1,Sigma_SI_Default);f3->SetParameter(2,(DM_Velocity_Aft_Colliding*1e3/3e8));f3->SetParameter(3,ATOM_KIND);
                double Random_Energy= f3->GetRandom();
                //cout << "Random_Energy: " << Random_Energy << endl;

                Energy_Lost_Total = Energy_Lost_Total + Random_Energy;
                DM_Energy_Aft_Colliding = (DM_Energy_Aft_Colliding - Random_Energy);
                DM_Velocity_Aft_Colliding = Velocity_DM(mx,DM_Energy_Aft_Colliding);
                //cout << "kkk: " << kkk << endl;
                //cout << "DM_Velocity_Aft_Colliding: " << DM_Velocity_Aft_Colliding << endl;
                Count = Count + 1;
                if(DM_Energy_Aft_Colliding<1e-2 or Random_Energy==0)
                {
                    //cout << "BREAK: " << endl;
                    break;
                }
                 
            }
        }
        else
        {
            //cout << "GHGH: " << endl;
            DM_Velocity_Aft_Colliding=Initial_Velocity;
        }
        //cout << "Count: " << Count << endl;
        //cout << "DM_Velocity_Aft_Colliding_Aft: " << DM_Velocity_Aft_Colliding << endl;
        RETURN_VALUE[0]=DM_Velocity_Aft_Colliding;
        RETURN_VALUE[1]=Energy_Lost_Total;
        RETURN_VALUE[2]=Energy_DM(mx,Initial_Velocity*1e3/3e8)-Energy_DM(mx,DM_Velocity_Aft_Colliding*1e3/3e8);
        RETURN_VALUE[3]=Count;
        return RETURN_VALUE;
    }
    else
    { RETURN_VALUE[0]=0; RETURN_VALUE[1]=0; RETURN_VALUE[2]=0;RETURN_VALUE[3]=0;
        return RETURN_VALUE;}
}



double total_Sigma(int Option=0, double Velocity=0, double Sigma_SI=0, double WIMP_mx =10, double A = AGe)//Velocity(m/s)
 {
    //cout << "Option: " << Option << endl;
    //cout << "Velocity: " << Velocity << endl;
    
    int reso_T=1000;double T[reso_T];double total_Sigma=0;
     double WIMP_max_T = 1000.0*max_recoil_A(WIMP_mx, 779.135*1000.0/2.99792458e8, A); //keV
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
        if(i==0) { dEx = T[0]; } else { dEx = T[i] - pEx; }
        if(Option==0 or (Option==1 and Velocity!=0))
            {
                total_Sigma = total_Sigma + (fdsigma_dT_keV(WIMP_mx, Sigma_SI, (Velocity*1e3/3e8), A, T[i])*dEx);
            }
        pEx = T[i];
    }
    return total_Sigma;
}
double total_C_AAAA(int Option=0, double Velocity=0, double Sigma_SI=0, double WIMP_mx =10, double A = AGe)//Velocity(m/s)
 {//Check for High-Mass Range above 10^10GeV
    //cout << "Option: " << Option << endl;
    //cout << "Velocity: " << Velocity << endl;
    
    int reso_T=1000;double T[reso_T];double total_Sigma=0;
     double WIMP_max_T = 1000.0*max_recoil_A(WIMP_mx, 779.135*1000.0/2.99792458e8, A); //keV
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
        if(i==0) { dEx = T[0]; } else { dEx = T[i] - pEx; }
        if(Option==0 or (Option==1 and Velocity!=0))
            {
                total_Sigma = total_Sigma + (AAAA_keV(WIMP_mx, Sigma_SI, (Velocity*1e3/3e8), A, T[i])*dEx);
            }
        pEx = T[i];
    }
    return total_Sigma;
}

double CDEX_Collision_Time_ATM(double Sigma_SI_Default, double PZ, double Velocity, double WIMP_mx,double CJPL_Length) //Velocity(km/s)
{
    //cout << "OK! Thanks!" << endl;
    //Theta(Between r and Z) Phi(Between r and X)
    double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
    double Lab_CARC[3]={0,0,0};// Original point X,Y,Z Function=>TMath::Pi()
    
    double Path_Lengths_for_atmosphere[19];
    double Density_of_Atmosphere_Layer[19];
    double Previous_Total=0;double N_Collision=0;
    
    for(int kkk=0; kkk<19; kkk++)
    {
        //cout << "kkk: " << kkk << endl;
        Density_of_Atmosphere_Layer[kkk] = (atm_table[kkk][4]+atm_table[kkk+1][4])/2;
    Path_Lengths_for_atmosphere[kkk] =Path_Length_Air(PZ,(atm_table[kkk+1][0]/1000),Path_Length_Earth_Rough(PZ),Previous_Total);
        //cout << "Path_Lengths_for_atmosphere_Bef[kkk]: " << Path_Lengths_for_atmosphere[kkk] << endl;
        Previous_Total = Previous_Total+ Path_Lengths_for_atmosphere[kkk];
        //cout << "Previous_Total : " << Previous_Total  << endl;
        if(Previous_Total <  CJPL_Length and PZ>=0) Path_Lengths_for_atmosphere[kkk]=0;
   else if(Previous_Total >= CJPL_Length and (Path_Lengths_for_atmosphere[kkk-1]==0 or kkk==0) and PZ>=0)Path_Lengths_for_atmosphere[kkk] = Previous_Total - CJPL_Length;
        //cout << "Density_of_Atmosphere_Layer[kkk] : " << Density_of_Atmosphere_Layer[kkk]  << endl;
        //cout << "Path_Lengths_for_atmosphere_Aft[kkk]: " << Path_Lengths_for_atmosphere[kkk] << endl;
        N_Collision = N_Collision + (kg_perm3_to_g_percm3(Density_of_Atmosphere_Layer[kkk])*Path_Lengths_for_atmosphere[kkk]*1e5)/(unified_atomic_mass_g*(0.8*AN+0.2*AO))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,0.8*AN+0.2*AO));
    }
    cout << "N_Collision_ATM: " << N_Collision << endl;
    return N_Collision;
}

double KS_Collision_Time_ATM(double Sigma_SI_Default, double PY, double PZ, double Velocity, double WIMP_Mass, double *Length_Components, double Earth_Path_Length) //Velocity(km/s)
{
    //cout << "OK! Thanks!" << endl;
    //Theta(Between r and Z) Phi(Between r and X)
    //cout << "PY: " << PY << "PZ: " << PZ << endl;
    double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
    double Lab_CARC[3]={0,0,0};// Original point X,Y,Z Function=>TMath::Pi()
    
    double Path_Lengths_for_atmosphere[19];
    double Density_of_Atmosphere_Layer[19];
    double Previous_Total_AIR=0;
    double N_Collision=0;
    double Previous_Total = ( Length_Components[1]+Length_Components[2]+Earth_Path_Length);

    for(int kkk=0; kkk<19; kkk++)
    {
        //cout << "AIR: " << kkk << endl;
        Density_of_Atmosphere_Layer[kkk] = (atm_table[kkk][4]+atm_table[kkk+1][4])/2;
        //cout << "Density_of_Atmosphere_Layer[kkk]: " << Density_of_Atmosphere_Layer[kkk] << endl;
        Path_Lengths_for_atmosphere[kkk] = KS_Air_Path_Length(PY,PZ,(atm_table[kkk+1][0]/1000),Previous_Total_AIR);
        //cout << "(atm_table[kkk+1][0]/1000): " << (atm_table[kkk+1][0]/1000) << endl;
        Previous_Total_AIR = Previous_Total_AIR+ Path_Lengths_for_atmosphere[kkk];
        if(kkk==0) Path_Lengths_for_atmosphere[kkk] = Path_Lengths_for_atmosphere[kkk] - Previous_Total;
        //cout << "Path_Lengths_for_atmosphere[kkk]: " << Path_Lengths_for_atmosphere[kkk] << "km " << endl;
        //cout << " Previous_Total_AIR: " <<  Previous_Total_AIR << endl;
        N_Collision = N_Collision + (kg_perm3_to_g_percm3(Density_of_Atmosphere_Layer[kkk])*Path_Lengths_for_atmosphere[kkk]*1e5)/(unified_atomic_mass_g*(0.8*AN+0.2*AO))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,0.8*AN+0.2*AO));
    }
    //cout << "N_Collision_AIR: " << N_Collision << endl;
    
    return N_Collision;
}
double Collision_Time_ATM(double Sigma_SI_Default, double Try_Z_direction_Unit, double Velocity, double WIMP_mx) //Velocity(km/s)
{
    //cout << "OK! Thanks!" << endl;
    //Theta(Between r and Z) Phi(Between r and X)
    double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
    double Lab_CARC[3]={0,0,0};// Original point X,Y,Z Function=>TMath::Pi()
    
    double Path_Lengths_for_atmosphere[19];
    double Density_of_Atmosphere_Layer[19];
    double Previous_Total=0;double N_Collision=0;
    
    for(int kkk=0; kkk<19; kkk++)
    {
        //cout << "kkk: " << kkk << endl;
        Density_of_Atmosphere_Layer[kkk] = (atm_table[kkk][4]+atm_table[kkk+1][4])/2;
        //cout << "Density_of_Atmosphere_Layer[kkk] : " << Density_of_Atmosphere_Layer[kkk]  << endl;
        Path_Lengths_for_atmosphere[kkk] = Path_Length_Air(Try_Z_direction_Unit,(atm_table[kkk+1][0]/1000),Path_Length_Earth_Rough(Try_Z_direction_Unit),Previous_Total);
        //cout << "Path_Lengths_for_atmosphere[kkk]: " << Path_Lengths_for_atmosphere[kkk] << endl;
        Previous_Total = Previous_Total+ Path_Lengths_for_atmosphere[kkk];
        //cout << "Previous_Total : " << Previous_Total  << endl;
        N_Collision = N_Collision + (kg_perm3_to_g_percm3(Density_of_Atmosphere_Layer[kkk])*Path_Lengths_for_atmosphere[kkk]*1000*100)/(unified_atomic_mass_g*(0.8*AN+0.2*AO))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,0.8*AN+0.2*AO));
        
    }
    cout << "N_Collision_ATM: " << N_Collision << endl;
    return N_Collision;
}

double CDEX_Collision_Time_EARTH(double Sigma_SI_Default,double PX, double PY, double PZ, double Velocity, double WIMP_mx, double CJPL_Length) //Velocity(km/s)
{
    //cout << "OK! Thanks!" << endl;
    //Theta(Between r and Z) Phi(Between r and X)
    double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
    double Lab_CARC[3]={0,0,0};// Original point X,Y,Z Function=>TMath::Pi()
    
    double Path_Lengths_for_earth[13];
    double Density_of_earth_Layer[13];
    double Previous_Total=0;double N_Collision=0;double N_Collision_Mountain=0;
    //cout << "Weighted_Atomic_Number" << Weighted_Atomic_Number << endl;
    double Check_Length_times_Density=0;
    for(int kkk=0; kkk<13; kkk++)
    {
        //cout << "kkk: " << kkk << endl;
        Density_of_earth_Layer[kkk] = (earth_table[kkk+1][1]+earth_table[kkk][1])/2;
        double Length_Temp=Path_Length_Earth(PZ,earth_table[kkk+1][0],Previous_Total);
        
        if(PZ>0)Length_Temp=0;
        if(Length_Temp>0 and Length_Temp<=6372*2){Path_Lengths_for_earth[kkk] = Length_Temp;}
        else{Path_Lengths_for_earth[kkk]=0;}
        
        if(Previous_Total>=0 and Previous_Total<=6372*2)Previous_Total = Previous_Total+ Path_Lengths_for_earth[kkk];
        else{Previous_Total = Previous_Total+ 0;}
        
        Check_Length_times_Density = Check_Length_times_Density + Path_Lengths_for_earth[kkk]*1e5;
        N_Collision = N_Collision + ( Density_of_earth_Layer[kkk]*Path_Lengths_for_earth[kkk]*1e5)/(unified_atomic_mass_g*(Weighted_Atomic_Number))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,Weighted_Atomic_Number) );
        
    }
    N_Collision_Mountain = N_Collision_Mountain + (CJPL_Length*1e5*shield_atom_density)*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,shield_A));
    cout << "N_Collision_Mountain: " << N_Collision_Mountain << endl;
    cout << "N_Collision_Earth: " << N_Collision << endl;
    return N_Collision+N_Collision_Mountain;
}

//KS Power Plane Earth Plane
double *KS_Collision_Time_EARTH(double Sigma_SI_Default, double PY, double PZ, double Velocity, double WIMP_Mass, double *Length_Components) //Velocity(km/s)
{
    //Theta(Between r and Z) Phi(Between r and X)
    double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;

    double Path_Lengths_for_earth[13];
    double Density_of_earth_Layer[13];
    double Previous_Total=0;
    double N_Collision_Earth=0;double N_Collision_Cement=0;double N_Collision_Reactor_Wall=0;double N_Collision_Reactor_Water=0;
    double N_Collision_Shielding=0;
    //cout << "Weighted_Atomic_Number" << Weighted_Atomic_Number << endl;
    double Check_Length_times_Density=0;
    static double RETURN_VALUE[6];
    /*
    cout << "Length_Components[0]: " << Length_Components[0] << endl;
    cout << "Length_Components[1]: " << Length_Components[1] << endl;
    cout << "Length_Components[2]: " << Length_Components[2] << endl;
     */
    for(int kkk=0; kkk<13; kkk++)
    {
      //  cout << "Earth: " << kkk << endl;
        Density_of_earth_Layer[kkk] = (earth_table[kkk+1][1]+earth_table[kkk][1])/2;
       // cout << "Density_of_earth_Layer[kkk]: " << Density_of_earth_Layer[kkk] << endl;
        double Length_Temp=KS_Earth_Path_Length(PY,PZ,earth_table[kkk+1][0],Previous_Total);
       // cout << "earth_table[kkk+1][0]: " << earth_table[kkk+1][0] << endl;
        if(Length_Temp>0 and Length_Temp<=6372*2 and PZ<0)
            {
                Path_Lengths_for_earth[kkk] = Length_Temp;
            }
        else if(kkk==12 and Length_Components[0]!=0 and PZ>=0)
            {
                //cout << "kkk: " << kkk << endl;
                //cout << "Length_Temp: " << Length_Temp << endl;
                //cout << "Length_Components[0]: " << Length_Components[0] << endl;
                Path_Lengths_for_earth[kkk] = Length_Temp-Length_Components[0];
            }
        else{Path_Lengths_for_earth[kkk]=0;}
        //cout << "Path_Lengths_for_earth[kkk]: " << Path_Lengths_for_earth[kkk] << endl;
        if(Previous_Total>=0 and Previous_Total<=6372*2)Previous_Total = Previous_Total+ Path_Lengths_for_earth[kkk];
        else{Previous_Total = Previous_Total+ 0;}
       // cout << "Previous_Total: " << Previous_Total  << endl;
    //Earth
    N_Collision_Earth = N_Collision_Earth + ( Density_of_earth_Layer[kkk]*Path_Lengths_for_earth[kkk]*1e5)/(unified_atomic_mass_g*(Weighted_Atomic_Number))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,Weighted_Atomic_Number));
    }
    //cout << "Earth_Collision: " <<  N_Collision_Earth << endl;
    //Cement
    N_Collision_Cement =  ((Density_of_Cement)*(Length_Components[1]*1e5))/(unified_atomic_mass_g*(Weighted_Atomic_Number_Cement))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,Weighted_Atomic_Number_Cement));
    //cout << "Cement_Collision: " <<  N_Collision_Cement << endl;
    //Reactor_Wall
    N_Collision_Reactor_Wall =
        (Density_of_Cement*((Length_Components[2]-Length_Components[3])*1e5))/(unified_atomic_mass_g*(Weighted_Atomic_Number_Cement))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,Weighted_Atomic_Number_Cement));
    //cout << "Reactor_Collision_Wall: " <<  N_Collision_Reactor_Wall << endl;
    //Reactor_Water
    N_Collision_Reactor_Water =
        (1*(Length_Components[3]*1e5))/(unified_atomic_mass_g*(AH2O))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,AH2O));
    //cout << "Reactor_Collision_Water: " <<  N_Collision_Reactor_Water << endl;//H2O density is 1g/cc
    //Shielding
    N_Collision_Shielding =
    (11.34*(15))/(unified_atomic_mass_g*(APb))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,APb))
    +(7.86*(5))/(unified_atomic_mass_g*(AFe))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,AFe))
    +(2.34*(25))/(unified_atomic_mass_g*(AB))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,AB))
    +(8.96*(5))/(unified_atomic_mass_g*(ACu))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,ACu));
    //cout << "N_Collision_Shielding: " <<  N_Collision_Shielding << endl;

    RETURN_VALUE[0]=Previous_Total;RETURN_VALUE[1]=N_Collision_Earth;RETURN_VALUE[2]=N_Collision_Cement;RETURN_VALUE[3]=N_Collision_Reactor_Wall;
    RETURN_VALUE[4]=N_Collision_Reactor_Water;RETURN_VALUE[5]=N_Collision_Shielding;
    return RETURN_VALUE;
}


double Collision_Time_EARTH(double Sigma_SI_Default, double Try_Z_direction_Unit, double Velocity, double WIMP_mx) //Velocity(km/s)
{
    //cout << "OK! Thanks!" << endl;
    //Theta(Between r and Z) Phi(Between r and X)
    double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
    double Lab_CARC[3]={0,0,0};// Original point X,Y,Z Function=>TMath::Pi()
    
    double Path_Lengths_for_earth[13];
    double Density_of_earth_Layer[13];
    double Previous_Total=0;double N_Collision=0;double N_Collision_Cement=0;
    //cout << "Weighted_Atomic_Number" << Weighted_Atomic_Number << endl;
    double Check_Length_times_Density=0;
    for(int kkk=0; kkk<13; kkk++)
    {
        //cout << "kkk: " << kkk << endl;
        Density_of_earth_Layer[kkk] = (earth_table[kkk+1][1]+earth_table[kkk][1])/2;
        double Length_Temp=Path_Length_Earth(Try_Z_direction_Unit,earth_table[kkk+1][0],Previous_Total);
        
        if(Try_Z_direction_Unit>0)Length_Temp=0;
        if(Length_Temp>0 and Length_Temp<=6372*2){Path_Lengths_for_earth[kkk] = Length_Temp;}
        else{Path_Lengths_for_earth[kkk]=0;}
        
        if(Previous_Total>=0 and Previous_Total<=6372*2)Previous_Total = Previous_Total+ Path_Lengths_for_earth[kkk];
        else{Previous_Total = Previous_Total+ 0;}
        
        Check_Length_times_Density = Check_Length_times_Density + Path_Lengths_for_earth[kkk]*1e5;
        N_Collision = N_Collision + ( Density_of_earth_Layer[kkk]*Path_Lengths_for_earth[kkk]*1e5)/(unified_atomic_mass_g*(Weighted_Atomic_Number))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,Weighted_Atomic_Number));
        
    }
    cout << "N_Collision_Earth: " << N_Collision << endl;
    return N_Collision;
}

double Mean_Free_Path_total(double *Lamda, int Number_of_material)
{
    double Lamda_total=0;
    for(int i=0; i<Number_of_material; i++)
    {
        //cout << "Lamda: " << Lamda[i] << endl;
        Lamda_total = Lamda_total + Lamda[i];
    }
    //cout << "Lamda_total: " << Lamda_total << endl;
    return (1/Lamda_total);
}

double Mean_Free_Path_of_Earth(int Option, double Velocity,double Sigma_SI)
{
    double Mean_Free_Path_array;//Total_Sigma*Mass_density
    double Weighted_Atomic_Number=(AFe*Fe_Percentage+AO*O_Percentage+ASi*Si_Percentage);
  Mean_Free_Path_array= 1/ ( total_Sigma(Option,Velocity,Sigma_SI,10,Weighted_Atomic_Number) * Mass_Density_of_Earth_New(5.67381,1,(unified_atomic_mass_g*(Weighted_Atomic_Number))) );
    
    return Mean_Free_Path_array;
    
}

  
//keV
//T:keV
double rate_scale_QF(double T)
{
  double dQF;
  dQF = 0.4343*(qf1+2*qf2*log10(T)+3*qf3*pow(log10(T),2)+4*qf4*pow(log10(T),3));
  
  return 1.0/(QF(T)+dQF);
}

double inv_TQF(double T_QF)
{
  TF1 *f_TQF = new TF1("f_TQF","TQF(x)",1e-9,1e9);

  if(T_QF<1e-9) { return f_TQF->GetX(1e-9); }
  else if(T_QF>1e9) { return f_TQF->GetX(1e9); }
  else { return f_TQF->GetX(T_QF); }
}

double insert(double x, double x0, double y0, double x1, double y1)
{
  double y  = y1 - (x1-x)*(y1-y0)/(x1-x0);
  return y;
}

 double Length_for_asking_the_collision(double Collision_Time, double WIMP_mx, double Velocity, double Sigma_SI_Default, double Density, double A_on_the_path)//km
 {//The example of count: (11.34*(15))/(unified_atomic_mass_g*(APb))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,APb))
     cout << "Collision_Time: " << Collision_Time << endl;
     cout << "WIMP_mx: " << WIMP_mx << endl;
     cout << "Velocity: " << Velocity << endl;
     cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
     cout << "Density: " << Density << endl;
     cout << "A_on_the_path:" << A_on_the_path << endl;
     
     double First_Term = ((Density)/(unified_atomic_mass_g*(A_on_the_path))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,A_on_the_path)));
     cout << "total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,A_on_the_path): " << total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,A_on_the_path) <<endl;
     cout << "First_Term: " << First_Term << endl;
     cout << "(Collision_Time/First_Term)*1e-5: " << (Collision_Time/First_Term)*1e-5 << endl;
     return (Collision_Time/First_Term)*1e-5;//km
 }
double Averaged_Collision_Time_with_Length(double Length, double WIMP_mx, double Velocity, double Sigma_SI_Default, double Density, double A_on_the_path)//km
{//The example of count: (11.34*(15))/(unified_atomic_mass_g*(APb))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,APb))
    double Time = ((Density)*Length*1e5/(unified_atomic_mass_g*(A_on_the_path))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,A_on_the_path)));
    //cout << "total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,A_on_the_path): " << total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_mx,A_on_the_path) <<endl;
    //cout << "(Collision_Time/First_Term)*1e-5: " << (Collision_Time/First_Term)*1e-5 << endl;
    return Time;//Collision_Time
}
 //===============================================Modified_Method_no_angular_consideration===============================================
 double *KS_Collision_Time_ATM_Aft_velocity(int index, double Sigma_SI_Default, double PY, double PZ, double Velocity, double WIMP_Mass, double *Length_Components, double Earth_Path_Length) //Velocity(km/s)
 {
     static double RETURN_VALUE[4];
     static int TEMP_INDEX=-1;
     //cout << "OK! Thanks!" << endl;
     //Theta(Between r and Z) Phi(Between r and X)
     //cout << "PY: " << PY << "PZ: " << PZ << endl;
     double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
     double Lab_CARC[3]={0,0,0};// Original point X,Y,Z Function=>TMath::Pi()
     
     double Path_Lengths_for_atmosphere[19];
     double Density_of_Atmosphere_Layer[19];
     double Previous_Total_AIR=0;
     double Previous_Total = ( Length_Components[1]+Length_Components[2]+Earth_Path_Length);
     //cout << "Earth_Path_Length: " << Earth_Path_Length  << endl;
     //cout << "Previous_Total: " << Previous_Total << endl;

     //Calculate the lengths for different layers
     for(int kkk=0; kkk<19; kkk++)
     {
         //cout << "AIR: " << kkk << endl;
         Density_of_Atmosphere_Layer[kkk] = (atm_table[kkk][4]+atm_table[kkk+1][4])/2;
         //cout << "Density_of_Atmosphere_Layer[kkk]: " << Density_of_Atmosphere_Layer[kkk] << endl;
         Path_Lengths_for_atmosphere[kkk] = KS_Air_Path_Length(PY,PZ,(atm_table[kkk+1][0]/1000),Previous_Total_AIR);
         //cout << "(atm_table[kkk+1][0]/1000): " << (atm_table[kkk+1][0]/1000) << endl;
         Previous_Total_AIR = Previous_Total_AIR+ Path_Lengths_for_atmosphere[kkk];
         if(kkk==0) Path_Lengths_for_atmosphere[kkk] = Path_Lengths_for_atmosphere[kkk] - Previous_Total;
         //cout << "Path_Lengths_for_atmosphere[kkk]: " << Path_Lengths_for_atmosphere[kkk] << "km " << endl;
         //cout << " Previous_Total_AIR: " <<  Previous_Total_AIR << endl;
     }
     

     //Step-to-Step(STS) collision
     double STS_Length[19];double STS_Density[19];
     //Check the lengths on the path
     for(int kkk=0; kkk<19; kkk++)
     {
         int OUTtoIN=18-kkk;
         if(Path_Lengths_for_atmosphere[OUTtoIN]>0){STS_Length[kkk]=(Path_Lengths_for_atmosphere[OUTtoIN]);}
         STS_Density[kkk]=kg_perm3_to_g_percm3(Density_of_Atmosphere_Layer[OUTtoIN]);
         //cout << "kkk: " << kkk << endl;
         //cout << "STS_Length[kkk]: " << STS_Length[kkk] << endl;
         //cout << "STS_Density[kkk]: " << STS_Density[kkk] << endl;
     }
     
     double DM_Velocity_Aft_Colliding=Velocity;
     //cout << "DM_Velocity_Aft_Colliding_Air: " << DM_Velocity_Aft_Colliding << endl;
     double Collision_TIME_Check=0;double Actual_collision=0;
     //double *V_Aft_Collision_AIR = Velocity_Aft_collision(AC,DM_mx,Sigma_SI,Velocity[kkk],2);
     double Temp_Length=0;
     double Path_Length=0;
     double Lamda_for_Average=0.001;
     
     cout << "Ene_Vel_Air_Initial: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
     cout << "Velocity_Air_Initial: " << DM_Velocity_Aft_Colliding << endl;

     static int Air_Threshold=0;
     cout << "index: " << index << endl;
     //cout << "================Check_The_turning_Point(Air)================" << endl;
     if(TEMP_INDEX!=index)
     {
         for(int EEE=0; EEE<5; EEE++)//i is earth and 2 is air
         {
             int Check   = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI_Default,779.135,2)[3];
             if(Air_Threshold<Check)Air_Threshold=Check;
         }
         cout << "Air_Threshold: " << Air_Threshold << endl;
     }
     TEMP_INDEX = index;
     //cout << "=======================================================" << endl;

     //===========================
     int Check_Threshold=0;
     for(int kkk=0; kkk<19; kkk++)
     {
         double Segment_Test = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,STS_Density[kkk],14.99);
         double Expectation_Test =  STS_Length[kkk]/Segment_Test*(Possion_GetRandom(2,Lamda_for_Average));
         //cout << "STS_Length[kkk]/Segment_Test: " << STS_Length[kkk]/Segment_Test << endl;

         if(Expectation_Test>Air_Threshold)
         {
             cout << " Expectation_Test:" << Expectation_Test << ">" << "Air_Threshold: " << Air_Threshold << endl;
             Check_Threshold  = Check_Threshold + 1;
             DM_Velocity_Aft_Colliding=1e-5;
         }
     }
     
     if(Check_Threshold==0)
     {
         /*
         for(int kkk=0; kkk<19; kkk++)
         {
             Path_Length= Path_Length + STS_Length[kkk];
             
             while(Temp_Length<Path_Length)
             {
                 //Every time it passes though a length for colliding once
                 double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,STS_Density[kkk],14.99);
                 //cout << "Segment_Number: " << STS_Length[kkk]/Segment << endl;
                 Collision_TIME_Check = Collision_TIME_Check + 1;
                 Temp_Length = Temp_Length + Segment;
                 //See whether it actually collides once or not
                 if(Possion_GetRandom(0,Lamda_for_Average)>=1)
                 {
                      //cout << "Air" << endl;
                      //cout << "Before_Vel: " <<DM_Velocity_Aft_Colliding << endl;
                     //cout << "kkk: " << kkk << endl;
                     Actual_collision = Actual_collision + 1;
                     double *V_Aft_Collision_AIR = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,2);
                     DM_Velocity_Aft_Colliding=V_Aft_Collision_AIR[0];
                     //cout << "Aft_Vel: " << DM_Velocity_Aft_Colliding << endl;
                     if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)break;
                 }
             }
         }
          */
         for(int kkk=0; kkk<19; kkk++)
         {
             Path_Length= Path_Length + STS_Length[kkk];
             
             while(Temp_Length<Path_Length)
             {
                 //Every time it passes though a length for colliding once
                 double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,STS_Density[kkk],14.99);
                 //cout << "Segment_Number: " << STS_Length[kkk]/Segment << endl;
                 int Times = Possion_GetRandom_Full(Lamda_for_Average);
                 Temp_Length = Temp_Length + Segment*Times;
                 
                 //Check the collision Time
                 if(Temp_Length<Path_Length)
                 {
                     Collision_TIME_Check = Collision_TIME_Check + (Times-1);
                     Actual_collision = Actual_collision + 1;
                     
                     double *V_Aft_Collision_AIR = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,2);
                     DM_Velocity_Aft_Colliding=V_Aft_Collision_AIR[0];
                     if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)break;
                 }
                 else if(Temp_Length>=Path_Length)
                 {
                     Temp_Length=Path_Length;
                 }
             }
         }

     }
     //=======================
     cout << "Ene_Vel_Air_Final: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
     cout << "Velocity_Air_Final: " << DM_Velocity_Aft_Colliding << endl;
     //cout << "Collision_TIME_Check: " <<Collision_TIME_Check << endl;
     //cout << " Actual_collision: " << Actual_collision << endl;

     RETURN_VALUE[0]=DM_Velocity_Aft_Colliding;RETURN_VALUE[1]=Collision_TIME_Check;RETURN_VALUE[2]=Actual_collision;RETURN_VALUE[3]=(Actual_collision/Collision_TIME_Check);
     return RETURN_VALUE;
 }



double *KS_Collision_Time_EARTH_Aft_velocity(int index, double Sigma_SI_Default, double PY, double PZ, double Velocity, double WIMP_Mass, double *Length_Components) //Velocity(km/s)
{
    static double RETURN_VALUE_1[4];
    static int TEMP_INDEX=-1;
    //cout << "OK! Thanks!" << endl;
    //Theta(Between r and Z) Phi(Between r and X)
    double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
    //cout << "PY: " << PY << "PZ: " << PZ << endl;

    double Path_Lengths_for_earth[13];
    double Density_of_earth_Layer[13];
    double Previous_Total=0;
    double N_Collision_Earth=0;double N_Collision_Cement=0;double N_Collision_Reactor_Wall=0;double N_Collision_Reactor_Water=0;
    double N_Collision_Shielding=0;
    //cout << "Weighted_Atomic_Number" << Weighted_Atomic_Number << endl;
    double Check_Length_times_Density=0;
    static double RETURN_VALUE[6];
    double Length_Temp_1=0;
    //Count the lengths for those layers
    for(int kkk=0; kkk<13; kkk++)
    {
      //  cout << "Earth: " << kkk << endl;
        Density_of_earth_Layer[kkk] = (earth_table[kkk+1][1]+earth_table[kkk][1])/2;
        cout << "Density_of_earth_Layer[kkk]: " << Density_of_earth_Layer[kkk] << endl;
        double Length_Temp  =KS_Earth_Path_Length_for_STS(0,PY,PZ,earth_table[kkk+1][0],Previous_Total);
               Length_Temp_1=KS_Earth_Path_Length_for_STS(1,PY,PZ,earth_table[13][0],Previous_Total);

       // cout << "earth_table[kkk+1][0]: " << earth_table[kkk+1][0] << endl;
        if(Length_Temp>0 and Length_Temp<=6372*2 and PZ<0)
            {
                Path_Lengths_for_earth[kkk] = Length_Temp;
            }
        else if(kkk==12 and Length_Components[0]!=0 and PZ>=0)
            {
                Path_Lengths_for_earth[kkk] = Length_Temp-Length_Components[0];
            }
        else{Path_Lengths_for_earth[kkk]=0;}
        //cout << "kkk: " << kkk << endl;
        //cout << "Path_Lengths_for_earth[kkk]_Original: " << Path_Lengths_for_earth[kkk] << endl;
        //cout << "Density_of_earth_Layer[kkk]_Original: " << Density_of_earth_Layer[kkk] << endl;

        if(Previous_Total>=0 and Previous_Total<=6372*2)Previous_Total = Previous_Total+ Path_Lengths_for_earth[kkk];
        else{Previous_Total = Previous_Total+ 0;}
        //Earth
    }
    //=========================================================
    //Calculate the velocity step-by-step
    double STS_Length[26];double STS_Density[26];

    if(PZ<0)
    {
        //Check the lengths on the path
        for(int kkk=0; kkk<13; kkk++)//Half distance
        {
            int OUTtoIN=12-kkk;
            //cout << "Path_Lengths_for_earth[kkk]: " << Path_Lengths_for_earth[OUTtoIN] << endl;
            if(Path_Lengths_for_earth[OUTtoIN]>0){STS_Length[kkk]=(Path_Lengths_for_earth[OUTtoIN]/2);}
            if(Path_Lengths_for_earth[OUTtoIN]<=0){STS_Length[kkk]=0;}
            STS_Density[kkk]=Density_of_earth_Layer[OUTtoIN];
        }
        for(int kkk=0; kkk<13; kkk++)//Half distance
        {
            int OUTtoIN=12-kkk;
            STS_Length[kkk+13] =STS_Length[OUTtoIN];
            STS_Density[kkk+13]=STS_Density[OUTtoIN];
        }
        STS_Length[25] = STS_Length[25] - Length_Temp_1;//the extra one for more accurate calculation
    }
    if(PZ>=0)
    {
        //Check the lengths on the path
        for(int kkk=0; kkk<25; kkk++)
        {
            STS_Length[kkk]=0;
            STS_Density[kkk]=0;
        }
        STS_Length[25] =Path_Lengths_for_earth[12];
        STS_Density[25]=Density_of_earth_Layer[12];
    }

    
    for(int LLL=0; LLL<26; LLL++)
    {
        //cout << "STS_Length[LLL]:" << STS_Length[LLL] << endl;
        //cout << "STS_Density[LLL]:" << STS_Density[LLL] << endl;
    }
    
    double DM_Velocity_Aft_Colliding=Velocity;
    
    cout << "Ene_Vel_Earth_Initial: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
    cout << "Velocity_Earth_Initial: " << DM_Velocity_Aft_Colliding << endl;

    double Collision_TIME_Check=0;double Actual_collision=0;
    //double *V_Aft_Collision_AIR = Velocity_Aft_collision(AC,DM_mx,Sigma_SI,Velocity[kkk],2);
    double Temp_Length=0;
    double Path_Length=0;
    double Lamda_for_Average=0.001;

    //cout << "================Check_The_turning_Point(Earth)================" << endl;
    static int Earth_Threshold=0;static int Concrete_Threshold=0;static int Water_Threshold=0;static int Shielding_Threshold=0;
    if(TEMP_INDEX!=index)
    {
        for(int EEE=0; EEE<5; EEE++)//i is earth and 2 is air
        {
            int Check   = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI_Default,779.135,1)[3];
            if(Earth_Threshold<Check)Earth_Threshold=Check;
            int Check1  = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI_Default,779.135,3)[3];
            if(Concrete_Threshold<Check1)Concrete_Threshold=Check1;
            int Check2  = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI_Default,779.135,5)[3];
            if(Water_Threshold<Check2)Water_Threshold=Check2;
            int Check3  = Velocity_Aft_collision(5000,WIMP_Mass,Sigma_SI_Default,779.135,6)[3];
            if(Shielding_Threshold<Check3)Shielding_Threshold=Check3;
        }
        cout << "Earth_Threshold: " << Earth_Threshold << endl;
        cout << "Concrete_Threshold: " << Concrete_Threshold << endl;
        cout << "Water_Threshold: " << Water_Threshold << endl;
        cout << "Shielding_Threshold: " << Shielding_Threshold << endl;

    }
    TEMP_INDEX=index;
    //cout << "=======================================================" << endl;
    //====================================Building_Factor=====================================//
    double Cement_Length=Length_Components[1];//km
    double Reactor_Wall_Length=(Length_Components[2]-Length_Components[3]);//km
    double Reactor_Water_Length=Length_Components[3];//km
    
    double Building_Length[3]={Cement_Length,Reactor_Wall_Length,Reactor_Water_Length};//km
    double Building_Density[3]={Density_of_Cement,Density_of_Cement,1};
    double Building_atom_number[3]={Weighted_Atomic_Number_Cement,Weighted_Atomic_Number_Cement,AH2O};
    int Building_Element_Codes[3]={3,3,5};
    int Threshold_for_Building[3]={Concrete_Threshold,Concrete_Threshold,Water_Threshold};
    //====================================Shielding_Factor=====================================//
    double Shielding_Length[4]={15*1e-5,5*1e-5,25*1e-5,5*1e-5};//km
    double Shielding_Density[4]={11.34,7.86,2.34,8.96};
    double Shielding_atom_number[4]={APb,AFe,AB,ACu};

    //======================================================================Check at first======================================================================//
    int Check_Threshold=0;
    for(int kkk=0; kkk<26; kkk++)
    {
        double Segment_Test = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,STS_Density[kkk],Weighted_Atomic_Number);
        double Expectation_Test =  STS_Length[kkk]/Segment_Test*(Possion_GetRandom(2,Lamda_for_Average));
        //cout << "STS_Length[kkk]/Segment_Test: " << STS_Length[kkk]/Segment_Test << endl;
        if(Expectation_Test>=Earth_Threshold)
        {
            cout << " Expectation_Test:" << Expectation_Test << ">" << "Earth_Threshold: " << Earth_Threshold << endl;
            Check_Threshold = Check_Threshold + 1;
            DM_Velocity_Aft_Colliding=1e-5;
        }
    }
    for(int kkk=0; kkk<3; kkk++)
    {
        double Segment_Test = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,Building_Density[kkk],Building_atom_number[kkk]);
        double Expectation_Test =  Building_Length[kkk]/Segment_Test*(Possion_GetRandom(2,Lamda_for_Average));
        //cout << "Building_Length[kkk]/Segment_Test: " << Building_Length[kkk]/Segment_Test << endl;

        if(Expectation_Test>=Threshold_for_Building[kkk])
        {
            cout << " Expectation_Test:" << Expectation_Test << ">" << "Threshold_for_Building: " << Threshold_for_Building[kkk] << endl;
            Check_Threshold = Check_Threshold + 1;
            DM_Velocity_Aft_Colliding=1e-5;
        }
    }
    for(int kkk=0; kkk<4; kkk++)
    {
        double Segment_Test = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,Shielding_Density[kkk],Shielding_atom_number[kkk]);
        double Expectation_Test =  Shielding_Length[kkk]/Segment_Test*(Possion_GetRandom(2,Lamda_for_Average));
        //cout << "Shielding_Length[kkk]/Segment_Test: " << Shielding_Length[kkk]/Segment_Test << endl;

        if(Expectation_Test>=Shielding_Threshold)
        {
            cout << " Expectation_Test:" << Expectation_Test << ">" << "Shielding_Threshold: " << Shielding_Threshold << endl;
            Check_Threshold = Check_Threshold + 1;
            DM_Velocity_Aft_Colliding=1e-5;
        }
    }
    //======================================================================Earth======================================================================//
    if(Check_Threshold==0){
        cout << "DM_Velocity_Aft_Colliding1: " << DM_Velocity_Aft_Colliding << endl;

        for(int kkk=0; kkk<26; kkk++)
        {
            Path_Length= Path_Length + STS_Length[kkk];
            while(Temp_Length<Path_Length)
            {
                double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,STS_Density[kkk],Weighted_Atomic_Number);
                //cout << "Segment_Number: " << STS_Length[kkk]/Segment << endl;
                int Times = Possion_GetRandom_Full(Lamda_for_Average);
                Temp_Length = Temp_Length + Segment*Times;
                
                if(Temp_Length<Path_Length)
                {
                    //Check the collision Time
                    Collision_TIME_Check = Collision_TIME_Check + (Times-1);
                    Actual_collision = Actual_collision + 1;
                    
                    //cout << "Earth=>kkk: " << kkk << endl;
                    double *V_Aft_Collision_AIR = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,1);
                    DM_Velocity_Aft_Colliding=V_Aft_Collision_AIR[0];
                    //cout << "Aft_Vel: " << DM_Velocity_Aft_Colliding << endl;
                    if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)break;
                }
                else if(Temp_Length>=Path_Length)
                {
                    Temp_Length=Path_Length;
                }
            }
        }
        //======================================================================building======================================================================//
            cout << "DM_Velocity_Aft_Colliding2: " << DM_Velocity_Aft_Colliding << endl;

        Path_Length = 0; Temp_Length = 0;
        for(int kkk=0; kkk<3; kkk++)
        {
            
            Path_Length= Path_Length + Building_Length[kkk];
            while(Temp_Length<Path_Length)
            {
                double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,Building_Density[kkk],Building_atom_number[kkk]);
                //cout << "Segment_Number: " << STS_Length[kkk]/Segment << endl;
                int Times = Possion_GetRandom_Full(Lamda_for_Average);
                Temp_Length = Temp_Length + Segment*Times;
                
                if(Temp_Length<Path_Length)
                {
                    //Check the collision Time
                    Collision_TIME_Check = Collision_TIME_Check + (Times-1);
                    Actual_collision = Actual_collision + 1;
                    //cout << "Build=>kkk: " << kkk << endl;

                    double *V_Aft_Collision_AIR = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,Building_Element_Codes[kkk]);
                    DM_Velocity_Aft_Colliding=V_Aft_Collision_AIR[0];
                    //cout << "Aft_Vel: " << DM_Velocity_Aft_Colliding << endl;
                    if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)break;
                }
                else if(Temp_Length>=Path_Length)
                {
                    Temp_Length=Path_Length;
                }
            }
        }
         
        //cout << "Ene_Vel_Building: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
        //cout << "Velocity_Aft_Building: " << DM_Velocity_Aft_Colliding << endl;

        //======================================================================Shieling======================================================================//
                cout << "DM_Velocity_Aft_Colliding3: " << DM_Velocity_Aft_Colliding << endl;

        Path_Length = 0; Temp_Length = 0;
        

        for(int kkk=0; kkk<4; kkk++)
        {
            Path_Length= Path_Length + Shielding_Length[kkk];
            while(Temp_Length<Path_Length)
            {
                double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,Shielding_Density[kkk],Shielding_atom_number[kkk]);
                //cout << "Segment_Number: " << STS_Length[kkk]/Segment << endl;
                int Times = Possion_GetRandom_Full(Lamda_for_Average);
                Temp_Length = Temp_Length + Segment*Times;
                
                if(Temp_Length<Path_Length)
                {
                    //Check the collision Time
                    Collision_TIME_Check = Collision_TIME_Check + (Times-1);
                    Actual_collision = Actual_collision + 1;
                    //cout << "Shielding=>kkk: " << kkk << endl;

                    double *V_Aft_Collision_AIR = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,6);
                    DM_Velocity_Aft_Colliding=V_Aft_Collision_AIR[0];
                    //cout << "Aft_Vel: " << DM_Velocity_Aft_Colliding << endl;
                    if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)break;
                }
                else if(Temp_Length>=Path_Length)
                {
                    Temp_Length=Path_Length;
                }
            }
        }
    }
    //cout << "Average_Expectation_Original: " << (Shielding_Length[kkk]*Shielding_Density[kkk]*1e5)/(unified_atomic_mass_g*(Shielding_atom_number[kkk]))*(total_Sigma(1,Velocity,Sigma_SI_Default,WIMP_Mass,Shielding_atom_number[kkk])) << endl;

    //cout << "Ene_Vel_Shielding: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
    //cout << "Velocity_Aft_Shielding: " << DM_Velocity_Aft_Colliding << endl;

    RETURN_VALUE_1[0]=DM_Velocity_Aft_Colliding;RETURN_VALUE_1[1]=Collision_TIME_Check;RETURN_VALUE_1[2]=Actual_collision;RETURN_VALUE_1[3]=(Actual_collision/Collision_TIME_Check);

    cout << "Ene_Vel_Earth_Aft: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
    cout << "Velocity_Aft_Earth: " << DM_Velocity_Aft_Colliding << endl;
    //cout << "Collision_TIME_Check: " << Collision_TIME_Check << endl;
    //cout << "Actual_collision: " << Actual_collision << endl;
    //cout << "(Actual_collision/Collision_TIME_Check): " << (Actual_collision/Collision_TIME_Check) << endl;
    
    RETURN_VALUE_1[0]=DM_Velocity_Aft_Colliding;RETURN_VALUE_1[1]=Collision_TIME_Check;RETURN_VALUE_1[2]=Actual_collision;RETURN_VALUE_1[3]=(Actual_collision/Collision_TIME_Check);
    return RETURN_VALUE_1;
}




 double *KS_Collision_Time_ATM_Aft_velocity_with_angle(int Straight_or_scattered, int Event, int index, double Sigma_SI_Default, double PY, double PZ, double Velocity, double WIMP_Mass, double *Length_Components, double Earth_Path_Length, double *Direction) //Velocity(km/s)
 {
     //double Direction[3]   ={ 1e-50, 1e-50, PZ };
     double Direction_VT[3]={ 0, 0, 0 };

     static double RETURN_VALUE[20];
     static int TEMP_INDEX=-1;
     //cout << "OK! Thanks!" << endl;
     //Theta(Between r and Z) Phi(Between r and X)
     //cout << "PY: " << PY << "PZ: " << PZ << endl;
     double Path_Vector[3]; double Path_Vector_Unit[3];double Scaling=0;
     double Lab_CARC[3]={0,0,0};// Original point X,Y,Z Function=>TMath::Pi()
     
     double Path_Lengths_for_atmosphere[19];
     double Density_of_Atmosphere_Layer[19];
     double Previous_Total_AIR=0;
     double Previous_Total = ( Length_Components[1]+Length_Components[2]+Earth_Path_Length);
     //cout << "Earth_Path_Length: " << Earth_Path_Length  << endl;
     //cout << "Previous_Total: " << Previous_Total << endl;

     //Calculate the lengths for different layers
     for(int kkk=0; kkk<19; kkk++)
     {
         //cout << "AIR: " << kkk << endl;
         Density_of_Atmosphere_Layer[kkk] = (atm_table[kkk][4]+atm_table[kkk+1][4])/2;
         //cout << "Density_of_Atmosphere_Layer[kkk]: " << Density_of_Atmosphere_Layer[kkk] << endl;
         Path_Lengths_for_atmosphere[kkk] = KS_Air_Path_Length(PY,PZ,(atm_table[kkk+1][0]/1000),Previous_Total_AIR);
         //cout << "(atm_table[kkk+1][0]/1000): " << (atm_table[kkk+1][0]/1000) << endl;
         Previous_Total_AIR = Previous_Total_AIR+ Path_Lengths_for_atmosphere[kkk];
         if(kkk==0) Path_Lengths_for_atmosphere[kkk] = Path_Lengths_for_atmosphere[kkk] - Previous_Total;
         //cout << "Path_Lengths_for_atmosphere[kkk]: " << Path_Lengths_for_atmosphere[kkk] << "km " << endl;
         //cout << " Previous_Total_AIR: " <<  Previous_Total_AIR << endl;
     }
     

     //Step-to-Step(STS) collision
     double STS_Length[19];double STS_Density[19];
     //Check the lengths on the path
     for(int kkk=0; kkk<19; kkk++)
     {
         int OUTtoIN=18-kkk;
         if(Path_Lengths_for_atmosphere[OUTtoIN]>0){STS_Length[kkk]=(Path_Lengths_for_atmosphere[OUTtoIN]);}
         STS_Density[kkk]=kg_perm3_to_g_percm3(Density_of_Atmosphere_Layer[OUTtoIN]);
         //cout << "kkk: " << kkk << endl;
         //cout << "STS_Length[kkk]: " << STS_Length[kkk] << endl;
         //cout << "STS_Density[kkk]: " << STS_Density[kkk] << endl;
     }
     
     double DM_Velocity_Aft_Colliding=Velocity;
     //cout << "DM_Velocity_Aft_Colliding_Air: " << DM_Velocity_Aft_Colliding << endl;
     double Collision_TIME_Check=0;double Actual_collision=0;
     //double *V_Aft_Collision_AIR = Velocity_Aft_collision(AC,DM_mx,Sigma_SI,Velocity[kkk],2);
     double Temp_Length=0;
     double Path_Length=0;
     double Check_the_dropping_point=0;
     double Lamda_for_Average=0.001;
     
     cout << "Ene_Vel_Air_Initial: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
     cout << "Velocity_Air_Initial: " << DM_Velocity_Aft_Colliding << endl;

     static int Air_Threshold=0;
     cout << "index: " << index << endl;
     //cout << "================Check_The_turning_Point(Air)================" << endl;
     //cout << "=======================================================" << endl;

     //===========================
     int Check_Threshold=0;
     int Collision_Time=0;
     int Arrival_or_not=0;
     /*
     for(int kkk=0; kkk<19; kkk++)
     {
         double Segment_Test = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,STS_Density[kkk],14.99);
         double Expectation_Test =  STS_Length[kkk]/Segment_Test*(Possion_GetRandom(2,Lamda_for_Average));
         //cout << "STS_Length[kkk]/Segment_Test: " << STS_Length[kkk]/Segment_Test << endl;

         if(Expectation_Test>Air_Threshold)
         {
             Check_Threshold  = Check_Threshold + 1;
             DM_Velocity_Aft_Colliding=1e-5;
         }
     }
      */
     cout << "Check_threshold " << Check_Threshold << endl;
     //
     
     TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,600,400);
        TGraph2D *dt = new TGraph2D();

     std::setprecision(9);
     double Track_Point[3]={0,0,-6371-80};//Start from the bottom point

     
     double  Path_Original= 0;
     double  Ratio_of_Energy_Loss_to_Atom = 0;
     
     if(Check_Threshold>-1)
     {
         double Center_of_Earth[3]={0,0,0};
         double Track_Point_interim[3]={0,0,0};//Applied for checking the points after the next step
         double Track_Point_interim_half[3]={0,0,0};//Applied for checking the points after the next "half" step to check the crossing case
         double Length_Interval[20];double Density_Index[20];int Length_Index[20];
         
         Length_Interval[0] = abs(-6371-80);Density_Index[0] = 0;Length_Index[0] = 0;//The first element is the reference point
         
         for(int kkk=0; kkk<19; kkk++)
         {
             Path_Length= Path_Length + STS_Length[kkk]; //From the bottom line to sum up
             Length_Interval[kkk+1] = abs((-6371-80)+Path_Length);
             Length_Index[kkk+1]    = kkk+1;
             Density_Index[kkk+1]   = STS_Density[kkk];
             
         }
         Length_Interval[19] = abs(-6371);//The first element is the reference point

         for(int kkk=0; kkk<20; kkk++)
         {
             
             /*cout << "=================" << endl;
             cout << "Length_Interval: " << Length_Interval[kkk] << endl;
             cout << "Length_Index: " << Length_Index[kkk] <<endl;
             cout << "Density_Index: " << Density_Index[kkk] << endl;
             cout << "=================" << endl;*/
         }
         
         int Step=0;
         
         Path_Original= Scale_for_scattering_process(1,6371,Track_Point,Direction);


         while(Distant(Track_Point,Center_of_Earth)-6371>1e-10 and Distant(Track_Point,Center_of_Earth)-6451<=1e-10 and Path_Original!=0)
         {
             dt->SetPoint(Step,Track_Point[0],Track_Point[1],Track_Point[2]);
             
             cout << "=================Distant(Track_Point,Center_of_Earth) Started================ " << Distant(Track_Point,Center_of_Earth) << endl;
             cout << "=========Path_Original========: " << Path_Original << endl;

             cout << "Event: " << Event << endl;
             cout << "Direction[0]: "   << Direction[0] << endl;
             cout << "Direction[1]: "   << Direction[1] << endl;
             cout << "Direction[2]: "   << Direction[2] << endl;

             cout << "Track_Point_Now[0]: " << Track_Point[0] << endl;
             cout << "Track_Point_Now[1]: " << Track_Point[1] << endl;
             cout << "Track_Point_Now[2]: " << Track_Point[2] << endl;
             
             int Place_Now=0;
             for(int kkk=0; kkk<20;kkk++)
             {
        if( Distant(Track_Point,Center_of_Earth)-Length_Interval[kkk+1]>1e-2 and (Distant(Track_Point,Center_of_Earth))-Length_Interval[kkk]<=1e-2)
                {
                    cout << "Distant(Track_Point,Center_of_Earth)-Length_Interval[kkk+1]" << Distant(Track_Point,Center_of_Earth)-Length_Interval[kkk+1] << endl;
                    cout << "Length_Interval[kkk]: "   << Length_Interval[kkk] << endl;
                    cout << "Length_Interval[kkk+1]: " << Length_Interval[kkk+1] << endl;
                    cout << "Distant(Track_Point,Center_of_Earth): " << Distant(Track_Point,Center_of_Earth) << endl;

                    cout << Length_Interval[kkk] << ">=" << Distant(Track_Point,Center_of_Earth) << "> " << Length_Interval[kkk+1] << endl;
                    cout << "Starting_Point " << Distant(Track_Point,Center_of_Earth) << endl;

                    Place_Now= Length_Index[kkk+1];
                }
             }
             cout << "Place_now" << Place_Now << endl;

             double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,Density_Index[Place_Now],14.99);
             cout << "Lamda_for_Average: " << Lamda_for_Average << endl;
             cout << "WIMP_Mass: " << WIMP_Mass << endl;
             cout << "DM_Velocity_Aft_Colliding: " << DM_Velocity_Aft_Colliding << endl;
             cout << "Density_Index[Place_Now]: " << Density_Index[Place_Now] << endl;
             cout << "Weighted_Atomic_Number: " << Weighted_Atomic_Number << endl;
             cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
             cout << "Segment:  " << Segment << endl;
             if(Segment>=1e20) Segment=0;
             int Times = Possion_GetRandom_Full(Lamda_for_Average);
             double Sprint_Length = Segment*Times;
             cout << "Sprint_Length: " << Sprint_Length << "WIMP_Mass: " << WIMP_Mass << endl;

             //Toward Inner or Outer
             double Inner_Scaling=0;double Outer_Scaling=0;
             double Final_Length_Scaling=0;
             

             if( ((Length_Interval[Place_Now-1])-(Distant(Track_Point,Center_of_Earth)))> 1e-2 )
             {
                 cout << "YES1" << endl;
                 Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now]  ,Track_Point,Direction);
                 Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-1],Track_Point,Direction);
             }
             else if( ( (Length_Interval[Place_Now-1]) - (Distant(Track_Point,Center_of_Earth)) ) < 1e-2)
             {
                 cout << "YES2" << endl;
                 if(Place_Now==1){
                 Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now]  ,Track_Point,Direction);
                 Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-1],Track_Point,Direction);}
                 else{
                     cout << "Length_Interval[Place_Now]: " << Length_Interval[Place_Now] << endl;
                 Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now]  ,Track_Point,Direction);
                 Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-2],Track_Point,Direction);}
                 
             }
             
             if(Inner_Scaling>1e20)Inner_Scaling=0;
             if(Outer_Scaling>1e20)Outer_Scaling=0;


             Final_Length_Scaling = Selection_for_STS(Inner_Scaling,Outer_Scaling);
                         
             if((Inner_Scaling<0 and Outer_Scaling==0))
             {
                 cout << "CASE5" << endl;
                 double Special_Length=0;
                 if(Place_Now==1)Special_Length = abs( (Length_Interval[Place_Now-1]) - (Distant(Track_Point,Center_of_Earth)) );
                 else{Special_Length = abs( (Length_Interval[Place_Now-2]) - (Distant(Track_Point,Center_of_Earth)) );}
                 Final_Length_Scaling = Special_Length;
             }

             else if((Inner_Scaling==0 and Outer_Scaling<0))
             {
                 cout << "CASE6" << endl;
                 double Special_Length=0;
                 Special_Length = abs( (Length_Interval[Place_Now]) - (Distant(Track_Point,Center_of_Earth)) );
                 Final_Length_Scaling = Special_Length;
             }

             
             cout << "Final_Length_Scaling: " << Final_Length_Scaling << endl;
             //====See the comparison length
             if(Sprint_Length>Final_Length_Scaling)
             {
                 Track_Point[0] = Track_Point[0]+(Final_Length_Scaling*Direction[0]);
                 Track_Point[1] = Track_Point[1]+(Final_Length_Scaling*Direction[1]);
                 Track_Point[2] = Track_Point[2]+(Final_Length_Scaling*Direction[2]);

                 Temp_Length = Temp_Length + Final_Length_Scaling;
                 cout << "GAGA1" << endl;

                 Step = Step + 1;
             }

             else if(Sprint_Length<=Final_Length_Scaling)
             {
                 Track_Point[0] = Track_Point[0]+(Sprint_Length*Direction[0]);
                 Track_Point[1] = Track_Point[1]+(Sprint_Length*Direction[1]);
                 Track_Point[2] = Track_Point[2]+(Sprint_Length*Direction[2]);
                 
                 Temp_Length = Temp_Length + Sprint_Length;

                 cout << "GAGA2" << endl;

                 cout << "Bef_Vel: " << DM_Velocity_Aft_Colliding << endl;

                 double *V_Aft_Collision_AIR = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,2);
                 double Energy_Loss_to_Atom   = Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) - Energy_DM(WIMP_Mass,V_Aft_Collision_AIR[0]*1e3/3e8);
                 double Ratio_of_Energy_Loss_to_Atom = Energy_Loss_to_Atom/Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8);
                 cout << "Ratio_of_Energy_Loss_to_Atom: " << Ratio_of_Energy_Loss_to_Atom << endl;
                 DM_Velocity_Aft_Colliding=V_Aft_Collision_AIR[0];
                 cout << "Aft_Vel: " << DM_Velocity_Aft_Colliding << endl;
                 
                 if(Straight_or_scattered==1){
                 double *Direction_aft = Aft_scatterd_Direction(2,WIMP_Mass,DM_Velocity_Aft_Colliding,Direction,Ratio_of_Energy_Loss_to_Atom);
                 Direction[0] = Direction_aft[0];Direction[1] = Direction_aft[1];Direction[2] = Direction_aft[2];}
                 
                 if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)
                 {
                     cout << "Energy Threshold" << endl;
                     break;
                 }
                 Step = Step + 1;
                 Collision_Time = Collision_Time + 1;
             }
             

             if( (Distant_Criteria(Track_Point,Center_of_Earth)-6371 <=1e-10) and Step!=0)
             {
                 cout << "In La~" << endl;
                 Arrival_or_not = 1;
                 cout << "Temp_Length: " << Temp_Length << endl;
                 break;
             }
             else if(  (6451 - (Distant_Criteria(Track_Point,Center_of_Earth))) <=1e-10 and Step!=0)
             {
                 cout << "Out La~" << endl;
                 Arrival_or_not = 0;
                 cout << "Temp_Length: " << Temp_Length << endl;
                 break;
             }

         }

         cout << "=========Temp_Length========: " << Temp_Length << endl;
         cout << "=========Path_Original========: " << Path_Original << endl;
         cout << "=================Distant(Track_Point,Center_of_Earth) ended================ " << Distant(Track_Point,Center_of_Earth) << endl;

         cout << "Track_Point_Now[0]: " << Track_Point[0] << endl;
         cout << "Track_Point_Now[1]: " << Track_Point[1] << endl;
         cout << "Track_Point_Now[2]: " << Track_Point[2] << endl;

         cout << "Step: " << Step << endl;
         if(Path_Original<=0) cout << "OK THEN THE EVENT IS OUT OF OUR CATEGORY" << endl;

         //cout << "=========Temp_Length========: " << Temp_Length << endl;

     }
     /*
     if(Temp_Length!=0)
     {
     gStyle->SetPalette(1);
     dt->SetMarkerStyle(20);
     dt->Draw("pcolLINE");
     }
      */
     //=======================
     //cout << "Ene_Vel_Air_Final: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
     //cout << "Velocity_Air_Final: " << DM_Velocity_Aft_Colliding << endl;
     //cout << "Collision_TIME_Check: " <<Collision_TIME_Check << endl;
     //cout << " Actual_collision: " << Actual_collision << endl;
     if(Check_Threshold==0)
     {
         if(Path_Original!=0){RETURN_VALUE[1]=(Temp_Length/Path_Original);}
         else{RETURN_VALUE[1]=0;}
     RETURN_VALUE[0]=DM_Velocity_Aft_Colliding;RETURN_VALUE[2]=Collision_Time;
     RETURN_VALUE[3]=Track_Point[0];RETURN_VALUE[4]=Track_Point[1];RETURN_VALUE[5]=Track_Point[2];RETURN_VALUE[6]=Arrival_or_not;
     RETURN_VALUE[7]=Direction[0];RETURN_VALUE[8]=Direction[1];RETURN_VALUE[9]=Direction[2];
     RETURN_VALUE[10]=Path_Original;RETURN_VALUE[11]=Temp_Length;RETURN_VALUE[12]=Collision_Time;
     }

     if(Check_Threshold!=0)
     {
    RETURN_VALUE[0]=0;RETURN_VALUE[1]=0;
     }

     return RETURN_VALUE;
 }


double *KS_Collision_Time_Earth_Aft_velocity_with_angle(int Straight_or_scattered, int Event, int index, double Sigma_SI_Default, double Velocity, double WIMP_Mass, double *Direction, double *Intermediate_Point) //Velocity(km/s)
{
    //double Direction[3]   ={ 1e-50, 1e-50, PZ };
    double Direction_VT[3]={ 0, 0, 0 };

    static double RETURN_VALUE[20];
    static int TEMP_INDEX=-1;
    
    double DM_Velocity_Aft_Colliding=Velocity;
    //cout << "DM_Velocity_Aft_Colliding_Air: " << DM_Velocity_Aft_Colliding << endl;
    double Collision_TIME_Check=0;double Actual_collision=0;
    //double *V_Aft_Collision_AIR = Velocity_Aft_collision(AC,DM_mx,Sigma_SI,Velocity[kkk],2);
    double Temp_Length=0;
    double Path_Length=0;
    double Check_the_dropping_point=0;
    double Lamda_for_Average=0.001;
    
    cout << "Ene_Vel_Air_Initial: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
    cout << "Velocity_Air_Initial: " << DM_Velocity_Aft_Colliding << endl;

    static int Earth_Threshold=0;

    //===========================
    int Check_Threshold=0;
    int Collision_Time=0;
    int Arrival_or_not = 0;
    TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,600,400);
       TGraph2D *dt = new TGraph2D();
     
    double Track_Point[3]={Intermediate_Point[0],Intermediate_Point[1],Intermediate_Point[2]};//Start from the bottom point

    double  Path_Original=0;
    if(Check_Threshold>-1)
    {
        double Center_of_Earth[3]={0,0,0};
        double Track_Point_interim[3]={0,0,0};//Applied for checking the points after the next step
        double Track_Point_interim_half[3]={0,0,0};//Applied for checking the points after the next "half" step to check the crossing case
        double Length_Interval[14];double Density_Index[14];int Length_Index[14];
        
        Length_Interval[0] = abs(-6371);Density_Index[0] = 0;Length_Index[0] = 0;//The first element is the reference point
        
        for(int kkk=0; kkk<13; kkk++)
        {
            Path_Length= Path_Length + (earth_table[13-kkk][0]-earth_table[12-kkk][0]); //From the bottom line to sum up
            Length_Interval[kkk+1] = abs((-6371)+Path_Length);
            Length_Index[kkk+1]    = kkk+1;
            Density_Index[kkk+1]   = (earth_table[13-kkk][1]+earth_table[12-kkk][1])*(0.5);
            
        }

        for(int kkk=0; kkk<14; kkk++)
        {
            
            cout << "=================" << endl;
            cout << "Length_Interval: " << Length_Interval[kkk] << endl;
            cout << "Length_Index: " << Length_Index[kkk] <<endl;
            cout << "Density_Index: " << Density_Index[kkk] << endl;
            cout << "=================" << endl;
        }
        
        int Step=0;
        
        Path_Original= Scale_for_scattering_process(1,6371,Track_Point,Direction);
        
        while(Distant(Track_Point,Center_of_Earth)-6371<1e-10)
        {
            dt->SetPoint(Step,Track_Point[0],Track_Point[1],Track_Point[2]);
            
            cout << "=================Distant(Track_Point,Center_of_Earth) Started================ " << Distant(Track_Point,Center_of_Earth) << endl;

            cout << "=========Path_Original========: " << Path_Original << endl;

            cout << "Event: " << Event << endl;
            cout << "Direction[0]: "   << Direction[0] << endl;
            cout << "Direction[1]: "   << Direction[1] << endl;
            cout << "Direction[2]: "   << Direction[2] << endl;

            cout << "Track_Point_Now[0]: " << Track_Point[0] << endl;
            cout << "Track_Point_Now[1]: " << Track_Point[1] << endl;
            cout << "Track_Point_Now[2]: " << Track_Point[2] << endl;
            
            int Place_Now=0;
            for(int kkk=0; kkk<14;kkk++)
            {
               if( Distant(Track_Point,Center_of_Earth)-Length_Interval[kkk+1]>1e-2 and (Distant(Track_Point,Center_of_Earth)-Length_Interval[kkk])<=1e-2)
               {
                   cout << "Distant(Track_Point,Center_of_Earth)-Length_Interval[kkk+1]" << Distant(Track_Point,Center_of_Earth)-Length_Interval[kkk+1] << endl;
                   cout << "Length_Interval[kkk]: "   << Length_Interval[kkk] << endl;
                   cout << "Length_Interval[kkk+1]: " << Length_Interval[kkk+1] << endl;
                   cout << "Distant(Track_Point,Center_of_Earth): " << Distant(Track_Point,Center_of_Earth) << endl;

                   cout << Length_Interval[kkk] << ">=" << Distant(Track_Point,Center_of_Earth) << "> " << Length_Interval[kkk+1] << endl;
                   cout << "Starting_Point " << Distant(Track_Point,Center_of_Earth) << endl;
                   Place_Now= Length_Index[kkk+1];
                   cout << "Now we are standing at//1: " << Place_Now << endl;
               }
            }
            cout << "Now we are standing at//2: " << Place_Now << endl;


            double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,Density_Index[Place_Now],Weighted_Atomic_Number);
            cout << "Lamda_for_Average: " << Lamda_for_Average << endl;
            cout << "WIMP_Mass: " << WIMP_Mass << endl;
            cout << "DM_Velocity_Aft_Colliding: " << DM_Velocity_Aft_Colliding << endl;
            cout << "Density_Index[Place_Now]: " << Density_Index[Place_Now] << endl;
            cout << "Weighted_Atomic_Number: " << Weighted_Atomic_Number << endl;
            cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
            cout << "Segment:  " << Segment << endl;
            if(Segment>=1e20) Segment=0;
            int Times = Possion_GetRandom_Full(Lamda_for_Average);

            double Sprint_Length = Segment*Times;
            cout << "Sprint_Length: " << Sprint_Length << "WIMP_Mass: " << WIMP_Mass << endl;
            //Toward Inner or Outer
            double Inner_Scaling=0;double Outer_Scaling=0;
            double Final_Length_Scaling=0;
            
            //ON THE SURFACE OF THE LATER BETWEEN DIFFERENT DENSITIES
            if( ((Length_Interval[Place_Now-1])-(Distant(Track_Point,Center_of_Earth)))> 1e-2 )
            {
                if(Place_Now<13){
                    cout << "OK3" << endl;
                Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now]  ,Track_Point,Direction);
                Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-1],Track_Point,Direction);}
                else{
                    cout << "OK4" << endl;
                    cout << "Length_Interval[Place_Now-1]: " << Length_Interval[Place_Now-1] << endl;
                Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-1]  ,Track_Point,Direction);
                Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-2]  ,Track_Point,Direction);}
            }
            else if( ( (Length_Interval[Place_Now-1]) - (Distant(Track_Point,Center_of_Earth)) ) < 1e-2)
            {
                if(Place_Now==1){
                    cout << "L1" << endl;
                Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now]  ,Track_Point,Direction);
                Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-1],Track_Point,Direction);}
                else if(Place_Now==13){
                    cout << "L1" << endl;
                Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-1]  ,Track_Point,Direction);
                Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-2]  ,Track_Point,Direction);}
                else{
                    cout << "L3" << endl;
                Inner_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now]  ,Track_Point,Direction);
                Outer_Scaling = Scale_for_scattering_process(1,Length_Interval[Place_Now-2],Track_Point,Direction);}
                
            }
            
            //DECIDE TOWARD INNER OR OUTER
            cout << "Inner_Scaling: " << Inner_Scaling << endl;
            cout << "Outer_Scaling: " << Outer_Scaling << endl;
            
            if(Inner_Scaling>1e20)Inner_Scaling=0;
            if(Outer_Scaling>1e20)Outer_Scaling=0;


            Final_Length_Scaling = Selection_for_STS(Inner_Scaling,Outer_Scaling);
                        
            if((Inner_Scaling<0 and Outer_Scaling==0))
            {
                cout << "CASE5" << endl;
                double Special_Length=0;
                if(Place_Now==1)Special_Length = abs( (Length_Interval[Place_Now-1]) - (Distant(Track_Point,Center_of_Earth)) );
                else{Special_Length = abs( (Length_Interval[Place_Now-2]) - (Distant(Track_Point,Center_of_Earth)) );}
                Final_Length_Scaling = Special_Length;
            }

            else if((Inner_Scaling==0 and Outer_Scaling<0))
            {
                cout << "CASE6" << endl;
                double Special_Length=0;
                if(Place_Now==13)Special_Length = abs( (Length_Interval[Place_Now-1]) - (Distant(Track_Point,Center_of_Earth)) );
                else{Special_Length = abs( (Length_Interval[Place_Now]) - (Distant(Track_Point,Center_of_Earth)) );}
                Final_Length_Scaling = Special_Length;
            }

            cout << "Final_Length_Scaling: " << Final_Length_Scaling << endl;

            if(Sprint_Length>Final_Length_Scaling)
            {
                Track_Point[0] = Track_Point[0]+(Final_Length_Scaling*Direction[0]);
                Track_Point[1] = Track_Point[1]+(Final_Length_Scaling*Direction[1]);
                Track_Point[2] = Track_Point[2]+(Final_Length_Scaling*Direction[2]);
                
                Temp_Length = Temp_Length + Final_Length_Scaling;
                
                cout << "CASE1" << endl;
                cout << "Final_Length_Scaling: " << Final_Length_Scaling << endl;
                cout << "Temp_Length: " << Temp_Length << endl;

                Step = Step + 1;
            }

            else if(Sprint_Length<=Final_Length_Scaling)
            {
                Track_Point[0] = Track_Point[0]+(Sprint_Length*Direction[0]);
                Track_Point[1] = Track_Point[1]+(Sprint_Length*Direction[1]);
                Track_Point[2] = Track_Point[2]+(Sprint_Length*Direction[2]);
                
                Temp_Length = Temp_Length + Sprint_Length;

                cout << "CASE2" << endl;
                cout << "Sprint: " << Sprint_Length << endl;
                cout << "Temp_Length: " << Temp_Length << endl;

                double *V_Aft_Collision_Earth = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,1);
                double Energy_Loss_to_Atom   = Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) - Energy_DM(WIMP_Mass,V_Aft_Collision_Earth[0]*1e3/3e8);
                double Ratio_of_Energy_Loss_to_Atom = Energy_Loss_to_Atom/Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8);
                cout << "Ratio_of_Energy_Loss_to_Atom: " << Ratio_of_Energy_Loss_to_Atom << endl;
                DM_Velocity_Aft_Colliding=V_Aft_Collision_Earth[0];
                cout << "Aft_Vel: " << DM_Velocity_Aft_Colliding << endl;
                
                if(Straight_or_scattered==1){
                double *Direction_aft = Aft_scatterd_Direction(1,WIMP_Mass,DM_Velocity_Aft_Colliding,Direction,Ratio_of_Energy_Loss_to_Atom);
                Direction[0] = Direction_aft[0];Direction[1] = Direction_aft[1];Direction[2] = Direction_aft[2];}
                
                if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)
                {
                    cout << "Energy Threshold" << endl;
                    break;
                }
                Step = Step + 1;
                Collision_Time = Collision_Time + 1;
            }
            
            cout << "///Distant(Track_Point,Center_of_Earth): " << (Distant(Track_Point,Center_of_Earth))<< endl;

            if( (6371 -Distant(Track_Point,Center_of_Earth)) <1e-10 and Step!=0)
            {
                cout << "(Distant(Track_Point,Center_of_Earth)-6371: " << (Distant(Track_Point,Center_of_Earth)-6371) << endl;
                cout << "Got you on earth!" << endl;
                cout << "Step: " << Step << endl;
                cout << "Temp_Length: " << Temp_Length << endl;
                Arrival_or_not = 1;
                break;
            }
            cout << "Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8): " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;

        }

        cout << "Track_Point[0]: " << Track_Point[0] << endl;
        cout << "Track_Point[1]: " << Track_Point[1] << endl;
        cout << "Track_Point[2]: " << Track_Point[2] << endl;
        
        cout << "=========Temp_Length========: " << Temp_Length << endl;
        cout << "=========Path_Original========: " << Path_Original << endl;
        
        cout << "=========Temp_Length-Path_Original=========: " << Temp_Length-Path_Original << endl;

        cout << "Distant(Track_Point,Center_of_Earth): " << Distant(Track_Point,Center_of_Earth) << endl;

        cout << "Step: " << Step << endl;
        cout << "Collision_Time: " << Collision_Time << endl;

        //cout << "=========Temp_Length========: " << Temp_Length << endl;
        
    }
    if(Temp_Length!=0)
    {
    gStyle->SetPalette(1);
    dt->SetMarkerStyle(20);
    dt->Draw("pcolLINE");
    }
    //=======================
    //cout << "Ene_Vel_Air_Final: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
    //cout << "Velocity_Air_Final: " << DM_Velocity_Aft_Colliding << endl;
    //cout << "Collision_TIME_Check: " <<Collision_TIME_Check << endl;
    //cout << " Actual_collision: " << Actual_collision << endl;
    if(Check_Threshold==0)
    {
        if(Path_Original!=0){RETURN_VALUE[1]=(Temp_Length/Path_Original);}
        else{RETURN_VALUE[1]=0;}
    RETURN_VALUE[0]=DM_Velocity_Aft_Colliding;RETURN_VALUE[2]=Collision_Time;
    RETURN_VALUE[3]=Arrival_or_not;
    RETURN_VALUE[4]=Path_Original;RETURN_VALUE[5]=Temp_Length;RETURN_VALUE[6]=Collision_Time;
    }
    
    if(Check_Threshold!=0)
    {
   RETURN_VALUE[0]=0;RETURN_VALUE[1]=0;RETURN_VALUE[2]=Earth_Threshold;
    }

    return RETURN_VALUE;
}

double *KS_Collision_Time_Building_Aft_velocity_with_angle(int Straight_or_scattered, double Sigma_SI_Default, double Velocity, double WIMP_Mass, double *Direction, double Length_of_Material, double Density_of_Material) //Velocity(km/s) Density(g/cm^3)
{
    //double Direction[3]   ={ 1e-50, 1e-50, PZ };
    double Direction_VT[3]={0,0,0};//Used for the calculation
    double Track_Point[3]={0,0,0};//Run from the bottom
    double LOF_km = Length_of_Material*1e-3;
    static double RETURN_VALUE[20];//Return the value back
    
    double DM_Velocity_Aft_Colliding=Velocity;double Lamda_for_Average=0.001;

    //===========================
    int Collision_Time=0;
    int Arrival_or_not = 0;
    TCanvas *c1 = new TCanvas("c1","Graph2D example",0,0,600,400);
       TGraph2D *dt = new TGraph2D();
     
    double Temp_Length=0;
    double Path_Original=10;
    int Check_Threshold=0;
    
    if(Check_Threshold>-1)
    {
        int Step=0;
        while( (Track_Point[2]-LOF_km)<1e-10)
        {
            dt->SetPoint(Step,Track_Point[0],Track_Point[1],Track_Point[2]);
            
            cout << "Direction[0]: "   << Direction[0] << endl;
            cout << "Direction[1]: "   << Direction[1] << endl;
            cout << "Direction[2]: "   << Direction[2] << endl;

            cout << "Track_Point_Now[0]: " << Track_Point[0] << endl;
            cout << "Track_Point_Now[1]: " << Track_Point[1] << endl;
            cout << "Track_Point_Now[2]: " << Track_Point[2] << endl;
            
            double Segment = Length_for_asking_the_collision(Lamda_for_Average,WIMP_Mass,DM_Velocity_Aft_Colliding,Sigma_SI_Default,Density_of_Material,Weighted_Atomic_Number);
            /*
            cout << "Lamda_for_Average: " << Lamda_for_Average << endl;
            cout << "WIMP_Mass: " << WIMP_Mass << endl;
            cout << "DM_Velocity_Aft_Colliding: " << DM_Velocity_Aft_Colliding << endl;
            cout << "Density_Index[Place_Now]: " << Density_Index[Place_Now] << endl;
            cout << "Weighted_Atomic_Number: " << Weighted_Atomic_Number << endl;
            cout << "Sigma_SI_Default: " << Sigma_SI_Default << endl;
            cout << "Segment:  " << Segment << endl;
             */
            if(Segment>=1e20) Segment=0;
            int Times = Possion_GetRandom_Full(Lamda_for_Average);
            double Sprint_Length = Segment*Times;
            cout << "Sprint_Length: " << Sprint_Length << "WIMP_Mass: " << WIMP_Mass << endl;
            
            double Critera_of_surpassing = Track_Point[2]+Sprint_Length*Direction[2];
            
            if(Critera_of_surpassing>LOF_km)
            {
                double Final_Length_Scaling = (1 - Track_Point[2]) / Direction[2] ;//Final_Path
                
                Track_Point[0] = Track_Point[0]+(Final_Length_Scaling*Direction[0]);
                Track_Point[1] = Track_Point[1]+(Final_Length_Scaling*Direction[1]);
                Track_Point[2] = Track_Point[2]+(Final_Length_Scaling*Direction[2]);
                
                Temp_Length = Temp_Length + Final_Length_Scaling;
                //======Print to check=====//
                cout << "CASE1" << endl;
                cout << "Final_Length_Scaling: " << Final_Length_Scaling << endl;
                cout << "Temp_Length: " << Temp_Length << endl;
                Step = Step + 1;
            }

            else if(Critera_of_surpassing<=LOF_km)
            {
                //In the material
                Track_Point[0] = Track_Point[0]+(Sprint_Length*Direction[0]);
                Track_Point[1] = Track_Point[1]+(Sprint_Length*Direction[1]);
                Track_Point[2] = Track_Point[2]+(Sprint_Length*Direction[2]);
                
                Temp_Length = Temp_Length + Sprint_Length;

                cout << "CASE2" << endl;
                cout << "Sprint: " << Sprint_Length << endl;
                cout << "Temp_Length: " << Temp_Length << endl;

                double *V_Aft_Collision_Earth = Velocity_Aft_collision(1,WIMP_Mass,Sigma_SI_Default,DM_Velocity_Aft_Colliding,1);
                double Energy_Loss_to_Atom   = Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) - Energy_DM(WIMP_Mass,V_Aft_Collision_Earth[0]*1e3/3e8);
                double Ratio_of_Energy_Loss_to_Atom = Energy_Loss_to_Atom/Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8);
                cout << "Ratio_of_Energy_Loss_to_Atom: " << Ratio_of_Energy_Loss_to_Atom << endl;
                DM_Velocity_Aft_Colliding=V_Aft_Collision_Earth[0];
                cout << "Aft_Vel: " << DM_Velocity_Aft_Colliding << endl;
                
                if(Straight_or_scattered==1){
                double *Direction_aft = Aft_scatterd_Direction(1,WIMP_Mass,DM_Velocity_Aft_Colliding,Direction,Ratio_of_Energy_Loss_to_Atom);
                Direction[0] = Direction_aft[0];Direction[1] = Direction_aft[1];Direction[2] = Direction_aft[2];}
                
                if(Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8)<0.01)
                {
                    cout << "Energy Threshold" << endl;
                    break;
                }
                Step = Step + 1;
                Collision_Time = Collision_Time + 1;
            }
            if( (Track_Point[2]-0.) <1e-10 and Step!=0)
            {
                cout << "No!Back_Scattering: " << Step << endl;
                cout << "Temp_Length: " << Temp_Length << endl;
                Arrival_or_not = 0;
                break;
            }
            if( (Track_Point[2]-10.) <1e-10)
            {
                cout << "Step: " << Step << endl;
                cout << "Temp_Length: " << Temp_Length << endl;
                Arrival_or_not = 1;
                break;
            }

        }

        cout << "Final_Track_Point[0]: " << Track_Point[0] << endl;
        cout << "Final_Track_Point[1]: " << Track_Point[1] << endl;
        cout << "Track_Track_Point[2]: " << Track_Point[2] << endl;
        
        
        cout << "=========Temp_Length-10=========: " << Temp_Length-10 << endl;


        cout << "Step: " << Step << endl;
        cout << "Collision_Time: " << Collision_Time << endl;

        
    }
    if(Temp_Length!=0)
    {
    gStyle->SetPalette(1);
    dt->SetMarkerStyle(20);
    dt->Draw("pcolLINE");
    }
    //=======================
    //cout << "Ene_Vel_Air_Final: " << Energy_DM(WIMP_Mass,DM_Velocity_Aft_Colliding*1e3/3e8) << endl;
    //cout << "Velocity_Air_Final: " << DM_Velocity_Aft_Colliding << endl;
    //cout << "Collision_TIME_Check: " <<Collision_TIME_Check << endl;
    //cout << " Actual_collision: " << Actual_collision << endl;
    if(Check_Threshold==0)
    {
        if(Path_Original!=0){RETURN_VALUE[1]=(Temp_Length/Path_Original);}
        else{RETURN_VALUE[1]=0;}
    RETURN_VALUE[0]=DM_Velocity_Aft_Colliding;RETURN_VALUE[2]=Collision_Time;
    RETURN_VALUE[3]=Arrival_or_not;
    RETURN_VALUE[4]=10;RETURN_VALUE[5]=Temp_Length;RETURN_VALUE[6]=Collision_Time;
    }
    
    if(Check_Threshold!=0)
    {
   RETURN_VALUE[0]=0;RETURN_VALUE[1]=0;RETURN_VALUE[2]=0;
    }

    return RETURN_VALUE;
}



