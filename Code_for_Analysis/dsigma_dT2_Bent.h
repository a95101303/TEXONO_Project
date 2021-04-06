//TMath
const double AT_mx = 931.494013*1e3; // keV/c^2   AT_mx: 931.494013MeV/c^2
const double Degree_to_Radian = 0.0174532925; //Degree*Degree_to_Radian=Radian

double Mass_to_Nucleon(double DM_mx, int Earth_or_Air)
{
    if(Earth_or_Air==1)
    {return (DM_mx)/(35.);}
    if(Earth_or_Air==2)
    {return (DM_mx)/(15.);}
}


double Square_Sum(double *Vector1,double *Vector2)//Momentum
{return Vector1[0]*Vector2[0]+Vector1[1]*Vector2[1]+Vector1[2]*Vector2[2];}

double Distant(double *The_first_point, double *The_second_point)
{
    setprecision(9);
    double X_Difference=(The_first_point[0]-The_second_point[0]);
    double Y_Difference=(The_first_point[1]-The_second_point[1]);
    double Z_Difference=(The_first_point[2]-The_second_point[2]);
    double Total_Difference = sqrt(TMath::Power(X_Difference,2)+TMath::Power(Y_Difference,2)+TMath::Power(Z_Difference,2));
    //cout << "Total_Difference" << Total_Difference << endl;
    return Total_Difference;
}
double Distant_Criteria(double *The_first_point, double *The_second_point)
{
    setprecision(9);
    double X_Difference=(The_first_point[0]-The_second_point[0]);
    double Y_Difference=(The_first_point[1]-The_second_point[1]);
    double Z_Difference=(The_first_point[2]-The_second_point[2]);
    double Total_Difference = sqrt(TMath::Power(X_Difference,2)+TMath::Power(Y_Difference,2)+TMath::Power(Z_Difference,2));
    //cout << "Total_Difference" << Total_Difference << endl;
    return Total_Difference;
}

double Phi_Rest_Frame(int Earth_or_Air, double DM_mx, double Theta)
{
    double tanPhi = TMath::Sin(Theta*Degree_to_Radian)/(TMath::Cos(Theta*Degree_to_Radian)* Mass_to_Nucleon(DM_mx,Earth_or_Air));
    double Phi = TMath::ATan(tanPhi)/Degree_to_Radian;
    return Phi;
}

//Lab-frame
//DM_mx:GeV/c^2
//IE_DM,FE_Nu,FE_DM:KeV
//Scattering Angle
double *Scattering_Angle(double DM_mx, double IE_DM, double FE_Nu, double FE_DM)//IM: Initial Energy, FM: Final Energy
{
    cout << "Initial_Energy: " << IE_DM << endl;
    cout << "Collision_Loss: " << FE_Nu << endl;
    cout << "Final_Energy: "   << FE_DM << endl;
    
    static double RETURN_VALUE[2];

    double Product_of_IMFM_DM = ((2*AT_mx*FE_Nu)-(2*DM_mx*1e6*IE_DM)-(2*DM_mx*1e6*FE_DM)) / (-2);
    cout << "Product_of_IMFM_DM: " << Product_of_IMFM_DM  << endl;
    double Angle_Value = Product_of_IMFM_DM/  ( sqrt(2*DM_mx*1e6*IE_DM)*sqrt(2*DM_mx*1e6*FE_DM) );
    cout << "sqrt(2*DM_mx*1e6*IE_DM)*sqrt(2*DM_mx*1e6*FE_DM): " << sqrt(2*DM_mx*1e6*IE_DM)*sqrt(2*DM_mx*1e6*FE_DM) << endl;
    double Degree_of_scattering = TMath::ACos(Angle_Value)/Degree_to_Radian;
    cout << "Angle_Value: " << Angle_Value << endl;
    RETURN_VALUE[0]=Product_of_IMFM_DM;RETURN_VALUE[1]=Degree_of_scattering;
    cout << "Degree_of_scattering: " << Degree_of_scattering << endl;
    cout << "Product_of_IMFM_DM: " << Product_of_IMFM_DM << endl;
    
    return RETURN_VALUE;
}

//=======First choose Pz
double PZ_PY_PX_Vector(int Option, double DM_mx, double Product_of_IMFM_DM, double *IM_DM_V, double FE_DM, double Determined_Z, double Determined_Y)
{
    //=================================Option=0:Determine Pz=================================
    double a=(Product_of_IMFM_DM/IM_DM_V[0]);
    double b=(IM_DM_V[1]/IM_DM_V[0]);
    double c=(IM_DM_V[2]/IM_DM_V[0]);
    double e=2*sqrt(b*b+1);
    
    double A=((2*b*c/e)*(2*b*c/e))-(c*c+1);
    double B=(2*a*c)-(2*(2*b*c/e)*(2*b*a/e));
    double C=((2*b*a/e)*(2*b*a/e))+(2*DM_mx*1e6*FE_DM)-(a*a);
    
    double DETERMINE1= sqrt((B*B-4*A*C)/(4*A*A))-(B/(2*A));
    double DETERMINE2=-sqrt((B*B-4*A*C)/(4*A*A))-(B/(2*A));
    double CRITERIA  = sqrt(2*DM_mx*1e6*FE_DM);
    
    //==========================
    double Positive_Z=0;double Negative_Z=0;double Final_Z=0;

    if(CRITERIA>DETERMINE1)
    {
         TF1 *f1 = new TF1("f1","1",DETERMINE1,CRITERIA);
         Positive_Z = f1->GetRandom();
    }
    if(-CRITERIA<DETERMINE2)
    {
         TF1 *f2 = new TF1("f2","1",DETERMINE2,-CRITERIA);
         Negative_Z = f2->GetRandom();
    }
    //==========================
    if(abs(Positive_Z)<CRITERIA and abs(Negative_Z)<CRITERIA)
    {
        TF1 *f3 = new TF1("f3","1",-1,1);
        double Random=f3->GetRandom();
        if(Random>=0){Final_Z=Positive_Z;}
        if(Random< 0){Final_Z=Negative_Z;}
    }
    else if(abs(Positive_Z)<CRITERIA and Negative_Z==0){Final_Z=Positive_Z;}
    else if(Positive_Z==0 and abs(Negative_Z)<CRITERIA){Final_Z=Negative_Z;}
    else if(Positive_Z==0 and Negative_Z==0){cout << "Weird!" << endl;}

    
    //=================================Option=1:Determine Py=================================
    cout << "Determined_Z: " << Determined_Z << endl;
    double Positive_Term1 =  sqrt(A*(Determined_Z)*(Determined_Z)+B*(Determined_Z)+C);
    double Negative_Term1 = -sqrt(A*(Determined_Z)*(Determined_Z)+B*(Determined_Z)+C);
    double Term2 = (2*b*c/e)*(Determined_Z)-(2*a*b/e);
    
    cout << "Positive_Term1: " << Positive_Term1 << endl;
    cout << "Negative_Term1: " << Negative_Term1 << endl;
    cout << "Term2: "          << Term2 << endl;

    double Positive_Y     = (2/e)*( Positive_Term1- Term2 );
    double Negative_Y     = (2/e)*( Negative_Term1- Term2 );
    double Final_Y = 0;
    
    cout << "Positive_Y: " << Positive_Y << endl;
    cout << "Negative_Y: " << Negative_Y << endl;

    if(Positive_Y>CRITERIA and Negative_Y>-CRITERIA)
    {
        TF1 *f4 = new TF1("f4","1",-1,1);
        double Random=f4->GetRandom();
        if(Random>=0){Final_Y=Positive_Y;}
        if(Random< 0){Final_Y=Negative_Y;}
    }
    else if(abs(Positive_Y)<=CRITERIA and abs(Negative_Y)>CRITERIA){Final_Y=Positive_Y;}
    else if(abs(Positive_Y)>CRITERIA and abs(Negative_Y)<=CRITERIA){Final_Y=Negative_Y;}
    else if(abs(Positive_Y)>=CRITERIA and abs(Negative_Y)>=CRITERIA){cout << "Weird!" << endl;}

    //=================================Option=1:Determine Px=================================
    double Final_X = a - (b*Final_Y) - (c*Final_Z);
    
    if(Option==0)return Final_Z;
    if(Option==1)return Final_Y;
    if(Option==2)return Final_X;
}
