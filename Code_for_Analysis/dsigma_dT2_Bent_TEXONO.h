//All values are expressed in cm, need to be changed to km when using. Thanks 1e-2*1e-3=1e-5
//Thickness(TNK)
double Cu_TKN=5;//OFHC Cu
double B_TKN=25;//Boron
double Fe_TKN=5;//Steel
double Pb_TKN=15;//Lead

double Ge_Layer[3] ={80,100,75};//Ge_Detector_Outmost_Layer
double Cu_Layer[3]={Ge_Layer[0]+Cu_TKN*2,Ge_Layer[1]+Cu_TKN*2,Ge_Layer[2]+Cu_TKN*2};//OFHC Cu
double B_Layer[3] ={Cu_Layer[0]+B_TKN*2,Cu_Layer[1]+B_TKN*2,Cu_Layer[2]+B_TKN*2};//Boron
double Fe_Layer[3]={B_Layer[0]+Fe_TKN*2,B_Layer[1]+Fe_TKN*2,B_Layer[2]+Fe_TKN*2};//Steel
double Pb_Layer[3]={Fe_Layer[0]+Pb_TKN*2,Fe_Layer[1]+Pb_TKN*2,Fe_Layer[2]+Pb_TKN*2};//Steel

double *Starting_Point()
{
    static double SPA[3];//Starting Point Array
    TF1 *X_Detector = new TF1("X_Detector","X_Detector",-Ge_Layer[0]*0.5,Ge_Layer[0]*0.5);
    TF1 *Y_Detector = new TF1("Y_Detector","Y_Detector",-Ge_Layer[1]*0.5,Ge_Layer[1]*0.5);
    TF1 *Z_Detector = new TF1("Z_Detector","Z_Detector",-Ge_Layer[2]*0.5,Ge_Layer[2]*0.5);

    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);
    //Choose the plane which the WIMP runs from
    TF1 *X_Y_Z_Choose = new TF1("X_Y_Z_Choose","X_Y_Z_Choose",0,3);
    int Plane_Chosen= X_Y_Z_Choose->GetRandom();

    TF1 *Two_Sides    = new TF1("Two_Sides","Two_Sides",0,2);
    int Side_confirmed = Two_Sides->GetRandom();

    if(Random_Value==0)//On X-Plane
    {
        if(Side_confirmed==0){SPA[0]=-Ge_Layer[0]*0.5;}
        if(Side_confirmed==1){SPA[0]= Ge_Layer[0]*0.5;}
        SPA[1]=Y_Detector->GetRandom();SPA[2]=Z_Detector->GetRandom();
        return(SPA);
    }
    if(Random_Value==1)//On Y-Plane
    {
        if(Side_confirmed==0){SPA[1]=-Ge_Layer[1]*0.5;}
        if(Side_confirmed==1){SPA[1]= Ge_Layer[1]*0.5;}
        SPA[0]=X_Detector->GetRandom();SPA[2]=Z_Detector->GetRandom();
        return(SPA);
    }
    if(Random_Value==2)//On Z-Plane
    {
        if(Side_confirmed==0){SPA[2]=-Ge_Layer[2]*0.5;}
        if(Side_confirmed==1){SPA[2]= Ge_Layer[2]*0.5;}
        SPA[0]=X_Detector->GetRandom();SPA[1]=Y_Detector->GetRandom();
        return(SPA);
    }
}

