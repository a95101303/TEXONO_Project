//All values are expressed in m, need to be changed to km when using. Thanks 1e-2*1e-3=1e-5
//Thickness(TNK)
double Cu_TKN=0.05;//OFHC Cu
double B_TKN=0.25;//Boron
double Fe_TKN=0.05;//Steel
double Pb_TKN=0.15;//Lead

double Ge_Layer[3] ={0.8,1,0.75};//Layer0, Ge_Detector_Outmost_Layer
double Cu_Layer[3]={Ge_Layer[0]+Cu_TKN*2,Ge_Layer[1]+Cu_TKN*2,Ge_Layer[2]+Cu_TKN*2};//Layer1, OFHC Cu
double B_Layer[3] ={Cu_Layer[0]+B_TKN*2,Cu_Layer[1]+B_TKN*2,Cu_Layer[2]+B_TKN*2};//Layer2, Boron
double Fe_Layer[3]={B_Layer[0]+Fe_TKN*2,B_Layer[1]+Fe_TKN*2,B_Layer[2]+Fe_TKN*2};//Layer3,Steel
double Pb_Layer[3]={Fe_Layer[0]+Pb_TKN*2,Fe_Layer[1]+Pb_TKN*2,Fe_Layer[2]+Pb_TKN*2};//Layer4, Steel

double H_OCT_KS    =70.;//m, High_Outer_Ceiling_Top_KS: 70m
double H_ICT_KS    =60.;//m, High_inner_Ceiling_Top_KS: 60m
double R_OW_KS     =57.;//m, Radius_Outer_Wall_KS: 57m
double R_IW_KS     =56.;//m, Radius_Inner_Wall_KS: 56m
double R_OWALL_RE  =28.;//m, Radius_Outer_Wall_Reactor: 28m
double R_IWATER_RE =27.;//m, Radius_Inner_Water_Reactor: 27m

//Starting point from the surface of the detector
double *Starting_Point()
{
    static double Return_Value[5];//Starting Point Array
    double SPA[3];//Starting_Points
    double VR[2];//Vector_Restriction[0]:Coordinate Vector_Restriction[1]:Plus or Minus
    
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
        if(Side_confirmed==0){SPA[0]=-Ge_Layer[0]*0.5;VR[0]=1;VR[1]=-1;}
        if(Side_confirmed==1){SPA[0]= Ge_Layer[0]*0.5;VR[0]=1;VR[1]= 1;}
        SPA[1]=Y_Detector->GetRandom();SPA[2]=Z_Detector->GetRandom();
    }
    if(Random_Value==1)//On Y-Plane
    {
        if(Side_confirmed==0){SPA[1]=-Ge_Layer[1]*0.5;VR[0]=2;VR[1]=-1;}
        if(Side_confirmed==1){SPA[1]= Ge_Layer[1]*0.5;VR[0]=2;VR[1]= 1;}
        SPA[0]=X_Detector->GetRandom();SPA[2]=Z_Detector->GetRandom();
    }
    if(Random_Value==2)//On Z-Plane
    {
        if(Side_confirmed==0){SPA[2]=-Ge_Layer[2]*0.5;VR[0]=3;VR[1]=-1;}
        if(Side_confirmed==1){SPA[2]= Ge_Layer[2]*0.5;VR[0]=3;VR[1]= 1;}
        SPA[0]=X_Detector->GetRandom();SPA[1]=Y_Detector->GetRandom();
    }
    Return_Value[0]=SPA[0];Return_Value[1]=SPA[1];Return_Value[2]=SPA[2];
    Return_Value[3]= VR[0];Return_Value[4]= VR[1];
    return(Return_Value);
}

//Now where the particle is in the shielding
double *Where_is_it(double *Point)
{
    static double Return_Value[2];
    int layer_Number=0;int On_Surface=0;
    
    //First Layer
    if( (abs(Point[0])<=Ge_Layer[0]*0.5) and (abs(Point[1])<=Ge_Layer[1]*0.5) and (abs(Point[1])<=Ge_Layer[2]*0.5) )
    {
        layer_Number=0;
    }
    //Second Layer
    if( (abs(Point[0])<=Cu_Layer[0]*0.5) and (abs(Point[1])<=Cu_Layer[1]*0.5) and (abs(Point[1])<=Cu_Layer[2]*0.5) )
    {
        if( (abs(Point[0])>Ge_Layer[0]*0.5) or (abs(Point[1])>Ge_Layer[1]*0.5) or (abs(Point[1])>Ge_Layer[2]*0.5) )
        {
            layer_Number=1;
        }
    }
        
    //Third Layer
    if( (abs(Point[0])<=B_Layer[0]*0.5) and (abs(Point[1])<=B_Layer[1]*0.5) and (abs(Point[1])<=B_Layer[2]*0.5) )
    {
        if( (abs(Point[0])>Cu_Layer[0]*0.5) or (abs(Point[1])>Cu_Layer[1]*0.5) or (abs(Point[1])>Cu_Layer[2]*0.5) )
        {
            layer_Number=2;
        }
    }

    //Fourth Layer
    if( (abs(Point[0])<=Fe_Layer[0]*0.5) and (abs(Point[1])<=Fe_Layer[1]*0.5) and (abs(Point[1])<=Fe_Layer[2]*0.5) )
    {
        if( (abs(Point[0])>B_Layer[0]*0.5) or (abs(Point[1])>B_Layer[1]*0.5) or (abs(Point[1])>B_Layer[2]*0.5) )
        {
            layer_Number=3;
        }
    }

    //Fourth Layer
    if( (abs(Point[0])<=Pb_Layer[0]*0.5) and (abs(Point[1])<=Pb_Layer[1]*0.5) and (abs(Point[1])<=Pb_Layer[2]*0.5) )
    {
        if( (abs(Point[0])>Fe_Layer[0]*0.5) or (abs(Point[1])>Fe_Layer[1]*0.5) or (abs(Point[1])>Fe_Layer[2]*0.5) )
        {
            layer_Number=4;
        }
    }

    //Out of the detector
    else{layer_Number=5;}
    
    return(layer_Number);
}

int In_or_Out(double *Point, double *Direction)
{
    double Criteria = Point[0]*Direction[0]+Point[1]*Direction[1]+Point[2]*Direction[2];
    int In_or_Out_N=0;
    if(Criteria>=0)In_or_Out_N=1;
    else{In_or_Out_N=1;}
    
    return(In_or_Out_N);
}

double Criteria_in_Shielding(int Layer, double *Direction)
{
    
}
