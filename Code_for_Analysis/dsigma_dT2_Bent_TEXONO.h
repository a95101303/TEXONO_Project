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

double The_Material_Layer[5][3]=
{
    Ge_Layer[0],Ge_Layer[1],Ge_Layer[2],
    Cu_Layer[0],Cu_Layer[1],Cu_Layer[2],
    B_Layer[0],B_Layer[1],B_Layer[2],
    Fe_Layer[0],Fe_Layer[1],Fe_Layer[2],
    Pb_Layer[0],Pb_Layer[1],Pb_Layer[2],
}

double Vector_for_axis[6][3]=
{
    1,0,0,
   -1,0,0,
    0,1,0,
    0,-1,0,
    0,0,1,
    0,0,-1,
}

double H_OCT_KS    =70.;//m, High_Outer_Ceiling_Top_KS: 70m
double H_ICT_KS    =60.;//m, High_inner_Ceiling_Top_KS: 60m
double R_OW_KS     =57.;//m, Radius_Outer_Wall_KS: 57m
double R_IW_KS     =56.;//m, Radius_Inner_Wall_KS: 56m
double R_OWALL_RE  =28.;//m, Radius_Outer_Wall_Reactor: 28m
double R_IWATER_RE =27.;//m, Radius_Inner_Water_Reactor: 27m

//Starting Position from the surface of the detector
double *Starting_Position()
{
    static double Return_Value[5];//Starting Position Array
    double SPA[3];//Starting_Positions
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
int *Where_is_it(double *Position)
{
    static double Return_Value[2];
    int layer_Number=0;int On_Surface=0;
    
    //First Layer
    if( (abs(Position[0])<=Ge_Layer[0]*0.5) and (abs(Position[1])<=Ge_Layer[1]*0.5) and (abs(Position[1])<=Ge_Layer[2]*0.5) )
    {
        layer_Number=0;
        if(  Position[0] ==Ge_Layer[0]*0.5) )On_Surface=1;  if(  Position[0] ==-Ge_Layer[0]*0.5) )On_Surface=2;
        if(  Position[1] ==Ge_Layer[1]*0.5) )On_Surface=3;  if(  Position[1] ==-Ge_Layer[1]*0.5) )On_Surface=4;
        if(  Position[2] ==Ge_Layer[2]*0.5) )On_Surface=5;  if(  Position[2] ==-Ge_Layer[2]*0.5) )On_Surface=6;

    }
    //Second Layer
    if( (abs(Position[0])<=Cu_Layer[0]*0.5) and (abs(Position[1])<=Cu_Layer[1]*0.5) and (abs(Position[1])<=Cu_Layer[2]*0.5) )
    {
        if( (abs(Position[0])>Ge_Layer[0]*0.5) or (abs(Position[1])>Ge_Layer[1]*0.5) or (abs(Position[1])>Ge_Layer[2]*0.5) )
        {
            layer_Number=1;
            if(  Position[0] ==Cu_Layer[0]*0.5) )On_Surface=1;  if(  Position[0] ==-Cu_Layer[0]*0.5) )On_Surface=2;
            if(  Position[1] ==Cu_Layer[1]*0.5) )On_Surface=3;  if(  Position[1] ==-Cu_Layer[1]*0.5) )On_Surface=4;
            if(  Position[2] ==Cu_Layer[2]*0.5) )On_Surface=5;  if(  Position[2] ==-Cu_Layer[2]*0.5) )On_Surface=6;
        }
    }
        
    //Third Layer
    if( (abs(Position[0])<=B_Layer[0]*0.5) and (abs(Position[1])<=B_Layer[1]*0.5) and (abs(Position[1])<=B_Layer[2]*0.5) )
    {
        if( (abs(Position[0])>Cu_Layer[0]*0.5) or (abs(Position[1])>Cu_Layer[1]*0.5) or (abs(Position[1])>Cu_Layer[2]*0.5) )
        {
            layer_Number=2;
            if(  Position[0] ==B_Layer[0]*0.5) )On_Surface=1;  if(  Position[0] ==-B_Layer[0]*0.5) )On_Surface=2;
            if(  Position[1] ==B_Layer[1]*0.5) )On_Surface=3;  if(  Position[1] ==-B_Layer[1]*0.5) )On_Surface=4;
            if(  Position[2] ==B_Layer[2]*0.5) )On_Surface=5;  if(  Position[2] ==-B_Layer[2]*0.5) )On_Surface=6;
        }
    }

    //Fourth Layer
    if( (abs(Position[0])<=Fe_Layer[0]*0.5) and (abs(Position[1])<=Fe_Layer[1]*0.5) and (abs(Position[1])<=Fe_Layer[2]*0.5) )
    {
        if( (abs(Position[0])>B_Layer[0]*0.5) or (abs(Position[1])>B_Layer[1]*0.5) or (abs(Position[1])>B_Layer[2]*0.5) )
        {
            layer_Number=3;
            if(  Position[0] ==Fe_Layer[0]*0.5) )On_Surface=1;  if(  Position[0] ==-Fe_Layer[0]*0.5) )On_Surface=2;
            if(  Position[1] ==Fe_Layer[1]*0.5) )On_Surface=3;  if(  Position[1] ==-Fe_Layer[1]*0.5) )On_Surface=4;
            if(  Position[2] ==Fe_Layer[2]*0.5) )On_Surface=5;  if(  Position[2] ==-Fe_Layer[2]*0.5) )On_Surface=6;
        }
    }

    //Fourth Layer
    if( (abs(Position[0])<=Pb_Layer[0]*0.5) and (abs(Position[1])<=Pb_Layer[1]*0.5) and (abs(Position[1])<=Pb_Layer[2]*0.5) )
    {
        if( (abs(Position[0])>Fe_Layer[0]*0.5) or (abs(Position[1])>Fe_Layer[1]*0.5) or (abs(Position[1])>Fe_Layer[2]*0.5) )
        {
            layer_Number=4;
            if(  Position[0] ==Pb_Layer[0]*0.5) )On_Surface=1;  if(  Position[0] ==-Pb_Layer[0]*0.5) )On_Surface=2;
            if(  Position[1] ==Pb_Layer[1]*0.5) )On_Surface=3;  if(  Position[1] ==-Pb_Layer[1]*0.5) )On_Surface=4;
            if(  Position[2] ==Pb_Layer[2]*0.5) )On_Surface=5;  if(  Position[2] ==-Pb_Layer[2]*0.5) )On_Surface=6;
        }
    }

    if( (abs(Position[0])>Pb_Layer[0]*0.5) or (abs(Position[1])>Pb_Layer[1]*0.5) or (abs(Position[1])>Pb_Layer[2]*0.5) )
    {
        layer_Number=5;
    }
    //Out of the detector
    
    Return_Value[0]=layer_Number;Return_Value[1]=On_Surface;
    return(Return_Value);
}

int DR_In_or_Out_on_surface(int IOS, double *DR)//Index_of_surface(IOS),Direction(DR)
{
    double Criteria = DR[0]*Vector_for_axis[IOS-1][0]+DR[1]*Vector_for_axis[IOS-1][1]+DR[2]*Vector_for_axis[IOS-1][2];
    int In_or_Out_N=0;
    if(Criteria>=0)In_or_Out_N=1;
    else{In_or_Out_N=0;}
    
    return(In_or_Out_N);
}

int DR_In_or_Out_off_surface(double *POS, double *DR)
{
    double Criteria = POS[0]*DR[0]+POS[1]*DR[1]+POS[2]*DR[2];
    int In_or_Out_N=0;
    if(Criteria>=0)In_or_Out_N=1;
    else{In_or_Out_N=0;}
    
    return(In_or_Out_N);
}

double Length_to_six_planes(int Layer, double *POS, double *DR)//Position(POS),
{
    double POS_Aft[3];
    double ET_XYZ[3][2];//Extended_Times_XYZ
    vector<double> Find_smallest_Vector;
    // Discover the length to six places
    for(int XYZ=0; XYZ<3; XYZ++)//The Axis
    {
        for(int PorN=0; PorN<2; PorN++)//Positive Or Negative
        {
            double ET = (The_Material_Layer[Layer][XYZ]-POS[XYZ])/(DR[XYZ]);//Extension_times
            for(int POS_After_Axis=0; POS_After_Axis<3; POS_After_Axis++){POS_Aft[POS_After_Axis] = POS[POS_After_Axis] + ET*DR[POS_After_Axis];}
            double *Fun_Layer_Aft= Where_is_it(POS_Aft);
            if(Fun_Layer_Aft[0]==Layer){ET_XYZ[XYZ][PorN]=ET;}else{ET_XYZ[XYZ][PorN]=0;}
            if(Fun_Layer_Aft[0]==Layer and ET_XYZ[XYZ][PorN]>0){Find_smallest_Vector.push_back(ET_XYZ[XYZ][PorN]);}
        }
    }
    
    std::sort(Find_smallest_Vector.begin(), Find_smallest_Vector.end());
    //The smallest one
    double The_final_scaling = Find_smallest_Vector[0];
    Find_smallest_Vector.clear();
    
    return The_final_scaling;

}
double Criteria_in_Shielding(double *POS, double *DR)//Position(POS),Direction(DR)
{
    double RSF;//Return_Scaling_Factor

    int *A_Position       = Where_is_it(POS);
    int  A_Layer          = A_Position[0];
    int  A_Surface_index  = A_Position[1];
    int  DR_IN_OR_OUT = 0;
    
    //In the shielding
    if(A_Layer<5)
    {
        if(A_Surface_index!=0)//On the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_on_surface(A_Surface_index,DR);
            if(DR_IN_OR_OUT==1)RSF=Length_to_six_planes(A_Layer+1,POS,DR);
            if(DR_IN_OR_OUT==0)RSF=Length_to_six_planes(A_Layer,POS,DR);
        }
        if(A_Surface_index==0)//Off the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_off_surface(POS,DR);
            if(DR_IN_OR_OUT==1)RSF = Length_to_six_planes(A_Layer,POS,DR);
            if(DR_IN_OR_OUT==0)
            {
               RSF = std::min(Length_to_six_planes(A_Layer,POS,DR),Length_to_six_planes(A_Layer-1,POS,DR))
            }
        }
    }//if(A_Layer<5)
    
    //Out of the shielding
    if(A_Layer==5)
    {
        
    }//if(A_Layer==5)
}





