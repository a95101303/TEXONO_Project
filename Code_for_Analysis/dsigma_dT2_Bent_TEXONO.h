//All values are expressed in m, need to be changed to km when using. Thanks 1e-2*1e-3=1e-5
double Density_for_Material[6]={Density_of_Cement,1,11.34,7.86,2.34,8.96};// Cement,H2O,Pb,Fe,B,Cu  g/cm^3
string ATOM_number_for_Material[6]={"Cement","AH2O","APb","AFe","AB","ACu"};
double Index_Atom[6]={Weighted_Atomic_Number_Cement,AH2O,APb,AFe,AB,ACu};

double ANIS[5] ={AGe,ACu,AB,AFe,APb};//Atomic Numbers Inside the Shielding
double DIS[5] ={5.323,8.96,2.34,7.86,11.34};//Densities Inside the Shielding(g/cm^3)
double ANOFS[3]={Weighted_Atomic_Number_Cement,AH2O,15};//Atomic Numbers Out Of the Shielding
double DOFS[3] ={Density_of_Cement,1,kg_perm3_to_g_percm3(1.225)};//Densities Out Of the Shielding(g/cm^3)

double Atomic_All[2][5]=
{
    ANIS[0],ANIS[1],ANIS[2],ANIS[3],ANIS[4],
    ANOFS[0],ANOFS[1],ANOFS[2],0,0,
};
double Density_All[2][5]=
{
    DIS[0],DIS[1],DIS[2],DIS[3],DIS[4],
    DOFS[0],DOFS[1],DOFS[2],0,0,
};
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
};

double Vector_for_axis[6][3]=
{
    1,0,0,
   -1,0,0,
    0,1,0,
    0,-1,0,
    0,0,1,
    0,0,-1,
};

double RP = 0.5;//Reference_Point
double H_OCT_KS    =70.-RP;//m, High_Outer_Ceiling_Top_KS: 70m
double H_ICT_KS    =60.-RP;//m, High_inner_Ceiling_Top_KS: 60m

double R_OW_KS     =57.;//m, Radius_Outer_Wall_KS: 57m
double R_IW_KS     =56.;//m, Radius_Inner_Wall_KS: 56m

double H_RE        =50.-RP;//m, High_Reactor: 50m
double R_OWALL_RE  =28.;//m, Radius_Outer_Wall_Reactor: 28m
double R_IWATER_RE =27.;//m, Radius_Inner_Water_Reactor: 27m

double The_smallest_in_a_vector(vector<double> Vector)
{
    double Return_the_smallest_one=0;
    if(Vector.size()>0)
    {
        std::sort(Vector.begin(), Vector.end());
        for(int kkk=0; kkk<Vector.size(); kkk++)
        {
            cout << "kkk: " << kkk << endl;
            cout << "Vector: " << Vector[kkk] << endl;
        }
        Return_the_smallest_one = Vector[0];
    }
    else{Return_the_smallest_one=0;}
    
    return Return_the_smallest_one;
}

//Starting Position from the surface of the detector
double *Starting_Position()
{
    static double Return_Value[5];//Starting Position Array
    double SPA[3];//Starting_Positions
    double VR[2];//Vector_Restriction[0]:Coordinate Vector_Restriction[1]:Plus or Minus
    
    TF1 *X_Detector = new TF1("X_Detector","1",-Ge_Layer[0]*0.5,Ge_Layer[0]*0.5);
    TF1 *Y_Detector = new TF1("Y_Detector","1",-Ge_Layer[1]*0.5,Ge_Layer[1]*0.5);
    TF1 *Z_Detector = new TF1("Z_Detector","1",-Ge_Layer[2]*0.5,Ge_Layer[2]*0.5);

    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);
    //Choose the plane where the WIMP starts running
    TF1 *X_Y_Z_Choose = new TF1("X_Y_Z_Choose","1",0,3);
    int Plane_Chosen= X_Y_Z_Choose->GetRandom();//Give out 0,1,2;
    cout << "Plane_Chosen: " << Plane_Chosen << endl;
    //Two sides for WIMPs
    TF1 *Two_Sides    = new TF1("Two_Sides","1",0,2);
    int Side_confirmed = Two_Sides->GetRandom();//Give out 0,1;
    double Confirmed_Side[2]={1,-1};
    //Fix the point on a certain plane
    cout << "Confirmed_Side[Side_confirmed]: " << Confirmed_Side[Side_confirmed] << endl;
    SPA[Plane_Chosen]= Confirmed_Side[Side_confirmed]*Ge_Layer[Plane_Chosen]*0.5;VR[0]=Plane_Chosen;VR[1]=Confirmed_Side[Side_confirmed];
    //Decide the rest of the coordination
    double X_Ran=X_Detector->GetRandom();
    double Y_Ran=Y_Detector->GetRandom();
    double Z_Ran=Z_Detector->GetRandom();

    if(abs(SPA[0]-0)<1e-5){SPA[0]=X_Ran;}
    if(abs(SPA[1]-0)<1e-5){SPA[1]=Y_Ran;}
    if(abs(SPA[2]-0)<1e-5){SPA[2]=Z_Ran;}

    Return_Value[0]=SPA[0];Return_Value[1]=SPA[1];Return_Value[2]=SPA[2];
    Return_Value[3]= VR[0];Return_Value[4]= VR[1];
    cout << "SPA[0]: " << SPA[0] << endl;
    cout << "SPA[1]: " << SPA[1] << endl;
    cout << "SPA[2]: " << SPA[2] << endl;

    return(Return_Value);
}

double Check_DR(int *ROP, double *DR)//Restriction_on_Direction
{
    double Sign=0;
    for(int kkk=0; kkk<3; kkk++)
    {
        if(ROP[0]==kkk)
        {
            if(DR[kkk]*ROP[1]>0)Sign= 1;
            if(DR[kkk]*ROP[1]<0)Sign=-1;
        }
    }
    return Sign;
}
//Now where the particle is in the shielding
int *Where_is_it_in_shielding(double *Position)
{
    static int Return_Value[3];
    int layer_Number=0;int On_Surface=0;
    int Yes_Or_No=0;
    //First Layer
    if( abs(Position[0])<=Ge_Layer[0]*0.5 and abs(Position[1])<=Ge_Layer[1]*0.5 and abs(Position[2])<=Ge_Layer[2]*0.5 )
    {
        layer_Number=0;
        if(  Position[0] ==Ge_Layer[0]*0.5  )On_Surface=1;  if(  Position[0] ==-Ge_Layer[0]*0.5  )On_Surface=2;
        if(  Position[1] ==Ge_Layer[1]*0.5  )On_Surface=3;  if(  Position[1] ==-Ge_Layer[1]*0.5  )On_Surface=4;
        if(  Position[2] ==Ge_Layer[2]*0.5  )On_Surface=5;  if(  Position[2] ==-Ge_Layer[2]*0.5  )On_Surface=6;

    }
    //Second to fourth layer
    for(int Layer=1; Layer<5; Layer++)
    {
        if( (abs(Position[0])<=The_Material_Layer[Layer][0]*0.5) and (abs(Position[1])<=The_Material_Layer[Layer][1]*0.5) and (abs(Position[2])<=The_Material_Layer[Layer][2]*0.5) )
        {
            if( (abs(Position[0])>The_Material_Layer[Layer-1][0]*0.5) or (abs(Position[1])>The_Material_Layer[Layer-1][1]*0.5) or (abs(Position[2])>The_Material_Layer[Layer-1][2]*0.5) )
            {
                layer_Number=Layer;
                if(  Position[0] ==The_Material_Layer[Layer][0]*0.5  )On_Surface=1;  if(  Position[0] ==-The_Material_Layer[Layer][0]*0.5  )On_Surface=2;
                if(  Position[1] ==The_Material_Layer[Layer][1]*0.5  )On_Surface=3;  if(  Position[1] ==-The_Material_Layer[Layer][1]*0.5  )On_Surface=4;
                if(  Position[2] ==The_Material_Layer[Layer][2]*0.5  )On_Surface=5;  if(  Position[2] ==-The_Material_Layer[Layer][2]*0.5  )On_Surface=6;
            }
        }
    }
    //Out of the fourth layer
    if(  abs(Position[0])>Pb_Layer[0]*0.5 or  abs(Position[1])>Pb_Layer[1]*0.5  or  abs(Position[2])>Pb_Layer[2]*0.5  )
    {
        layer_Number=5;
    }
    //Out of the detector
    
    if(layer_Number<4 or (layer_Number==4 and On_Surface==0) ) Yes_Or_No=1;
    else if( (layer_Number==4 and On_Surface>0) ) Yes_Or_No=2;
    else{Yes_Or_No=0;}
    Return_Value[0]=layer_Number;Return_Value[1]=On_Surface;Return_Value[2]=Yes_Or_No;
    return(Return_Value);
}


double Length_to_six_planes(int Layer, double *POS, double *DR)//Position(POS),
{
    double POS_Aft[3];
    double ET_XYZ[3][2];//Extended_Times_XYZ
    vector<double> Find_smallest_Vector;
    double Confirmed_Side[2]={1,-1};
    double The_final_scaling=0;
    // Discover the length to six places
    for(int XYZ=0; XYZ<3; XYZ++)//The Axis
    {
        for(int PorN=0; PorN<2; PorN++)//Positive Or Negative
        {
            double ET = (Confirmed_Side[PorN]*The_Material_Layer[Layer][XYZ]*0.5-POS[XYZ])/(DR[XYZ]);//Extension_times
            cout << "ET: " << ET << endl;
            for(int POS_After_Axis=0; POS_After_Axis<3; POS_After_Axis++){POS_Aft[POS_After_Axis] = POS[POS_After_Axis] + ET*DR[POS_After_Axis];}
            int *Fun_Layer_Aft= Where_is_it_in_shielding(POS_Aft);
            cout << "Layer: " << Layer << endl;
            cout << "Fun_Layer_Aft[0]: " << Fun_Layer_Aft[0] << endl;
            if(Fun_Layer_Aft[0]==Layer){ET_XYZ[XYZ][PorN]=ET;}else{ET_XYZ[XYZ][PorN]=0;}
            if(Fun_Layer_Aft[0]==Layer and ET>0){Find_smallest_Vector.push_back(ET);cout << "ET1: " << ET << endl;}
        }
    }
    
    The_final_scaling = The_smallest_in_a_vector(Find_smallest_Vector);
    cout << "The_final_scaling: " << The_final_scaling << endl;
    //The smallest one
    Find_smallest_Vector.clear();
    
    return The_final_scaling;

}

int DR_In_or_Out_on_surface(int IOS, double *DR)//Index_of_surface(IOS),Direction(DR)
{//For the shielding
    double Criteria = DR[0]*Vector_for_axis[IOS-1][0]+DR[1]*Vector_for_axis[IOS-1][1]+DR[2]*Vector_for_axis[IOS-1][2];
    int In_or_Out_N=0;
    if(Criteria>=0)In_or_Out_N=1;
    else{In_or_Out_N=0;}
    
    return(In_or_Out_N);
}

int DR_In_or_Out_off_surface(double *POS, double *DR)
{//For the shielding
    double Criteria = POS[0]*DR[0]+POS[1]*DR[1]+POS[2]*DR[2];
    int In_or_Out_N=0;
    if(Criteria>=0)In_or_Out_N=1;
    else{In_or_Out_N=0;}
    
    return(In_or_Out_N);
}

double *CID(double *POS, double *DR)//Criteria_In_Shielding(CID),Position(POS),Direction(DR)
{
    static double Return_Factor[2];//Return_Scaling_Factor, Layer the WIMP runs in
    double RSF;
    int *A_Position       = Where_is_it_in_shielding(POS);
    int  A_Layer          = A_Position[0];
    int  A_Surface_index  = A_Position[1];
    int  DR_IN_OR_OUT = 0; double Layer_run=0;
    vector<double> Return_the_smallest_one;

    //In the shielding
    if(  A_Layer<4 )
    {
        if(A_Surface_index!=0)//On the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_on_surface(A_Surface_index,DR);
            if(DR_IN_OR_OUT==1)
            {
                cout << "L" << A_Layer << "OS OutDR" << endl;
                RSF           = Length_to_six_planes(A_Layer+1,POS,DR);
                Layer_run =  A_Layer+1;
            }
            if( (A_Layer>0 and DR_IN_OR_OUT==0) )
            {
               cout << "L" << A_Layer << "OS InDR" << endl;
               double A1=Length_to_six_planes(A_Layer,POS,DR);
               double A2=Length_to_six_planes(A_Layer-1,POS,DR);
               if(A1>0) Return_the_smallest_one.push_back(A1);if(A2>0) Return_the_smallest_one.push_back(A2);
               RSF = The_smallest_in_a_vector(Return_the_smallest_one);
               Layer_run =  A_Layer;
            }
        }
        if(A_Surface_index==0)//Off the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_off_surface(POS,DR);
            Layer_run =  A_Layer;
            if(DR_IN_OR_OUT==1)
            {
                cout << "L" << A_Layer << "FS OutDR" << endl;
                RSF = Length_to_six_planes(A_Layer,POS,DR);
            }
            if(DR_IN_OR_OUT==0)
            {
                cout << "L" << A_Layer << "FS InDR" << endl;
                double A1=Length_to_six_planes(A_Layer,POS,DR);
                double A2=Length_to_six_planes(A_Layer-1,POS,DR);
                if(A1>0) Return_the_smallest_one.push_back(A1);if(A2>0) Return_the_smallest_one.push_back(A2);
                RSF = The_smallest_in_a_vector(Return_the_smallest_one);
            }
        }
    }//if(A_Layer<5)
    
    if( A_Layer==4 )
    {
        Layer_run =  A_Layer;
        if(A_Surface_index!=0)//On the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_on_surface(A_Surface_index,DR);
            if(DR_IN_OR_OUT==0)
            {
                double A1=Length_to_six_planes(A_Layer,POS,DR);
                double A2=Length_to_six_planes(A_Layer-1,POS,DR);
                if(A1>0) Return_the_smallest_one.push_back(A1);if(A2>0) Return_the_smallest_one.push_back(A2);
                RSF = The_smallest_in_a_vector(Return_the_smallest_one);
            }
        }
        if(A_Surface_index==0)//Off the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_off_surface(POS,DR);
            if(DR_IN_OR_OUT==1)
            {
                RSF = Length_to_six_planes(A_Layer,POS,DR);
            }
            if(DR_IN_OR_OUT==0)
            {
                double A1=Length_to_six_planes(A_Layer,POS,DR);
                double A2=Length_to_six_planes(A_Layer-1,POS,DR);
                if(A1>0) Return_the_smallest_one.push_back(A1);if(A2>0) Return_the_smallest_one.push_back(A2);
                RSF = The_smallest_in_a_vector(Return_the_smallest_one);
            }
        }

    }
    cout << "RSF: " << RSF << endl;
    cout << "Layer_run: " << Layer_run << endl;
    Return_Factor[0]=RSF;Return_Factor[1]=Layer_run;
    return Return_Factor;
}






double SFLE(double a, double b, double c) // Solution_for_Linear_Equation[a(X^2) + b(X) + c = 0]
{//Check OK!
    double Denominator=0; //Down
    double Numerator1=0;double Numerator2=0;double Numerator_Final=0; //Up
    vector<double> Find_smallest_Vector;
    
    Denominator = 2*a;
    Numerator1 = -b + sqrt(b*b-4*a*c);
    if(Numerator1>0)Find_smallest_Vector.push_back(Numerator1);
    Numerator2 = -b - sqrt(b*b-4*a*c);
    if(Numerator2>0)Find_smallest_Vector.push_back(Numerator2);
    
    Numerator_Final = The_smallest_in_a_vector(Find_smallest_Vector);

    return (Numerator_Final/Denominator);
}

double Scaling_to_others_XY(double *POS, double *DR, double R)//Position(POS),Direction(DR),Radius(R)
{
    return SFLE( (DR[0]*DR[0]+DR[1]*DR[1]),(2*POS[0]*DR[0]+2*(POS[1]-30)*DR[1]),(POS[0]*POS[0]+(POS[1]-30)*(POS[1]-30)-R*R) );
}
double Scaling_to_others_Z(double *POS, double *DR, double H)//Position(POS),Direction(DR),Radius(R)
{
    return (H-POS[2])/DR[2];
}
double Radius_for_ceiling(double *POS, double *DR, double S)//Position(POS),Direction(DR),Scaling(S),Radius(R)
{
    double X = POS[0]+DR[0]*S;
    double Y = (POS[1]-30)+DR[1]*S;
    return sqrt(X*X+Y*Y);
}
double RTC(double *POS)//Radius_To_Center(RTC) on XY plane
{
    double X = (POS[0]-0);
    double Y = (POS[1]-30);
    return sqrt(X*X+Y*Y);
}

int *Where_is_it_out_of_shielding(double *POS)
{
    static int Return_Value[2];
    int Component=0;int On_Surface=0;
    
    if( (RTC(POS)>R_IWATER_RE) and (RTC(POS)<=R_OWALL_RE) and (POS[2]<=H_RE) )//Outer Wall of Reactor
    {
        Component=1;
        if(POS[2]==H_RE)On_Surface=2;if(RTC(POS)==R_OWALL_RE)On_Surface=22;
    }
    if( (RTC(POS)<=R_IWATER_RE) and (POS[2]<=H_RE) )//Inner Wall of Reactor
    {
        Component=2;
        if(POS[2]==H_RE)On_Surface=2;if(RTC(POS)==R_IWATER_RE)On_Surface=22;
    }
    if( (RTC(POS)<=R_OW_KS) and (RTC(POS)>=R_IW_KS) and (POS[2]<=H_ICT_KS) )//Side Wall of KS
    {
        Component=3;//Wall of KS
        if(RTC(POS)==R_IW_KS)On_Surface=1;if(RTC(POS)==R_OW_KS)On_Surface=2;
    }
    if( (RTC(POS)<=R_IW_KS) and (POS[2]>=H_ICT_KS) and (POS[2]<=H_OCT_KS) )//Upper ceiling of KS
    {
        Component=4;//Wall of KS
        if(POS[2]==H_ICT_KS)On_Surface=1;if(POS[2]==H_OCT_KS)On_Surface=2;
    }
    if( (RTC(POS)>R_IW_KS) and (RTC(POS)<=R_OW_KS) and (POS[2]>H_ICT_KS) and (POS[2]<=H_OCT_KS))//Corner of the wall
    {
        Component=5;
        if(RTC(POS)==R_OW_KS)On_Surface=2;if(POS[2]==H_OCT_KS)On_Surface=22;
    }
    Return_Value[0]=Component;Return_Value[1]=On_Surface;
    
    return(Return_Value);
}

double DOLOUD(double *Component, double *POS, double *DR)//Decision_On_Length_Out_Of_Detector(DOLOUD)
{
    vector<double> Find_smallest_Vector;
    
    if(Component[0]==1){//Floor
        double Down_floor_Scaling     = Scaling_to_others_Z(POS,DR,-0.5);
        if(Radius_for_ceiling(POS,DR,Down_floor_Scaling)<=R_IW_KS)Find_smallest_Vector.push_back(Down_floor_Scaling);
    }
    if(Component[1]==1){//Shielding
        double Shielding_Scaling      = Length_to_six_planes(4,POS,DR);
        if(Shielding_Scaling>=0)Find_smallest_Vector.push_back(Shielding_Scaling);
    }
    if(Component[2]==1){//Outer_Reactor
        double Outer_RE_Scaling = Scaling_to_others_XY(POS,DR,R_OWALL_RE);
        if(POS[2]+Outer_RE_Scaling*DR[2]<=H_RE)Find_smallest_Vector.push_back(Outer_RE_Scaling);
    }
    if(Component[3]==1){//Inner_Reactor
        double Inner_RE_Scaling = Scaling_to_others_XY(POS,DR,R_IWATER_RE);
        if(POS[2]+Inner_RE_Scaling*DR[2]<=H_RE)Find_smallest_Vector.push_back(Inner_RE_Scaling);
    }
    if(Component[4]==1){//Inner_Wall
        double Inner_Wall_Scaling = Scaling_to_others_XY(POS,DR,R_IW_KS);
        if(POS[2]+Inner_Wall_Scaling*DR[2]<=H_ICT_KS)Find_smallest_Vector.push_back(Inner_Wall_Scaling);
    }
    if(Component[5]==1){//Outer_Wall
        double Outer_Wall_Scaling = Scaling_to_others_XY(POS,DR,R_OW_KS);
        if(POS[2]+Outer_Wall_Scaling*DR[2]<=H_OCT_KS)Find_smallest_Vector.push_back(Outer_Wall_Scaling);
    }
    if(Component[6]==1){//Low_ceiling
        double Low_ceiling_Scaling = Scaling_to_others_Z(POS,DR,H_ICT_KS);
        if(Radius_for_ceiling(POS,DR,Low_ceiling_Scaling)<=R_IW_KS)Find_smallest_Vector.push_back(Low_ceiling_Scaling);
    }
    if(Component[7]==1){//Top_ceiling
        double Top_ceiling_Scaling   = Scaling_to_others_Z(POS,DR,H_OCT_KS);
        if(Radius_for_ceiling(POS,DR,Top_ceiling_Scaling)<=R_OW_KS)Find_smallest_Vector.push_back(Top_ceiling_Scaling);
    }
    if(Component[8]==1){//Top_reactor_Outer
        double Top_Outer_Reactor_Scaling = Scaling_to_others_Z(POS,DR,H_RE);
        if(Radius_for_ceiling(POS,DR,Top_Outer_Reactor_Scaling)<=R_OWALL_RE and Radius_for_ceiling(POS,DR,Top_Outer_Reactor_Scaling)>R_IWATER_RE)Find_smallest_Vector.push_back(Top_Outer_Reactor_Scaling);
    }
    if(Component[9]==1){//Top_reactor_Inner
        double Top_Inner_Reactor_Scaling = Scaling_to_others_Z(POS,DR,H_RE);
        if(Radius_for_ceiling(POS,DR,Top_Inner_Reactor_Scaling)<=R_IWATER_RE)Find_smallest_Vector.push_back(Top_Inner_Reactor_Scaling);
    }
    double Final_choice = The_smallest_in_a_vector(Find_smallest_Vector);

    return Final_choice;
}
    
//For the component out of the shielding
int DIOOTC(double *POS, double *DR)//DR_In_or_Out_to_Center(DIOOTC)
{
    double Criteria = POS[0]*DR[0]+(POS[1]-30)*DR[1];
    int In_or_Out_N=0;
    if(Criteria>=0)In_or_Out_N=1;
    else{In_or_Out_N=0;}
    
    return In_or_Out_N ;
}
    
double *COFD(double *POS, double *DR)//Criteria_out_of_Shielding(COFD),Position(POS),Direction(DR)
{
    static double RETURN_VALUE[2];
    int *A_Position  = Where_is_it_out_of_shielding(POS);
    int     A_Layer     = A_Position[0];
    int A_Surface_index = A_Position[1];

    int   *B_Position  = Where_is_it_in_shielding(POS);
    int     Compoent    = B_Position[0];
    int B_Surface_index = B_Position[1];
    
    double Scaling_Length;double Layer_run;
    
    if(A_Layer==4 and A_Surface_index!=0)
    {
        double Arrival_of_point[10]={1,1,1,0,1,0,1,0,0,0};
        Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
        Layer_run = 2;
    }
    if(Compoent==1)//Outer Wall of Reactor
    {
        if(A_Surface_index==2)//High
        {
            if(DR[2]>=0)
            {
                double Arrival_of_point[10]={0,0,0,0,1,0,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 2;
            }
            if(DR[2]<0)
            {
                double Arrival_of_point[10]={1,0,1,1,0,0,0,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 0;
            }
        }

        if(A_Surface_index==22)//Side
        {
            if(DIOOTC(POS,DR)>=0)
            {
                double Arrival_of_point[10]={1,1,0,0,1,0,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 2;
            }
        }
        if( (A_Surface_index==22 and DIOOTC(POS,DR)<0) or (A_Surface_index==0) )
        {
                double Arrival_of_point[10]={1,0,1,1,0,0,0,0,1,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 0;
        }
    }
    if(Compoent==2)//Inner Wall of Reactor
    {
        if(A_Surface_index==2)//High
        {
            if(DR[2]>=0)
            {
                double Arrival_of_point[10]={0,0,0,0,1,0,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 2;
            }
            if(DR[2]<0)
            {
                double Arrival_of_point[10]={1,0,0,1,0,0,0,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 1;
            }
        }
        if(A_Surface_index==22)//Side
        {
            if(DIOOTC(POS,DR)>=0)
            {
                double Arrival_of_point[10]={1,0,1,1,0,0,0,0,1,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 2;
            }
        }
        if( (A_Surface_index==22 and DIOOTC(POS,DR)<0) or (A_Surface_index==0) )
        {
                double Arrival_of_point[10]={1,0,0,1,0,0,0,0,0,1};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 1;
        }
    }
    if(Compoent==3)//
    {
        if( (A_Surface_index==1 and DIOOTC(POS,DR)<0) )
        {
                double Arrival_of_point[10]={1,1,1,0,0,0,1,0,1,1};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 2;
        }
        if( (A_Surface_index==1 and DIOOTC(POS,DR)>=0) )
        {
                double Arrival_of_point[10]={1,0,0,0,0,1,0,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 0;
        }
        if( (A_Surface_index==2 and DIOOTC(POS,DR)<0) )
        {
                double Arrival_of_point[10]={1,0,0,0,1,0,0,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 0;
        }
        if( A_Surface_index==0 )
        {
                double Arrival_of_point[10]={1,0,0,0,1,1,0,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 0;
        }
    }
    if(Compoent==4)//
    {
        if( (A_Surface_index==1 and DR[2]<0) )
        {
                double Arrival_of_point[10]={1,1,1,0,1,0,0,0,1,1};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 2;
        }
        if( A_Surface_index==0 )
        {
                double Arrival_of_point[10]={0,0,0,0,0,1,1,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 0;
        }
        if( (A_Surface_index==2) and DR[2]<0 )
        {
                double Arrival_of_point[10]={0,0,0,0,0,1,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                Layer_run = 0;
        }
    }
    if(Compoent==5)//
    {
        Layer_run = 0;
        if( (A_Surface_index==2) and DIOOTC(POS,DR)<0 )
        {
            double Arrival_of_point[10]={0,0,0,0,1,0,1,1,0,0};
            Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
        }
        if( (A_Surface_index==22) and DR[2]<0 )
        {
            double Arrival_of_point[10]={0,0,0,0,1,1,1,0,0,0};
            Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
        }
        if( A_Surface_index==0 )
        {
            double Arrival_of_point[10]={0,0,0,0,1,1,1,1,0,0};
            Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
        }
    }

    RETURN_VALUE[0]=Scaling_Length;RETURN_VALUE[1]=Layer_run;
    return (RETURN_VALUE);
}

double *PAP(double S, double *POS_Int, double *DR)//Pos_Aft_Position
{
    static double POS_Aft[3];
    POS_Aft[0] = POS_Int[0]+(S*DR[0]);
    POS_Aft[1] = POS_Int[1]+(S*DR[1]);
    POS_Aft[2] = POS_Int[2]+(S*DR[2]);
    return POS_Aft;
}

double *VAC(double Mx, double Sigma_SI, double V, double AN)//Velocity_Aft_Collision(VAC),Mx(Mass of WIMP),Sigma_SI, V(Velocity), Atomic Number(AN)
{
    static double Return_Value[2];
    cout << "V: " << V << endl;
    double *V_Aft = Velocity_Aft_collision_Bent(1,Mx,Sigma_SI,V,AN);
    double Energy_Loss_to_Atom   = Energy_DM(Mx,V*1e3/3e8) - Energy_DM(Mx,V_Aft[0]*1e3/3e8);
    double Ratio_of_Energy_Loss_to_Atom = (Energy_Loss_to_Atom)/(Energy_DM(Mx,V*1e3/3e8));
    double DM_Velocity_Aft_Colliding=V_Aft[0];
    cout << "DM_Velocity_Aft_Colliding: " << DM_Velocity_Aft_Colliding << endl;
    Return_Value[0]=DM_Velocity_Aft_Colliding;Return_Value[1]=Ratio_of_Energy_Loss_to_Atom;
    
    return Return_Value;
}
double *DRAC(double AN, double Mx, double V_Aft, double *DR, double ROELTA)//Angle_Aft_Collision,Direction,Ratio_of_Energy_Loss_to_Atom(ROELTA)
{
    static double DR_Aft[3];
    double *Direction_aft = Aft_scatterd_Direction(3,AN,Mx,V_Aft,DR,ROELTA);
    DR_Aft[0] = Direction_aft[0];DR_Aft[1] = Direction_aft[1];DR_Aft[2] = Direction_aft[2];
    return DR_Aft;
}

double *NP(int IorO, double *POS_Int, double *DR, double Mx, double V, double Sigma_SI, int SorT)//InorOut_Shieing(IorO),Next_Point_In_Shielding(NPIS)
{
    static double Return_Value[7];
    
    double LFA=0.001;
    double V_aft=V;double ROELTA=0;
    
    double *SLU;//Scaling_Length_Used,The layer a WIMP runs in
    cout << "IorO : " << IorO << endl;
    if(IorO==0){cout << "InsideOfShielding" << endl;SLU= CID(POS_Int,DR);}
    if(IorO==1){cout << "OutOfShieling" << endl; SLU= COFD(POS_Int,DR);}
    int Layer_run = SLU[1];
    double Segment = 1e3*Length_for_asking_the_collision(LFA,Mx,V_aft,Sigma_SI,Density_All[IorO][Layer_run],Atomic_All[IorO][Layer_run]);//Atomic Number Of Material(ANOM)
    cout << "SLU[0]: " << SLU[0] << endl;
    cout << "Segment: " << Segment << endl;
    if(Segment<=SLU[0])
    {
        POS_Int = PAP(Segment,POS_Int,DR);
        double *VAC_end = VAC(Mx,Sigma_SI,V_aft,Atomic_All[IorO][Layer_run]);//
        V_aft = VAC_end[0]; ROELTA = VAC_end[1];
        if(SorT==1)DR = DRAC(Atomic_All[IorO][Layer_run],Mx,V_aft,DR,ROELTA);
    }
    if(Segment>SLU[0])  POS_Int = PAP(SLU[0],POS_Int,DR);
    
    Return_Value[0]=POS_Int[0];Return_Value[1]=POS_Int[1];Return_Value[2]=POS_Int[2];
    Return_Value[3]=DR[0];Return_Value[4]=DR[1];Return_Value[5]=DR[2];
    Return_Value[6]=V_aft;
    cout << "V_aft_Func: " << V_aft << endl;
    return Return_Value;
}
double *KS_Real_N_With_Angle(int SorT, double Sigma_SI, double V_Int, double Mx, double *DR, double *POS_Int) //Mx(Mass of WIMP),Velocity(km/s) Density(g/cm^3)
{//Straight_or_scattered(SorT)
    static double RETURN_VALUE[20];//Return the value back
    double Direction_VT[3]={0,0,0};//Used for the calculation
    double V_aft=V_Int;double LFA=0.001;//Lamda_for_Average(LFA)
    double Segment; //Ratio_of_Energy_Loss_to_Atom(ROELTA)
    //===========================
    while( RTC(POS_Int)<=R_OW_KS and POS_Int[0]<=H_OCT_KS and Energy_DM(Mx,V_aft*1e3/3e8)>=0.01)
    {
        cout << "V_aft: " << V_aft << endl;
        cout << "Energy_DM: " << Energy_DM(Mx,V_aft*1e3/3e8);
        
        int *IOS_Position = Where_is_it_in_shielding(POS_Int);//Inside_Of_Shieing_Position
        cout << "IOS_Position[0]: " << IOS_Position[0] << endl;
        cout << "IOS_Position[1]: " << IOS_Position[1] << endl;
        cout << "IOS_Position[2]: " << IOS_Position[2] << endl;

        if(IOS_Position[2]==1)
        {
            cout << "Inside the shielding1" << endl;
            double *NP_1 = NP(0,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
        }
        else if(IOS_Position[2]==2 and DR_In_or_Out_on_surface(IOS_Position[1],DR)==0)
        {
            cout << "Inside the shielding2" << endl;
            double *NP_1 = NP(0,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
        }
        else if(IOS_Position[2]==2 and DR_In_or_Out_on_surface(IOS_Position[1],DR)==1)
        {
            cout << "Inside the shielding3" << endl;
            double *NP_1 = NP(1,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
        }
        else
        {
            cout << "Out the shielding" << endl;
            double *NP_1 = NP(1,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
        }
        cout << "V_aft: " << V_aft << endl;
        cout << "POS_Int[0]: " << POS_Int[0] << endl;cout << "POS_Int[1]: " << POS_Int[1] << endl;cout << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "RTC(POS_Int): " << RTC(POS_Int) << ">? " << "R_OW_KS" << R_OW_KS << endl;
        cout << "Energy_DM(Mx,V_aft*1e3/3e8): " << Energy_DM(Mx,V_aft*1e3/3e8) << endl;
    }
    RETURN_VALUE[0]=V_aft;
    if(RTC(POS_Int)>=R_OW_KS) cout << "Radius>Outer_Radius" << endl;
    if(POS_Int[0]>=H_OCT_KS) cout << "High>Cement" << endl;
    if(Energy_DM(Mx,V_aft*1e3/3e8)<=0.01) cout << "Energy<threshold" << endl;


    return RETURN_VALUE;
}


