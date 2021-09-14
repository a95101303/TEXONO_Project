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
        cout << "Vector[0]: " << Vector[0] << endl;
    }
    else{Return_the_smallest_one=0;}
    
    return Return_the_smallest_one;
}

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
double Ge_Layer_Half[3] ={0.4,0.5,0.375};//Layer0, Ge_Detector_Outmost_Layer

double Cu_Layer[3]={Ge_Layer[0]+Cu_TKN*2,Ge_Layer[1]+Cu_TKN*2,Ge_Layer[2]+Cu_TKN*2};//Layer1, OFHC Cu
double B_Layer[3] ={Cu_Layer[0]+B_TKN*2,Cu_Layer[1]+B_TKN*2,Cu_Layer[2]+B_TKN*2};//Layer2, Boron
double Fe_Layer[3]={B_Layer[0]+Fe_TKN*2,B_Layer[1]+Fe_TKN*2,B_Layer[2]+Fe_TKN*2};//Layer3,Steel
double Pb_Layer[3]={Fe_Layer[0]+Pb_TKN*2,Fe_Layer[1]+Pb_TKN*2,Fe_Layer[2]+Pb_TKN*2};//Layer4, Steel

double TCZ= Ge_Layer[2]*0.5+Cu_TKN+B_TKN+Fe_TKN+Pb_TKN;//Thickness_Constrain_Z(TCZ)

double *Starting_Point_with_DR(double *DR)
{
    static double Return_Value[3];
    
    double Nominator[5]={Ge_Layer_Half[0],-Ge_Layer_Half[0],Ge_Layer_Half[1],-Ge_Layer_Half[1],Ge_Layer_Half[2]};//X,Y and Z Plane
    int Denominator[5]={0,0,1,1,2};//Direction
    vector<double> Six_Values;

    for(int kkk=0; kkk<5; kkk++)
    {
        double Ratio=(Nominator[kkk]/DR[Denominator[kkk]]);
        if(Ratio>0)
        {
            Six_Values.push_back(Ratio);
        }
    }
    double Scaling = The_smallest_in_a_vector(Six_Values);
    Return_Value[0]=DR[0]*Scaling;Return_Value[1]=DR[1]*Scaling;Return_Value[2]=DR[2]*Scaling;
    
    return Return_Value;
}


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


int DBT(double A1, double A2)
{
    double A3 = A1-A2;
    if(abs(A3)<1e-10)return 1;
    if(abs(A3)>1e-10)return 0;
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
    if( abs(Position[0])-Ge_Layer[0]*0.5<=1e-10 and abs(Position[1])-Ge_Layer[1]*0.5<=1e-10 and abs(Position[2])-Ge_Layer[2]*0.5<=1e-10 )
    {
        layer_Number=0;
        if(  DBT(Position[0],Ge_Layer[0]*0.5)==1  )On_Surface=1;  if(  DBT(Position[0],-Ge_Layer[0]*0.5)==1  )On_Surface=2;
        if(  DBT(Position[1],Ge_Layer[1]*0.5)==1  )On_Surface=3;  if(  DBT(Position[1],-Ge_Layer[1]*0.5)==1  )On_Surface=4;
        if(  DBT(Position[2],Ge_Layer[2]*0.5)==1  )On_Surface=5;  if(  DBT(Position[2],-Ge_Layer[2]*0.5)==1  )On_Surface=6;

    }
    //Second to fourth layer
    for(int Layer=1; Layer<5; Layer++)
    {
        if( (abs(Position[0])-The_Material_Layer[Layer][0]*0.5<1e-10) and (abs(Position[1])-The_Material_Layer[Layer][1]*0.5<1e-10) and (abs(Position[2])-The_Material_Layer[Layer][2]*0.5<1e-10) )
        {
            if( (abs(Position[0])-The_Material_Layer[Layer-1][0]*0.5>1e-10) or (abs(Position[1])-The_Material_Layer[Layer-1][1]*0.5>1e-10) or (abs(Position[2])-The_Material_Layer[Layer-1][2]*0.5>1e-10) )
            {
                /*
                cout << "The_Material_Layer[Layer][0]: " << The_Material_Layer[Layer][0]*0.5 << endl;
                cout << "The_Material_Layer[Layer-1][0]*0.5: " << The_Material_Layer[Layer-1][0]*0.5 << endl;
                cout << "The_Material_Layer[Layer][1]: " << The_Material_Layer[Layer][1]*0.5 << endl;
                cout << "The_Material_Layer[Layer-1][1]*0.5: " << The_Material_Layer[Layer-1][1]*0.5 << endl;
                cout << "The_Material_Layer[Layer][2]: " << The_Material_Layer[Layer][2]*0.5 << endl;
                cout << "The_Material_Layer[Layer-1][2]*0.5: " << The_Material_Layer[Layer-1][2]*0.5 << endl;
                 */
                layer_Number=Layer;
                if(  DBT(Position[0],The_Material_Layer[Layer][0]*0.5)==1  )On_Surface=1;  if(  DBT(Position[0],-The_Material_Layer[Layer][0]*0.5)==1  )On_Surface=2;
                if(  DBT(Position[1],The_Material_Layer[Layer][1]*0.5)==1  )On_Surface=3;  if(  DBT(Position[1],-The_Material_Layer[Layer][1]*0.5)==1  )On_Surface=4;
                if(  DBT(Position[2],The_Material_Layer[Layer][2]*0.5)==1  )On_Surface=5;  if(  DBT(Position[2],-The_Material_Layer[Layer][2]*0.5)==1  )On_Surface=6;
            }
        }
    }
    //Out of the fourth layer
    if(  abs(Position[0])-Pb_Layer[0]*0.5>1e-10 or  abs(Position[1])-Pb_Layer[1]*0.5>1e-10  or  abs(Position[2])-Pb_Layer[2]*0.5>1e-10 )
    {
        /*
        cout << "Pb_Layer[0]*0.5: " << Pb_Layer[0]*0.5 << endl;
        cout << "Pb_Layer[1]*0.5: " << Pb_Layer[1]*0.5 << endl;
        cout << "Pb_Layer[2]*0.5: " << Pb_Layer[2]*0.5 << endl;
        cout << "abs(Position[0])-Pb_Layer[0]*0.5: " << abs(Position[0])-Pb_Layer[0]*0.5 << endl;
        cout << "abs(Position[1])-Pb_Layer[1]*0.5: " << abs(Position[1])-Pb_Layer[1]*0.5 << endl;
        cout << "abs(Position[2])-Pb_Layer[2]*0.5: " << abs(Position[2])-Pb_Layer[2]*0.5 << endl;
         */
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
            for(int POS_After_Axis=0; POS_After_Axis<3; POS_After_Axis++){POS_Aft[POS_After_Axis] = POS[POS_After_Axis] + ET*DR[POS_After_Axis];}
            int *Fun_Layer_Aft= Where_is_it_in_shielding(POS_Aft);
            /*
            cout << "==================================" << endl;
            cout << "XYZ: " << XYZ << "PorN: " << PorN << endl;
            cout << "Confirmed_Side[PorN]*The_Material_Layer[Layer][XYZ]*0.5: " << Confirmed_Side[PorN]*The_Material_Layer[Layer][XYZ]*0.5 << endl;
            cout << "DR[0]: " << DR[0] << "DR[1]: " << DR[1] << "DR[2]: " << DR[2] << endl;
            cout << "POS_Aft[0]: " << POS_Aft[0] << endl;cout << "POS_Aft[1]: " << POS_Aft[1] << endl;cout << "POS_Aft[2]: " << POS_Aft[2] << endl;
            cout << "ET: " << ET << endl;
            cout << "Layer: " << Layer << endl;
            cout << "Fun_Layer_Aft[0]: " << Fun_Layer_Aft[0] << endl;
            cout << "==================================" << endl;
             */
            if(Fun_Layer_Aft[0]==Layer){ET_XYZ[XYZ][PorN]=ET;}else{ET_XYZ[XYZ][PorN]=0;}
            if(Fun_Layer_Aft[0]==Layer and ET>0){Find_smallest_Vector.push_back(ET);cout << "ET1: " << ET << endl;}
        }
    }
    
    The_final_scaling = The_smallest_in_a_vector(Find_smallest_Vector);
    cout << "Length_to_six_planes(The_final_scaling): " << The_final_scaling << endl;
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

double FSBTLS(int A_Layer, double *POS, double *DR)//Find_Scaling_Between_Two_Layers(FSBTLS)
{
    double RSF;
    vector<double> Return_the_smallest_one;
    double A1=Length_to_six_planes(A_Layer,POS,DR);
    double A2=Length_to_six_planes(A_Layer-1,POS,DR);
    if(A1>0) Return_the_smallest_one.push_back(A1);if(A2>0) Return_the_smallest_one.push_back(A2);
    RSF = The_smallest_in_a_vector(Return_the_smallest_one);
    
    return RSF;
}
double *CID(double *POS, double *DR)//Criteria_In_Shielding(CID),Position(POS),Direction(DR)
{
    static double Return_Factor[2];//Return_Scaling_Factor, Layer the WIMP runs in
    double RSF;
    int *A_Position       = Where_is_it_in_shielding(POS);
    int  A_Layer          = A_Position[0];
    cout << "A_Layer: " << A_Layer  << endl;
    int  A_Surface_index  = A_Position[1];
    cout << "A_Surface_index: " << A_Surface_index << endl;
    int  DR_IN_OR_OUT = 0; double Layer_run=0;
    vector<double> Return_the_smallest_one;

    //In the shielding
    if(  A_Layer<4 )
    {
        if(A_Surface_index!=0)//On the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_on_surface(A_Surface_index,DR);
            cout << "DR_IN_OR_OUT: " << DR_IN_OR_OUT << endl;
            if(DR_IN_OR_OUT==1)
            {
                cout << "L" << A_Layer << "OS OutDR" << endl;
                RSF           = Length_to_six_planes(A_Layer+1,POS,DR);
                Layer_run =  A_Layer+1;
            }
            if( (A_Layer>0 and DR_IN_OR_OUT==0) )
            {
               cout << "L" << A_Layer << "OS InDR" << endl;
               RSF = FSBTLS(A_Layer,POS,DR);
               Layer_run =  A_Layer;
            }
        }
        if(A_Surface_index==0)//Off the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_off_surface(POS,DR);
            Layer_run =  A_Layer;
            RSF = FSBTLS(A_Layer,POS,DR);
            /*
            if(DR_IN_OR_OUT==1)
            {
                cout << "L" << A_Layer << "FS OutDR" << endl;
                RSF = FSBTLS(A_Layer,POS,DR);
            }
            if(DR_IN_OR_OUT==0)
            {
                cout << "L" << A_Layer << "FS InDR" << endl;
                RSF = FSBTLS(A_Layer,POS,DR);
            }
             */
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
                RSF = FSBTLS(A_Layer,POS,DR);
            }
        }
        if(A_Surface_index==0)//Off the surface
        {
            DR_IN_OR_OUT = DR_In_or_Out_off_surface(POS,DR);
            if(DR_IN_OR_OUT==1)
            {
                RSF = FSBTLS(A_Layer,POS,DR);
            }
            if(DR_IN_OR_OUT==0)
            {
                RSF = FSBTLS(A_Layer,POS,DR);
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
    //cout << "Numerator1: " << Numerator1 << endl;
    if(Numerator1>1e-10)Find_smallest_Vector.push_back(Numerator1);
    Numerator2 = -b - sqrt(b*b-4*a*c);
    //cout << "Numerator2: " << Numerator2 << endl;
    if(Numerator2>1e-10)Find_smallest_Vector.push_back(Numerator2);
    
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
double RFC(double *POS, double *DR, double S)//Radius_for_Ceiling(RFC),Position(POS),Direction(DR),Scaling(S),Radius(R)
{
    double X = POS[0]+DR[0]*S;
    double Y = (POS[1]-30)+DR[1]*S;
    double Z = POS[2]+DR[2]*S;
    cout << "X: " << X << "Y: " << Y << "Z: " << Z << endl;
    return sqrt(X*X+Y*Y);
}
double RTC(double *POS)//Radius_To_Center(RTC) on XY plane
{
    double X = (POS[0]-0);
    double Y = (POS[1]-30);
    return sqrt(X*X+Y*Y);
}

int Yes_123(double A1, double A2)
{
    //cout << "A1-A2: " << (A1-A2) << endl;
    if( (A1-A2)>-1e-10 )return 1;

    else{return 0;}
}
int *Where_is_it_out_of_shielding(double *POS)
{
    static int Return_Value[2];
    int Component=0;int On_Surface=0;
    cout << "Where_is_it_out_of_shielding_START" << endl;
    cout << "R_IWATER_RE - RTC(POS)" << R_IWATER_RE - RTC(POS) << endl;
    cout << "R_OWALL_RE-RTC(POS)"    << R_OWALL_RE - RTC(POS) << endl;
    cout << "R_OW_KS-RTC(POS)" << R_OW_KS-RTC(POS) << endl;
    cout << "R_IW_KS-RTC(POS)" << R_IW_KS-RTC(POS) << endl;
    cout << "H_RE-POS[2]" << H_RE-POS[2] << endl;
    cout << "H_ICT_KS-POS[2]" << H_ICT_KS-POS[2] << endl;
    cout << "H_OCT_KS-POS[2]" << H_OCT_KS-POS[2] << endl;
    Yes_123(H_RE,POS[2]);
    if( Yes_123(RTC(POS),R_IWATER_RE)==1 and Yes_123(R_OWALL_RE,RTC(POS))==1 and Yes_123(H_RE,POS[2])==1 )//Outer Wall of Reactor
    {
        cout << "Outer_Wall_of_Reactor" << endl;
        Component=1;
        if(DBT(POS[2],H_RE)==1)On_Surface=2;if(DBT(RTC(POS),R_OWALL_RE)==1)On_Surface=22;
    }
    if( Yes_123(R_IWATER_RE,RTC(POS))==1 and Yes_123(H_RE,POS[2])==1 )//Inner Wall of Reactor
    {
        cout << "Inner_Wall_of_Reactor" << endl;
        Component=2;
        if(DBT(POS[2],H_RE)==1)On_Surface=2;if(DBT(RTC(POS),R_IWATER_RE)==1)On_Surface=22;
    }
    if( Yes_123(R_OW_KS,RTC(POS))==1 and Yes_123(RTC(POS),R_IW_KS)==1 and Yes_123(H_ICT_KS,POS[2])==1 )//Side Wall of KS
    {
        cout << "Side Wall of KS" << endl;
        Component=3;//Wall of KS
        if(DBT(RTC(POS),R_IW_KS)==1)On_Surface=1;if(DBT(RTC(POS),R_OW_KS)==1)On_Surface=2;
    }
    if( Yes_123(R_IW_KS,RTC(POS))==1 and Yes_123(POS[2],H_ICT_KS)==1 and Yes_123(H_OCT_KS,POS[2])==1 )//Upper ceiling of KS
    {
        cout << "Upper ceiling of KS" << endl;
        Component=4;//Wall of KS
        if(DBT(POS[2],H_ICT_KS)==1)On_Surface=1;if(DBT(POS[2],H_OCT_KS)==1)On_Surface=2;
    }
    if( Yes_123(RTC(POS),R_IW_KS)==1 and Yes_123(R_OW_KS,RTC(POS))==1 and Yes_123(POS[2],H_ICT_KS)==1 and Yes_123(H_OCT_KS,POS[2])==1)//Corner of the wall
    {
        cout << "Corner of the wall of KS" << endl;
        Component=5;
        if(DBT(RTC(POS),R_OW_KS)==1)On_Surface=2;if(DBT(POS[2],H_OCT_KS)==1)On_Surface=22;
    }

    Return_Value[0]=Component;Return_Value[1]=On_Surface;
    cout << "On_Surface" << On_Surface << endl;
    cout << "Where_is_it_out_of_shielding_END" << endl;

    return(Return_Value);
}

double DOLOUD(double *Component, double *POS, double *DR)//Decision_On_Length_Out_Of_Detector(DOLOUD)
{
    vector<double> Find_smallest_Vector;
    cout << "DOLOUD_START " << endl;
    if(Component[0]==1){//Floor
        cout << "Floor" << endl;
        double DFS    = Scaling_to_others_Z(POS,DR,-1.375);//Down_floor_Scaling(DFS)
        cout << "Down_floor_Scaling: " << DFS  << endl;
        if(Yes_123(R_IW_KS,RFC(POS,DR,DFS))==1 and DFS>1e-10)
        {
            cout << "DFS_Added" << endl;
            Find_smallest_Vector.push_back(DFS);
        }
    }
    if(Component[1]==1){//Shielding
        cout << "Shielding" << endl;
        double SS     = Length_to_six_planes(4,POS,DR);//Shielding_Scaling(SS)
        cout << "Shielding_Scaling: " << SS << endl;
        if(SS>1e-10)
        {
            cout << "SS_Added" << endl;
            Find_smallest_Vector.push_back(SS);
        }
    }
    if(Component[2]==1){//Outer_Reactor
        cout << "Outer_Reactor" << endl;
        double ORS = Scaling_to_others_XY(POS,DR,R_OWALL_RE);//Outer_RE_Scaling(ORS)
        cout << "Outer_RE_Scaling: " << ORS << endl;
        if(Yes_123(H_RE,POS[2]+ORS*DR[2])==1 and ORS>1e-10)
        {
            cout << "ORS_Added" << endl;
            Find_smallest_Vector.push_back(ORS);
        }
    }
    if(Component[3]==1){//Inner_Reactor
        cout << "Inner_Reactor" << endl;
        double IRS = Scaling_to_others_XY(POS,DR,R_IWATER_RE);//Inner_RE_Scaling(IRS)
        cout << "Inner_RE_Scaling: " << IRS << endl;
        if(Yes_123(H_RE,POS[2]+IRS*DR[2])==1 and IRS>1e-10)
        {
            cout << "IRS_Added" << endl;
            Find_smallest_Vector.push_back(IRS);
        }
    }
    if(Component[4]==1){//Inner_Wall
        cout << "Inner_Wall" << endl;
        double IWS = Scaling_to_others_XY(POS,DR,R_IW_KS);//Inner_Wall_Scaling(IWS)
        cout << "Inner_Wall_Scaling: " << IWS << endl;
        if(Yes_123(H_ICT_KS,POS[2]+IWS*DR[2])==1 and IWS>1e-10)
        {
            cout << "IWS_Added" << endl;
            Find_smallest_Vector.push_back(IWS);
        }
    }
    if(Component[5]==1){//Outer_Wall
        cout << "Outer_Wall" << endl;
        double OWS = Scaling_to_others_XY(POS,DR,R_OW_KS);//Outer_Wall_Scaling(OWS)
        cout << "Outer_Wall_Scaling: " << OWS << endl;
        if(Yes_123(H_OCT_KS,POS[2]+OWS*DR[2])==1 and OWS>1e-10)
        {
            cout << "OWS_Added" << endl;
            Find_smallest_Vector.push_back(OWS);
        }
    }
    if(Component[6]==1){//Low_ceiling
        cout << "Low_ceiling" << endl;
        double LCS = Scaling_to_others_Z(POS,DR,H_ICT_KS);//Low_Ceiling_Scaling(LCS)
        cout << "Low_ceiling_Scaling: " << LCS << endl;
        if(Yes_123(R_IW_KS,RFC(POS,DR,LCS))==1 and LCS>1e-10)
        {
            cout << "LCS_Added" << endl;
            Find_smallest_Vector.push_back(LCS);
        }
    }
    if(Component[7]==1){//Top_ceiling
        cout << "Top_ceiling" << endl;
        double TCS   = Scaling_to_others_Z(POS,DR,H_OCT_KS);//Top_Ceiling_Scaling(TCS)
        cout << "Top_ceiling_Scaling: " << TCS << endl;
        if(Yes_123(R_OW_KS,RFC(POS,DR,TCS))==1 and TCS>1e-10)
        {
            cout << "TCS_Added" << endl;
            Find_smallest_Vector.push_back(TCS);
        }
    }
    if(Component[8]==1){//Top_reactor_Outer
        cout << "Top_reactor_Outer" << endl;
        double TORS = Scaling_to_others_Z(POS,DR,H_RE);//Top_Outer_Reactor_Scaling(TORS)
        cout << "Top_Outer_Reactor_Scaling: " << TORS << endl;
        cout << "RFC(POS,DR,TORS): " << RFC(POS,DR,TORS) << endl;
        if(Yes_123(R_OWALL_RE,RFC(POS,DR,TORS))==1 and Yes_123(RFC(POS,DR,TORS),R_IWATER_RE)==1 and TORS>1e-10)
        {
            cout << "TORS_Added" << endl;
            Find_smallest_Vector.push_back(TORS);
        }
    }
    if(Component[9]==1){//Top_reactor_Inner
        cout << "Top_reactor_Inner" << endl;
        double TIRS = Scaling_to_others_Z(POS,DR,H_RE);//Top_Inner_Reactor_Scaling(TIRS)
        cout << "Top_Inner_Reactor_Scaling: " << TIRS << endl;
        cout << "RFC(POS,DR,TIRS): " << RFC(POS,DR,TIRS) << endl;
        if(Yes_123(R_IWATER_RE,RFC(POS,DR,TIRS))==1 and TIRS>1e-10)
        {
            cout << "TIRS_Added" << endl;
            Find_smallest_Vector.push_back(TIRS);
        }
    }
    cout << "DOLOUD_END " << endl;
    double Final_choice = The_smallest_in_a_vector(Find_smallest_Vector);

    return Final_choice;
}
    
//For the component out of the shielding
int DIOOTC(double *POS, double *DR)//DR_In_or_Out_to_Center(DIOOTC)
{
    double Criteria = POS[0]*DR[0]+(POS[1]-30)*DR[1];
    cout << "DIOOTC_Criteria" << Criteria << endl;
    int In_or_Out_N=0;
    if(Criteria>=0)In_or_Out_N=1;
    else{In_or_Out_N=0;}
    
    return In_or_Out_N ;
}
    
double *COFD(double *POS, double *DR)//Criteria_out_of_Shielding(COFD),Position(POS),Direction(DR)
{
    cout << "COFD_STATLE!!" << endl;
    static double RETURN_VALUE[3];
    int *A_Position  = Where_is_it_in_shielding(POS);
    int     A_Layer     = A_Position[0];
    int A_Surface_index = A_Position[1];
    cout << "A_Layer: " << A_Layer << endl;
    cout << "A_Surface_index: " << A_Surface_index << endl;

    int   *B_Position  = Where_is_it_out_of_shielding(POS);
    int     Compoent    = B_Position[0];
    int B_Surface_index = B_Position[1];
    cout << "Compoent: " << Compoent << endl;
    cout << "B_Surface_index: " << B_Surface_index << endl;
    
    double Scaling_Length;double Layer_run;
    cout << "END1" << endl;
    if(A_Layer==4 and A_Surface_index!=0)
    {
        double Arrival_of_point[10]={1,1,1,0,1,0,1,0,0,0};
        Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
        cout << "A_Layer==4 and A_Surface_index!=0: " << endl;
        cout << "Scaling_Length: " << Scaling_Length << endl;
        Layer_run = 2;
    }
    if(Compoent==0)//Air
    {
        double Arrival_of_point[10]={1,1,1,0,1,0,1,0,1,1};
        Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
        cout << "Compoent==0" << endl;
        cout << "Scaling_Length: " << Scaling_Length << endl;
        Layer_run = 2;
    }
    if(Compoent==1)//Outer Wall of Reactor
    {
        cout << "Compoent==1" << endl;
        if(B_Surface_index==2)//High
        {
            if(DR[2]>=0)
            {
                double Arrival_of_point[10]={0,0,0,0,1,0,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                cout << "B_Surface_index==2, DR[2]>=0" << endl;
                cout << "Scaling_Length: " << Scaling_Length << endl;
                Layer_run = 2;
            }
            if(DR[2]<0)
            {
                double Arrival_of_point[10]={1,0,1,1,0,0,0,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                cout << "B_Surface_index==2, DR[2]<0" << endl;
                cout << "Scaling_Length: " << Scaling_Length << endl;
                Layer_run = 0;
            }
        }

        if(B_Surface_index==22)//Side
        {
            if(DIOOTC(POS,DR)==1)
            {
                double Arrival_of_point[10]={1,1,0,0,1,0,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                cout << "B_Surface_index==22, DIOOTC(POS,DR)==1" << endl;
                cout << "Scaling_Length: " << Scaling_Length << endl;
                Layer_run = 2;
            }
        }
        if( (B_Surface_index==22 and DIOOTC(POS,DR)==0) or (B_Surface_index==0) )
        {
                double Arrival_of_point[10]={1,0,1,1,0,0,0,0,1,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "DIOOTC(POS,DR): " << DIOOTC(POS,DR) << endl;
            cout << "B_Surface_index==22, DIOOTC(POS,DR)==0 or (B_Surface_index==0)" << endl;

                Layer_run = 0;
        }
    }
    if(Compoent==2)//Inner Wall of Reactor
    {
        cout << "Compoent==2" << endl;
        if(B_Surface_index==2)//High
        {
            if(DR[2]>=0)
            {
                double Arrival_of_point[10]={0,0,0,0,1,0,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                cout << "B_Surface_index==2, DR[2]>=0" << endl;
                cout << "Scaling_Length: " << Scaling_Length << endl;
                Layer_run = 2;
            }
            if(DR[2]<0)
            {
                double Arrival_of_point[10]={1,0,0,1,0,0,0,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                cout << "B_Surface_index==2, DR[2]<0" << endl;
                cout << "Scaling_Length: " << Scaling_Length << endl;
                Layer_run = 1;
            }
        }
        if(B_Surface_index==22)//Side
        {
            if(DIOOTC(POS,DR)==1)
            {
                double Arrival_of_point[10]={1,0,1,1,0,0,0,0,1,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
                cout << "A_Surface_index==22, DIOOTC(POS,DR)==0" << endl;
                cout << "Scaling_Length: " << Scaling_Length << endl;
                Layer_run = 2;
            }
        }
        if( (B_Surface_index==22 and DIOOTC(POS,DR)==0) or (B_Surface_index==0) )
        {
                double Arrival_of_point[10]={1,0,0,1,0,0,0,0,0,1};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "B_Surface_index==22, DIOOTC(POS,DR)==0 or (B_Surface_index==0) " << endl;
            cout << "Scaling_Length: " << Scaling_Length << endl;
                Layer_run = 1;
        }
    }
    if(Compoent==3)//
    {
        cout << "Compoent==3" << endl;
        if( (B_Surface_index==1 and DIOOTC(POS,DR)==0) )
        {
                double Arrival_of_point[10]={1,1,1,0,1,0,1,0,1,1};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "(B_Surface_index==1 and DIOOTC(POS,DR)==0) " << endl;
                Layer_run = 2;
        }
        if( (B_Surface_index==1 and DIOOTC(POS,DR)==1) )
        {
                double Arrival_of_point[10]={1,0,0,0,0,1,0,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "(B_Surface_index==1 and DIOOTC(POS,DR)==1)" << endl;
                Layer_run = 0;
        }
        if( (B_Surface_index==2 and DIOOTC(POS,DR)==0) )
        {
                double Arrival_of_point[10]={1,0,0,0,1,0,0,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "(B_Surface_index==2 and DIOOTC(POS,DR)==0)" << endl;
                Layer_run = 0;
        }
        if( B_Surface_index==0 )
        {
                double Arrival_of_point[10]={1,0,0,0,1,1,0,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "B_Surface_index==0 " << endl;

                Layer_run = 0;
        }
    }
    if(Compoent==4)//
    {
        cout << "Compoent==4" << endl;
        if( (B_Surface_index==1 and DR[2]<0) )
        {
                double Arrival_of_point[10]={1,1,1,0,1,0,0,0,1,1};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "(B_Surface_index==1 and DR[2]<0)" << endl;
                Layer_run = 2;
        }
        if( B_Surface_index==0 or (B_Surface_index==1 and DR[2]>=0) )
        {
                double Arrival_of_point[10]={0,0,0,0,0,1,1,1,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "B_Surface_index==0" << endl;
                Layer_run = 0;
        }
        if( (B_Surface_index==2) and DR[2]<0 )
        {
                double Arrival_of_point[10]={0,0,0,0,0,1,1,0,0,0};
                Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "(B_Surface_index==2) and DR[2]<0 " << endl;
                Layer_run = 0;
        }
    }
    if(Compoent==5)//
    {
        cout << "Compoent==5" << endl;
        Layer_run = 0;
        if( (B_Surface_index==2) and DIOOTC(POS,DR)==0 )
        {
            double Arrival_of_point[10]={0,0,0,0,1,0,1,1,0,0};
            Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "(B_Surface_index==2) and DIOOTC(POS,DR)==0" << endl;

        }
        if( (B_Surface_index==22) and DR[2]<0 )
        {
            double Arrival_of_point[10]={0,0,0,0,1,1,1,0,0,0};
            Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "Scaling_Length: " << Scaling_Length << endl;
            cout << "(B_Surface_index==22) and DR[2]<0" << endl;
        }
        if( (B_Surface_index==0) or ( (B_Surface_index==2) and DIOOTC(POS,DR)==1 ) )
        {
            double Arrival_of_point[10]={0,0,0,0,1,1,1,1,0,0};
            Scaling_Length = DOLOUD(Arrival_of_point,POS,DR);
            cout << "B_Surface_index==0" << endl;
            cout << "Scaling_Length: " << Scaling_Length << endl;
        }
    }
    cout << "END2" << endl;

    RETURN_VALUE[0]=Scaling_Length;RETURN_VALUE[1]=Layer_run;RETURN_VALUE[2]=Compoent;
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
    static double Return_Value[10];
    
    double LFA=0.001;
    double V_aft=V;double ROELTA=0;
    int Collision_Time=0;
    double *SLU;//Scaling_Length_Used,The layer a WIMP runs in
    cout << "IorO : " << IorO << endl;
    if(IorO==0){cout << "InsideOfShielding" << endl;SLU= CID(POS_Int,DR);}
    if(IorO==1){cout << "OutOfShieling" << endl; SLU= COFD(POS_Int,DR);}
    int Layer_run = SLU[1];//Inside the shielding
    int Compoent  = SLU[2];//Out of the shielding
    cout << "Compoent_Check: " << Compoent << endl;
    cout << "SLU[0]:Going Length: " << SLU[0] << endl;

    int Times = Possion_GetRandom_Full(LFA);
    double Segment = 1e3*Length_for_asking_the_collision(LFA,Mx,V_aft,Sigma_SI,Density_All[IorO][Layer_run],Atomic_All[IorO][Layer_run]);//Atomic Number Of Material(ANOM)
    cout << "IorO: " << IorO << "Layer_run: " << Layer_run << endl;
    cout << "Atomic_All[IorO][Layer_run]: " << Atomic_All[IorO][Layer_run] << endl;
    cout << "Segment: " << Segment << endl;
    double Sprint_GO = Segment*Times;
    if(Sprint_GO<=SLU[0])
    {
        cout << "Segment<=SLU[0]" << endl;
        POS_Int = PAP(Sprint_GO,POS_Int,DR);
        double *VAC_end = VAC(Mx,Sigma_SI,V_aft,Atomic_All[IorO][Layer_run]);//
        cout << "Atomic_All[IorO][Layer_run]: " << Atomic_All[IorO][Layer_run] << endl;
        V_aft = VAC_end[0]; ROELTA = VAC_end[1];
        if(SorT==1)DR = DRAC(Atomic_All[IorO][Layer_run],Mx,V_aft,DR,ROELTA);
        Collision_Time=1;
    }
    if(Sprint_GO>SLU[0])
    {
        cout << "Segment>SLU[0]" << endl;
        POS_Int = PAP(SLU[0],POS_Int,DR);
    }
    Return_Value[0]=POS_Int[0];Return_Value[1]=POS_Int[1];Return_Value[2]=POS_Int[2];
    Return_Value[3]=DR[0];Return_Value[4]=DR[1];Return_Value[5]=DR[2];
    Return_Value[6]=V_aft;Return_Value[7]=Collision_Time;Return_Value[8]=Layer_run;Return_Value[9]=Compoent;
    
    cout << "V_aft_Func: " << V_aft << endl;
    cout << "Layer_run: " << Layer_run << endl;
    return Return_Value;
}
//KS-Angle
double *KS_Real_N_With_Angle(int SorT, double Sigma_SI, double V_Int, double Mx, double *DR, double *POS_Int) //Mx(Mass of WIMP),Velocity(km/s) Density(g/cm^3)
{//Straight_or_scattered(SorT)
    static double RETURN_VALUE[20];//Return the value back
    double Collision_Time1[5]={0,0,0,0,0};
    double Collision_Time2[6]={0,0,0,0,0,0};

    double Direction_VT[3]={0,0,0};//Used for the calculation
    double V_aft=V_Int;double LFA=0.001;//Lamda_for_Average(LFA)
    double Segment; //Ratio_of_Energy_Loss_to_Atom(ROELTA)
    int Step=0;int Collision_Time=0;
    //===========================
    while( (R_OW_KS-RTC(POS_Int))>1e-10 and (H_OCT_KS-POS_Int[2])>1e-10 and POS_Int[2]+(TCZ)>1e-10 and Energy_DM(Mx,V_aft*1e3/3e8)>=0.01)
    {
        cout << "V_aft: " << V_aft << endl;
        cout << "Energy_DM: " << Energy_DM(Mx,V_aft*1e3/3e8);
        
        int *IOS_Position = Where_is_it_in_shielding(POS_Int);//Inside_Of_Shieing_Position
        cout << "DR[0]: " << DR[0] << "DR[1]:" << DR[1] << "DR[2]:" << DR[2] << endl;
        cout << "IOS_Position[0]: " << IOS_Position[0] << endl;//Layer
        cout << "IOS_Position[1]: " << IOS_Position[1] << endl;//On which surface
        cout << "IOS_Position[2]: " << IOS_Position[2] << endl;//
        if(IOS_Position[0]==0 and Step!=0)
        {
            RETURN_VALUE[2]=0;
            cout << "OKOKOKOK!!!!!" << endl;
            break;
        }
        
        if(IOS_Position[2]==1)
        {
            cout << "Inside the shielding1" << endl;
            double *NP_1 = NP(0,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
            if(int(NP_1[7])==1)
            {
                cout << "int(NP_1[8]): " << int(NP_1[8]) << endl;
                Collision_Time1[int(NP_1[8])] = Collision_Time1[int(NP_1[8])] + 1;
            }
            Step = Step + 1;
        }
        else if(IOS_Position[2]==2 and DR_In_or_Out_on_surface(IOS_Position[1],DR)==0)
        {
            cout << "Inside the shielding2" << endl;
            double *NP_1 = NP(0,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
            if(int(NP_1[7])==1)
            {
                cout << "int(NP_1[8]): " << int(NP_1[8]) << endl;
                Collision_Time1[int(NP_1[8])] = Collision_Time1[int(NP_1[8])] + 1;
            }
            Step = Step + 1;
        }
        else if(IOS_Position[2]==2 and DR_In_or_Out_on_surface(IOS_Position[1],DR)==1)
        {
            cout << "Inside the shielding3" << endl;
            double *NP_1 = NP(1,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
            if(int(NP_1[7])==1)
            {
                cout << "int(NP_1[8]): " << int(NP_1[8]) << endl;
                Collision_Time1[int(NP_1[8])] = Collision_Time1[int(NP_1[8])] + 1;
            }
            Step = Step + 1;
        }
        else
        {
            cout << "Out the shielding" << endl;
            double *NP_1 = NP(1,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
            if(int(NP_1[7])==1)
            {
                cout << "int(NP_1[9]): " << int(NP_1[9]) << endl;
                Collision_Time2[int(NP_1[9])] = Collision_Time2[int(NP_1[9])] + 1;
            }
            Step = Step + 1;
        }
        cout << "V_aft: " << V_aft << endl;
        cout << "POS_Int[0]: " << POS_Int[0] << endl;cout << "POS_Int[1]: " << POS_Int[1] << endl;cout << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "RTC(POS_Int): " << RTC(POS_Int) << ">? " << "R_OW_KS" << R_OW_KS << endl;
        cout << "RTC(POS_Int)-R_OW_KS: " << RTC(POS_Int)-R_OW_KS << endl;
        cout << "========" << endl;
        cout << "R_OW_KS-RTC(POS_Int): " << R_OW_KS-RTC(POS_Int) << endl;
        cout << "H_OCT_KS-POS_Int[2]: " << H_OCT_KS-POS_Int[2] << endl;
        cout << "POS_Int[2]+(TCZ): " << POS_Int[2]+(TCZ) << endl;
        cout << "Energy_DM(Mx,V_aft*1e3/3e8): " << Energy_DM(Mx,V_aft*1e3/3e8) << endl;
    }
    
    for(int kkk=0; kkk<5; kkk++) Collision_Time = Collision_Time + Collision_Time1[kkk];
    RETURN_VALUE[0]=V_aft;RETURN_VALUE[1]=Collision_Time;
    cout << "Collision_Time: " << Collision_Time << endl;
    if(V_aft>100) cout << "YoMAN Cool!" << endl;
    if((R_OW_KS-RTC(POS_Int))<1e-10)
    {
        if(POS_Int[2]<13.)
        {
            cout << "Through Earth: " << endl;
            RETURN_VALUE[2]=0;
        }
        else
        {
        RETURN_VALUE[2]=1;
        cout << "Radius>Outer_Radius" << endl;
        }
    }
    if(H_OCT_KS-POS_Int[2]<1e-10)
    {
        RETURN_VALUE[2]=1;
        cout << "High>Cement" << endl;
    }
    if(POS_Int[2]+(TCZ)<1e-10)
    {
        RETURN_VALUE[2]=0;
        cout << "Z at the ground" << endl;
    }
    if(Energy_DM(Mx,V_aft*1e3/3e8)<0.01)
    {
        RETURN_VALUE[2]=0;
        cout << "Energy<threshold" << endl;
    }

    for(int kkk=1; kkk<5; kkk++){cout << "Collision_Time1: " << Collision_Time1[kkk] << endl;}
    for(int jjj=1; jjj<6; jjj++){cout << "Collision_Time2: " << Collision_Time2[jjj] << endl;}


    return RETURN_VALUE;
}

double SqrtN2(double *A)
{
    return sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
}

double SqrtN2_ABC(double A, double B, double C)
{
    //cout << "sqrt(A*A+B*B+C*C): " << sqrt(A*A+B*B+C*C) << endl;
    return sqrt(A*A+B*B+C*C);
}

void Check_Fun(double *POS_Int, double *DR, double S)
{
    //cout << "Lcjpl_NonNL(POS_Int[0]+S*DR[0],POS_Int[1]+S*DR[1],POS_Int[2]+S*DR[2]): " << Lcjpl_NonNL(POS_Int[0]+S*DR[0],POS_Int[1]+S*DR[1],POS_Int[2]+S*DR[2]) << endl;
    //cout << "SqrtN2_ABC([POS_Int[0]+S*DR[0],POS_Int[1]+S*DR[1],POS_Int[2]+S*DR[2]): " << SqrtN2_ABC(POS_Int[0]+S*DR[0],POS_Int[1]+S*DR[1],POS_Int[2]+S*DR[2]) << endl;
                                                                                                     
     return 0;
}

double *DAC_to_LEG(double *POS_Int, double *DR, double PL)//Direction_Aft_Collision_to_Length(PAC),Predicted_Length(PL)
{
    static double POS_Aft[5];
    //POS_Aft[3]=Count or not
    
    POS_Aft[0] = POS_Int[0]+PL*DR[0];
    POS_Aft[1] = POS_Int[1]+PL*DR[1];
    POS_Aft[2] = POS_Int[2]+PL*DR[2];

    double ML  = Lcjpl_NonNL(POS_Aft[0],POS_Aft[1],POS_Aft[2]);//Mountain_Length(ML)
    double LML = SqrtN2(POS_Aft);//Later_Moving_Length
    
    double STU; //Boundary Or Not
    
    if(0-POS_Aft[2]>=1e-10)
    {
        cout << "Yes1 " << endl;
        cout << "Below the mountain^>^" << endl;
        POS_Aft[3]=-1;//Below the mountain, counted as an event from the bottom earth and will be blocked at the boundary
    }
    else if(LML<ML)
    {
        cout << "Inside the mountain^>^ " << endl;
        POS_Aft[3]=0;//Inside the mountain and this event has an interaction
    }
    else if(LML>ML)
    {
        cout << "Above the mountain^>^ " << endl;
        POS_Aft[3]=1;//Go out of the mountain and counted as an event
    }
            
    return POS_Aft;
}

double *NP_CDEX(double *POS_Int, double *DR, double Mx, double V, double Sigma_SI, int SorT)//Next_Point_In_Shielding(NPIS)
{
    static double Return_Value[8];
    
    double LFA=0.001;
    double V_aft=V;double ROELTA=0;
    int Collision_Time=0;

    int Times = Possion_GetRandom_Full(LFA);
    //double Segment = 1e3*Length_for_asking_the_collision(LFA,Mx,V_aft,Sigma_SI,1.81,Weighted_Atomic_Number);//Atomic Number Of Material(ANOM)
    double Segment = Length_for_asking_the_collision(LFA,Mx,V_aft,Sigma_SI,1.81,Weighted_Atomic_Number);//km, Atomic Number Of Material(ANOM)

    double Sprint_GO = Segment*Times;
    cout << "Sprint_GO: " << Sprint_GO << endl;
    double *SLU=DAC_to_LEG(POS_Int,DR,Sprint_GO);//Scaling_Length_Used,The layer a WIMP runs in
    cout << "Below(-1), In(0) or Out(1) of the mountain: " << SLU[3] << endl;//Below(-1), In(0) or Out(1) of the mountain

    if( int(SLU[3])==0 )
    {
        POS_Int = PAP(Sprint_GO,POS_Int,DR);
        double *VAC_end = VAC(Mx,Sigma_SI,V_aft,Weighted_Atomic_Number);//
        V_aft = VAC_end[0]; ROELTA = VAC_end[1];
        if(SorT==1)DR = DRAC(Weighted_Atomic_Number,Mx,V_aft,DR,ROELTA);
        Collision_Time=1;
        
    }

    if( int(SLU[3])==-1 )
    {
        double Scaling = (0-POS_Int[2])/DR[2];
        if(Sprint_GO<Scaling)
        {
            cout << "Sprint_GO<Scaling " << endl;
            POS_Int = PAP(Sprint_GO,POS_Int,DR);
            double *VAC_end = VAC(Mx,Sigma_SI,V_aft,Weighted_Atomic_Number);//
            V_aft = VAC_end[0]; ROELTA = VAC_end[1];
            if(SorT==1)DR = DRAC(Weighted_Atomic_Number,Mx,V_aft,DR,ROELTA);
            Collision_Time=1;
        }
        if(Sprint_GO>=Scaling)
        {
            cout << "Sprint_GO>=Scaling" << endl;
            POS_Int = PAP(Sprint_GO,POS_Int,DR);
        }
    }
    if( int(SLU[3])==1)
    {
        cout << "Out_of_the_boundary" << endl;
        POS_Int = PAP(Sprint_GO,POS_Int,DR);
    }


    Return_Value[0]=POS_Int[0];Return_Value[1]=POS_Int[1];Return_Value[2]=POS_Int[2];
    Return_Value[3]=DR[0];Return_Value[4]=DR[1];Return_Value[5]=DR[2];
    Return_Value[6]=V_aft;Return_Value[7]=Collision_Time;
    cout << "V_aft_Func: " << V_aft << endl;
     
    return Return_Value;
}

int Out_or_Not(double *POS)
{
    int Judge=0;
    double Length = Lcjpl_NonNL(POS[0],POS[1],POS[2]);//m
    double NF     = SqrtN2(POS);//m
    if(NF<Length){Judge=1;}
    else{Judge=0;}
    
    return Judge;
}

double *CDEX_Real_N_With_Angle(int SorT, double Sigma_SI, double V_Int, double Mx, double *DR, double *POS_Int) //Mx(Mass of WIMP),Velocity(km/s) Density(g/cm^3)
{//Straight_or_scattered(SorT)
    static double RETURN_VALUE[20];//Return the value back
    double Direction_VT[3]={0,0,0};//Used for the calculation
    double V_aft=V_Int;double LFA=0.001;//Lamda_for_Average(LFA)
    int Step=0;int Collision_Time=0;
    //===========================
    while( Out_or_Not(POS_Int)==1 and 0-POS_Int[2]<1e-10 and Energy_DM(Mx,V_aft*1e3/3e8)>=0.01)
    {
        cout << "V_aft: " << V_aft << endl;
        cout << "Energy_DM: " << Energy_DM(Mx,V_aft*1e3/3e8);
        cout << "DR[0]: " << DR[0] << "DR[1]: " << DR[1] << "DR[2]: " << DR[2] << endl;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        double *NP_1 = NP_CDEX(POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
        POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
        if(NP_1[7]>0)Collision_Time = Collision_Time + 1;
    }
    RETURN_VALUE[0]=V_aft;RETURN_VALUE[1]=Collision_Time;
    cout << "Collision_Time: " << Collision_Time << endl;
    cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;

    if(Out_or_Not(POS_Int)==0)
    {
        RETURN_VALUE[2]=1;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "Outside" << endl;
    }
    if(0-POS_Int[2]>=1e-10)
    {
        RETURN_VALUE[2]=0;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "Underground" << endl;
    }
    if(Energy_DM(Mx,V_aft*1e3/3e8)<0.01)
    {
        RETURN_VALUE[2]=0;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "Energy<threshold" << endl;
    }
    return RETURN_VALUE;
}
double Radius_To_Center_of_Earth(double *POS)
{
    return sqrt(POS[0]*POS[0]+POS[1]*POS[1]+POS[2]*POS[2]);
}
int In_or_Out_E(double *DR, double *POS)
{
    int In_or_Out=0;//1 for in and 0 for out
    double A=1;
    double B=2*(DR[0]*POS[0]+DR[1]*POS[1]+DR[2]*POS[2]);
    double C=POS[0]*POS[0]+POS[1]*POS[1]+POS[2]*POS[2]-6371*6371;
    
    if(B*B-4*A*C>0)In_or_Out=1;
    else In_or_Out=0;
    
    return In_or_Out;
}

double *DAC_to_LEG_30MWE(double *POS_Int, double *DR, double PL)//Direction_Aft_Collision_to_Length(PAC),Predicted_Length(PL)
{
    static double POS_Aft[5];
    //POS_Aft[3]=Count or not
    
    POS_Aft[0] = POS_Int[0]+PL*DR[0];
    POS_Aft[1] = POS_Int[1]+PL*DR[1];
    POS_Aft[2] = POS_Int[2]+PL*DR[2];
    
    double STU; //Boundary Or Not
    cout << "POS_Aft[0]: " << POS_Aft[0] << "POS_Aft[1]: " << POS_Aft[1] << "POS_Aft[2]: " << POS_Aft[2] << endl;
    cout << "Radius_To_Center_of_Earth(POS_Aft): " << Radius_To_Center_of_Earth(POS_Aft) << endl;
    cout << "Radius_To_Center_of_Earth(POS_Aft)-6371.: "   << Radius_To_Center_of_Earth(POS_Aft)-6371. << endl;
    cout << "Radius_To_Center_of_Earth(POS_Aft)-6371.01: " << Radius_To_Center_of_Earth(POS_Aft)-6371.01 << endl;

    int Check_on_Earth = In_or_Out_E(DR,POS_Int);
    if(Radius_To_Center_of_Earth(POS_Aft)-6371.>=1E-10 and Radius_To_Center_of_Earth(POS_Aft)-6371.01<=1E-10)
    {
        cout << "Yes1 " << endl;
        cout << "in the shielding^>^" << endl;
        POS_Aft[3]=1;//in the shielding
    }
    else if(Check_on_Earth==1)
    {
        cout << "Inside the earth>^ " << endl;
        cout << "Check_on_Earth: " << Check_on_Earth << endl;
        POS_Aft[3]=2;//Inside the mountain and this event has an interaction
    }
    else if(Check_on_Earth==0)
    {
        cout << "out of the shielding^>^ " << endl;
        cout << "Check_on_Earth: " << Check_on_Earth << endl;
        POS_Aft[3]=0;//Go out of the mountain and counted as an event
    }
            
    return POS_Aft;
}

double *NP_30MWE(double *POS_Int, double *DR, double Mx, double V, double Sigma_SI, int SorT)//Next_Point_In_Shielding(NPIS)
{
    static double Return_Value[8];
    
    double LFA=0.001;
    double V_aft=V;double ROELTA=0;
    int Collision_Time=0;

    int Times = Possion_GetRandom_Full(LFA);
    //double Segment = 1e3*Length_for_asking_the_collision(LFA,Mx,V_aft,Sigma_SI,1.81,Weighted_Atomic_Number);//Atomic Number Of Material(ANOM)
    double Segment = Length_for_asking_the_collision(LFA,Mx,V_aft,Sigma_SI,1.81,Weighted_Atomic_Number);//km, Atomic Number Of Material(ANOM)

    double Sprint_GO = Segment*Times;
    cout << "Sprint_GO: " << Sprint_GO << endl;
    double *SLU=DAC_to_LEG_30MWE(POS_Int,DR,Sprint_GO);//Scaling_Length_Used,The layer a WIMP runs in
    cout << "in(-1), Inside the earth(0) or Out of the shielding(1): " << SLU[3] << endl;//Below(-1), In(0) or Out(1) of the mountain

    if( int(SLU[3])==1 )//In the shielding
    {
        POS_Int = PAP(Sprint_GO,POS_Int,DR);
        double *VAC_end = VAC(Mx,Sigma_SI,V_aft,Weighted_Atomic_Number);//
        V_aft = VAC_end[0]; ROELTA = VAC_end[1];
        if(SorT==1)DR = DRAC(Weighted_Atomic_Number,Mx,V_aft,DR,ROELTA);
        Collision_Time=1;
        
    }

    if( int(SLU[3])==2 )//Inside the earth
    {
        POS_Int = PAP(Sprint_GO,POS_Int,DR);
    }
    if( int(SLU[3])==0)//Out of the shielding(Above the air)
    {
        cout << "Out_of_the_boundary" << endl;
        POS_Int = PAP(Sprint_GO,POS_Int,DR);
    }


    Return_Value[0]=POS_Int[0];Return_Value[1]=POS_Int[1];Return_Value[2]=POS_Int[2];
    Return_Value[3]=DR[0];Return_Value[4]=DR[1];Return_Value[5]=DR[2];
    Return_Value[6]=V_aft;Return_Value[7]=Collision_Time;
    cout << "V_aft_Func: " << V_aft << endl;
     
    return Return_Value;
}

double *KS_Real_N_With_Angle_30MWE(int SorT, double Sigma_SI, double V_Int, double Mx, double *DR) //Mx(Mass of WIMP),Velocity(km/s) Density(g/cm^3)
{//Straight_or_scattered(SorT)
    static double RETURN_VALUE[20];//Return the value back
    double Direction_VT[3]={0,0,0};//Used for the calculation
    double V_aft=V_Int;double LFA=0.001;//Lamda_for_Average(LFA)
    int Step=0;int Collision_Time=0;
    double POS_Int[3]={0,0,-(6371+0.01)};
    double R_i = 6371.+0.01;//Starting Point
    double R_f = 6371.;//Starting Point
    //===========================
    while( Radius_To_Center_of_Earth(POS_Int)-6371.>1e-10 and Radius_To_Center_of_Earth(POS_Int)-(6371.+0.01) <1e-10 and Energy_DM(Mx,V_aft*1e3/3e8)>=0.01)
    {
        cout << "V_aft: " << V_aft << endl;
        cout << "Energy_DM: " << Energy_DM(Mx,V_aft*1e3/3e8);
        cout << "DR[0]: " << DR[0] << "DR[1]: " << DR[1] << "DR[2]: " << DR[2] << endl;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        double *NP_1 = NP_30MWE(POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
        POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
        if(NP_1[7]>0)Collision_Time = Collision_Time + 1;
    }
    RETURN_VALUE[0]=V_aft;RETURN_VALUE[1]=Collision_Time;
    cout << "Collision_Time: " << Collision_Time << endl;
    cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;

    if(Radius_To_Center_of_Earth(POS_Int)- 6371.<=1e-10)
    {
        RETURN_VALUE[2]=1;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "Inside the earth" << endl;
    }
    if(Radius_To_Center_of_Earth(POS_Int)-(6371.+0.01)>=1e-10)
    {
        RETURN_VALUE[2]=0;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "Out of the shielding" << endl;
    }
    if(Energy_DM(Mx,V_aft*1e3/3e8)<0.01)
    {
        RETURN_VALUE[2]=0;
        cout << "POS_Int[0]: " << POS_Int[0] << "POS_Int[1]: " << POS_Int[1] << "POS_Int[2]: " << POS_Int[2] << endl;
        cout << "Energy<threshold" << endl;
    }
    return RETURN_VALUE;
}

/*
TF1 *ROLTML = new TF1("ROLTML","Lcjpl_NonNL([0]+x*[1],[2]+x*[3],[4]+x*[5])/SqrtN2_ABC([0]+x*[1],[2]+x*[3],[4]+x*[5])",0,PL);//Ratio_Of_Location_to_Mountain_Length(RLTML)
ROLTML->SetParameter(0,POS_Int[0]);
ROLTML->SetParameter(1,DR[0]);
ROLTML->SetParameter(2,POS_Int[1]);
ROLTML->SetParameter(3,DR[1]);
ROLTML->SetParameter(4,POS_Int[2]);
ROLTML->SetParameter(5,DR[2]);
STU= ROLTML->GetX(1,0,PL);//Scaling_To_Use(STU)
Check_Fun(POS_Int,DR,STU);
cout << "DAC_to_LEG->GetX(1,0,S); " << STU << endl;
 */
