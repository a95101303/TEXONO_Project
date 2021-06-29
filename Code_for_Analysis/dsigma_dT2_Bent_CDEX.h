//All values are expressed in m, need to be changed to km when using. Thanks 1e-2*1e-3=1e-5

//CDEX-Angle
double *CDEX_Real_N_With_Angle(int SorT, double Sigma_SI, double V_Int, double Mx, double *DR, double *POS_Int) //Mx(Mass of WIMP),Velocity(km/s) Density(g/cm^3)
{//Straight_or_scattered(SorT)
    static double RETURN_VALUE[20];//Return the value back
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
        cout << "IOS_Position[0]: " << IOS_Position[0] << endl;
        cout << "IOS_Position[1]: " << IOS_Position[1] << endl;
        cout << "IOS_Position[2]: " << IOS_Position[2] << endl;
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
            if(NP_1[7]==1)Collision_Time = Collision_Time + 1;
            Step = Step + 1;
        }
        else if(IOS_Position[2]==2 and DR_In_or_Out_on_surface(IOS_Position[1],DR)==0)
        {
            cout << "Inside the shielding2" << endl;
            double *NP_1 = NP(0,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
            if(NP_1[7]==1)Collision_Time = Collision_Time + 1;
            Step = Step + 1;
        }
        else if(IOS_Position[2]==2 and DR_In_or_Out_on_surface(IOS_Position[1],DR)==1)
        {
            cout << "Inside the shielding3" << endl;
            double *NP_1 = NP(1,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
            if(NP_1[7]==1)Collision_Time = Collision_Time + 1;
            Step = Step + 1;
        }
        else
        {
            cout << "Out the shielding" << endl;
            double *NP_1 = NP(1,POS_Int,DR,Mx,V_aft,Sigma_SI,SorT);
            POS_Int[0]=NP_1[0];POS_Int[1]=NP_1[1];POS_Int[2]=NP_1[2];DR[0]=NP_1[3];DR[1]=NP_1[4];DR[2]=NP_1[5];V_aft=NP_1[6];
            if(NP_1[7]==1)Collision_Time = Collision_Time + 1;
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
    RETURN_VALUE[0]=V_aft;RETURN_VALUE[1]=Collision_Time;
    cout << "Collision_Time: " << Collision_Time << endl;
    if((R_OW_KS-RTC(POS_Int))<1e-10)
    {
        RETURN_VALUE[2]=1;
        cout << "Radius>Outer_Radius" << endl;
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


    return RETURN_VALUE;
}


