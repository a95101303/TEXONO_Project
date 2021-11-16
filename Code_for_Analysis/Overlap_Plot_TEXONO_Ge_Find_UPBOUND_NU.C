#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
//#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "velocity_distribution_2000_Ave_ER.h"

#include "cpkkd_calculation_New.h"

double DM_E(double Mx, double v)//Mx(GeV/c^2),v(km/s)
{
    return 0.5*(Mx*1e6)*(v*1e3/3e8)*(v*1e3/3e8);
}
void Overlap_Plot_TEXONO_Ge_Find_UPBOUND_NU()
{
    vector<double> WIMP_mx_Array ={2.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.12,0.09,0.07,0.05};
    //Read the file of DCS for different masses
    vector<double> Collision_Time_Array;
    vector<double> Energy_Loss_Per_Collision;
    vector<double> Cross_Section_Set;
    
    for(int kkk=0; kkk<WIMP_mx_Array.size(); kkk++)
    {
        double Energy_Loss_all_Event=0;
        double Event_Number=0;
        cout << "WIMP_mx_Array[kkk]: " << WIMP_mx_Array[kkk] << endl;
        while(Event_Number<50)
        {
            cout << "Event_Number: " << Event_Number << endl;
            double *NU_Collision = Velocity_Aft_collision(1000,WIMP_mx_Array[kkk],1e-31,779,1);
            Energy_Loss_all_Event = Energy_Loss_all_Event + NU_Collision[4];
            Event_Number = Event_Number + 1;
        }
        cout << "Energy_DM(WIMP_mx_Array[kkk],779*1e3/3e8): " << Energy_DM(WIMP_mx_Array[kkk],779*1e3/3e8) << endl;
        double Ener_Diff = Energy_DM(WIMP_mx_Array[kkk],779*1e3/3e8)-0.01;
        Energy_Loss_all_Event = (1./(Ener_Diff))*(Energy_Loss_all_Event/Event_Number);
        cout << "Energy_Loss_all_Event: " << Energy_Loss_all_Event << endl;
        Energy_Loss_Per_Collision.push_back(Energy_Loss_all_Event);
    }
    
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< WIMP_mx_Array.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Energy_Loss_Per_Collision[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;

    /*
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Cross_Section_Set[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Collision_Time_Array[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
    cout << "========================================================" << endl;
    for(int kkk=0; kkk< Cross_Section_Set.size(); kkk++)
    {
        cout << WIMP_mx_Array[kkk] << "," << Energy_Loss_Per_Collision[kkk] << "," << endl;
    }
    cout << "========================================================" << endl;
     */
    //==================================================//

}//End_Main
    
     


