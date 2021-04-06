#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "Hist_SetLimit_Plot_v2_Extract_Peak.h"
#include "B_L_Henke_data_PE_f1_f2.h"
#include "dsigma_dT2.h"
#include "velocity_distribution_2000_Ave.h"
#include "cpkkd_calculation_New.h"

void PRL_Check_For_Heavy_DM_Died()//CDEX:Threshold=160eV, TEXONO:Threshold=200eV
{
    
    double DM_Mass=1e10; double Sigma_SI=5e-23;
    //1km Rock with no shielding, highest speed with the upper bound
    double Cross_Section = total_C_AAAA(1,779,Sigma_SI,DM_Mass,34);
    double Cross_Section1 = total_C_AAAA(1,779,Sigma_SI,DM_Mass,14);
    //Counts through the earth ( Air is so small )
    double COUNT_EARTH=1.8*0.01*1e5/(unified_atomic_mass_g*(Weighted_Atomic_Number)) * Cross_Section ;
    double COUNT_AIR=0.4*1e-3*40*1e5/(unified_atomic_mass_g*(14)) * Cross_Section1 ;
    cout << "COUNT_EARTH: " << COUNT_EARTH << endl;
    cout << "COUNT_AIR: "   << COUNT_AIR << endl;
    //Mean free Path = N atom amount*total Cross-section
    //double Mean_Free_Path_MAJORANA = N_atom_Ge_1kg*(29.7) * total_C_AAAA(1,750,1e-27,1e10,72.64)*1e-9;
    double Mean_Free_Path_TEXONO = 1/( Number_density_Ge * total_C_AAAA(1,180,1e-32,DM_Mass,72.64)) ;
    cout << "Mean_Free_Path_TEXONO: " << Mean_Free_Path_TEXONO << endl;
    //double Mean_Free_Path_XENON1T = N_atom_Xe_1kg*(1.5) * total_C_AAAA(1,750,1e-27,1e10,72.64);
    //double Mean_Free_Path_XENON1T =  1/( Number_density_Xe * total_C_AAAA(1,220,2e-29,1e11,131.64)) ;//V OK!
    //cout << "Mean_Free_Path_XENON1T: " << Mean_Free_Path_XENON1T << endl;
    double Mean_Free_Path_XENON1T_3 =  1/( Number_density_Xe * total_C_AAAA(1,220,2e-29,1e11,131.64)) ;//V OK!
    cout << "Mean_Free_Path_XENON1T_3: " << Mean_Free_Path_XENON1T_3 << endl;
    //Energy_DM(1e10,799*1e3/3e8);
    //max_recoil_A_keV(1e10,799*1e3/3e8,34)
    double Total_Collision_Time=0;
    
    for(int Initial=0; Initial<779; Initial++)
    {
        double Energy_Diff = Energy_DM(DM_Mass,(779-Initial)*1e3/3e8)-Energy_DM(DM_Mass,(778-Initial)*1e3/3e8);
        //cout << "Energy_Diff: " << Energy_Diff << endl;
        double Recoil_Energy = max_recoil_A_keV(DM_Mass,(779-Initial)*1e3/3e8,14);
        double Count_Real = Energy_Diff/(Recoil_Energy*0.5);
        Total_Collision_Time = Total_Collision_Time + Count_Real;
        /*
        if(Total_Collision_Time> COUNT_EARTH+COUNT_AIR){cout << "Final Velocity: " << (779-Initial) << endl;
            break;
         
        }*/
    }
     
    cout << "Total_Collision_Time: " << Total_Collision_Time << endl;
    
                                                                                                                            
}
    

