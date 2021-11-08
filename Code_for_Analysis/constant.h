#include "dm_parameters.h"
//Planck Constant
const double h_eV_s=4.135667662*1e-15;//(eV*s)
const double h_J_s =6.626069934*1e-34;//(J*s)
// Unit convert
const double Fine_Structure=1.0/137.0;
// 1 cm^-1 = 1.97327053323587611e-11 MeV
const double cm1_MeV = 1.97327053323587611e-11;
// 1 fm^-1 = 1.97327053323587611e2 MeV
const double fm1_MeV = cm1_MeV*1e13;
// 1 cm = 1.0/cm1_MeV (MeV^-1)
const double cm_MeV1 = 1.0/cm1_MeV;
// 1 fm = 1.0/fm1_MeV (MeV^-1)
const double fm_MeV1 = 1.0/fm1_MeV;

// a barns = 10^-24 cm^2
// 1 barns = 10^-28 m^2 = 1e-28*((100.0*(1.0/cm1_MeV))^2)
const double barn_cm = 1e-24; 
const double barn_MeV1 = 1.0e-28*pow((100.0*cm_MeV1),2);

const double PI = 3.14159265358979323846;

// Physics constants
const double kms1_to_cmday1 = 100000.0*86400.0; // km/s to cm/day
const double kms1_to_c = 1000.0/2.99792458e8; // km/s to natural unit
const double MeV1_to_keV1 = 1.0/1000.0; // per MeV to per keV
//
const double Me = 0.51099906; // MeV
const double Mp = 938.272046; // MeV
const double Alpha = 7.29735307964481886e-03;
const double GF = 1.16639e-11; // MeV^-2
const double sin2_theta_w = 0.23867;
const double gV = 2.0*sin2_theta_w+0.5;
const double gA = -0.5; // for anti-nu, =0.5 for nu

// Classical electron radius (m)
const double r_e_m  = 2.8179402894e-15;
// Classical electron radius (cm)
const double r_e_cm = 2.8179402894e-13;
// Classical electron radius (MeV^-1)
const double r_e = r_e_m*100.0*cm_MeV1;

// Initial parameter for R2 and MM
const double MM = 1e-10;
const double R2 = 1e-32/pow(cm1_MeV,2);

const double ASi = 28.0855;
const double AGe = 72.64;
const double AXe = 131.293 ;
const double AAL2O3 = 19.4;
// const double AGe = 72.64;
// const double AGe = 72.59;
const double ZGe = 32;
const double ZXe = 54;
const double NGe = AGe-ZGe;
const double NXe = AXe-ZXe;

const double unified_atomic_mass_g = 1.66053873e-24; // g
const double unified_atomic_mass_MeV = 931.494013; // MeV
const double unified_atomic_mass_GeV = 0.931494013; // GeV

const double atom_mass_Ge_g = unified_atomic_mass_g*AGe; // g
const double atom_mass_Xe_g = unified_atomic_mass_g*AXe; // g
const double atom_mass_Si_g = unified_atomic_mass_g*ASi; // g

const double Si_Number_density = 2.33/(atom_mass_Si_g);
const double atom_mass_Ge_MeV = unified_atomic_mass_MeV*AGe; // MeV
const double MNGe = atom_mass_Ge_MeV;
// number of atom in 1kg of Ge
const double N_atom_Ge_1kg = 1000.0/atom_mass_Ge_g;
const double N_atom_Xe_1kg = 1000.0/atom_mass_Xe_g;

const double avogadro_number = 6.02214129e+23;
const double N_atom_1kg_Ge_Electron = (1000.0*avogadro_number)/AGe;

const double DM_density_ER = 0.4; //GeV cm^{-3}
const double Electron_number_ER = N_atom_1kg_Ge_Electron; //kg^{-1}

// mass_density in g/cm^3
const double mass_density_Ge = 5.323;
const double mass_density_Xe = 5.761*1e-3;
// number of atom per cm^3
const double Number_density_Ge = mass_density_Ge/atom_mass_Ge_g;
const double Number_density_Xe = mass_density_Xe/atom_mass_Xe_g;

// number of electron per cm^3
const double electron_density_Ge = ZGe*mass_density_Ge/atom_mass_Ge_g;

// cm^2/g = 0.00829393 barns/atom in Ge (Henke)
const double cm2g_barn = barn_cm/(unified_atomic_mass_g*AGe);
// const double cm2g_barn = 0.00829393; // Henke
// const double cm2g_barn = 0.008297; // Storm

// barns/atom = 120.57 cm^2/g in Ge (Henke)
const double barn_cm2g = 1.0/cm2g_barn;
// const double barn_cm2g = 120.57; // Henke
// const double barn_cm2g = 120.5249114; // Storm

// cross-section(cm^2/g) = f2_factor*f2/E(keV)
// f2_factor = 579.51; // Henke

const double NORMALIZE_kg_keV_day = N_atom_Ge_1kg*0.001*86400;

// Unit convert
//const double unified_atomic_mass_g = 1.66053873e-24; // g
//const double unified_atomic_mass_MeV = 931.494013; // MeV
//ATM
const double AO = 15.99 ;
const double AN = 14.0;
//Cement
const double AFe = 55.854;
const double ACa = 40.078;
const double AAl = 26.982;
const double AH2O = 18;
//Shielding
const double APb = 207.2;
const double AB  = 10.811;
const double ACu = 63.546;
//*const double AFe = 55.854;
//
const double atom_mass_Fe_g = unified_atomic_mass_g*AFe; // g
const double atom_mass_O_g = unified_atomic_mass_g*AO; // g
const double atom_mass_N_g = unified_atomic_mass_g*AN; // g


// mass_density in g/cm^3
const double mass_density_Fe = 7.874 ;
const double mass_density_O = 1.429/1000;
const double mass_density_Si = 1.141;

// Percentage of the material of the earth
const double Fe_Percentage = 0.321/(0.321+0.301+0.151);
const double O_Percentage  = 0.301/(0.321+0.301+0.151);
const double Si_Percentage = 0.151/(0.321+0.301+0.151);

//Percentage of the material of ATM
const double O_Percentage_ATM = 0.2;
const double N_Percentage_ATM = 0.8;

//Percentage of the material of Cement
const double H2O_Percentage=0.25;
const double Others=0.75;
const double Density_of_Cement=2.8;

const double shield_A = 16.0*0.6+28.0*0.4; // average of O and Si
const double shield_density = 2.6; // g cm^-3
const double shield_atom_density = shield_density/(unified_atomic_mass_g*shield_A); // number of atom per cm^3
const double shield_depth = 2.4*1000.0*100.0; // 2.4 km in unit: cm
// 6,371

const double Earth_Radius = 6371; //km

///////////////////////////////////////////////////
// atmosphere
///////////////////////////////////////////////////
//atm_table[i][0] : hight above sea level, m
//atm_table[i][4] : density, kg/m^3
double atm_table[20][6] = {
    0, 15.00, 9.807, 10.13, 1.225, 1.789,
    1000, 8.50, 9.804, 8.988, 1.112, 1.758,
    2000, 2.00, 9.801, 7.950, 1.007, 1.726,
    3000, -4.49, 9.797, 7.012, 0.9093, 1.694,
    4000, -10.98, 9.794, 6.166, 0.8194, 1.661,
    5000, -17.47, 9.791, 5.405, 0.7364, 1.628,
    6000, -23.96, 9.788, 4.722, 0.6601, 1.595,
    7000, -30.45, 9.785, 4.111, 0.5900, 1.561,
    8000, -36.94, 9.782, 3.565, 0.5258, 1.527,
    9000, -43.42, 9.779, 3.080, 0.4671, 1.493,
    10000, -49.90, 9.776, 2.650, 0.4135, 1.458,
    15000, -56.50, 9.761, 1.211, 0.1948, 1.422,
    20000, -56.50, 9.745, 0.5529, 0.08891, 1.422,
    25000, -51.60, 9.730, 0.2549, 0.04008, 1.448,
    30000, -46.64, 9.715, 0.1197, 0.01841, 1.475,
    40000, -22.80, 9.684, 0.0287, 0.003996, 1.601,
    50000, -2.5, 9.654, 0.007978, 0.001027, 1.70,
    60000, -26.13, 9.624, 0.002196, 0.0003097, 1.584,
    70000, -53.57, 9.594, 0.00052, 0.00008283, 1.438,
    80000, -74.51, 9.564, 0.00011, 0.00001846, 1.321
};
///////////////////////////////////////////////////
// earth
///////////////////////////////////////////////////
//earth_table[i][0] : distance from center, km
//earth_table[i][5] : average density from this layer, calculated from earth_table[i][1-4], g/cm^3
/*
double earth_table[14][6] = {
    0     , 13.090, 0.0,  0.0   , 0.0, 0.0,
    1221.5, 12.760, 0.0, -8.8381, 0.0, 12.893569,
    3480.0, 9.9, -1.2638, -3.6426, -5.5281, 10.900696,
    3630.0, 5.51, -6.4761, 5.5283, -3.0807, 5.528405,
    5600.0, 7.9565, -6.4761, 5.5283, -3.0807, 4.912992,
    5701.0, 4.44, -6.4761, 5.5283, -3.0807, 4.411899,
    5771.0, 5.3197, -1.4836, 0.0, 0.0, 3.983938,
    5971.0, 11.2494, -8.0298, 0.0, 0.0, 3.848353,
    6151.0, 7.1089, -3.8045, 0.0, 0.0, 3.488987,
    6291.0, 2.6910, 0.6924, 0.0, 0.0, 3.367155,
    6346.6, 2.6910, 0.6924, 0.0, 0.0, 3.377736,
    6356.0, 2.9, 0.0, 0.0, 0.0, 2.900000,
    6368.0, 2.6, 0.0, 0.0, 0.0, 2.600000,
    6371.0, 1.02, 0.0, 0.0, 0.0, 1.020000
}*/

    
double earth_table[14][6] = {
    0     , 13.0885, 0.0,  0.0   , 0.0, 0.0,
    1221.5, 12.760, 0.0, -8.8381, 0.0, 12.893569,
    3480.0, 9.9, -1.2638, -3.6426, -5.5281, 10.900696,
    3630.0, 5.50, -6.4761, 5.5283, -3.0807, 5.528405,
    5600.0, 4.44, -6.4761, 5.5283, -3.0807, 4.912992,
    5701.0, 4.38, -6.4761, 5.5283, -3.0807, 4.411899,
    5771.0, 3.97, -1.4836, 0.0, 0.0, 3.983938,
    5971.0, 3.72, -8.0298, 0.0, 0.0, 3.848353,
    6151.0, 3.43, -3.8045, 0.0, 0.0, 3.488987,
    6291.0, 3.37, 0.6924, 0.0, 0.0, 3.367155,
    6346.6, 3.38, 0.6924, 0.0, 0.0, 3.377736,
    6356.0, 2.9, 0.0, 0.0, 0.0, 2.900000,
    6368.0, 2.6, 0.0, 0.0, 0.0, 2.600000,
    6371.0, 1.02, 0.0, 0.0, 0.0, 1.020000
};
double eV_to_keV = 1e-3;
//Energy Level possibility Constant(Integral)
double N_1_System[1]={1.8e-5};//1s
double N_2_System[2]={1.3e-4,7.3e-4};//2s,2p
double N_3_System[3]={5.5e-4,2.4e-3,2.8e-2};//3s,3p,3d
double N_4_System[2]={6.1e-4,2.6e-2};//4s,4p
//Energy Level Constant(eV)
double N_1_Energy[1]={1.1e+4*eV_to_keV};//1s
double N_2_Energy[2]={1.4e+3*eV_to_keV,1.2e+3*eV_to_keV};//2s,2p
double N_3_Energy[3]={1.7e+2*eV_to_keV,1.2e+2*eV_to_keV,3.5e+1*eV_to_keV};//3s,3p,3d
double N_4_Energy[2]={1.5e+1*eV_to_keV,6.5e+0*eV_to_keV};//4s,4p

