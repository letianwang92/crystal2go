
#define PI (3.14159265)
#define kb (1.38e-23) // Boltzmann Const J/K
#define kbe (8.62e-5)  // Boltzmann Const eV/K
#define Evd 1.08 // eV


//crystal
#define Tm 1685.0
#define rhos 2328.0 // density 
#define Lca 8.345e-20  // latent heat per a atom 
#define Lc 1.787e6  //latent heat  J/kg
#define st_cl 0.34 //0.34 //////////////////// Surface tension N/m, crystal-liquid
#define fCon 0.3341325 // 77deg, (2-3*cos(a)+(cos(a))^3)/4  // Contact angle cSi on SiOx , Im et al.

#define st_al 0.18 // amorphous-liquid
#define st_ac 0.06 // amorphous - crystal

//Amorphous  /////// Especially for SPC process, steady state parameters are used but needed to be reconsidered later.
#define Ta 1460.0 //1348.0
#define rhoa 2200.0 // density Kg/m^3
#define Laa 6.174e-20  // latent heat per an atom 
#define La 1.323e6  //latent heat J/kg
#define Ena 5.3 // Activation Energy for SPC-nucleation, eV
#define Ina 1.7e44 // Prefactor for SPC-nucleation, s^-1 m^-3
#define Ega 3.1 // Activation energy for SPC-growth, eV
#define Iga 2.1e7 // Prefactor for SPC-growth, m/s
