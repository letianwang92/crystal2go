
#define PI 3.14159265

// Crystal Properties
#define Tm 1685.0 		// Melting Temperature of cSi
#define rhos 2328.0 // density

// Amorphous Properties
#define Ta 1460.0 //1348.0 // Melting T of aSi  , 1420 in other references
#define ka 1.806 // J/smK
#define rhoa 2200.0 // density Kg/m^3

// Oxide Properties // FO=1
#define rhox 2203.0
#define Cpx 703.0
//#define kx 1.38

// Metal Properties (Copper) // FO=2
#define rhom 8940.0
#define Cpm 385.0
#define km 401.0

// Metal Properties (Copper) // FO=99 , Pseudo values
#define rhov 1.0
#define Cpv 1.0
#define kv 1.0e-9

// Liquid Properties
#define Cpl 860.7 

// optical properties
// lambda(nm) = 1240/E(eV)
// 532 nm -> 2.33 eV
#define Eg 3.652  // direct bandgap of cSi, eV
#define Eg2 3.648 
//#define RF 0.34 // reflectance
//#define AB 1.02e6 // 1.037e7 (aSi) // 1.02e6 (cSi) // Absorption coefficient (1/m)

#define AB_aSi 7.0e6 // absorption coefficient of aSi 
#define RF_aSi 0.44 // reflectance of aSi , Webber et al APL 43(7) 1983

#define OptN_liq 3.0 // n of liquid Si
#define OptK_liq 4.84 // k of liquid Si 
//#define RF_liq 	// reflectivity of liquid silicon
//#define AB_liq 	// Absorption coefficient

