/* Revision Note
Aug 31 2012, optical absorption of aSi
Sept 4 2012, correct absorption & transmission functions.
 
 
 
*/

#include "mex.h"
#include "matrix.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include "sConst1.h"
#include "hNode.hpp"

double OptN_cSi(double lambda, double T)
{
    // Temperature dependent refractive index n
	// lambda(nm) = 1240/E(eV)
	double E;
	E = 1240.0/lambda; 
	T = T - 273.0 ; // Celcius
	return ( sqrt(4.386-0.00343*T+(99.14+0.062*T)/(Eg*Eg-E*E)) );
}

double OptK_cSi(double lambda, double T)
{
	// Temperature dependent refractive index n
	// lambda(nm) = 1240/E(eV)
	double E;
	E = 1240.0/lambda; 
	T = T - 273.0 ; // Celcius
	
	double K0, T0;
	K0 = -0.0805 + exp(-3.1893+7.946/(Eg2*Eg2-E*E)); 
	T0 = 369.9 - exp(-12.92+5.509*E);
	
	return (K0*exp(T/T0));
}

double fAbsorb(double lambda, double k)
{	// unit of lambda: nm
	return (4.0*PI*k/(lambda*1.0e-9)); 
}

double fReflect(double n, double k)
{
	return ( ((n-1.0)*(n-1.0)+k*k)/((n+1.0)*(n+1.0)+k*k) );
}

void Laser_surf(double lambda, double W, double heating, hNode *pnode)
{
	double rf, ab ;
	// double lambda = 532.0 ; // Wavelength of the laser
	//double absorb;
	double absorb_c, absorb_l, absorb_a;
	double rf_cSi, ab_cSi; 
	double rf_lSi, ab_lSi; 
	double rf_aSi, ab_aSi; 
	
	if ((W>0.0) | (heating>0.0) ){
	
		switch (pnode->m)
		{
			case 0: // Silicon
				if (lambda == 248.0)	{ // properties for excimer laser
					rf_cSi = 0.66;
					ab_cSi = 1.8e8;
					rf_lSi = 0.685;	
					ab_lSi = 1.56e8;
					rf_aSi = 0.55; 
					ab_aSi = 1.5e8; 
				}
				else {
					rf_cSi = fReflect(OptN_cSi(lambda, pnode->T), OptK_cSi(lambda, pnode->T)) ;
					ab_cSi = fAbsorb(lambda, OptK_cSi(lambda, pnode->T)) ;
					rf_lSi = fReflect(OptN_liq, OptK_liq) ;
					ab_lSi = fAbsorb(lambda, OptK_liq) ;
					rf_aSi = RF_aSi ; 
					ab_aSi = AB_aSi ; 
					//mexPrintf("ref= %e\n", rf_cSi);
				}

				
				absorb_c = (1.0-rf_cSi)*ab_cSi*exp(ab_cSi*(-(pnode->dz)));
				absorb_l = (1.0-rf_lSi)*ab_lSi*exp(ab_lSi*(-(pnode->dz)));
				absorb_a = (1.0-rf_aSi)*ab_aSi*exp(ab_aSi*(-(pnode->dz))); 
				pnode->Q_laser_1 = W*((pnode->c)*absorb_c+(pnode->a)*absorb_a+(1.0-(pnode->c)-(pnode->a))*absorb_l); // W/m^3
				pnode->Tr_laser_1 = W*((pnode->c)*(1.0-rf_cSi)+(pnode->a)*(1.0-rf_aSi)+(1-(pnode->c)-(pnode->a))*(1.0-rf_lSi))
					-(pnode->Q_laser_1)*(pnode->dz); // W/m^2 
				pnode->Q_laser_2 = heating*((pnode->c)*absorb_c + (pnode->a)*absorb_a + (1.0 - (pnode->c) - (pnode->a))*absorb_l); // W/m^3
				pnode->Tr_laser_2 = heating*((pnode->c)*(1.0 - rf_cSi) + (pnode->a)*(1.0 - rf_aSi) + (1 - (pnode->c) - (pnode->a))*(1.0 - rf_lSi))
					- (pnode->Q_laser_2)*(pnode->dz); // W/m^2 
				pnode->IsReflected = true;
			break;
			
			case 1: // Silicon Oxide
				pnode->Q_laser_1 = 0.0;
				pnode->Tr_laser_1 = W;
				pnode->Q_laser_2 = 0.0;
				pnode->Tr_laser_2 = heating;
				//pnode->IsReflected = false;
			break; 
			
			case 2: // Metal
				////
				pnode->IsReflected = true;
			break;
			
			case 3: // Virtual aSi
				////
				pnode->IsReflected = true;
			break;
			
			case 99: // Vacuum
				pnode->Q_laser_1 = 0.0;
				pnode->Tr_laser_1 = W;
				pnode->Q_laser_2 = 0.0;
				pnode->Tr_laser_2 = heating;
			break; 
		}
		
	}
	else 
	{
		pnode->Q_laser_1 = 0.0; 
		pnode->Tr_laser_1 = 0.0; 
		pnode->Q_laser_2 = 0.0;
		pnode->Tr_laser_2 = 0.0;
	}
	
	return ;
}

void Laser_bulk(double lambda, double W, double heating, hNode *pnode)
{	
	double ab ;
	// double lambda = 532.0 ; // Wavelength of the laser
	//double absorb;
	double absorb_c, absorb_l, absorb_a;
	double ab_cSi; 
	double ab_lSi; 
	double ab_aSi; 
	
	if ((W>0.0)|(heating>0.0)) {
		
		switch (pnode->m) 
		{
			case 0: // Silicon
				if (lambda == 248.0)	{
					ab_cSi = 1.8e8;
					ab_lSi = 1.56e8;
					ab_aSi = 1.5e8;
				}
				else {
					ab_cSi = fAbsorb(lambda, OptK_cSi(lambda, pnode->T)) ;
					ab_lSi = fAbsorb(lambda, OptK_liq) ;
					ab_aSi = AB_aSi; 
				}
				
				
				absorb_c = ab_cSi*exp(ab_cSi*(-(pnode->dz)));
				absorb_l = ab_lSi*exp(ab_lSi*(-(pnode->dz)));
				absorb_a = ab_aSi*exp(ab_aSi*(-(pnode->dz)));
				pnode->Q_laser_1 = W*((pnode->c)*absorb_c+(pnode->a)*absorb_a+(1.0-(pnode->c)-(pnode->a))*absorb_l); // W/m^3
				pnode->Tr_laser_1 = W-(pnode->Q_laser_1)*(pnode->dz); // W/m^2 
				pnode->Q_laser_2 = heating*((pnode->c)*absorb_c + (pnode->a)*absorb_a + (1.0 - (pnode->c) - (pnode->a))*absorb_l); // W/m^3
				pnode->Tr_laser_2 = heating - (pnode->Q_laser_2)*(pnode->dz); // W/m^2 
			break;
			
			case 1: // Silicon Oxide
				pnode->Q_laser_1 = 0.0;
				pnode->Tr_laser_1 = W;
				pnode->Q_laser_2 = 0.0;
				pnode->Tr_laser_2 = heating;
			break; 
			
			case 2: // Metal
			
			break;
			
			case 3: // Virtual aSi
			
			break;
			
			case 99: // Vacuum
				pnode->Q_laser_1 = 0.0;
				pnode->Tr_laser_1 = W;
				pnode->Q_laser_2 = 0.0;
				pnode->Tr_laser_2 = heating;
			break; 
		}
					
	}
	else 
	{
		pnode->Q_laser_1 = 0.0;
		pnode->Tr_laser_1 = 0.0;
		pnode->Q_laser_2 = 0.0;
		pnode->Tr_laser_2 = 0.0;
	}
	
	return ;
}


	