#include <math.h>
#include <stdlib.h>
#include "sConst1.h"
#include "fProp1.hpp"

// Thermal properties of silicon
double rhol(double T) {
	// Density of Liquid
	return (2533.0-0.152*(T-Tm)) ;
}

double Cpcr(double T) {
	// Cp of crystal
	return ((0.844+1.18e-4*T-1.55e4*pow(T,-2))*1.0e3) ;
}

double Cpa(double T) {
	// Cp of amorphous
	return (Cpcr(T)+4200.0*(-0.00191+0.0409*T/1685.0)) ;
}

double kcr(double T) {
	// heat conductivity of crystal
  return (1*((1.521e5/pow(T,1.226))+((8.97e2/sqrt(T))-(1.521e5/pow(T,1.226)))*(T>1200.0))) ;
  //return 3;
}

double kl(double T) {
	// heat conductivity of liquid
	return (0.010*(T-Tm)+56.0) ;
}



// Thermal properites of silicon oxide
double kx(double T) {
	// heat conductivity of SiOx (W/mk)
	return (4.639e-10 *(T*T*T) - 2.308e-6*(T*T) + 3.9063e-3*T); 
}

//FO = 0 : Silicon (can change phase)
//FO = 1 : oxide capping (SiOx)
//FO = 2 : Metal
//FO = 3 : Virtual amorphous ( phase change like crystal, heat transfer like aSi
//		Need to change this to value -1.
//FO = 99 : Vacuum ( zero thermal conductivity, zero mass)
		
double fRho(double Fa, double Fc, double T, int Fo) {
	// density
	if (Fo == 1) {
		return rhox;}
	else if (Fo == 2) {
		return rhom;}
	else if (Fo == 3) {
		return (Fa*rhoa+Fc*rhoa+(1.0-Fa-Fc)*rhol(T)); }
	else if (Fo == 99) {
		return rhov; }
	else {		
		return (Fa*rhoa+Fc*rhos+(1.0-Fa-Fc)*rhol(T)) ;
	}
}

double fCp(double Fa, double Fc, double T, int Fo) {
	if (Fo == 1) {
		return Cpx; }
	else if (Fo == 2) {
		return Cpm; }
	else if (Fo == 3) {
		return (Fa*Cpa(T)+Fc*Cpa(T)+(1.0-Fa-Fc)*Cpl); }
	else if (Fo == 99) {
		return Cpv; }
	else {
		return (Fa*Cpa(T)+Fc*Cpcr(T)+(1.0-Fa-Fc)*Cpl) ; 
	}
}

double fK(double Fa, double Fc, double T, int Fo) {
	// Thermal conductivity
	if (Fo == 1) {
		return kx(T); }
	else if (Fo == 2) {
		return km; }
	else if (Fo == 3) {
		return (Fa*ka+Fc*ka+(1.0-Fa-Fc)*kl(T)); }
	else if (Fo == 99) {
		return kv; }
	else {
		return (Fa*ka+Fc*kcr(T)+(1.0-Fa-Fc)*kl(T)) ;	
	}	
}

