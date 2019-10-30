#include "mex.h"
#include "matrix.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>
#include "sConst2.h"
#include "cNode.hpp"

int fUpdateType1(cNode fc222, cNode fc122, cNode fc322,
        cNode fc212,cNode fc232,cNode fc221,cNode fc223, double tol)
{
    //const double tol = 1e-5; 
    if (((fc222.c+fc222.a)>tol)&&((fc222.c+fc222.a)<(1.0-tol))) {
        return 1;
    }
    else{
        return (int) ((fc122.c+fc122.a)<tol) || ((fc322.c+fc322.a)<tol) || ((fc212.c+fc212.a)<tol) || 
                ((fc232.c+fc232.a)<tol) || ((fc221.c+fc221.a)<tol) || ((fc223.c+fc223.a)<tol) ;
    }
}


// Nucleation functions

// Crystalline Si

double Grc(double T)	// Growth rate 
{
    return (1.8e5*exp(-Evd/kbe/T)*(1-exp(-Lca/(kb*T)*(Tm-T)/Tm))) ;
}  	
	
double Vgm_c(double T)  // melting speed of crystalline silicon
{
    return fabs(1.8e5*exp(-Evd/kbe/T)*(1-exp(-Lca/(kb*T)*(Tm-T)/Tm))) ;
}  

double D(double T)
{
    return (2.3e-8*exp(0.4/1800.0/kbe)*exp(-0.4/kbe/T));
}  

double I2(double T)
{
    const double vl = 1.8423e-29; // atomic volume in liquid phase
    return (8*D(T)/pow(vl,4.0/3.0)*pow(PI/6.0,1.0/6.0)*sqrt(st_cl/(6*PI*kb*T)));
}  

double Gc(double T)
{
    return (16*PI*pow(st_cl,3.0)*pow(Tm,2.0)/(3*pow(Lc,2.0)*pow(rhos,2.0)*pow((T-Tm),2.0)));
}  

double Nrc(double T)	// Homogeneous nucleation rate
{
    return (1.0e39*exp(-Gc(T)/(kb*T))); //(I2(T)*exp(-Gc(T)/(kb*T)));
}  

double Nrc_h(double T)
{ // Heterogeneous nucleation rate
    return (1.0e27*exp(-Gc(T)*fCon/(kb*T))); 
} 

double Rc(double T)	// Critical radius for nucleation
{
    return (2.0*st_cl*Tm/(Lc*rhos*fabs(T-Tm)));
} 

double Nrc_hEx(double T)
{ // Heterogeneous nucleation rate on aSi surf
	double fConEx, fAngEx ;
	fAngEx = (st_al-st_ac)/st_cl ;
	fConEx = (2.0-3.0*(fAngEx)+pow((fAngEx),3.0))/4.0; 
    return (1.0e32*exp(-Gc(T)*fConEx/(kb*T)));  //(1.0e39*Rc(T)*
}  

// aSi nucleation

double Vgm_a(double T)  // melting speed of amorphous silicon
{
    return fabs(1.1e6*exp(-Evd/kbe/T)*(1-exp(-Laa/(kb*T)*(Ta-T)/Ta))) ;
} 

double Nra(double T)
{
    return (Ina*exp(-Ena/(kbe*T)));
} 

double Gra(double T)
{
    return (Iga*exp(-Ega/(kbe*T)));
} 


//

void getID26(int *pID, int *pTY, double *psigma,
cNode fc122, cNode fc322, cNode fc212, cNode fc232,
cNode fc221, cNode fc223, cNode fc112, cNode fc132, cNode fc312, 
cNode fc332, cNode fc121, cNode fc123, cNode fc321, cNode fc323, 
cNode fc211, cNode fc213, cNode fc231, cNode fc233, cNode fc111, 
cNode fc131, cNode fc311, cNode fc331, cNode fc113, cNode fc133, 
cNode fc313, cNode fc333)
{
	double tol = 1.0-1.0e-5; 
	*pID = 0;
	*pTY = 1;
	*psigma = 0.0;
	
	*pID=((int)(((fc122.c+fc122.a)>tol)&&(fc122.c>0.0)))*fc122.id;
    if (*pID != 0) { *pTY=4; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc322.c+fc322.a)>tol)&&(fc322.c>0.0)))*fc322.id;
    if (*pID != 0) { *pTY=4; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc212.c+fc212.a)>tol)&&(fc212.c>0.0)))*fc212.id;
    if (*pID != 0) { *pTY=4; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc232.c+fc232.a)>tol)&&(fc232.c>0.0)))*fc232.id;
    if (*pID != 0) { *pTY=4; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc221.c+fc221.a)>tol)&&(fc221.c>0.0)))*fc221.id;
    if (*pID != 0) { *pTY=4; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc223.c+fc223.a)>tol)&&(fc223.c>0.0)))*fc223.id;
    if (*pID != 0) { *pTY=4; *psigma = fc122.sigma; return; }
    
    *pID=((int)(((fc112.c+fc112.a)>tol)&&(fc112.c>0.0)))*fc112.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc132.c+fc132.a)>tol)&&(fc132.c>0.0)))*fc132.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc312.c+fc312.a)>tol)&&(fc312.c>0.0)))*fc312.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc332.c+fc332.a)>tol)&&(fc332.c>0.0)))*fc332.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc121.c+fc121.a)>tol)&&(fc121.c>0.0)))*fc121.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc123.c+fc123.a)>tol)&&(fc123.c>0.0)))*fc123.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc321.c+fc321.a)>tol)&&(fc321.c>0.0)))*fc321.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc323.c+fc323.a)>tol)&&(fc323.c>0.0)))*fc323.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc211.c+fc211.a)>tol)&&(fc211.c>0.0)))*fc211.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc213.c+fc213.a)>tol)&&(fc213.c>0.0)))*fc213.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; } 
    *pID=((int)(((fc231.c+fc231.a)>tol)&&(fc231.c>0.0)))*fc231.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc233.c+fc233.a)>tol)&&(fc233.c>0.0)))*fc233.id;
    if (*pID != 0) { *pTY=5; *psigma = fc122.sigma; return; }
    
    *pID=((int)(((fc111.c+fc111.a)>tol)&&(fc111.c>0.0)))*fc111.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc131.c+fc131.a)>tol)&&(fc131.c>0.0)))*fc131.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc311.c+fc311.a)>tol)&&(fc311.c>0.0)))*fc311.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc331.c+fc331.a)>tol)&&(fc331.c>0.0)))*fc331.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc113.c+fc113.a)>tol)&&(fc113.c>0.0)))*fc113.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc133.c+fc133.a)>tol)&&(fc133.c>0.0)))*fc133.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc313.c+fc313.a)>tol)&&(fc313.c>0.0)))*fc313.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
    *pID=((int)(((fc333.c+fc333.a)>tol)&&(fc333.c>0.0)))*fc333.id;
    if (*pID != 0) { *pTY=6; *psigma = fc122.sigma; return; }
	
	return;
}

void fMinPosFC(int *n, int *m, int *indices, cNode *fc, int size)
{
	int i;
	double temp = 1.0; 
	*n = 0; 
	*m = 0; 
	
	for (i=0;i<size;i++) {
		if ((fc[i].c>0.0)&&(fc[i].ty>-1)&&(fc[i].m==0)) {
			*n = *n + 1; 
			indices[*n-1] = i; 
			if (fc[i].c <= temp){
				temp = fc[i].c ;
				*m = i; 
			}
		}
	}
	
	return;
}

void fMinPosFC_a(int *n, int *m, int *indices, cNode *fc, int size)
{
	int i;
	double temp = 1.0; 
	*n = 0; 
	*m = 0; 
	
	for (i=0;i<size;i++) {
		if ((fc[i].a>0.0)&&(fc[i].ty>-1)&&(fc[i].m==0)) {
			*n = *n + 1; 
			indices[*n-1] = i; 
			if (fc[i].a <= temp){
				temp = fc[i].a ;
				*m = i; 
			}
		}
	}
	
	return;
}

void fMaxPosFC(int *n, int *m, int *indices, cNode *fc, int size)
{

	int i;
	double temp = 0.0; 
	*n = 0; 
	*m = 0; 
	
	for (i=0;i<size;i++) {
		if ((fc[i].c<1.0)&&(fc[i].ty>-1)&&(fc[i].m==0)) {
			*n = *n + 1; 
			indices[*n-1] = i; 
			if (fc[i].c >= temp){
				temp = fc[i].c ;
				*m = i; 
			}
		}
	}
	
	return;
}

void fSpread(double exFC, double exFA, int id, double sigma, 
cNode *fc122, cNode *fc322, cNode *fc212, cNode *fc232,
cNode *fc221, cNode *fc223, double tol)
{
	cNode afc[6]; 
	int const no = 6; 
	int n, n2, i, m; 
	int indices[6];
	bool eflag1, eflag2; 
	double mi, ma, dc, dc2;
	
	afc[0] = *fc122; 
	afc[1] = *fc322;
	afc[2] = *fc212;
	afc[3] = *fc232;
	afc[4] = *fc221;
	afc[5] = *fc223;
	
	n=0; m=0; dc=0.0; 
	
	if (exFC<-tol) {
		fMinPosFC(&n, &m, indices, afc, no);
		if (n>0) {
			mi = afc[m].c ;
			if ((mi+exFC/n)>=0) {
				dc = exFC/n; 
				eflag1 = true;
			}
			else {
				dc = -mi;
				eflag1 = false;
			}
			
			for (i=0;i<n;i++) {
				afc[indices[i]].c += dc;
				afc[indices[i]].ty = 1; 
				if ((afc[indices[i]].c+afc[indices[i]].a) < tol) {
				// full melt
					afc[indices[i]].id = 1;
					afc[indices[i]].ty = 0;
					afc[indices[i]].sigma = 0.0;
				}
			}
		}
		else {
			eflag1 = true;
		}
		
	}
	else if (exFC>tol) {
		fMaxPosFC(&n, &m, indices, afc, no);
		if (n>0) {
			ma = afc[m].c ;
			if ((ma+exFC/n)<=1) {
				dc = exFC/n; 
				eflag1 = true;
			}
			else {
				dc = 1-ma;
				eflag1 = false;
			}
			
			for (i=0;i<n;i++) {
				if (afc[indices[i]].c < (tol)) {
					afc[indices[i]].id = id;  // id transfer 
					afc[indices[i]].sigma = sigma;  // sigma transfer 
					afc[indices[i]].ty = 4;
				}
				afc[indices[i]].c += dc;
				
				if (afc[indices[i]].c > (1-tol)) {
					afc[indices[i]].ty = 0;
				}
			}
		}
		else {
			eflag1 = true;
		}
		
	}
	else {
		eflag1 = true;
	}
	
	if (exFA<-tol) {
		fMinPosFC_a(&n2, &m, indices, afc, no);
		if (n2>0) {
			mi = afc[m].a ;
			if ((mi+exFA/n2)>=0) {
				dc2 = exFA/n2; 
				eflag2 = true;
			}
			else {
				dc2 = -mi;
				eflag2 = false;
			}
			
			for (i=0;i<n2;i++) {
				afc[indices[i]].a += dc2;
				afc[indices[i]].ty = 1; 
				if ((afc[indices[i]].c+afc[indices[i]].a) < tol) {
					afc[indices[i]].id = 1;
					afc[indices[i]].ty = 0;
					afc[indices[i]].sigma = 0.0;
				}
			}
		}
		else {
			eflag2 = true;
		}
	}
    else {
		eflag2 = true;
	} 
	
	 (*fc122).c = afc[0].c; (*fc122).a = afc[0].a; (*fc122).id = afc[0].id; (*fc122).ty = afc[0].ty; (*fc122).sigma = afc[0].sigma;
	 (*fc322).c = afc[1].c; (*fc322).a = afc[1].a; (*fc322).id = afc[1].id; (*fc322).ty = afc[1].ty; (*fc322).sigma = afc[1].sigma;
	 (*fc212).c = afc[2].c; (*fc212).a = afc[2].a; (*fc212).id = afc[2].id; (*fc212).ty = afc[2].ty; (*fc212).sigma = afc[2].sigma;
	 (*fc232).c = afc[3].c; (*fc232).a = afc[3].a; (*fc232).id = afc[3].id; (*fc232).ty = afc[3].ty; (*fc232).sigma = afc[3].sigma;
	 (*fc221).c = afc[4].c; (*fc221).a = afc[4].a; (*fc221).id = afc[4].id; (*fc221).ty = afc[4].ty; (*fc221).sigma = afc[4].sigma;
	 (*fc223).c = afc[5].c; (*fc223).a = afc[5].a; (*fc223).id = afc[5].id; (*fc223).ty = afc[5].ty; (*fc223).sigma = afc[5].sigma;

	if ((eflag1 == false)||(eflag2 == false)) {
		// recursion
		exFC = exFC - dc*(double)n;
		exFA = exFA - dc2*(double)n2;
		fSpread(exFC, exFA, id, sigma, fc122, fc322, fc212, fc232, fc221, fc223, tol);
	}
	
	return;	
} 
