// Sapphire ver 2.0.1

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include "dSFMT.hpp"
#include "sConst2.h"
#include "pFunc.hpp"
#include "cNode.hpp"

int cNode::nID = 0;
int cNode::nx = 0 ;
int cNode::ny = 0;
int cNode::dnz1 = 0;
int cNode::fsb = 0; 
dsfmt_t cNode::dsfmt;

cNode::cNode()
{
	// Pseudo Cells
	c = 0.0; a = 0.0; l = 1.0;
	id = 0; ty = -1; m = 99; // initialized to vacuum mat.
	exFC = 0.0 ; exFA = 0.0 ;
	sigma = 0.0;
}

void cNode::init(double c_, double a_, int id_, int ty_, int m_, double sigma_, double T_, double dh_, int x_, int y_, int z_)
{
	c = c_ ;
	a = a_ ;
	l = 1.0 - c - a; 
	
	id = id_ ;
	ty = ty_;
	m = m_; 
	sigma = sigma_;
	T = T_ ;
	
	dh = dh_ ;
	x = x_; y = y_; z = z_;
}

void cNode::init2(cNode *fc122_, cNode *fc322_, cNode *fc212_, cNode *fc232_, cNode *fc221_, cNode *fc223_, 
cNode *fc112_, cNode *fc132_, cNode *fc312_, cNode *fc332_, cNode *fc121_, cNode *fc123_, cNode *fc321_, 
cNode *fc323_, cNode *fc211_, cNode *fc213_, cNode *fc231_, cNode *fc233_, cNode *fc111_, cNode *fc131_, 
cNode *fc311_, cNode *fc331_, cNode *fc113_, cNode *fc133_, cNode *fc313_, cNode *fc333_)
{
	fc122 = fc122_; fc322 = fc322_; fc212 = fc212_; fc232 = fc232_; fc221 = fc221_; fc223 = fc223_;
	fc112 = fc112_; fc132 = fc132_; fc312 = fc312_; fc332 = fc332_; fc121 = fc121_; fc123 = fc123_; fc321 = fc321_;
	fc323 = fc323_; fc211 = fc211_; fc213 = fc213_; fc231 = fc231_; fc233 = fc233_; fc111 = fc111_; fc131 = fc131_;
	fc311 = fc311_; fc331 = fc331_; fc113 = fc113_; fc133 = fc133_; fc313 = fc313_; fc333 = fc333_; 
}

void cNode::init3(int nID_, dsfmt_t dsfmt_)
{
	nID = nID_;
	dsfmt = dsfmt_;
}

bool cNode::ready()
{
	if ( (((c>0.0)&&(T>Tm)) || ((a>0.0)&&(T>Ta)) || (((c+a)<1.0) && (T<Tm))) && (m==0) )
	{
		return true;
	}
	else {
	return false;
	}						
}
	
int cNode::IsNearOx() 
{	// related to phase change at the interface, need to clarify the interface characteristics of MAT 1 and MAT 2. 
	return ((int)(fc122->m==1)+(int)(fc322->m==1)+(int)(fc212->m==1)+(int)(fc232->m==1)+
		(int)(fc221->m==1)+(int)(fc223->m==1));
}

void cNode::melt(double dt,double tol1)
{
		
	//Melting of aSi 
	
	if ((T>Ta) && (a>0.0)) {
		if ((ty!=1)&&(ty!=7)){ // definition of type 7: hetero nucleus on aSi surface
			ty = fUpdateType1(*this,*fc122,*fc322,*fc212,*fc232,*fc221,*fc223,tol1) ;
		}                          
		else {
			a = a - Vgm_a(T)*dt/dh; 
			if (a<0.0) {
				if (c<tol1) { //fully-melt
				id=1;	// id=1 for liquid silicon
				ty = 0;
				sigma = 0.0;
				}
				exFA = a; //negative , excessive melting
				a = 0.0; 
			}
		}
	}
	
	
	// Melting of cSi
	if ((T>Tm) && (c>0.0)) {
		if (ty!=1){
			ty = fUpdateType1(*this,*fc122,*fc322,*fc212,*fc232,*fc221,*fc223,tol1) ;
		}  
								
		else {
			c = c - Vgm_c(T)*dt/dh; // melting!
			if (c<0.0) {
				if (a<tol1) {
				id=1; // full-liquid
				ty = 0;
				sigma = 0.0;
				}
				exFC = c; //negative 
				c = 0.0; 
			}
		}
	}

}

void cNode::solidify(double dt,double tol1, int ppp,int ppp2,int ppp3)
{
	double prob, prob_he;
	double rnc, rnac; // random number
	double vgc; //growth speed
	int nox; // number of near oxide cells for het-nucleation
	int *pID = (int*) mxMalloc(sizeof(int)); 
	int *pTY = (int*) mxMalloc(sizeof(int));
	double *psigma = (double*) mxMalloc(sizeof(double));
	
	// Liquid crystallization (in cell exist liquid)
	if ((T<Tm) && ((c+a)<1.0) )
    {
		vgc = Grc(T);
		if ((c+a)<tol1)
        {  // fully liquid, ready for nucleation
			// heterogeneous nucleation on Oxide
			prob_he = 0.0; 
			nox = IsNearOx();
			
			if (nox>0)
            {
				
				prob_he = 1.0-exp(-(((double)nox)*dh*dh)*Nrc_h(T)*dt);
			}
			
			rnc = dsfmt_genrand_close_open(&dsfmt); 
			if (rnc < prob_he)
            { // Heterogeneous nucleation
				id = nID; 
				ty = 3; 
				c = Rc(T)/(dh) ; 
				mexPrintf("Het-Nucleation! %d at %f K in the position of (%d %d %d)\n", nID, T, ppp,ppp2,ppp3);
				mexPrintf("where %e < %e \n", rnc, prob_he); 
				nID++; // ID increment 
			}
			//
			
			prob = 1.0-exp(-pow(dh,3) * Nrc(T)*dt) ;
			if (rnc < prob)
            { // Homogeneous nucleation
				id = nID; 
				ty = 2; 
				c = Rc(T)/(dh/2) ; 
				mexPrintf("Ho-Nucleation! %d at %f K in the position of (%d %d %d)\n", nID, T, ppp, ppp2, ppp3);
				mexPrintf("where %e < %e \n", rnc, prob); 
				nID++; 
			}
			else
            {
				// Search Adjacent Cells
				getID26(pID, pTY, psigma, *fc122, *fc322, *fc212, *fc232,
					*fc221, *fc223, *fc112, *fc132, *fc312, *fc332, *fc121, *fc123,
					*fc321, *fc323, *fc211, *fc213, *fc231, *fc233, *fc111, *fc131,
					*fc311, *fc331, *fc113, *fc133, *fc313, *fc333) ;	
				
				if ((*pTY)>3)
                {
					
					id = *pID; // transfer of the adjacent cell information
					ty = *pTY;
					sigma = *psigma;
					
					//mexPrintf(" transfer %d  type %d \n",*pID,*pTY);
					//mexPrintf(" sigma  %f \n", sigma);
					
					switch (ty)
                    {
						case 4:
						c = c + (1+sigma)*vgc*dt/dh; 
						break;
						
						case 5:
						c = c + (1+sigma)*vgc*dt/(dh*sqrt(2.0));
						break; 
						
						case 6:
						c = c + (1+sigma)*vgc*dt/(dh*sqrt(3.0)) ;
						break;
					}
                }
            }
			
		}
		
		else if ((a>0.0) && (a<1.0) && (c<tol1)) 
		{ // mixture of aSi and lSi , heterogeneous nucleation at aSi - lSi interface
			prob = 1.0-exp(-(dh*dh) * Nrc_hEx(T)*dt) ;
			
			// random#(0~1) generation, double type
			rnac = dsfmt_genrand_close_open(&dsfmt);
			
			if (rnac < prob) { 
				// Disabled on Sept 4 2012
				/*
				id = nID; 
				ty = 7;  // definition of type 7: hetero nucleus on aSi surface
				c = Rc(T)/(dh/2.0) ;  
				mexPrintf("l-a Nucleation! %d \n",nID); 
// 					                     	mexPrintf("FC: %f , FA: %f \n", c, a) ;
				nID++; 
				 */
			}
 					              /*   	else {
 			                             	// Search Adjacent Cells
 			                            	getID26(ptempID, ptempTY, fc122, fc322, fc212, fc232,
 			                            	fc221, fc223, fc112, fc132, fc312, fc332, fc121, fc123,
 			                            	fc321, fc323, fc211, fc213, fc231, fc233, fc111, fc131,
 			                            	fc311, fc331, fc113, fc133, fc313, fc333) ;	
 		                            	 	
 			                            	if ((*ptempTY)>3) {
 				                            	
 				                            	id = *ptempID ; 
 				                            	ty = *ptempTY;
 				                            	
 			                            		switch (*ptempTY) {
 				                            		case 4:
 				                            		c = c + vgc*dt/dh; 
 				                            		break;
 				                            		
 				                            		case 5:
 				                            		c = c + vgc*dt/(dh*sqrt(2.0));
 				                            		break; 
 				                            		
 				                            		case 6:
 				                            		c = c + vgc*dt/(dh*sqrt(3.0)) ;
 				                            		break;
 			                            		}
 		                            		}
 	                            		}
										*/
		}
			
		else if ((c>0.0)&&(ty==1)) {
			ty = 3; // growth from partial melt
			c = c +(1+sigma)*vgc*dt/dh;
		}
		
		else if ((c>0.0)&&(ty==2)) { // growth of nucleus from homogeneous nucleaetion
			c = c + (1+sigma)*vgc*dt/(dh/2);
		}
		
		else if ((c>0.0)&&(ty==3)) {
			c = c + (1+sigma)*vgc*dt/dh;
		}
		
		else if ((c>0.0)&&(ty==4)) {
			/* debug
			if ((id==1)) {
			mexPrintf("id %d and ty %d is growing. \n", (int)id, (int)ty);
			} */
			c = c + (1+sigma)*vgc*dt/dh;
		}
		
		else if ((c>0.0)&&(ty==5)) {
			c = c + (1+sigma)*vgc*dt/(sqrt(2.0)*dh);
		}
		
		else if ((c>0.0)&&(ty==6)) {
			c = c + (1+sigma)*vgc*dt/(sqrt(3.0)*dh);
		}
		else if ((c>0.0)&&(ty==7)) {
			// cSi nuc on aSi type.
			c = c + (1+sigma)*vgc*dt/(dh/2);
		}
		
		if ((c+a)>1.0) {  
			exFC = c+a - 1.0 ; 
			c = 1.0 ; 
			ty = 0 ;
		}		                             	
	}
	
	mxFree(pID);
	mxFree(pTY);
	mxFree(psigma);
}

void cNode::relax(double dt,double tol1)
{
	if ((fabs(exFC)>tol1)||(fabs(exFA)>tol1))	
	{
		fSpread(exFC,exFA,id,sigma,fc122,fc322,fc212,fc232,fc221,fc223, tol1);
		//mexPrintf("exFA = %f \n", exFA); 
	}
	
	exFC = 0.0;
	exFA = 0.0;
	
	return;
}


