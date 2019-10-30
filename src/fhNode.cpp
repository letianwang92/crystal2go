#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include "sConst2.h"
#include "fProp1.hpp"
#include "hNode.hpp"

hNode::hNode()
{
	c = 0.0; a = 0.0;
	m = 0; 
	//T= 0.0; Told = 0.0; 
	Q_laser_1 = 0.0; Tr_laser_1 = 0.0; 
	Q_laser_2 = 0.0; Tr_laser_2 = 0.0;
	IsReflected = false;
}


void hNode::init(double c_, double a_, int m_, double T_, double Told_,
	double dx_, double dy_, double dz_, int x_, int y_, int z_) 
{
	c = c_ ;
	a = a_ ;
	m = m_; 
	T = T_ ;
	Told = Told_;
	// Q_phase= 0.0; 
	
	dx = dx_ ; dy = dy_; dz= dz_;
	x = x_; y = y_; z = z_;

}


double hNode::get_rho()
{	
	return fRho(a, c, Told, m);
}

double hNode::get_Cp()
{
	return fCp(a, c, Told, m);
}

double hNode::get_k()
{
	return fK(a, c, Told, m);
}

double hNode::get_alpha(double dt)
{
	if (m != 99) {
		return (dt/(get_rho()*get_Cp()));
	}
	else {
		return 0.0; 
	}
}