#include "dSFMT.hpp"

class cNode
{
public:
	double c, a, l; // l is redundant.
	int id, ty;  
	int m; // m (0):silicon, (1):oxide, (2):Metal, (3):Virtual amorphous-phase change like crystal, heat transfer like aSi     
	double T; 
	double exFC, exFA; 
	double dh;
	double sigma;
	int x, y, z;
	
	cNode *fc122, *fc322, *fc212, *fc232, *fc221, *fc223;
	cNode *fc112, *fc132, *fc312, *fc332, *fc121, *fc123, *fc321, 
		*fc323, *fc211, *fc213, *fc231, *fc233;
	cNode *fc111, *fc131, *fc311, *fc331, *fc113, *fc133, *fc313, *fc333;
	
	static int nID; 
	static dsfmt_t dsfmt; 
	static int nx, ny, dnz1, fsb; 
	
	cNode(); 
	void init(double c, double a, int id, int ty, int m, double sigma, double T, double dh, int x_, int y_, int z_);
	void init2(cNode *fc122_, cNode *fc322_, cNode *fc212_, cNode *fc232_, cNode *fc221_, cNode *fc223_, 
		cNode *fc112_, cNode *fc132_, cNode *fc312_, cNode *fc332_, cNode *fc121_, cNode *fc123_, cNode *fc321_, 
		cNode *fc323_, cNode *fc211_, cNode *fc213_, cNode *fc231_, cNode *fc233_, cNode *fc111_, cNode *fc131_, 
		cNode *fc311_, cNode *fc331_, cNode *fc113_, cNode *fc133_, cNode *fc313_, cNode *fc333_);
	static void init3(int nID, dsfmt_t dsfmt);
	bool ready();
	int IsNearOx();
	
	void melt(double dt, double tol1);
	void solidify(double dt, double tol1,int ppp,int ppp2,int ppp3);
	void relax(double dt, double tol1) ; 

};
