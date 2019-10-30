
int fUpdateType1(class cNode fc222,class cNode fc122,class cNode fc322,
       class cNode fc212,class cNode fc232,class cNode fc221,class cNode fc223, double tol);



// Nucleation funcitons

// Crystalline Si

double Grc(double T);
	
double Vgm_c(double T);  // melting speed of crystalline silicon

double D(double T) ;

double I2(double T) ;

double Gc(double T) ;

double Nrc(double T) ;

double Nrc_h(double T) ;

double Rc(double T) ;

double Nrc_hEx(double T) ;

// aSi nucleation

double Vgm_a(double T)  ;

double Nra(double T) ;

double Gra(double T) ;


//

void getID26(int *ptempID, int *ptempTY, double *psigma,
class cNode fc122, class cNode fc322, class cNode fc212,class cNode fc232,
class cNode fc221, class cNode fc223, class cNode fc112, class cNode fc132, class cNode fc312, 
class cNode fc332, class cNode fc121, class cNode fc123, class cNode fc321, class cNode fc323, 
class cNode fc211, class cNode fc213, class cNode fc231, class cNode fc233, class cNode fc111, 
class cNode fc131, class cNode fc311, class cNode fc331, class cNode fc113, class cNode fc133, 
class cNode fc313, class cNode fc333) ;

void fMinPosFC(int *n, int *m, int *indices, class cNode fc, int size) ;

void fMinPosFC_a(int *n, int *m, int *indices, class cNode fc, int size) ;

void fMaxPosFC(int *n, int *m, int *indices, class cNode fc, int size) ;

void fSpread(double exFC, double exFA, int id, double sigma,
class cNode *fc122, class cNode *fc322, class cNode *fc212,class cNode *fc232,
class cNode *fc221,class cNode *fc223, double tol) ;
