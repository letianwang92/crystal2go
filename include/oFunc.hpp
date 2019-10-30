double OptN_cSi(double lambda, double T);

double OptK_cSi(double lambda, double T);

double fAbsorb(double lambda, double k);

double fReflect(double n, double k);

//double LH1(double x, double y, double z, double t, double W, double tau, double Fa, double Fc, int Fo) ;
//
//double LH2(double z, double t, double W, double tau) ;
//
//double LH3(double z, double t, double W, double tau) ;
//
//double LH4(double z, double t, double W, double tau) ;

void Laser_surf(double lambda, double W, double heating, class hNode *pnode);

void Laser_bulk(double lambda, double W, double heating, class hNode *pnode);

