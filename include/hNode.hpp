class hNode
{
public:
	double c, a; // l is redundant.
	int m; // m (0):silicon, (1):oxide, (2):Metal, (3):Virtual amorphous-phase change like crystal, heat transfer like aSi     
	double T, Told; 
	//double Q_phase, 
	double Q_laser_1, Tr_laser_1, Q_laser_2, Tr_laser_2;
	bool IsReflected;
	double dx, dy, dz;
	int x,y,z;
	
	//hNode *fc122, *fc322, *fc212, *fc232, *fc221, *fc223;
	
	hNode(); 
	//void init(double c_, double a_, int m_, double T_, double Told_, double q, 
	//	double dx_, double dy_, double dz_); 
	void init(double c_, double a_, int m_, double T_, double Told_, 
		double dx_, double dy_, double dz_, int x_, int y_, int z_); 
	//void init2(hNode *fc122_, hNode *fc322_, hNode *fc212_, hNode *fc232_, hNode *fc221_, hNode *fc223_);

	double get_rho();
	double get_Cp();
	double get_k();
	double get_alpha(double dt);
	
};
