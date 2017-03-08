#ifndef STRUCTS_H
#define STRUCTS_H

struct RSM_struct {
  double beta[6][7];
  double xmin[6][7];
  double xrange[6][7];
  double mean[6];
  double sd[6];
  double lamz[6];
  double lamws[6];
  double x[6][1000][7];
  double w[6][1000];
};

struct config_struct {
  int EarthGrav;
  int DensityModel;
  int Use2Body;
  int UseSphHarmon;
  int GDO;
  int UseSRP;
  int UseDrag;
  int DynamicMSIS;
  int ReadSW;
  int UseBCscale;
  int UseHFCd;
  int RSM;
  int UseSTimeAvg;
  int UseSun;
  int UseMoon;
  int DensityTimeLookup;
  int RungeKutta;
  int StateEnsembles;
  int Keplerian;
};

struct constants_struct {
  int DRAGFLAG;
  int SPHFLAG;
  int DOGF;
  int ATMDENMODEL;
  int EARTHGRAVMODEL;
  int SRPFLAG;
  int CDFLAG;
  int BCFLAG;
  int TAVGFLAG;
  int MOONFLAG;
  int SUNFLAG;
  int RKORDER;
  double JGM3RE;
  double JGM3RP;
  double JGM3GME;
  double JGM3OME;
  double EGM96GME;
  double EGM96RE;
  double MOONGM;
  double SUNGM;
  double EARTHECC;
  int GEOMETRY;
  int ADSMODEL;
  int GSIMODEL;
};

struct densities_struct {
  int Nfiles;
  double etlist[500];
  int Nvert;
  double AltKm[429];
  int Nlat;
  double LatDeg[37];
  int Nlon;
  double LonDeg[73];
  double Rho[2][73][37][429];
  double Temperature[2][73][37][429];
  double ndenvec[6][2][73][37][429];
  char files[500][100];
};

struct initCd_struct {
  int geometry;
  int gsi_model;
  int ads_model;
  double pitch;
  double yaw;
  double radius;
  double length;
  double height;
  double width;
  double surf_mass;
  double sat_mass;
  double proj_area;
  double constant_Cd;
  double SRP_Cr;
};

struct msis_struct {
  double F107;
  double F107A;
  double ap_index[9];
};

struct sat_struct {
  double pos[3];
  double vel[3];
  double kep[6];
  double rho;
  double temp;
  double nden[6];
  double Cd;
  double Cr;
  double BC;
  double Adrag;
  double Asrp;
  double m;
  double satacc[6][3];
  double phi;
  double theta;
  double shadowflag;
};

#endif  // STRUCTS_H
