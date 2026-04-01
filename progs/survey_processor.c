#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#define DEGS(DD,MM,SS)  (DD + (MM + SS/60.0e0)/60.0e0)
#define DEG2RAD (M_PI/180.0e0)
#define _GNU_SOURCE
#define NMAX 200
#define ORIGIN_LATITUDE  (39.296917937e0) 
#define ORIGIN_LONGITUDE (-112.908732522e0) 
#define ORIGIN_HEIGHT    (1370.046e0)

typedef struct POSITION
{
  char name[20];
  double eta, xi;
  double latitude, longitude, height;
  double east, north, up;
  double azimuth, elevation, distance;
}POSITION;

typedef struct MEASUREMENT
{
  int    set;
  int    count;
  int    status;
  char   base_name[20];
  double base_height;
  char   target_name[20];
  double target_height;
  int    azi_dd, azi_mm, azi_ss;
  char   elev_sign;
  int    elev_dd, elev_mm, elev_ss;
  double S, H, V;
  char notes[128];
  char name[20];
  double baseline_azimuth;
  POSITION base, target;
}MEASUREMENT;

int main( void);
int siteid( char name[]);
int mirid( char name[]);
int stationid( char name[]);
int print_measurements(int N, MEASUREMENT meas[], int ref);
int process_measurements(int N, MEASUREMENT meas[], int ref);
int process_measurement( MEASUREMENT *meas, int ref);
int set_base( MEASUREMENT *meas, int ref);
int set_target( MEASUREMENT *meas);
void printf_ddmmss( double degrees);
void printf2_ddmmss( double degrees);
int read_measurements( MEASUREMENT meas[]);
void read_monuments( void);
int set_base_position( POSITION *base, int ref, int site);
double compute_azimuth(POSITION *base, POSITION *target);
void xyz2enu(double eta, double xi, 
             double latitude, double longitude, 
             double dx, double dy, double dz,
             double *east, double *north, double *up);
void enu2xyz(double eta, double xi, 
             double latitude, double longitude, 
             double east, double north, double up,
             double *dx, double *dy, double *dz);
void xyz2geo(double x, double y, double z,
             double *latitude, double *longitude, double *height);
void geo2xyz(double latitude, double longitude, double height, 
	     double *x, double *y, double *z);
void geo2enu(double eta, double xi,
	     double lat1, double longi1, double height1,
	     double lat2, double longi2, double height2,
	     double *east, double *north, double *up);
void design_distances(int is, int im, int jm, double *S, double *H, double *V);


static POSITION TA_base[5][8];
static POSITION TA_site[4][15];
static char *refname[5] = {"OPUS","CSRS","AusPos","APPS","ITRF00"};
static char *stationname[4] = {"MD","BR","LR","CLF"};

int main( void)
{
  int i, N = 0;
  MEASUREMENT meas[NMAX];

  N = read_measurements( meas);
  for(i=0;i<5;i++)
    {
      process_measurements(N, meas, i);
      print_measurements(N, meas, i);
    }
  return 0;
}


int print_measurements(int N, MEASUREMENT meas[], int ref)
{
  static int FIRSTPASS = 1;
  int i,j, im, jm, is, js;
  double S, H, V;
  double S0, H0, V0;
  double an[4][15];
  double east, north, up;

  if(FIRSTPASS==1)
    {
#ifdef DISTANCES
      printf("name1,name2,ref,latitude,longitude,alt,dS,dH,dV\n");
#endif

#ifdef SUMMARY
      printf("name,ref,latitude,longitude,alt,east,north,up\n");
#endif
      FIRSTPASS = 0;
    }

#ifdef DISTANCES
  for(i=0;i<N;i++)meas[i].count = 0;

  for(i=0;i<N;i++)
    {
      im = mirid(meas[i].target_name);
      if((im<0)||(im>14))continue;
      is = stationid(meas[i].base_name);
      if((is<0)||(is>2))continue;
  
      if(meas[i].status==0)
	{
	  for(j=i+1;j<N;j++)
	    {
	      jm = mirid(meas[j].target_name);
	      if((jm<0)||(jm>14))continue;
	      js = stationid(meas[j].base_name);
	      if((js<0)||(js>2))continue;
	      if(meas[j].status==0)
		{
		  S  = pow(meas[j].target.east-meas[i].target.east,2.0e0);
		  S += pow(meas[j].target.north-meas[i].target.north,2.0e0);
		  H  = S;
		  V = meas[j].target.up - meas[i].target.up;
		  S += pow(V,2.0e0);
		  S  = sqrt(S);
		  H  = sqrt(H);

		  if(S<1000.0e0)
		    {
		      design_distances(is, im, jm, &S0, &H0, &V0);
		      
		      if(is==0)
			{
			  if(fabs(S-S0)*100.0e0>15.0e0)continue;
			  if(fabs(H-H0)*100.0e0>15.0e0)continue;
			  if(fabs(V-V0)*100.0e0>5.0e0)continue;
			}
		      else
			{
			  if(fabs(S-S0)*100.0e0>5.0e0)continue;
			  if(fabs(H-H0)*100.0e0>5.0e0)continue;
			  if(fabs(V-V0)*100.0e0>5.0e0)continue;
			}

		      meas[i].count++;
		      meas[j].count++;
		    }
		}
	    }
	}
    }
  
  for(is=0;is<3;is++)
    {
      for(im=0;im<15;im++)
	{
	  an[is][im] = 0.0e0;
	  TA_site[is][im].latitude  = 0.0e0;
	  TA_site[is][im].longitude = 0.0e0;
	  TA_site[is][im].height    = 0.0e0;
	}
    }
  
  for(i=0;i<N;i++)
    {
      im = mirid(meas[i].target_name);
      if((im<0)||(im>15))continue;
      is = stationid(meas[i].base_name);
      if((is<0)||(is>2))continue;
      if(meas[i].count<13)continue;
      if((is>0)&&(im>11))break;

      if(meas[i].status==0)
	{
	  an[is][im]                += 1.0e0;
	  TA_site[is][im].latitude  += meas[i].target.latitude;
	  TA_site[is][im].longitude += meas[i].target.longitude;
	  TA_site[is][im].height    += meas[i].target.height;
	}
    }

  printf("{|border=\"1\" align=\"center\" style=\"text-align:center;\"!Site !! Latitude [deg] !! Longitude [deg] !! height [m] !! East [m] !! North [m] !! Up [m] !! Processor\n");
  
  printf("|-\n");
  printf("| \'\'\'ORIGIN_CLF\'\'\' || ");
  printf2_ddmmss(ORIGIN_LATITUDE);
  printf(" || ");
  printf2_ddmmss(ORIGIN_LONGITUDE);
  printf(" || %.3lf", ORIGIN_HEIGHT);
  printf(" || %.3lf", 0.0e0);
  printf(" || %.3lf", 0.0e0);
  printf(" || %.3lf", 0.0e0);
  printf(" || %s\n", refname[ref]);
  printf("|}\n\n");

  for(is=0;is<3;is++)
    {
      printf("{|border=\"1\" align=\"center\" style=\"text-align:center;\"!Site !! Latitude [deg] !! Longitude [deg] !! height [m] !! East [m] !! North [m] !! Up [m] !! Processor\n");
      for(im=0;im<15;im++)
	{
	  if((is>0)&&(im>11))break;
	  if(an[is][im]<0.5e0)continue;

	  TA_site[is][im].latitude /= an[is][im];
	  TA_site[is][im].longitude /= an[is][im];
	  TA_site[is][im].height /= an[is][im];

	  printf("|-\n");
	  printf("| \'\'\'%s_M%.2d\'\'\' || ", stationname[is], im);
	  printf2_ddmmss(TA_site[is][im].latitude);
	  printf(" || ");
	  printf2_ddmmss(TA_site[is][im].longitude);
	  printf(" || %.3lf", TA_site[is][im].height);

	  geo2enu(0.0e0, 0.0e0, 
		  ORIGIN_LATITUDE*DEG2RAD,-fabs(ORIGIN_LONGITUDE*DEG2RAD),ORIGIN_HEIGHT,
		  (TA_site[is][im].latitude)*DEG2RAD,
		  -fabs((TA_site[is][im].longitude)*DEG2RAD),
		  TA_site[is][im].height,
		  &east, &north,  &up);

	  printf(" || %.3lf", east);
	  printf(" || %.3lf", north);
	  printf(" || %.3lf", up);

	  printf(" || %s\n", refname[ref]);
	}
      printf("|}\n\n");
    }

#endif

#ifdef SUMMARY
  for(i=0;i<N;i++)
    {
      if(meas[i].status==0)
	{
	  printf("%s_", meas[i].base_name);
	  printf("%s,", meas[i].target_name);
	  printf("%s,", refname[ref]);
	  printf("%.9lf,",meas[i].target.latitude);
	  printf("%.9lf,",meas[i].target.longitude);
	  printf("%.3lf,",meas[i].target.height);
	  printf("%.3lf,",meas[i].target.east);
	  printf("%.3lf,",meas[i].target.north);
	  printf("%.3lf\n",meas[i].target.up);
	}
    }
#endif

#ifdef FULL
  for(i=0;i<N;i++)
    {
      if(meas[i].status==0)
	{
	  printf("%d,", i+1);
	  printf("%s,", meas[i].base_name);
	  printf("%.3lf,", meas[i].base_height);
	  printf("%s,", meas[i].target_name);
	  printf("%.3lf,", meas[i].target_height);
	  printf("%d,%.2d,%.2d,", meas[i].azi_dd, meas[i].azi_mm, meas[i].azi_ss);
	  printf("%c,%d,%.2d,%.2d,", meas[i].elev_sign, meas[i].elev_dd, meas[i].elev_mm, meas[i].elev_ss);
	  printf("%.3lf,", meas[i].S);
	  printf("%.3lf,", meas[i].H);
	  printf("%.3lf,", meas[i].V);
	  printf_ddmmss(meas[i].base.latitude);
	  printf_ddmmss(meas[i].base.longitude);
	  printf("%.3lf,", meas[i].base.height);
	  printf_ddmmss(meas[i].base.azimuth);
	  printf_ddmmss(meas[i].target.latitude);
	  printf_ddmmss(meas[i].target.longitude);
	  printf("%.3lf,", meas[i].target.height);
	  printf("%.3lf,", meas[i].target.east);
	  printf("%.3lf,", meas[i].target.north);
	  printf("%.3lf,", meas[i].target.up);
	  printf_ddmmss(meas[i].target.azimuth);
	  printf("%.3lf,", meas[i].H - sqrt(pow(meas[i].target.east,2.0e0)+pow(meas[i].target.north,2.0e0)));
	  printf("%s,", meas[i].notes);
	  printf("%s\n", refname[ref]);
	}
    }
#endif

return 0;
}

int process_measurements(int N, MEASUREMENT meas[], int ref)
{
  double an, sum;
  double baseline_azimuth;
  int i;
  int j, n0, set;

  // Process measurements

  set = 0;
  for(i=0;i<N;i++)
    {
      for(n0=i;n0<N;n0++) // Find end of this measurement set
	{
	  if(siteid(meas[i].base_name)!=siteid(meas[n0].base_name))break;
	  if(fabs(meas[i].base_height-meas[n0].base_height)>0.001e0)break;
	}

      // Scan to the end of this measurement set for baseline measurements

      an = sum = 0.0e0;
      for(j=i;j<n0;j++)
	{
	  if(siteid(meas[j].base_name)>=0)
	    {
	      if(siteid(meas[j].target_name)>=0)
		{
		  meas[j].baseline_azimuth = (-1.0e0); // Flag baseline 
		  if(process_measurement( &meas[j], ref) == 0)
		    {
		      sum += meas[j].target.azimuth;
		      an += 1.0e0;
		      if(fabs(meas[j].target.azimuth-sum/an)>1.0e0/60.0e0)
			{
			  fprintf(stderr,"Invalid baseline measurement - ABORTED!\7\7\n");
			  exit(-1);
			}
		    }
		}
	    }
	}

      // Set average baseline measurement

      for(j=i;j<n0;j++)
	{
	  if(an>0.5e0) 
	    meas[j].baseline_azimuth = sum/an;
	  else
	    meas[j].baseline_azimuth = 0.0e0;
	  meas[j].set = set;
	}

      // Now process the set of measurements
      
      for(;i<n0;i++)process_measurement( &meas[i], ref);
      i--;
      set++;
    }

  return 0;
}

void printf_ddmmss( double degrees)
{
  int dd, mm;
  double ss;
  double degs;

  degs = fabs(degrees);

  dd = (int)floor(degs);
  degs -= (double)dd;
  degs *= 60.0e0;
  mm = (int)floor(degs);
  degs -= (double)mm;
  ss = 60.0e0*degs;

  if(degrees/degs< (-0.5e0))
    printf("-,%d,%.2d,%.5lf,", dd, mm, ss);
  else
    printf("+,%d,%.2d,%.5lf,", dd, mm, ss);

  return;
}

void printf2_ddmmss( double degrees)
{
  int dd, mm;
  double ss;
  double degs;

  degs = fabs(degrees);

  dd = (int)floor(degs);
  degs -= (double)dd;
  degs *= 60.0e0;
  mm = (int)floor(degs);
  degs -= (double)mm;
  ss = 60.0e0*degs;

  printf("%d %.2d %.5lf", dd, mm, ss);

  return;
}

int set_target(MEASUREMENT *meas)
{
  double azi, ele;

  meas->target.east = meas->target.north = meas->target.up = 0.0e0;
  meas->target.latitude = meas->target.longitude = meas->target.height = 0.0e0;
  meas->target.azimuth = meas->target.elevation = meas->target.distance = 0.0e0;

  meas->target.distance  = meas->S; // Convert to meters

  if(strcasestr(meas->notes,"cube") == NULL)
    {
      meas->target.distance += 30.0e0/1000.0e0; // Add 30 mm offset if not a cube reflector
    }

  meas->target.azimuth   = DEGS(meas->azi_dd, meas->azi_mm, meas->azi_ss);
  meas->target.elevation = DEGS(meas->elev_dd, meas->elev_mm, meas->elev_ss);
  if(meas->elev_sign == '-')meas->target.elevation = -fabs(meas->target.elevation);

  // Correct the target azimuth to true value CCW from north

  if(meas->baseline_azimuth >= 0.0e0)
    {
      meas->target.azimuth -= meas->baseline_azimuth; // Subract the measured baseline azimuth
      meas->target.azimuth += meas->base.azimuth;     // Add the official baseline azimuth
      if(meas->target.azimuth < 0.0e0)meas->target.azimuth = 360.0e0 - meas->target.azimuth;
    }

  // Compute the east, north, up coordinates of the target

  azi = meas->target.azimuth*DEG2RAD;
  ele = meas->target.elevation*DEG2RAD;

  meas->target.east  = cos(ele)*sin(azi)*meas->target.distance;
  meas->target.north = cos(ele)*cos(azi)*meas->target.distance;
  meas->target.up    = sin(ele)*meas->target.distance;

  // Correct the relative target height for instrument and target heights

  meas->target.up += meas->base_height;
  meas->target.up -= meas->target_height;
  
  // Compute the real target distance and elevation from the base
  
  meas->target.distance  = sqrt(pow(meas->target.east,2.0e0) + 
				pow(meas->target.north,2.0e0) + 
				pow(meas->target.up,2.0e0));

  meas->target.elevation = atan(meas->target.up/meas->target.distance)/DEG2RAD; 

  return 0;
}

int set_base(MEASUREMENT *meas, int ref)
{
  char *name;
  int index, num;
  int site = -1;


  meas->base.eta = meas->base.xi = 0.0e0;
  meas->base.east = meas->base.north = meas->base.up = 0.0e0;
  meas->base.latitude = meas->base.longitude = meas->base.height = 0.0e0;
  meas->base.azimuth = meas->base.elevation = meas->base.distance = 0.0e0;

  if(set_base_position(&(meas->base), ref, siteid(meas->base_name)) != 0)return (-1);

  return 0;
}

int process_measurement(MEASUREMENT *meas, int ref)
{
  double X0, Y0, Z0;
  double dX, dY, dZ;
  double X, Y, Z;
  double phi, lambda, h;

  meas->status = (-1);

  if(set_base(meas, ref) != 0){meas->status = -1; return -1;}
  if(set_target(meas) != 0){meas->status = -1; return -1;}

  geo2xyz(meas->base.latitude*DEG2RAD, -fabs(meas->base.longitude)*DEG2RAD, meas->base.height, &X0, &Y0, &Z0);
  enu2xyz(meas->base.eta*DEG2RAD/3600.0e0, meas->base.xi*DEG2RAD/3600.0e0,
	  meas->base.latitude*DEG2RAD, -fabs(meas->base.longitude)*DEG2RAD,
	  meas->target.east, meas->target.north, meas->target.up, 
	  &dX, &dY, &dZ);

  X = X0 + dX;
  Y = Y0 + dY;
  Z = Z0 + dZ;

  xyz2geo(X, Y, Z, &phi, &lambda, &h);

  meas->target.latitude  = phi/DEG2RAD;
  meas->target.longitude = -fabs(lambda)/DEG2RAD;
  meas->target.height    = h;

  geo2enu(0.0e0, 0.0e0, 
	  ORIGIN_LATITUDE*DEG2RAD,-fabs(ORIGIN_LONGITUDE*DEG2RAD),ORIGIN_HEIGHT,
	  (meas->target.latitude)*DEG2RAD,-fabs((meas->target.longitude)*DEG2RAD),meas->target.height,
	  &(meas->target.east), &(meas->target.north),  &(meas->target.up));

  meas->status = 0;

  return meas->status;
}

int read_measurements( MEASUREMENT meas[])
{
  FILE *fp;
  char buffer[1024], *token;
  int N = 0;

  if((fp=fopen("Survey_Record.csv","r"))==NULL)exit(-1);

  while(fgets(buffer,1023,fp)!=(char *)NULL)
    {
      if((token = strtok(buffer,","))==(char *)NULL)continue;
      strncpy(meas[N].base_name, token, 19);

      meas[N].base_height = 0.0e0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].base_height = atof(token);
	  if(meas[N].base_height < 0.0e0)continue;  // Can't be too small.
	  if(meas[N].base_height > 3.0e0)continue;  // Or too large
	}

      if((token = strtok(NULL,","))==(char *)NULL)continue;
      strncpy(meas[N].target_name, token, 19);

      meas[N].target_height = 0.0e0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].target_height = atof(token);
	  if(fabs(meas[N].target_height) > 3.0e0)continue;
	}

      meas[N].azi_dd = 0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].azi_dd = atoi(token);
	  if((meas[N].azi_dd<0)||(meas[N].azi_dd>359))continue;
	}

      meas[N].azi_mm = 0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].azi_mm = atoi(token);
	  if((meas[N].azi_mm<0)||(meas[N].azi_mm>59))continue;
	}

      meas[N].azi_ss = 0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].azi_ss = atoi(token);
	  if((meas[N].azi_ss<0)||(meas[N].azi_ss>59))continue;
	}

      meas[N].elev_sign = ' ';
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].elev_sign = *token;
	  if((meas[N].elev_sign != '-')&&(meas[N].elev_sign != '+'))continue;
	}

      meas[N].elev_dd = 0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].elev_dd = atoi(token);
	  if((meas[N].elev_dd<0)||(meas[N].elev_dd>359))continue;
	}

      meas[N].elev_mm = 0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].elev_mm = atoi(token);
	  if((meas[N].elev_mm<0)||(meas[N].elev_mm>59))continue;
	}

      meas[N].elev_ss = 0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].elev_ss = atoi(token);
	  if((meas[N].elev_ss<0)||(meas[N].elev_ss>59))continue;
	}

      meas[N].S = 0.0e0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].S = atof(token)*12.0e0*2.54e0/100.0e0;
	}

      meas[N].H = 0.0e0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].H = atof(token)*12.0e0*2.54e0/100.0e0;
	}

      meas[N].V = 0.0e0;
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  meas[N].V = atof(token)*12.0e0*2.54e0/100.0e0;
	}

      meas[N].notes[0] = '\0';
      if((token = strtok(NULL,","))!=(char *)NULL)
	{
	  strncpy(meas[N].notes, token, 127);
	  meas[N].notes[strlen(meas[N].notes)-1] = '\0';
	}
      
      if((meas[N].S+meas[N].H+meas[N].V)>0.0e0)N++;
    }
     
  fclose(fp);

  return N;
}

void read_monuments(void)
{
  FILE *fp;
  char buffer[1024], *token;
  int dd, mm;
  double ss, hh;
  int site, ref;

  if((fp=fopen("Survey_monuments.csv","r"))==NULL)exit(-1);

  while(fgets(buffer,1023,fp)!=(char *)NULL)
    {
      if((token = strtok(buffer,","))==(char *)NULL)continue;

      site = siteid( token);
      if((site<0)||(site>7))continue;

      if((token = strtok(NULL,","))==(char *)NULL)continue;
      for(ref=0;ref<5;ref++)if(strstr(token,refname[ref])!=(char *)NULL)break;
      if(ref>4)continue;

      dd = 0;
      if((token = strtok(NULL,","))==(char *)NULL)continue;
      dd = atoi(token);
      if((dd<0)||(dd>89))continue;

      mm = 0;
      if((token = strtok(NULL,","))==(char *)NULL)continue;
      mm = atoi(token);
      if((mm<0)||(mm>59))continue;

      ss = 0.0e0;
      if((token = strtok(NULL,","))==(char *)NULL)continue;
      ss = atof(token);
      if((ss<0.0e0)||(ss>=60.0e0))continue;

      if((token = strtok(NULL,","))==(char *)NULL)continue;

      TA_base[ref][site].latitude = (double)dd+((double)mm+ss/60.0e0)/60.0e0;

      dd = 0;
      if((token = strtok(NULL,","))==(char *)NULL)continue;
      dd = atoi(token);
      if((dd<0)||(dd>179))continue;

      mm = 0;
      if((token = strtok(NULL,","))==(char *)NULL)continue;
      mm = atoi(token);
      if((mm<0)||(mm>59))continue;

      ss = 0.0e0;
      if((token = strtok(NULL,","))==(char *)NULL)continue;
      ss = atof(token);
      if((ss<0.0e0)||(ss>=60.0e0))continue;

      if((token = strtok(NULL,","))==(char *)NULL)continue;

      TA_base[ref][site].longitude = -fabs((double)dd+((double)mm+ss/60.0e0)/60.0e0);

      hh = 0.0e0;
      if((token = strtok(NULL,","))==(char *)NULL)continue;
      hh = atof(token);
      if((hh<0.0e0)||(hh>=6000.0e0))continue;

      TA_base[ref][site].height = hh;

    }
     
  fclose(fp);

  return;
}


double compute_azimuth(POSITION *base, POSITION *target)
{
  double X0, Y0, Z0;
  double X1, Y1, Z1;
  double east, north, up;
  double azimuth;

  geo2xyz(base->latitude*DEG2RAD, -fabs(base->longitude)*DEG2RAD, base->height, &X0, &Y0, &Z0);
  geo2xyz(target->latitude*DEG2RAD, -fabs(target->longitude)*DEG2RAD, target->height, &X1, &Y1, &Z1);
  xyz2enu(base->eta*DEG2RAD/3600.0e0, base->xi*DEG2RAD/3600.0e0,
	  base->latitude*DEG2RAD, -fabs(base->longitude)*DEG2RAD,
	  X1-X0, Y1-Y0, Z1-Z0, 
	  &east, &north, &up);

  azimuth = atan2(east, north)/DEG2RAD;
  if(azimuth<0.0e0)azimuth += 360.0e0;

  return azimuth;
}

int set_base_position( POSITION *base, int ref, int site)
{
  static int FIRSTPASS = 1;
  int i, j;

  if((ref<0)||(ref>4))return (-1);
  if((site<0)||(site>7))return (-1);

  if(FIRSTPASS == 1)
    {
      // Load monument data from file

      read_monuments();

      for(j=0;j<5;j++)
	{
	  // BR1 
	  
	  TA_base[j][0].eta       = -1.00e0;
	  TA_base[j][0].xi        = -2.05e0;
	  
	  // BR2 
	  
	  TA_base[j][1].eta       = -1.02e0;
	  TA_base[j][1].xi        = -2.05e0;
	  
	  // LR1 
	  
	  TA_base[j][2].eta       = -0.08e0;
	  TA_base[j][2].xi        = -3.51e0;
	  
	  // LR2 
	  
	  TA_base[j][3].eta       = -0.07e0;
	  TA_base[j][3].xi        = -3.51e0;
	  
	  // MD1 
	  
	  TA_base[j][4].eta       = +1.77e0;
	  TA_base[j][4].xi        = -4.29e0;
	  
	  // MD2 
	  
	  TA_base[j][5].eta       = +1.77e0;
	  TA_base[j][5].xi        = -4.26e0;
	  
	  // CLF1 
	  
	  TA_base[j][6].eta       = -1.42e0;
	  TA_base[j][6].xi        = -1.88e0;
	  
	  // CLF2 
	  
	  TA_base[j][7].eta       = -1.44e0;
	  TA_base[j][7].xi        = -1.87e0;
	  
	  for(i=0;i<8;i++)
	    {
	      TA_base[j][i].azimuth = compute_azimuth(&TA_base[j][i], &TA_base[j][i + (i+1)%2 - i%2]);
	    }
	}
    }
  
  FIRSTPASS = 0;

  base->eta       = TA_base[ref][site].eta;
  base->xi        = TA_base[ref][site].xi;
  base->latitude  = TA_base[ref][site].latitude;
  base->longitude = -fabs(TA_base[ref][site].longitude);
  base->height    = TA_base[ref][site].height;
  base->azimuth   = TA_base[ref][site].azimuth;

  return 0;
}

// Radians and meters
void xyz2enu(double eta, double xi, 
             double latitude, double longitude, 
             double dx, double dy, double dz,
             double *east, double *north, double *up)
{
  double x, y, z;
  double dn, de, du, ds;

  longitude = -fabs(longitude);

  // Correction for the deflection of the vertical
  
  latitude   = asin(sin(latitude)/cos(eta)) + xi;
  longitude += asin(sin(eta)/cos(latitude));

  // rotate by longitude

  x = dx*cos(longitude) + dy*sin(longitude);
  y = dy*cos(longitude) - dx*sin(longitude);
  z = dz;

  // rotate by latitutude

  ds = x*cos(M_PI/2.0e0-latitude) - z*sin(M_PI/2.0e0-latitude);
  de = y;
  du = z*cos(M_PI/2.0e0-latitude) + x*sin(M_PI/2.0e0-latitude);

  dn = -ds;
  
  *east  = de;
  *north = dn;
  *up    = du;
  
  return;
}

// Radians and meters
void enu2xyz(double eta, double xi, 
             double latitude, double longitude, 
             double east, double north, double up,
             double *dx, double *dy, double *dz)
{
  double dn, de, du, ds;
  double s, x, y, z;
  
  longitude = -fabs(longitude);

 /* Correction for deflection of the vertical */
      
  latitude   = asin(sin(latitude)/cos(eta)) + xi;
  longitude += asin(sin(eta)/cos(latitude));

  dn = north;
  de = east;
  du = up;
  ds = -north;

  // rotate by s = (LAT-90 degrees)

  s  = latitude - M_PI/2.0e0;

  x  = ds*cos(s) - du*sin(s);
  y  = de;
  z  = du*cos(s) + ds*sin(s);

  // rotate by s = -(LON degrees)

  s  = -longitude;
 
  *dx = x*cos(s) + y*sin(s);
  *dy = y*cos(s) - x*sin(s);
  *dz = z;
  
  return;
}

void xyz2geo(double x, double y, double z,
             double *lat, double *longi, double *height)
{
  // WGS84 parameters

  double a = 6378137.0e0;
  double one_f = 298.257223563e0;
  double e2;
  double p, b, theta, S3, C3, N, ep2;

  /* Eccentricity squared */

  e2 = 2.0E0/one_f-1.0E0/one_f/one_f;

  /* p is the magnitude of the vector in the x/y plane */

  p  = sqrt(x*x+y*y);
  b  = a-(a/one_f);
  ep2 = (a*a - b*b)/(b*b);

  /* The source of the following calculations can be found at 
     www.colorado.edu/geography/gcraft/notes/datum/gif/11hxyz.gif */

  theta = atan(z*a/(p*b)); 
  S3 = pow(sin(theta), 3.0e0);
  C3 = pow(cos(theta), 3.0e0);

  *lat = atan((z+ep2*b*S3)/(p-e2*a*C3)); //latitude calculation
  *longi = atan2(y,x);                  //longitude calculation

  N = a/sqrt(1.0e0-e2*sin(*lat)*sin(*lat)); //radius of curvature in prime vertical
  *height = p/cos(*lat)-N;               //height calculation

  return;
}

void geo2xyz(double lat, double longi, double height, 
	     double *x, double *y, double *z)
{
  // WGS84 parameters

  double a = 6378137.0e0;
  double one_f = 298.257223563e0;
  double e2;
  double C, S, w, en;
  
  e2 = 2.0/one_f-1.0/one_f/one_f;  

  C  = cos(lat);
  S  = sin(lat);
  w  = sqrt(1.0E0-e2*S*S);
  en = a/w;
  
  *x = (en+height)*C*cos(longi);
  *y = (en+height)*C*sin(longi);
  *z = (en-en*e2+height)*S;
  
  return;
}

int siteid( char name[])
{
  if(strcasecmp(name, "BR1")==0)return 0;
  if(strcasecmp(name, "BR2")==0)return 1;
  if(strcasecmp(name, "LR1")==0)return 2;
  if(strcasecmp(name, "LR2")==0)return 3;
  if(strcasecmp(name, "MD1")==0)return 4;
  if(strcasecmp(name, "MD2")==0)return 5;
  if(strcasecmp(name, "CLF1")==0)return 6;
  if(strcasecmp(name, "CLF2")==0)return 7;
  return (-1);
}

int stationid( char name[])
{
  if(strcasestr(name, "MD")!=(char *)NULL)return 0;
  if(strcasestr(name, "BR")!=(char *)NULL)return 1;
  if(strcasestr(name, "LR")!=(char *)NULL)return 2;
  if(strcasestr(name, "CLF")!=(char *)NULL)return 3;
  return (-1);
}

int mirid( char name[])
{
  if(strcasecmp(name, "M00")==0)return 0;
  if(strcasecmp(name, "M01")==0)return 1;
  if(strcasecmp(name, "M02")==0)return 2;
  if(strcasecmp(name, "M03")==0)return 3;
  if(strcasecmp(name, "M04")==0)return 4;
  if(strcasecmp(name, "M05")==0)return 5;
  if(strcasecmp(name, "M06")==0)return 6;
  if(strcasecmp(name, "M07")==0)return 7;
  if(strcasecmp(name, "M08")==0)return 8;
  if(strcasecmp(name, "M09")==0)return 9;
  if(strcasecmp(name, "M10")==0)return 10;
  if(strcasecmp(name, "M11")==0)return 11;
  if(strcasecmp(name, "M12")==0)return 12;
  if(strcasecmp(name, "M13")==0)return 13;
  if(strcasecmp(name, "M14")==0)return 14;
  return (-1);
}

void geo2enu(double eta, double xi,
	     double lat1, double longi1, double height1,
	     double lat2, double longi2, double height2,
	     double *east, double *north, double *up)
{
  double x1, y1, z1; // geocentric x, y, z coordinates of the e/n/u origin
  double x2, y2, z2; // geocentric x, y, z coordinates of your point
  double dx, dy, dz; /* the resultant non-geocentric x, y, z vector from the
			e/n/u origin to your point that you get by subtracting
			(x1, y1, z1) from (x2, y2, z2) */

  /* Computes the Earth centered Cartesian coordinates for the origin */

  geo2xyz(lat1, -fabs(longi1), height1, &x1, &y1, &z1);

  /* Computes the Earth centered Cartesian coordinates for site */

  geo2xyz(lat2, -fabs(longi2), height2, &x2, &y2, &z2);

  /* Computes the origin-centered dx, dy, dz vector from the origin to the
     GPS position */

  dx = x2 - x1; 
  dy = y2 - y1; 
  dz = z2 - z1; 

  /* Rotate the e/n/u origin-centered dx, dy, dz vector to east, north, up 
     coordinates at the e/n/u origin */

  xyz2enu(eta, xi, 
	  lat1, -fabs(longi1), 
	  dx, dy, dz, 
	  east, north, up);
  
  return;
}

// This routine is to check the measured distances between mirrors against
// the original design distances. This should help to root out any bad 
// measurements. Given mirror number i and mirror number j it returns the 
// design values for S (Slant distance), H (Horizontal distance), and V (Vertical distance).
//
// This applies to LR and BR at this point in time.

void design_distances(int is, int im, int jm, double *S, double *H, double *V)
{
  static double X[12] = {-5.6984e0,-5.4371e0,-9.2736e0,-9.1059e0,2.2096e0,2.2674e0,
			 -2.2096e0,-2.2674e0,9.2736e0,9.1059e0,5.6984e0,5.4371e0};
  static double Y[12] = {11.6024e0,11.3411e0,9.0048e0,8.6756e0,12.7359e0,12.3710e0,
			 12.7359e0,12.3710e0,9.0048e0,8.6756e0,11.6024e0,11.3411e0};
  static double Z[12] = {1.5271e0,5.4980e0,1.5271e0,5.4980e0,1.5271e0,5.4980e0,
			 1.5271e0,5.4980e0,1.5271e0,5.4980e0,1.5271e0,5.4980e0};

  static double X0[15] = {0.0e0,20.5715e0,18.9705e0,15.1648e0,13.1114e0,8.5359e0,
			  6.1958e0,1.2192e0,-1.2192e0,-6.1958e0,-8.5359e0,
			  -13.1114e0,-15.1648e0,-18.9705e0,-20.5715e0};

  static double Y0[15] = {0.0e0,9.9214e0,8.0823e0,4.7972e0,3.4820e0,1.3988e0,
			  0.7136e0,0.0000e0,0.0000e0,0.7136e0,1.3988e0,
			  3.4820e0,4.7972e0,8.0823e0,9.9214e0};

  static double Z0[15] = {0.0e0,0.1065e0,0.0e0,0.1065e0,0.0e0,
			  0.1065e0,0.0e0,0.1065e0,0.0e0,0.1065e0,
			  0.0e0,0.1065e0,0.0e0,0.1065e0,0.0e0};

  static double W = 8.0e0*12.0e0*2.54e0/100.0e0;
  static double D = 11.0e0*12.0e0*2.54e0/100.0e0;
  static double R = 75.0e0*12.0e0*2.54e0/100.0e0;

  if(is==0)
    {
      *H = 1.013e0*sqrt(pow(X0[jm]-X0[im],2.0e0)+pow(Y0[jm]-Y0[im],2.0e0));
      *S = 1.013e0*sqrt(pow(X0[jm]-X0[im],2.0e0)+pow(Y0[jm]-Y0[im],2.0e0)+pow(Z0[jm]-Z0[im],2.0e0));
      *V = Z0[jm]-Z0[im];
    }
  else if((is==1)||(is==2))
    {
      *H = sqrt(pow(X[jm]-X[im],2.0e0)+pow(Y[jm]-Y[im],2.0e0));
      *S = sqrt(pow(X[jm]-X[im],2.0e0)+pow(Y[jm]-Y[im],2.0e0)+pow(Z[jm]-Z[im],2.0e0));
      *V = Z[jm]-Z[im];
    }
  else
    {
      *H = 0.0e0;
      *S = 0.0e0;
      *V = 0.0e0;
    }

  return;
}
