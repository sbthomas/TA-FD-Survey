#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv[]);

int main(int argc, char *argv[])
{
  FILE *pp;
  char *token, command[512], buffer[1025];
  char name0[128], name1[128];
  double latitude, longitude, height;
  int first = 1;
  double an;
  double latitude_sum, latitude_sum2, latitude_mean, latitude_sigma;
  double longitude_sum, longitude_sum2, longitude_mean, longitude_sigma;
  double height_sum, height_sum2, height_mean, height_sigma;

  if(argc!=2)
    {
      fprintf(stderr,"Usage: %s <.csv file>\n", argv[0]);
      exit(-1);
    }

  sprintf(command, "sort -r %s", argv[1]);

  pp = popen(command,"r");

  name0[0] = '\0';
  name1[0] = '\0';
    
  an = 0.0e0;
  while(fgets(buffer,1024,pp)!=NULL)
    {
      if((token=strtok(buffer,","))==NULL)continue;
      strcpy(name1,token);

      if((token=strtok(NULL,","))==NULL)continue;
      latitude = atof(token);
      if(fabs(fabs(latitude)-39.5e0) > 2.0e0)continue;

      if((token=strtok(NULL,","))==NULL)continue;
      longitude = atof(token);
      if(fabs(fabs(longitude)-112.5e0) > 2.0e0)continue;

      if((token=strtok(NULL,","))==NULL)continue;
      height = atof(token);
      if(fabs(fabs(height)-1600.0e0) > 500.0e0)continue;

      if(first==1)
	{
	  strcpy(name0,name1);
	  first = 0;
	}

      if(strcmp(name0,name1) != 0)
	{
	  if(an>=1.0e0)
	    {
	      latitude_mean  = latitude_sum/an;
	      if(an>1.0e0)
		{
		  latitude_sigma = ((latitude_sum2-an*pow(latitude_mean,2.0e0))/(an-1.0e0));
		  if(latitude_sigma<0.0e0)latitude_sigma = 0.0e0;
		  latitude_sigma = sqrt(latitude_sigma);
		}
	      else
		latitude_sigma = 0.0e0;

	      longitude_mean  = longitude_sum/an;
	      if(an>1.0e0)
		{
		  longitude_sigma = ((longitude_sum2-an*pow(longitude_mean,2.0e0))/(an-1.0e0));
		  if(longitude_sigma<0.0e0)longitude_sigma = 0.0e0;
		  longitude_sigma = sqrt(longitude_sigma);
		}
	      else
		longitude_sigma = 0.0e0;
	      
	      height_mean  = height_sum/an;
	      if(an>1.0e0)
		{
		  height_sigma = ((height_sum2-an*pow(height_mean,2.0e0))/(an-1.0e0));
		  if(height_sigma<0.0e0)height_sigma = 0.0e0;
		  height_sigma = sqrt(height_sigma);
		}
	      else
		height_sigma = 0.0e0;

	      printf("%s,", name0);
	      printf("%d,", (int)floor(an+0.5e0));
	      printf("%.9lf,", latitude_mean);
	      printf("%.9lf,", longitude_mean);
	      printf("%.3lf,", height_mean);
	      printf("%.9lf,", latitude_sigma);
	      printf("%.9lf,", longitude_sigma);
	      printf("%.3lf\n", height_sigma);

	      strcpy(name0,name1);
	    }

	  an = 1.0e0;
	  latitude_sum = latitude;
	  latitude_sum2 = latitude*latitude;
	  longitude_sum = longitude;
	  longitude_sum2 = longitude*longitude;
	  height_sum = height;
	  height_sum2 = height*height;
	}
      else
	{
	  an += 1.0e0;
	  latitude_sum += latitude;
	  latitude_sum2 += latitude*latitude;
	  longitude_sum += longitude;
	  longitude_sum2 += longitude*longitude;
	  height_sum += height;
	  height_sum2 += height*height;
	}

    }

  if(an>=1.0e0)
    {
      latitude_mean  = latitude_sum/an;
      if(an>1.0e0)
	{
	  latitude_sigma = ((latitude_sum2-an*pow(latitude_mean,2.0e0))/(an-1.0e0));
	  if(latitude_sigma<0.0e0)latitude_sigma = 0.0e0;
	  latitude_sigma = sqrt(latitude_sigma);
	}
      else
	latitude_sigma = 0.0e0;
      
      longitude_mean  = longitude_sum/an;
      if(an>1.0e0)
	{
	  longitude_sigma = ((longitude_sum2-an*pow(longitude_mean,2.0e0))/(an-1.0e0));
	  if(longitude_sigma<0.0e0)longitude_sigma = 0.0e0;
	  longitude_sigma = sqrt(longitude_sigma);
	}
      else
	longitude_sigma = 0.0e0;
      
      height_mean  = height_sum/an;
      if(an>1.0e0)
	{
	  height_sigma = ((height_sum2-an*pow(height_mean,2.0e0))/(an-1.0e0));
	  if(height_sigma<0.0e0)height_sigma = 0.0e0;
	  height_sigma = sqrt(height_sigma);
	}
      else
	height_sigma = 0.0e0;
      
      printf("%s,", name0);
      printf("%d,", (int)floor(an+0.5e0));
      printf("%.9lf,", latitude_mean);
      printf("%.9lf,", longitude_mean);
      printf("%.3lf,", height_mean);
      printf("%.9lf,", latitude_sigma);
      printf("%.9lf,", longitude_sigma);
      printf("%.3lf\n", height_sigma);

      strcpy(name0,name1);
    }
  
  return 0;
}
