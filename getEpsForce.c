/* Title: Getting Energy
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Sep 06 2021
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "tag.h"
#include "heights.h"

scalar f[];
double Ohd, mu1, eps1, We, pForce;

char filename[80], nameEnergy[80];

int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf(nameEnergy, "%s", arguments[2]);
  Ohd = atof(arguments[3]);
  We = atof(arguments[4]);

  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
  // fprintf(ferr, "Ohd %3.2e, We %g\n", Ohd, We);
  // return 1;

  // boundary conditions
  u.t[left] = f[]*dirichlet(0.) + (1-f[])*neumann(0.0);
  f[left] = 0.;
  u.n[right] = neumann(0.);
  p[right] = dirichlet(0.0);
  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.0);


  mu1 = Ohd/sqrt(We);
  f.prolongation = refine_bilinear;
  boundary((scalar *){f, u.x, u.y, p});

  // fprintf(ferr, "Ohd %3.2e, We %g\n", Ohd, We);
  // return 1;

  // tag all liquid parts starts
  scalar d[];
  double threshold = 1e-4;
  foreach(){
    d[] = (f[] > threshold);
  }
  int n = tag (d), size[n];
  for (int i = 0; i < n; i++){
    size[i] = 0;
  }
  foreach_leaf(){
    if (d[] > 0){
      size[((int) d[]) - 1]++;
    }
  }
  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n; i++){
    // fprintf(ferr, "%d %d\n",i, size[i]);
    if (size[i] > MaxSize){
      MaxSize = size[i];
      MainPhase = i+1;
    }
  }
  // tag all liquid parts ends

  scalar sf[];
  foreach()
    sf[] = (4.*f[] + 
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
  sf.prolongation = refine_bilinear;
  boundary ({sf});
  /*
  Do calculations start
  */
  eps1 = 0., pForce = 0.;

  face vector s[];
  s.x.i = -1;
  double yMax = -HUGE;
  double xMax = -HUGE;
  double vTip = 0., xTP = 0.;
  double uTip = 0., yTP = 0.;
  
  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase) {
      coord n1 = facet_normal (point, f, s);
      double alpha1 = plane_alpha (f[], n1);
      coord segment1[2];
      if (facets (n1, alpha1, segment1) == 2){
        double x1 = x + (segment1[0].x+segment1[1].x)*Delta/2.;
        double y1 = y + (segment1[0].y+segment1[1].y)*Delta/2.;
        if (y1 > yMax){
          yMax = y1;
          xTP = x1;
          vTip = interpolate (u.y, xTP, yMax);
        }
        if (y1 < 0.01){
          if (x1 > xMax){
            xMax = x1;
            yTP = y1;
            uTip = interpolate (u.x, xMax, yTP);
          }
        }
      }
    }
  }

  foreach (){
    if (d[] == MainPhase){
      double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
      double D22 = (u.y[]/max(y,1e-20));
      double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
      double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
      double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
      eps1 += (2*pi*y)*( 2*mu1*clamp(sf[], 0., 1.)*D2 )*sq(Delta);
    }
  }
  
  // calculate the force on the substrate
  double pdatum = 0, wt = 0;
  foreach_boundary(top){
    pdatum += 2*pi*y*p[]*(Delta);
    wt += 2*pi*y*(Delta);
  }
  if (wt >0){
    pdatum /= wt;
  }
  foreach_boundary(left){
    pForce += 2*pi*y*(p[]-pdatum)*(Delta);
  }

  boundary((scalar *){f, u.x, u.y, p});

  /*
  Do calculations end
  */

  FILE *fp;
  fp = fopen (nameEnergy, "a");
  restore (file = filename);

  if (t == 0){
    fprintf(ferr, "Ohd %3.2e\n", Ohd);
    fprintf(ferr, "t eps pforce betaMax vBeta Hmax vH\n");
    fprintf(fp, "t eps pforce betaMax vBeta Hmax vH\n");    
  }

  fprintf(ferr, "%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n", t, eps1, pForce, yMax, vTip, xMax, uTip);
  fprintf(fp, "%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n", t, eps1, pForce, yMax, vTip, xMax, uTip);
  fclose(fp);
}
