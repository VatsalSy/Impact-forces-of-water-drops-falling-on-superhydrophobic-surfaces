/* Title: getting Data from simulation snapshot
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[80];
int nx, ny, len;
double xmin, ymin, xmax, ymax, Deltax, Deltay, Ohd, Ohs;
scalar f[], D2c[], vel[];
scalar * list = NULL;

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]); ymin = atof(arguments[3]);
  xmax = atof(arguments[4]); ymax = atof(arguments[5]);
  ny = atoi(arguments[6]);
  Ohd = atof(arguments[7]); Ohs = atof(arguments[8]);

  list = list_add (list, D2c);
  list = list_add (list, vel);

  // boundary conditions
  u.t[left] = f[]*dirichlet(0.) + (1-f[])*neumann(0.0);
  f[left] = dirichlet(0.0);
  u.n[right] = neumann(0.);
  p[right] = dirichlet(0.0);
  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.0);


  /*
  Actual run and codes!
  */
  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});

  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/y);
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = 2*(clamp(f[], 0., 1.) * (Ohd-Ohs) + Ohs)*D2;

    if (D2c[] > 0.){
      D2c[] = log(D2c[])/log(10);
    } else {
      D2c[] = -10;
    }
    vel[] = sqrt(sq(u.x[])+sq(u.y[]));
  }
  boundary ((scalar *){D2c, vel});

  FILE * fp = ferr;
  Deltay = (double)(ymax-ymin)/(ny);
  // fprintf(ferr, "%g\n", Deltay);
  nx = (int)(xmax - xmin)/Deltay;
  // fprintf(ferr, "%d\n", nx);
  Deltax = (double)(xmax-xmin)/(nx);
  // fprintf(ferr, "%g\n", Deltax);
  len = list_len(list);
  // fprintf(ferr, "%d\n", len);
  double ** field = (double **) matrix_new (nx, ny+1, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*(i+1./2) + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*(j+1./2) + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  fclose (fp);
  matrix_free (field);
}
