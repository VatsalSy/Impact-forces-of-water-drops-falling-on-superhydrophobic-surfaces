/* Title: Getting Energy
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Dec 07 2020
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

trace
double interface_energy (scalar c){
  double se = 0.;
  foreach (reduction(+:se)){
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      coord p, n = interface_normal (point, c);
      double alpha = plane_alpha (c[], n);
      double len = line_length_center(n, alpha, &p);
      se += 2.*pi*( y + p.y*Delta )*(len*Delta); // 2*pi*\int_l (r_c)dl
    }
  }
  return se;
}

scalar f[];
double ke1, ke, se, gpe, gpe1, rho1, rho2, Rhor, Ohd, mu1, mu2, Ohs, eps1, eps, Bo, vcm, zcm, wt, We;

char filename[80], nameEnergy[80];



int main(int a, char const *arguments[]) {
  sprintf (filename, "%s", arguments[1]);
  sprintf(nameEnergy, "%s", arguments[2]);
  Rhor = atof(arguments[3]);
  Ohd = atof(arguments[4]);
  Ohs = atof(arguments[5]);
  Bo = atof(arguments[6]);
  We = atof(arguments[7]);
  // fprintf(ferr, "Rhor %g, Ohd %3.2e, Ohs %3.2e, Bo %g\n", Rhor, Ohd, Ohs, Bo);
  // return 1;
  // boundary conditions
  u.t[left] = f[]*dirichlet(0.) + (1-f[])*neumann(0.0);
  f[left] = 0.;
  u.n[right] = neumann(0.);
  p[right] = dirichlet(0.0);
  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.0);


  FILE *fp;
  fp = fopen (nameEnergy, "a");
  restore (file = filename);

  Bo /= We;
  rho1 = 1.0; mu1 = Ohd/sqrt(We);
  rho2 = Rhor; mu2 = Ohs/sqrt(We);

  f.prolongation = refine_bilinear;

  boundary((scalar *){f, u.x, u.y});

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
  ke1 = 0., ke = 0., gpe = 0., gpe1 = 0., se = 0., eps1 = 0., eps = 0., vcm = 0., zcm = 0., wt = 0.;

  foreach (){
    double rho = clamp(sf[], 0., 1.)*(rho1 - rho2) + rho2;
    ke1 += (2*pi*y)*(0.5*clamp(sf[], 0., 1.)*rho1*(sq(u.x[]) + sq(u.y[])))*sq(Delta); // 2*pi*\int_A(0.5*rho*v^2)*r*dr*dz
    ke += (2*pi*y)*(0.5*rho*(sq(u.x[]) + sq(u.y[])))*sq(Delta); // 2*pi*\int_A(0.5*rho*v^2)*r*dr*dz

    gpe1 += (2*pi*y)*(clamp(sf[], 0., 1.)*rho1*Bo*x)*sq(Delta);
    gpe += (2*pi*y)*(rho*Bo*x)*sq(Delta); // 2*pi*\int_A(rho*g*z)rdrdz

    zcm += (2*pi*y)*(rho1*clamp(sf[], 0., 1.)*x)*sq(Delta);
    vcm += (2*pi*y)*(rho1*clamp(sf[], 0., 1.)*u.x[])*sq(Delta);
    wt += (2*pi*y)*rho1*clamp(sf[], 0., 1.)*sq(Delta);

    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));

    double mu = clamp(sf[], 0., 1.)*(mu1 - mu2) + mu2;
    eps1 += (2*pi*y)*( 2*mu1*clamp(sf[], 0., 1.)*D2 )*sq(Delta);
    eps += (2*pi*y)*( 2*mu*D2 )*sq(Delta);
  }
  zcm /= wt; vcm /= wt;


  boundary((scalar *){f, u.x, u.y});

  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});
  sf.prolongation = fraction_refine;
  boundary ({sf});
  
  se = (interface_energy (f) - 4*pi)/We;

  double Zmin = 0., temp = 0.;
  temp = HUGE;
  foreach_boundary(bottom){
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      // fprintf(ferr, "%f\n", x);
      if (fabs(x) < temp){
        temp = fabs(x);
      }
      if (temp < 2.5*Delta){
        temp = 0.;
      }
    }
  }
  Zmin = temp;
  /*
  Do calculations end
  */
  if (t == 0){
    fprintf(ferr, "Rhor %g, Ohd %3.2e, Ohs %3.2e, Bo %g\n", Rhor, Ohd, Ohs, Bo);
    fprintf(ferr, "t ke ke1 gpe gpe1 se eps eps1 vcm zcm Zmin\n");
    fprintf(fp, "t ke ke1 gpe gpe1 se eps eps1 vcm zcm Zmin\n");    
  }

  fprintf(ferr, "%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n", t, ke, ke1, gpe, gpe1, se, eps, eps1, vcm, zcm, Zmin);
  fprintf(fp, "%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n", t, ke, ke1, gpe, gpe1, se, eps, eps1, vcm, zcm, Zmin);
}
