// Author: Radu Cimpeanu
// Date: 26/09/2022

#include "axi.h"                     // axisymmetric geometry
#include "navier-stokes/centered.h"  // solve NS equations
#define FILTERED                     // Smear density and viscosity jumps
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phaseDOD.h"
#include "tension.h"                 // include surface tension between phases
#include "vof.h"                     // solve using VoF method
#include "fractions.h"               // initially defined fraction of fluid phases
#include "view.h"                    // need to make the animations
#include "draw.h"                    // visualisation helper
#include "tag.h"                     // helps track droplet properties

// Dimensional quantities (to be passed as arguments in main below)
double rhoLiquid; // liquid phase density (kg/m^3)
double rhoGas;    // gas phase density (kg/m^3)

double muLiquid;  // liquid dynamic viscosity (kg/ms)
double muGas;     // gas dynamic viscosity(kg/ms)

double sig;       // surface tension (N/m)

double g_accel;   // gravitational acceleration (m/s^2)

double dRadius;   // drop radius (m)

double v_init;    // initial drop velocity (m/s)

// dimensionless specifications (key groupings defined in main below)
#define rho_ratio   (rhoGas/rhoLiquid) // density ratio
#define mu_ratio    (muGas/muLiquid)   // viscosity ratio 

#define poolHeight 10.0                // Pool height (in radii)
#define domainSize 20.0                // Computational box size (in radii)

face vector av[];

FILE * fp_stats;
FILE * fp_vol;
FILE * fp_droplets;

double ND_Weber;
double ND_Reynolds;
double ND_Froude;
double ND_Bond;
double ND_Ohnesorge;

double filmHeight;

int minLevel = 6;
int maxLevel; // = 11;

double tEnd;

// Bottom of Pool = LEFT of domain
// No slip, no permeability
u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

// Side of Pool = TOP of domain
// Currently no slip and no permeability
// May change to slip condition
u.n[top] = dirichlet(0.); // Impermeability
u.t[top] = neumann(0.); // Slip

// Above the Pool = RIGHT of domain
// Outflow
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

// Default for bottom is symmetry

int main(int argc, char * argv[]) {

  rhoLiquid = atof(argv[1]);  // prescribed liquid density 
  rhoGas = atof(argv[2]);     // prescribed gas density
  muLiquid = atof(argv[3]);   // prescribed liquid (dynamic) viscosity
  muGas = atof(argv[4]);      // prescribed gas (dynamic) viscosity
  sig = atof(argv[5]);        // prescribed surface tension coefficient
  g_accel = atof(argv[6]);    // prescribed gravitational acceleration
  dRadius = atof(argv[7]);    // prescribed drop radius
  v_init = atof(argv[8]);     // prescribed initial velocity
  tEnd = atof(argv[9]);       // prescribed simulation end time
  maxLevel = atof(argv[10]);  // prescribed maximum resolution level
  
  ND_Weber = (rhoLiquid*pow(v_init,2.0)*dRadius)/sig;
  ND_Reynolds = (rhoLiquid*v_init*dRadius)/muLiquid;
  ND_Froude = v_init/pow(dRadius*g_accel,0.5);
  ND_Bond = rhoLiquid*g_accel*pow(dRadius,2.0)/sig;
  ND_Ohnesorge = muLiquid/pow(rhoLiquid*sig*dRadius,0.5);
    
  init_grid(1 << 8);
  
  size(domainSize);                     
  origin(-0.5*domainSize, 0.0);

  // Create output folders
  mkdir("Slices", 0700);
  mkdir("Animations", 0700);
  mkdir("Interfaces", 0700);

  // Print dimensionless numbers for verification
  fprintf(stdout, "Reynolds number = %0.6f \n", ND_Reynolds); fflush(stdout);
  fprintf(stdout, "Weber number = %0.6f \n", ND_Weber); fflush(stdout);
  fprintf(stdout, "Froude number = %0.6f \n", ND_Froude); fflush(stdout);
  fprintf(stdout, "Bond number = %0.6f \n", ND_Bond); fflush(stdout);
  fprintf(stdout, "Ohnesorge number = %0.6f \n", ND_Ohnesorge); fflush(stdout);

  rho1 = 1.;
  rho2 = rho_ratio;
  
  mu1 = 1./ND_Reynolds;
  mu2 = mu_ratio*mu1;
  
  f1.sigma = 1./ND_Weber;
  f2.sigma = 1./ND_Weber;

  a = av;

  // Pointer of the file to save stats
  {
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");
  }

  // Pointer of the file to droplet count info
  {
    char name[200];
    sprintf(name, "logdroplets.dat");
    fp_droplets = fopen(name, "w");
  }

  DT = 1e-2;
  NITERMIN = 1; // default 1
  NITERMAX = 200; // default 100
  TOLERANCE = 1e-6; // default 1e-3
  
  run();

  fclose(fp_stats);
  fclose(fp_droplets);
}

event acceleration (i++) {
  foreach_face(x)  
    av.x[] -= 1./pow(ND_Froude,2.0);
  foreach_face(y)  
    av.y[] += 0.0;
}

scalar omega[], viewingfield[], mylevel[], velnorm[];

event init (t = 0.0) {

  filmHeight = -domainSize/2. + poolHeight;

  // Strong refinement around the interfacial regions
  refine (((sq(x - (filmHeight + 1.0 + 0.5)) + sq(y) < sq(1.0*1.05) && sq(x - (filmHeight + 1.0 + 0.5)) + sq(y) > sq(1.0*0.95)) || fabs(x - filmHeight) <= 0.005) && level < maxLevel);
  
  // Create active liquid phase as union between drop and film
  fraction (f1, sq(1.0) - sq(x - (filmHeight + 1.0 + 0.5)) - sq(y));
  fraction (f2, - x + filmHeight);
  
  // Initialise uniform velocity field inside droplet
  foreach()
  {
  	u.x[] = -1.0*f1[];
        u.y[] = 0.0;
        p[] = 0.0;
	omega[] = 0.0;
  }
}

event adapt (i++) {

  // Refine only with respect to interfacial shape(s) location and velocity component magnitude
  adapt_wavelet ((scalar *){f1, f2, u}, (double[]){1e-5, 1e-5, 1e-3, 1e-3}, maxLevel, minLevel);

}

event gfsview (t = 0.0; t += 0.1; t <= tEnd) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/DropImpact-%0.1f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);
}

event saveInterfaces (t += 0.1) {

    char nameInterfaces1[200];

    sprintf(nameInterfaces1,"Interfaces/interfaceDrop-%0.1f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    output_facets (f1, fp1);	
    fclose(fp1);

    char nameInterfaces2[200];

    sprintf(nameInterfaces2,"Interfaces/interfacePool-%0.1f.dat",t);

    FILE * fp2 = fopen(nameInterfaces2, "w");
    output_facets (f2, fp2);	
    fclose(fp2);
}

event extractPressureData (t = 0.0; t += 0.1) {

  char namePressureData[200];

  sprintf(namePressureData,"Interfaces/customPressureData-%0.1f.dat",t);

  FILE *fpPressureData = fopen(namePressureData, "w");
  
  // Extract pressure data in tailored contact region (for visualisation and analysis inside gas film)
  for(double x = -1.; x < 0.5; x += 0.001){
	for(double y = 0.; y < 1.5; y += 0.001){
	    fprintf (fpPressureData, "%g\n", interpolate (p, x, y));
	}
  }

  fclose(fpPressureData);

}

// Fluid volume metrics
event droplets (t += 0.01)
{
  scalar m[];
  foreach()
    m[] = f1[] > 1e-2;
  int n = tag (m);

  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f1[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f1[]*p.x;
    }

  #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif
  for (int j = 0; j < n; j++)
    fprintf (fp_droplets, "%d %g %d %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fp_droplets);
}


// Output animations
event movies (t += 0.05){

  char timestring[100];
  
  foreach(){
	omega[] = (u.y[1,0] - u.y[-1,0])/(2.*Delta) - (u.x[0,1] - u.x[0,-1])/(2.*Delta);
        velnorm[] = sqrt(sq(u.x[]) + sq(u.y[]));
  	viewingfield[] = 1.0 - f1[] - 0.5*f2[];
  	mylevel[] = level;
  }

  view(width=1900, height=1050, fov=7.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });
	
  clear();
  draw_vof("f1", lw=2);
  draw_vof("f2", lw=2);
  squares("viewingfield", map = cool_warm, min = -0.5, max = 2.5);
  mirror({0,1}) {
	draw_vof("f1", lw=2);
	draw_vof("f2", lw=2);		
	cells(lw=0.5);
	squares("mylevel", map = cool_warm, min = minLevel, max = maxLevel);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactSummary.mp4");

  view(width=1900, height=1050, fov=7.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });;
  clear();
  
  draw_vof("f1", lw=2);
  draw_vof("f2", lw=2);
  squares("u.x", map = cool_warm, min = -1., max = 0.5,);
  mirror({0,1}) {
	draw_vof("f1", lw=2);
	draw_vof("f2", lw=2);		
	squares("u.y", map = cool_warm, min = -0.5, max = 2.);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactVelocities.mp4");

  view(width=1900, height=1050, fov=7.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });
  clear();
  
  draw_vof("f1", lw=2);
  draw_vof("f2", lw=2);
  squares("omega", map = cool_warm, min = -3., max = 3.);
  mirror({0,1}) {
	draw_vof("f1", lw=2);
	draw_vof("f2", lw=2);	
	squares("p", map = cool_warm, min = -0.25, max = 4.0);
  } 

  sprintf(timestring, "t=%2.03f",t);
  draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
  
  save ("Animations/ImpactPVort.mp4");

}

event logstats (t += 0.01) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // Output i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}
