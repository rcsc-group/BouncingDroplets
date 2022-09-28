// Author: Vatsal Sanjay
// Link: https://github.com/VatsalSy/Lifting-a-sessile-drop/blob/master/CaseI/two-phaseDOD.h

#include "vof.h"

scalar f1[], f2[], *interfaces = {f1, f2};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  if (mu1 || mu2)
    mu = new face vector;
}

#ifndef rho
#define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
#define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

#ifdef FILTERED
scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};
#else
#define sf1 f1
#define sf2 f2
scalar *smearInterfaces = {sf1, sf2};
#endif

event properties (i++) {
  #ifdef FILTERED
    int counter1 = 0;
    for (scalar sf in smearInterfaces){
      counter1++;
      int counter2 = 0;
      for (scalar f in interfaces){
        counter2++;
        if (counter1 == counter2){
          // fprintf(ferr, "%s %s\n", sf.name, f.name);
        #if dimension <= 2
            foreach(){
              sf[] = (4.*f[] +
          	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
          	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
            }
        #else // dimension == 3
            foreach(){
              sf[] = (8.*f[] +
          	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
          	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
          		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
          		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
          	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
          	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
            }
        #endif
        }
      }
    }
    #endif
  #if TREE
    for (scalar sf in smearInterfaces){
      sf.prolongation = refine_bilinear;
      boundary ({sf});
    }
  #endif



  foreach_face() {
    double ff1 = (sf1[] + sf1[-1])/2.;
    double ff2 = (sf2[] + sf2[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff1+ff2);
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff1+ff2);
  }
  foreach()
    rhov[] = cm[]*rho(sf1[]+sf2[]);

#if TREE
  for (scalar sf in smearInterfaces){
    sf.prolongation = fraction_refine;
    boundary ({sf});
  }
#endif
}
