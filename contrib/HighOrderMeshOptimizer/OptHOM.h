#ifndef _OPTHOM_H_
#define _OPTHOM_H_

#include <string>
#include <math.h>



#ifdef HAVE_BFGS

#include "ap.h"



class Mesh;

class OptHOM
{

public:

  Mesh mesh;

  OptHOM(GEntity *gf, std::set<MVertex*> & toFix, int method);
  void getDistances(double &distMaxBND, double &distAvgBND, double &minJac, double &maxJac);
  int optimize(double lambda, double lambda2, double barrier, int pInt, int itMax);  // optimize one list of elements
  void updateMesh(const alglib::real_1d_array &x);
  void evalObjGrad(const alglib::real_1d_array &x, double &Obj, alglib::real_1d_array &gradObj);
  void printProgress(const alglib::real_1d_array &x, double Obj);

  double barrier;

private:

//  double lambda, lambda2, powM, powP, invLengthScaleSq;
  double lambda, lambda2, jacBar, bTerm, invLengthScaleSq;
  int iter, progressInterv;            // Current iteration, interval of iterations for reporting
  double initObj, initMaxD, initAvgD;  // Values for reporting

  inline double setBarrierTerm(double jacBarrier) { bTerm = jacBarrier/(1.-jacBarrier); }
  inline double compute_f(double v);
  inline double compute_f1(double v);
  bool addJacObjGrad(double &Obj, alglib::real_1d_array &gradObj);
  bool addDistObjGrad(double Fact, double Fact2, double &Obj, alglib::real_1d_array &gradObj);
  void OptimPass(alglib::real_1d_array &x, const alglib::real_1d_array &initGradObj, int itMax);

};



inline double OptHOM::compute_f(double v)
{
  if (v > jacBar) {
    const double l = log((1 + bTerm) * v - bTerm);
    const double m = (v - 1);
    return l * l + m * m;
  }
  else return 1.e300;
//  if (v < 1.) return pow(1.-v,powM);
//  if (v < 1.) return exp((long double)pow(1.-v,3));
//  else return pow(v-1.,powP);
}



inline double OptHOM::compute_f1(double v)
{
  if (v > jacBar) {
    const double veps = (1 + bTerm) * v - bTerm;
    const double m = 2 * (v - 1);
    return m + 2 * log(veps) / veps * (1 + bTerm);
  }
  else return -1.e300;
//  if (v < 1.) return -powM*pow(1.-v,powM-1.);
//  if (v < 1.) return -3.*pow(1.-v,2)*exp((long double)pow(1.-v,3));
//  else return powP*pow(v-1.,powP-1.);
}



#endif



#endif
