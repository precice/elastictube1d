#ifndef AD_H_
#define AD_H_

#include <vector>

namespace AD {

  struct dualReal{
    
    double u;
    double v;
    dualReal() : u(0), v(0) {}
    dualReal(double _u) : u(_u), v(0) {}
    dualReal(double _u, double _v) : u(_u), v(_v) {}
  };
  
  
  // allows to take the derivative in all (or multiple) directions with a single function evaluation.
  struct dReal_grad{
    
    double u;
    std::vector<double> v;
    dReal_grad()                                   : u(0),  v()  {v.resize(1); v[0] = 0.;}
    dReal_grad(int n)                              : u(0),  v()  {v.resize(n); for (int i=0; i < n; i++){ v[i] = 0.;}}
    dReal_grad(double _u)                          : u(_u), v()  {v.resize(1); v[0] = 0.;}
    dReal_grad(double _u, int n)                   : u(_u), v()  {v.resize(n); for (int i=0; i < n; i++){ v[i] = 0.;}}
    dReal_grad(double _u, double _v)               : u(_u), v(1) {v[0] = _v;}
    dReal_grad(double _u, std::vector<double>& _v) : u(_u), v()  {v.resize(_v.size()); for (std::size_t i=0; i < _v.size(); i++){ v[i] = _v[i];}}
  };
  
  // dReal_grad methods
  
  dReal_grad operator+(const dReal_grad& lhs, const dReal_grad& rhs);
  
  dReal_grad operator-(const dReal_grad& lhs, const dReal_grad& rhs);
  
  void operator-(dReal_grad& rhs);
  
  dReal_grad operator*(const dReal_grad& lhs, const dReal_grad& rhs);
  
  dReal_grad operator/(const dReal_grad& lhs, const dReal_grad& rhs);
  
  void sin(dReal_grad& x);
  
  void cos(dReal_grad& x);
  
  void exp(dReal_grad& x);
  
  void log(dReal_grad& x);
  
  void pow(dReal_grad& x, int e);
  
  void abs(dReal_grad& x);
  
  void sqrt(dReal_grad& x);
  
  // dualReal methods
  
  dualReal operator+(const dualReal& lhs, const dualReal& rhs);
  
  dualReal operator-(const dualReal& lhs, const dualReal& rhs);
  
  dualReal operator-(const dualReal& rhs);
  
  dualReal operator*(const dualReal& lhs, const dualReal& rhs);
  
  dualReal operator/(const dualReal& lhs, const dualReal& rhs);
  
  dualReal sin(const dualReal& x);
  
  dualReal cos(const dualReal& x);
  
  dualReal exp(const dualReal& x);
  
  dualReal log(const dualReal& x);
  
  dualReal pow(const dualReal& x, int e);
  
  dualReal abs(const dualReal& x);
  
  dualReal sqrt(const dualReal& x);
}
#endif
