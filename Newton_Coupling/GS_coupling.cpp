#include "fluid_nl.h"
#include "precice/SolverInterface.hpp"
#include <iostream>
#include <stdlib.h>
#include "../utils/AD.hpp"
#include "Eigen/Dense"
#include <assert.h>


using std::cout;
using std::endl;

using namespace AD;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrixRM;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector;

int main(int argc, char** argv)
{
  cout << "Starting Fluid Solver..." << endl;

  if (argc < 4) {
    cout << endl;
    cout << "Usage: " << argv[0] << " N tau kappa" << endl;
    cout << endl;
    cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
    cout << "tau:   Dimensionless time step size." << endl;
    cout << "kappa: Dimensionless structural stiffness." << endl;
    cout << "eps:   Threshold for relative convergence criterion. " << endl;
    return -1;
  }

  int N = atoi(argv[1]);
  double tau = atof(argv[2]);
  double kappa = atof(argv[3]);
  double eps = (argc > 4) ? atof(argv[4]) : 1e-6;

  std::cout << "N: " << N << " tau: " << tau << " kappa: " << kappa << std::endl;

  // init data
  vector velocity                      = vector::Zero(N+1);
  vector velocity_n                    = vector::Zero(N+1);
  vector pressure                      = vector::Zero(N+1);
  vector pressure_n                    = vector::Zero(N+1);
  vector crossSectionLength            = vector::Zero(N+1);
  vector crossSectionLength_n          = vector::Zero(N+1);
  vector p_k                           = vector::Zero(N+1);
  vector cSL_k                         = vector::Zero(N+1);
  vector v_k                           = vector::Zero(N+1);
  
  vector positions                     = vector::Zero(N+1);
  
  matrix J_F = matrix::Zero(N+1, N+1);
  matrix J_S = matrix::Zero(N+1, N+1);
  matrix J   = matrix::Zero(2*(N+1), 2*(N+1));
  
  dualReal *v, *v_prev, *p, *p_prev, *cSL, *cSL_prev; 
  
  v         = new dualReal[N+1];
  v_prev    = new dualReal[N+1];
  p         = new dualReal[N+1];
  p_prev    = new dualReal[N+1];
  cSL       = new dualReal[N+1];
  cSL_prev  = new dualReal[N+1];

  for (int i = 0; i < pressure.size(); i++) {
    velocity(i)             = 1.0 / (kappa * 1.0);
    velocity_n(i)           = 1.0 / (kappa * 1.0);
    v_k(i)                  = 1.0 / (kappa * 1.0);
    crossSectionLength(i)   = 1.0;
    crossSectionLength_n(i) = 1.0;    
    cSL_k(i)                = 1.0;
    positions(i)            = double(i);
  }
  
  // coupling iterations counter
  int iter = 1; 
  int max_iter = 100;
  std::vector<int> coupling_iterations;
 
  // TODO: initialize Data
  
  // coupling iteration
  for(int t=1; t<=100; t++) {
    
    iter = 1;
    
    while(true)
    {
      /**
      * fluid solver call
      * 
      * ------------------------------------------------------------------------------------------
      */
      
      // save pressure, crossSectionArea and velocity for automatic differentiation of fluid solver
      // needs N+1 calls to fluid solver with identical state variables
      // p_k, cSL_k, v_k from previous iteration
      
      // build Jacobian J_F of Fluid solver, N+1 function evaluations
      for(int r = 0; r < crossSectionLength.size(); r++)
      {
	// get new dual datatypes
	for (int i = 0; i <= N; i++) {
	  double ei = (i==r) ? 1. : 0.;
	cSL[i] = dualReal(cSL_k[i], ei); 
	cSL_prev[i] = dualReal(crossSectionLength_n[i]);
	p[i] = dualReal(p_k[i]);
	p_prev[i] = dualReal(pressure_n[i]);
	v[i] = dualReal(v_k[i]);
	v_prev[i] = dualReal(velocity_n[i]);
	}
	
	// p_old is not used for gamma = 0.0
	fluid_nl(cSL, cSL_prev,                               // crossSectionLength
		v, v_prev,                                   // velocity
		p, p_prev, p,                                // pressure
		t/100.,                                      // scaled time for inflow condition (sample sine curve)
		N, dualReal(kappa), dualReal(tau), 0.0);     // dimensionless parameters
	
	// write back
	for (int i = 0; i <= N; i++){
	crossSectionLength[i] = cSL[i].u;
	crossSectionLength_n[i] = cSL_prev[i].u;
	pressure[i] = p[i].u;
	pressure_n[i] = p_prev[i].u;
	velocity[i] = v[i].u;
	velocity_n[i] = v_prev[i].u;
	}
	
	// \frac{ \partial F }{ \partial xr } is rth column of Jacobian
	for (int i = 0; i <= N; i++) 
	  J_F(i,r) = p[i].v;
      }
      //vector ptil = pressure_n + J_F*(crossSectionLength - crossSectionLength_n);
      //vector diff = pressure - ptil;
      //Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
      //std::cout<<"\n norm(p-ptil): "<<diff.norm()<<"\n p := "<<pressure.format(CommaInitFmt)<<"\n ptil:= "<<ptil.format(CommaInitFmt)<<std::endl;
      
      
      /**
      * solid solver call
      * 
      * ------------------------------------------------------------------------------------------
      */
      
      // build Jacobian J_F of Fluid solver, N+1 function evaluations
      for(int r = 0; r < pressure.size(); r++)
      {
	// get new dual datatypes
	for (int i = 0; i <= N; i++) {
	  double ei = (i==r) ? 1. : 0.;
	  cSL[i] = dualReal(cSL_k[i]); 
	  p[i] = dualReal(pressure[i], ei);     // directional derivative
	}
	
	// structure solver
	for (int i = 0; i <= N; i++) {  
	  cSL[i] = dualReal(4.0) / ((dualReal(2.0) - p[i]) * (dualReal(2.0) - p[i]));
	}
	      
	// \frac{ \partial S }{ \partial xr } is rth column of Jacobian
	for (int i = 0; i <= N; i++){ 
	  J_S(i,r) = cSL[i].v;
	}
      }
      
      // write back
      for (int i = 0; i <= N; i++){
	crossSectionLength[i] = cSL[i].u;
      }
      
      
      /**
      * checking for convergence
      * 
      * ------------------------------------------------------------------------------------------
      */
      vector res_d = crossSectionLength - cSL_k;
      vector res_p = pressure - p_k;
      double norm_res_d = res_d.norm();
      double norm_d = crossSectionLength.norm();
      double norm_res_p = res_p.norm();
      double norm_p = pressure.norm();
      
      std::cout<<"relative convergence meassure, cross-sectional area: ||res_a|| = "<<norm_res_d<<", limit = "<<eps*norm_d<<", conv: "<<(norm_res_d <= eps * norm_d)<<std::endl;
      std::cout<<"relative convergence meassure, pressure area:        ||res_p|| = "<<norm_res_p<<", limit = "<<eps*norm_p<<", conv: "<<(norm_res_p <= eps * norm_p)<<std::endl;
      if(norm_res_d <= eps * norm_d && norm_res_p <= eps * norm_p)
      {
	// CONVERGENCE
	std::cout<<"\n  ## time step "<<t<<" converged within "<<iter<<" iterations, advancing in time ##"<<std::endl;
	coupling_iterations.push_back(iter);
	
	break;
      }
      
      // check for divergence
      if(iter >= max_iter || isnan(crossSectionLength.norm())){
	std::cout<<"\n  ## divergence, residual does not come down, relative p res: "<<norm_res_p/norm_p<<", relative a res: "<<norm_res_d/norm_d<<", iter: "<<iter<<" ##"<<std::endl;
	return -1;
      }
      
      
      /**
      * coupling using Newton
      * 
      * ------------------------------------------------------------------------------------------
      */
      
    
      // building system Jacobian 
      J = J_S*J_F - Eigen::MatrixXd::Identity(N+1, N+1);
      
      // solve linear system for Newton update: J \Delta x = rhs
      // right hand side: -( S o F (d) - d ) = -( d_til - d )
      vector delta_x = J.fullPivHouseholderQr().solve(-res_d);
      
      // update displacements (x)
      crossSectionLength = cSL_k + delta_x;
      
      // update variables
      p_k = pressure;
      v_k = velocity;
      cSL_k = crossSectionLength;
      iter++;      
      
      // also reset current state variables, but keep crossSectionLength
      // (if no subcycling is enabled, coupling iterations are slightly worse if we reset
      //  velocity and pressure. For the case (tau, kappa) = (0.01, 10) we get 6.15 its
      //  without reset and 6.20 iterations if velocity and pressure is reset)
      //velocity = velocity_n; 
      //pressure = pressure_n;
    }
    
    
    /**
      * coupling for current time step converged, advancing in time
      * 
      * ------------------------------------------------------------------------------------------
      */
    
    // store state variables from last time step (required in fluid_nl
    velocity_n = velocity;
    pressure_n = pressure;
    crossSectionLength_n = crossSectionLength;       
    
  }
  std::cout<<"\n --------------- finished ----------------"<<std::endl;
  std::cout  <<" time step | iterations | total iterations"<<std::endl;
  int total = 0;
  for(int i=0; i< coupling_iterations.size(); i++){
    total += coupling_iterations[i];
    printf(" %i  %i  %i \n", i+1, coupling_iterations[i], total);
  }
  std::cout<<"\n average coupling iterations: "<<total/100.<<std::endl;
  std::cout<<" -----------------------------------------"<<std::endl;
  cout << "\n\nRUN FINISHED, t=100, exiting..." << endl;
  
  
  delete v; 
  delete v_prev;
  delete p;
  delete p_prev;
  delete cSL;
  delete cSL_prev;

  return 0;
}
