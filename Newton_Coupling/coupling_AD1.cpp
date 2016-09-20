#include "fluid_nl.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
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
  
  int N = 0;
  double tau = 0., kappa = 0.;
  bool parallel_scheme = false;
  double eps = 1e-6;
  double scaling_F = 1e-7;
  
  
  for(int i = 1; i < argc; ++i){
    if(std::string(argv[i]) == "-vectorial"){
      parallel_scheme = true;
    }else if(std::string(argv[i]) == "-serial"){
      parallel_scheme = false;
    }else if(std::string(argv[i]) == "-kappa"){
      kappa = atof(argv[i+1]);
      ++i;
    }else if(std::string(argv[i]) == "-N"){
      N = atoi(argv[i+1]);
      ++i;
    }else if(std::string(argv[i]) == "-tau"){
      tau = atof(argv[i+1]);
      ++i;
    }else if(std::string(argv[i]) == "-eps"){
      eps = atof(argv[i+1]);
      ++i;
    }else if(std::string(argv[i]) == "-lambda"){
      scaling_F = atof(argv[i+1]);
      ++i;
    }        
  }
  
  
  std::ofstream fout;
  fout.setf (std::ios::showpoint);
  fout.precision (3);
  fout.open("convergence-history.txt");

  std::cout << "N: " << N << " tau: " << tau << " kappa: " << kappa << std::endl;
  fout<<std::left<<"\n -- 1d flexible tube, N = "<<N<<", tau = "<<tau<<", kappa = "<<kappa<<" --\n"<<std::endl;
  fout.precision(3);
  fout<<std::left<<std::setw(15)<<"time step"<<std::setw(15)<<"iterations"<<std::setw(10)<<"alpha"<<std::setw(10)<<"rho(a)"<<std::setw(10)<<"rho(p)";
  fout.precision(16);
  fout<<std::setw(30)<<"res(a)"<<std::setw(30)<<"res(p)"<<std::setw(30)<<"relative res(a)"<<std::setw(30)<<"relative res(p)"<<std::endl;

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
  matrix J;
  
  if(parallel_scheme)
    J = matrix::Zero(2*(N+1), 2*(N+1));
  else 
    J = matrix::Zero(N+1, N+1);
  
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
  
  vector res_d, res_p, res_x;
  double norm_res_d = 1, norm_res_p = 1, norm_d, norm_p = 1;
  double norm_res_d_n = 1, norm_res_p_n = 1;
  double norm_res_d_nn = 1, norm_res_p_nn = 1;
  double norm_res_d_relative_0 = 1.;
  double norm_res_p_relative_0 = 1.;
  double alpha = 0.;               // estimated convergence rate
  double rho_d = 0., rho_p = 0.;   // average convergence rate per time step 
 
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
	fluid_nl(cSL, cSL_prev,                              // crossSectionLength
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
	  
	  if(parallel_scheme)
	    p[i] = dualReal(p_k[i], ei);          // directional derivative
	  else
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
      res_d = crossSectionLength - cSL_k;
      res_p = pressure - p_k;
      norm_res_d = res_d.norm();
      norm_d = crossSectionLength.norm();
      norm_res_p = res_p.norm();
      norm_p = pressure.norm();
      
      // estimate convergence rate
      if (iter > 2){
	alpha = std::log(std::fabs(norm_res_d/norm_res_d_n))/std::log(std::fabs(norm_res_d_n/norm_res_d_nn));	
      }
      if(iter > 1){
	rho_d   = std::pow((norm_res_d / norm_d) / norm_res_d_relative_0, 1./(double)iter); 
	rho_p   = std::pow((norm_res_p / norm_p) / norm_res_p_relative_0, 1./(double)iter); 
      }else if(iter == 1){
	norm_res_d_relative_0 = norm_res_d / norm_d;
        norm_res_p_relative_0 = norm_res_p / norm_p;
      }
	
      
      std::cout<<"relative convergence meassure, cross-sectional area: ||res_a|| = "<<norm_res_d<<", limit = "<<eps*norm_d<<", conv: "<<(norm_res_d <= eps * norm_d)<<std::endl;
      std::cout<<"relative convergence meassure, pressure area:        ||res_p|| = "<<norm_res_p<<", limit = "<<eps*norm_p<<", conv: "<<(norm_res_p <= eps * norm_p)<<std::endl;
      std::cout<<"convergence speed, alpha = "<<alpha<<", rho(res_a) = "<<rho_d<<", rho(res_p) = "<<rho_p<<std::endl;
      fout.precision(3);
      fout<<std::left<<std::setw(15)<<t<<std::setw(15)<<iter<<std::setw(10)<<alpha<<std::setw(10)<<rho_d<<std::setw(10)<<rho_p;
      fout.precision(16);
      fout<<std::setw(30)<<norm_res_d<<std::setw(30)<<norm_res_p<<std::setw(30)<<norm_res_d/norm_d<<std::setw(30)<<norm_res_p/norm_p<<std::endl;
      if(norm_res_d <= eps * norm_d && norm_res_p <= eps * norm_p)
      {
	// CONVERGENCE
	std::cout<<"\n  ## time step "<<t<<" converged within "<<iter<<" iterations, advancing in time ##"<<std::endl;
	coupling_iterations.push_back(iter);
	fout<<"\n";
	rho_d = 0; rho_p = 0; alpha = 0;
	
	break;
      }
      
      // save norms from previous iteration
      norm_res_d_nn = norm_res_d_n;
      norm_res_p_nn = norm_res_p_n;
      norm_res_d_n = norm_res_d;
      norm_res_p_n = norm_res_p;      
      
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
      if(parallel_scheme){
	J.block(0,0,N+1,N+1)     = - Eigen::MatrixXd::Identity(N+1, N+1);
	J.block(0,N+1,N+1,N+1)   = J_S;
	J.block(N+1,0,N+1,N+1)   = scaling_F * J_F;
	J.block(N+1,N+1,N+1,N+1) = - scaling_F* Eigen::MatrixXd::Identity(N+1, N+1);
	res_x = Eigen::VectorXd::Zero(2*N+2);
	res_x.segment(0,N+1)   = res_d;
	res_x.segment(N+1,N+1) = scaling_F * res_p;	
      }else{
	J = J_S*J_F - Eigen::MatrixXd::Identity(N+1, N+1);
	res_x = res_d;
      }
      
      // solve linear system for Newton update: J \Delta x = rhs
      // right hand side: -( S o F (d) - d ) = -( d_til - d )
      vector delta_x = J.fullPivHouseholderQr().solve(-res_x);
      
      // update displacements (x)
      if(parallel_scheme){
        crossSectionLength = cSL_k + delta_x.segment(0,N+1);
	pressure           = p_k   + delta_x.segment(N+1,N+1)/scaling_F;
      }else{
        crossSectionLength = cSL_k + delta_x;
      }
      
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
    
        
    Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    std::stringstream ss;
    ss << "Jacobian_t"<<t<<".m";
    std::ofstream mats;
    mats.open(ss.str());
    mats<<"J_t"<<t<<" = ...\n"<<J.format(OctaveFmt)<<";"<<std::endl;
    mats.close();
    
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
