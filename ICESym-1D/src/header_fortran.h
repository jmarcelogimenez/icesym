#include "structs.h"

extern "C" {
  
  void state_initial_cylinder(int* icyl, dataCylinder* myData, double* atm, dataSim* globalData, 
			      double* state_ini, double* mass_c, double* twall);
  void state_initial_junction(dataJunction* myData, double* atm, dataSim* globalData, double* state_ini, 
			      int* type_end);
  void state_initial_tube(dataTube* myData, double* atm,dataSim* globalData, double* state_ini, double* xnod, 
			  double* Area, double* twall, double* curvature, double* dAreax);
  void state_initial_tank(dataTank* myData, double* atm,dataSim* globalData, double* state_ini, int* type_end, 
			  double* Cd_ports);
  
  void solve_cylinder(int* icyl,dataCylinder* myData, dataSim* globalData, double* state, 
		      double* new_state, double* mass_C, double* twall);
  void solve_tube(dataTube* myData, dataSim* globalData, double* state, double* new_state, 
		  double* xnod, double* hele, double* Area, double* twall, double* curvature, 
		  double* dAreax, int* itube);
  void solve_tank(dataTank* myData, dataSim* globalData, double* state, double* new_state, 
		  int* type_end, double* Area_tube, double* twall_tube, double* dAreax_tube, 
		  double* Cd_ports);
  void solve_junction(dataJunction* myData, dataSim* globalData, double* state, double* new_state, 
		      int* type_end, double* dx_tube, double* Area_tube, double* twall_tube, double* dAreax_tube);
  
}
