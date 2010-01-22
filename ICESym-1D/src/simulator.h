/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * sim-c
 * Copyright (C) Juan Marcelo Gimenez 2009 <jmarcelogimenez@gmail.com>
 * 
 * sim-c is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * sim-c is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "cylinder.h"
#include "junction.h"
#include "tank.h"
#include "tube.h"
#include "atmosphere.h"

using namespace std;

/**
	 \class Class to control the Simulation
 */
class Simulator
{
public:
	Simulator(double,double,int,vector<double>,vector<double>,vector<double>,int,int,int,int,int,int,int,double,int,double,double,int,int,double,
			  int,int,int,char*,char*,char*,char*,char*,vector<int>,vector<Cylinder>,vector<Tube>,vector<Junction>,vector<Tank>,int,vector<Atmosphere>,int, int);
	void solver();
	void solveStep();
	Simulator(){};
	void printData();

	void actualizeDt();
	double GetTimeStep();
	vector<double> GetOldState();
    vector<double> GetNewState();
	int GetNewStateSize();
	int GetOldStateSize();
protected:

private:
	void solverEngine();
	void solverNOEngine();
	void createDir(bool);
	void createHeaderFile();
	void makeStruct(dataSim &data);
	void undoStruct(dataSim &data){};
	void makeSimulation();
	void saveState();
	void readState();
	void saveHisto(bool first);
	void openFortranUnits();
	void closeFortranUnits();
	void compressFiles();
	int get_state;					/**< 0: inicializa con state_ini, 1: inicializa con Xn, 2: inicializa en fortran*/
    double dt;
	double tf;
	double t; 
	int nrpms;
	int irpm;
	vector<double> rpms;
	int iteration;
	double time;
    vector<double> Xn;              /**< Global State at time "n" */
    vector<double> Xn1;             /**< Global State at time "n+1" */
	vector<double> final_times;
	int ntubes;						/**< Number of tubes */
	int ncyl;						/**< Number of Cylinders */
	int ntank;						/**< Number of Tanks */
	int njunc;						/**< Number of Junctions */
	int natm;						/**< Number of Atmospheres */
	int iter_sim1d;					/**< */
	int nsave;						/**< Frequency to save state (number of steps) */
	int nappend;					/**< Frequency to save historial (number of steps) */
	double dtheta_rpm;
	bool inicia;
	double Courant;					/**< Courant Number */
	double ga;						/**< Ga number*/
	int viscous_flow;				/**< */
	int heat_flow;					/**< */
	double R_gas;					/**< Gas constant */
	int engine_type;     	   		/**< Int que identifica el tipo de geometria del motor 0: "alternative", 1: "Opposed-Piston"*/
	double theta;					/**< Actual angle */
	double omega; 					/**< Actual angular velocity */
	int nstroke;					/**< Number of strokes for cycle */
	int ncycles;					/**< Number of cycles for rpm */
	int icycle;						/**< Actual cycle */
	unsigned int iAtm;				/**< Indice de la atmosfera en el arreglo Global */
	vector<int> ig_order;			/**< Ignition Order of Cylinder's */
	double rho_atm;					/**< Reference density */
	double u_atm;					/**< Reference velocity */
	double p_atm;					/**< Reference pressure */
	double crank_angle;				/**< Crank angle */

	char* filein_state;				/**< File to up the state*/
	char* filesave_state;			/**< File to save state */
	char* filein_spd;		
	char* filesave_spd;
	char* folder_name;				/**< Folder name for actual simulation */
	char* folderRPM;				/**< Folder name for actual rpm */
	char* folderGral;				/**< Folder name for actual simulation */

	bool calc_engine_data;			/**< Indicates if PostProcess calculates engine performance characteristics (power, torque, IMEP, BMEP, SFC, etc)*/
	
	vector<Cylinder> cylinders;
	vector<Tube> tubes;
	vector<Junction> junctions;
	vector<Tank> tanks;
	vector<Atmosphere> atmospheres;

	static const double pi = 3.1415926;
};

extern "C"{
	void initialize_cylinders(int* ncyl);
}

#endif // _SOLVER_H_
