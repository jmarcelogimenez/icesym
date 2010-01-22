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

#ifndef _TANK_H_
#define _TANK_H_

#include "component.h"

using namespace std;

/**
 Clase que modela los Tanques
 */
class Tank: public Component 
{
public:
	Tank(unsigned int, unsigned int, unsigned int, int, vector<double>, vector<int>, char*, double, double, double, double, double, vector<int>, 
	     vector<double>, vector<int>, vector<int>,int);
	Tank(Tank* t);
	Tank(){};
	void makeStruct(dataTank &data);
	void undoStruct(dataTank &data);
	void calculate_state(double* atm, dataSim &globalData);
	void openFortranUnit(char* file);
	void closeFortranUnit();
	
	vector<int> int2tube;	   	/**< Mapping intake nodes - index of tube	 */
	vector<int> exh2tube;	   	/**< Mapping exhaust nodes - index of tube	 */
	vector<double*> Area_tube;  /**< Array of areas of tube's element connected */
	vector<double*> twall_tube;	/**< Array of temperatures of tube's element connected */
	vector<double*> dAreax_tube;/**< Array of area's differential of tube's element connected */
	bool extras;				/**< Save extra data for postprocess*/
	int nunit;					/**< Number of fortran unit file */
protected:

private:
	double Volume;				/**< Tank Volume */
	double mass;				/**< Mass in Tank */
	double h_film;				/**< Peculiar Coefficient of Heat Trasference */
	double Area_wall;			/**< Wall Area */
	double T_wall;				/**< Temperature of Wall */
	vector<int> type_end;		/**< External normal to tank (+1,-1) */
	vector<double> Cd_ports;	/**< Discharge coefficients */
	void calculate(dataSim &globalData);
};

extern "C"{
	void open_unit(int* nu, char* file, int* nfsize);
	void close_unit(int* nu);
}
#endif // _TANK_H_
