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

#ifndef _TUBE_H_
#define _TUBE_H_
#include "component.h"

using namespace std;

/**
   Clase que modela los Tubos
*/
class Tube: public Component 
{
 public:
	Tube(unsigned int, unsigned int, unsigned int, int, vector<double>, vector<int>, char*,
	     double, vector<double>, vector<double>, vector<double>, vector<double>, 
		 vector<double>, char*, unsigned int, char*, unsigned int);
	Tube(Tube* t);
	Tube(){};
	void makeStruct(dataTube &data);
	void undoStruct(dataTube &data);
	void calculate_state(double* atm, dataSim &globalData);
	
	char* tleft;					/**< Component type connected at left */
	unsigned int nleft;				/**< Component number connected at left in global array*/
	char* tright;					/**< Component type connected at right */
	unsigned int nright;			/**< Component number connected at right in global array*/
	double dt_max;					/**< Maximum time step */
	vector<double> hele;			/**< Size of each element*/
	vector<double> Area;			/**< Transversal area for each node */
	vector<double> twall;			/**< Temperature tube wall for each node */
	vector<double> dAreax;			/**< Differential area in each node */
	int itube;						/**< index in global array */
	vector<double> xnod;			/**< Coordinates for each node */
 protected:

 private:
	double longitud;				/**< Tube lenght */
	vector<double> curvature;	    /**< Curvature for each node */
	void calculate(dataSim &globalData);
	
};

#endif // _TUBE_H_
