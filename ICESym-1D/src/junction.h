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

#ifndef _JUNCTION_H_
#define _JUNCTION_H_
#include "component.h"

/**
	\class Junction Class
 */
class Junction: public Component 
{
public:
	Junction(unsigned int, unsigned int, unsigned int, int, vector<double>, vector<int> , char*, vector<int>, int, vector<int>);
	Junction(Junction* j);
	Junction(){};
	void makeStruct(dataJunction &data);
	void undoStruct(dataJunction &data);
	void calculate_state(double* atm, dataSim &globalData);
	vector<int> node2tube;  		/**< Mapping number of node (index) - connected tube (value) */
	vector<double*> dx_tube;		/**< Array of sizes of tube's element connected */
	vector<double*> Area_tube;		/**< Array of areas of tube's element connected */
	vector<double*> twall_tube;		/**< Array of temperatures of tube's element connected */
	vector<double*> dAreax_tube;	/**< Array of area's differential of tube's element connected */
protected:

private:
	vector<int> type_end;   		/**< External Normal of each connected tube (+1,-1) */
	int modelo_junc;				/**< Junction Model (0: basic) */
	void calculate(dataSim &globalData);
};

#endif // _JUNCTION_H_
