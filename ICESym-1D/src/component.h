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

#ifndef _COMPONENT_H_
#define _COMPONENT_H_
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "header_fortran.h"

using namespace std;

/**
 Componente genérico
 */
class Component
{
public:
	bool implicit;					/**< Resolution Type: implicit = true, explicit = false */
	unsigned int nnod;				/**< Number of own nodes */
    unsigned int ndof;				/**< Number of degress of freedom */
    unsigned int nnod_input;		/**< Number of Referency nodes */
	vector<double> state_ini;		/**< Array with the initial state, only used if get_state = 0 */
	vector<unsigned int> i_state;   /**< Indexs to global state in time "n"*/
	vector<double> state;			/**< Actual component state ((nnod+nnod_input)*ndof) */
	vector<double> new_state;	   	/**< Next component state (nnod*ndof) */
	vector<int> histo;				/**< Node's indexs to save state*/
	char* label;					/**< Identificator for this element */
	Component(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int type, vector<double> state_ini, vector<int> histo, char* label);
	void solver(dataSim &globalData, vector<double>&Xn, vector<double>&Xn1, bool is_tube);
	Component(){};
	Component(Component* c);
	void Make(vector<double> &Xn, unsigned int &iXn, int get_state);
	virtual void calculate_state(double* atm){}; 
	void saveHisto(bool first, char* file, int icycle, double crank_angle, double time,vector<double> &xnod);
protected:

private:
	vector<unsigned int>  i_new_state; /**< Indexs to global state in time "n+1" */
	void build(vector<double>&Xn); 
	virtual void calculate(dataSim &globalData){}; 
	void actualize(vector<double>&Xn1,vector<double>&Xn, bool is_tube); 
	vector<double> f(vector<double>&state, double t){vector<double> dx; return dx;};
	void F(vector<double>& state, vector<double>& new_state,double t){};
};


void strcopy(char* &r, char* s);
bool strcomp(char* r, char* s);
char* strconcat(char* r, char* s);
char* int2char(int n);
#endif // _COMPONENT_H_

