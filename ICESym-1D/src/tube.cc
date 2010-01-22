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

#include "tube.h"

/**
   	\brief Tube's Constructor
   	\param all: each Tube attribute 
*/
Tube::Tube(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, 
		   vector<double> state_ini, vector<int> histo, char* label, double longitud, vector<double> xnod, 
		   vector<double> Area, vector<double> twall, vector<double> curvature, 
		   vector<double> dAreax, char* tleft, unsigned int nleft, char* tright, unsigned int nright):
Component(nnod,ndof,nnod_input,implicit,state_ini,histo,label){
	this->longitud  = longitud;
	this->xnod		= xnod;
	this->Area		= Area;
	this->twall		= twall;
	this->curvature = curvature;
	strcopy(this->tleft,tleft);
	strcopy(this->tright,tright);
	this->nleft		= nleft;
	this->nright	= nright;
	this->dt_max	= 1000;
	this->hele.resize(this->nnod-1);
	for(unsigned int i=0;i<(nnod-1);i++){
		this->hele[i] = xnod[i+1]-xnod[i];
	}
	if(dAreax.size()==0){
		this->dAreax.resize(nnod);
		this->dAreax[0] = (Area[1]-Area[0])/(xnod[1]-xnod[0]);
		for(unsigned int i=1;i<(nnod-1);i++){
			this->dAreax[i] = (Area[i+1]-Area[i-1])/(xnod[i+1]+xnod[i-1]);
		}
		this->dAreax[nnod-1] = (Area[nnod-1]-Area[nnod-2])/(xnod[nnod-1]-xnod[nnod-2]);
	}
		
}

/**
	\brief Tube's Copy Constructor
*/
Tube::Tube(Tube* t):Component(t->nnod,t->ndof,t->nnod_input,t->implicit,t->state_ini,t->histo,t->label){
	this->longitud  = t->longitud;
	this->xnod		= t->xnod;
	this->hele		= t->hele;
	this->Area		= t->Area;
	this->twall		= t->twall;
	this->dAreax	= t->dAreax;
	this->curvature = t->curvature;
	strcopy(this->tleft,t->tleft);
	this->nleft		= t->nleft;
	strcopy(this->tright,t->tright);
	this->nright	= t->nright;
	this->dt_max	= t->dt_max; 
}


/**
	   \brief Makes the struct for send data to Fortran
*/
void Tube::makeStruct(dataTube &data){
	data.nnod		= this->nnod;
	data.ndof		= this->ndof;
	data.nnod_input	= this->nnod_input;
	data.longitud   = this->longitud;
	data.dt_max	    = this->dt_max;
}

/**
	\brief Actualizes the component data from the struct that received from Fortran
*/
void Tube::undoStruct(dataTube &data){
	//this->nnod		= data.nnod;
	//this->ndof		= data.ndof;
	//this->nnod_input	= data.nnod_input;
	//this->longitud	= data.longitud;
	this->dt_max	    = data.dt_max;
}

/**
	\brief Calcules the tube's new state calling Fortran
*/
void Tube::calculate(dataSim &globalData){
	dataTube myData;
	makeStruct(myData);
	solve_tube(&myData, &globalData, &(this->state[0]), &(this->new_state[0]), &(this->xnod[0]), 
			   &(this->hele[0]), &(this->Area[0]), &(this->twall[0]), &(this->curvature[0]), 
			   &(this->dAreax[0]), &(this->itube));
	undoStruct(myData);
}


/**
   \brief Calcules the tube's initial state calling Fortran
*/
void Tube::calculate_state(double* atm, dataSim &globalData){
	dataTube myData;
	makeStruct(myData);
	state_initial_tube(&myData, &atm[0], &globalData, &(this->state_ini[0]),&(this->xnod[0]), 
					   &(this->Area[0]), &(this->twall[0]), &(this->curvature[0]), &(this->dAreax[0]));
	undoStruct(myData);
}
