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

#include "junction.h"


/**
   	\brief Junction's Constructor
   	\param all: each Junction attribute 
*/
Junction::Junction(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int type, vector<double> state_ini, 
				   vector<int> histo, char* label, vector<int> type_end,int modelo_junc,vector<int> node2tube):
Component(nnod,ndof,nnod_input,implicit,state_ini,histo,label){
	
	this->type_end	  = type_end;
	this->modelo_junc = modelo_junc;
	this->dx_tube.resize(nnod-1);
	this->Area_tube.resize(nnod);
	this->twall_tube.resize(nnod);
	this->dAreax_tube.resize(nnod);
	this->node2tube	  = node2tube;
}


/**
	\brief Junction's Copy Constructor
*/
Junction::Junction(Junction* j):Component(j->nnod,j->ndof,j->nnod_input,j->implicit,j->state_ini,j->histo,j->label){
	this->type_end	  = j->type_end;
	this->modelo_junc = j->modelo_junc;
	this->dx_tube     = j->dx_tube;
	this->Area_tube	  = j->Area_tube;
	this->twall_tube  = j->twall_tube;
	this->dAreax_tube = j->dAreax_tube;
	this->node2tube	  = j->node2tube;
}

/**
	   \brief Makes the struct for send data to Fortran
*/
void Junction::makeStruct(dataJunction &data){
	data.nnod			= this->nnod;
	data.nnod_input		= this->nnod_input;
	data.ndof			= this->ndof;
	data.modelo_junc    = this->modelo_junc;
}

/**
	\brief Actualizes the component data from the struct that received from Fortran
*/
void Junction::undoStruct(dataJunction &data){
	//this->nnod		= data.nnod;
	//this->nnod_input	= data.nnod_input;
	//this->ndof		= data.ndof;
	//this->modelo_junc = data.modelo_junc;
}

/**
	\brief Calcules the cylinder's new state calling Fortran
*/
void Junction::calculate(dataSim &globalData){
	dataJunction myData;
	makeStruct(myData);
	vector<double> aux_dx_tube; 
	vector<double> aux_Area_tube; 
	vector<double> aux_twall_tube;
	vector<double> aux_dAreax_tube;
	for(unsigned int i=0;i<nnod;i++){
		aux_dx_tube.push_back(*Area_tube[i]);
		aux_Area_tube.push_back(*Area_tube[i]);
		aux_twall_tube.push_back(*twall_tube[i]);
		aux_dAreax_tube.push_back(*dAreax_tube[i]);
	}
	solve_junction(&myData, &globalData, &(this->state[0]), &(this->new_state[0]),&(this->type_end[0]), 
				   &(aux_dx_tube[0]), &(aux_Area_tube[0]), &(aux_twall_tube[0]), &(aux_dAreax_tube[0]));
	undoStruct(myData);
}

/**
   \brief Calcules the cylinder's initial state calling Fortran
*/
void Junction::calculate_state(double* atm, dataSim &globalData){
	dataJunction myData;
	makeStruct(myData);
	state_initial_junction(&myData, &atm[0], &globalData, &(this->state_ini[0]),&(this->type_end[0]));
	undoStruct(myData);
}
