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

#include "tank.h"

/**
   	\brief Tank's Constructor
   	\param all each Tank attribute 
*/
Tank::Tank(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, vector<double> state_ini, 
	        vector<int> histo, char* label, double Volume, double mass, double h_film, double Area_wall, double T_wall,
	 		vector<int> type_end,vector<double> Cd_ports, vector<int> int2tube, vector<int> exh2tube, int extras):
	Component(nnod,ndof,nnod_input,implicit,state_ini,histo,label){
		this->Volume	= Volume;
		this->mass		= mass;
		this->h_film	= h_film;
		this->Area_wall = Area_wall;
		this->T_wall	= T_wall;
		this->type_end  = type_end;
		this->Area_tube.resize(nnod_input);
		this->dAreax_tube.resize(nnod_input);
		this->twall_tube.resize(nnod_input);
		this->Cd_ports  = Cd_ports;
		this->int2tube  = int2tube;
		this->exh2tube  = exh2tube;
		this->extras	= extras;
		
}

/**
	\brief Tank's Copy Constructor
*/
Tank::Tank(Tank* t):Component(t->nnod,t->ndof,t->nnod_input,t->implicit,t->state_ini,t->histo, t->label){

		this->Volume	= t->Volume;
		this->mass		= t->mass;
		this->h_film	= t->h_film;
		this->Area_wall = t->Area_wall;
		this->T_wall	= t->T_wall;
		this->type_end  = t->type_end;
		this->Area_tube = t->Area_tube;
		this->twall_tube = t->twall_tube;
		this->dAreax_tube = t->dAreax_tube;
		this->Cd_ports  = t->Cd_ports;
		this->int2tube  = t->int2tube;
		this->exh2tube  = t->exh2tube;
		this->extras	= t->extras;
		
}


/**
	   \brief Makes the struct for send data to Fortran
*/
void Tank::makeStruct(dataTank &data){
		data.nnod		= this->nnod;
		data.nnod_input	= this->nnod_input;
		data.ndof		= this->ndof;
		data.Volume		= this->Volume;
		data.mass	    = this->mass;
		data.h_film		= this->h_film;
		data.Area_wall  = this->Area_wall;
		data.T_wall		= this->T_wall;
		data.nunit 		= this->nunit;
}

/**
	\brief Actualizes the component data from the struct that received from Fortran
*/
void Tank::undoStruct(dataTank &data){
		//this->nnod			= data.nnod;
		//this->nnod_input	= data.nnod_input;
		//this->ndof			= data.ndof;
		this->Volume		= data.Volume;
		this->mass			= data.mass;
		this->h_film		= data.h_film;
		this->Area_wall		= data.Area_wall;
		this->T_wall		= data.T_wall;
}

/**
	\brief Calcules the tank's new state calling Fortran
*/
void Tank::calculate(dataSim &globalData){
	dataTank myData;
	makeStruct(myData);
	vector<double> aux_Area_tube; 
	vector<double> aux_twall_tube;
	vector<double> aux_dAreax_tube;
	for(unsigned int i=0;i<this->nnod_input;i++){
		aux_Area_tube.push_back(*Area_tube[i]);
		aux_twall_tube.push_back(*twall_tube[i]);
		aux_dAreax_tube.push_back(*dAreax_tube[i]);
	}
	solve_tank(&myData, &globalData, &(this->state[0]), &(this->new_state[0]), &(this->type_end[0]), &(aux_Area_tube[0]),&(aux_twall_tube[0]),&(aux_dAreax_tube[0]), &(this->Cd_ports[0]));
	undoStruct(myData);
}

/**
   \brief Calcules the tank's initial state calling Fortran
*/
void Tank::calculate_state(double* atm, dataSim &globalData){
	dataTank myData;
	makeStruct(myData);
	state_initial_tank(&myData, &atm[0], &globalData, &(this->state_ini[0]), &(this->type_end[0]),&(this->Cd_ports[0]));
	undoStruct(myData);
}

/**
   \brief Open Fortran file units for save Extra data
*/
void Tank::openFortranUnit(char* file){
	int nfsize = strlen(file);
	open_unit(&(this->nunit),&file[0],&nfsize);
}

/**
   \brief Close Fortran file units
*/
void Tank::closeFortranUnit(){
	close_unit(&(this->nunit));
}
