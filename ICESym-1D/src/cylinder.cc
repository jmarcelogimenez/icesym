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

#include "cylinder.h"

/**
   	\brief Cylinder's Constructor
   	\param all each Cylinder attribute 
*/
Cylinder::Cylinder(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, vector<double> state_ini, vector<int> histo, char* label, double Bore, 
				   double crank_radius, double Vol_clearance, double rod_length, double head_chamber_area, double piston_area,
				   double theta_0, double delta_ca, vector<double> Twall, vector<double> prop, vector<double> U_crevice, vector<double> data_crevice, 
				   vector<double> mass_C, int model_ht, double factor_ht, int scavenge,char* scavenge_type, int type_ig, int full_implicit, fuel fuel_data, 
				   combustion combustion_data, injection injection_data,vector<valve> intake_valves, vector<valve> exhaust_valves, Scavenge scavenge_data, int extras,int species_model):
Component(nnod,ndof,nnod_input,implicit,state_ini,histo,label){
	this->fuel_data				= fuel_data;
	this->combustion_data		= combustion_data;
	this->injection_data		= injection_data;
	this->intake_valves			= intake_valves;
	this->exhaust_valves		= exhaust_valves;
	this->Bore					= Bore;
	this->crank_radius			= crank_radius;
	this->Vol_clearance			= Vol_clearance;
	this->rod_length			= rod_length;
	this->head_chamber_area		= head_chamber_area;
	this->piston_area			= piston_area;
	this->theta_0				= theta_0;
	this->delta_ca				= delta_ca;
	this->Twall					= Twall;
	this->prop					= prop;
	this->U_crevice				= U_crevice;
	this->data_crevice			= data_crevice;
	this->mass_C				= mass_C;
	this->model_ht				= model_ht;
	this->factor_ht				= factor_ht;
	this->scavenge				= scavenge;
	strcopy(this->scavenge_type,scavenge_type);
	this->scavenge_data			= scavenge_data;
	this->type_ig				= type_ig;
	this->full_implicit			= full_implicit;
	this->extras				= extras;
	this->species_model			= species_model;
	this->nvi					= this->intake_valves.size();
	this->nve					= this->exhaust_valves.size();
}

/**
	\brief Cylinder's Copy Constructor
*/
Cylinder::Cylinder(Cylinder* c):Component(c->nnod,c->ndof,c->nnod_input,c->implicit,c->state_ini, c->histo, c->label){
	this->fuel_data			= c->fuel_data;
	this->combustion_data   = c->combustion_data;
	this->injection_data	= c->injection_data;
	this->intake_valves		= c->intake_valves;
	this->exhaust_valves	= c->exhaust_valves;
	
	this->Bore				= c->Bore;
	this->crank_radius		= c->crank_radius;
	this->Vol_clearance		= c->Vol_clearance;
	this->rod_length		= c->rod_length;
	this->head_chamber_area = c->head_chamber_area;
	this->piston_area		= c->piston_area;
	this->theta_0			= c->theta_0;
	this->delta_ca			= c->delta_ca;
	this->Twall				= c->Twall;
	this->prop				= c->prop;
	this->U_crevice			= c->U_crevice;
	this->data_crevice		= c->data_crevice;
	this->mass_C			= c->mass_C;
	this->model_ht			= c->model_ht;
	this->factor_ht			= c->factor_ht;
	this->scavenge			= c->scavenge;
	strcopy(this->scavenge_type,c->scavenge_type);
	this->type_ig			= c->type_ig;
	this->full_implicit		= c->full_implicit;
	this->extras			= c->extras;
	this->species_model		= c->species_model;
	this->nvi			    = c->nvi;
	this->nve				= c->nve;
}

/**
	\brief Initialize Fortran's types for def_cylinder.f90 module
	\param icycle: actual cycle
*/
void Cylinder::initFortran(int icyl){
	this->icyl = icyl;
	int nvi = this->nvi;
	int nve = this->nve;
	initialize_valves(&icyl, &nvi, &nve);
	initialize_fuel(&icyl,&(this->fuel_data));
	initialize_scavenge(&icyl, &(this->scavenge_data));
	int l1,l2; 
	if(this->type_ig == 1){
		l1 = this->injection_data.mfdot_array.size();
		initialize_injection(&icyl, &(this->injection_data.pulse), &(this->injection_data.m_inj), 
							 &(this->injection_data.dtheta_inj), &(this->injection_data.T_fuel), 
							 &(this->injection_data.theta_inj_ini), &(this->injection_data.theta_id), 
							 &(this->injection_data.integral), &(this->injection_data.mfdot_array[0]),
							 &(this->injection_data.ignition_delay_model), &l1);
	}
	
	l1 = this->combustion_data.xbdot_array.size();
	initialize_combustion(&icyl, &(this->combustion_data.theta_ig_0), &(this->combustion_data.dtheta_comb), 
						  &(this->combustion_data.phi), &(this->combustion_data.phi_ig), &(this->combustion_data.a_wiebe),
						  &(this->combustion_data.m_wiebe), &(this->combustion_data.xbdot_array[0]),
						  &(this->combustion_data.combustion_model), &(this->combustion_data.start_comb), &l1);
	for(unsigned int ival=0;ival<this->intake_valves.size();ival++){
		l1 = this->intake_valves[ival].Cd.size();
		l2 = this->intake_valves[ival].Lv.size();
		//cout<<"sizes: "<<l1<<", "<<l2<<endl;
		int ivalMod = ival+1;
		initialize_intake_valves(&icyl, &ivalMod,&(this->intake_valves[ival].Nval), &(this->intake_valves[ival].type_dat), 
								 &(this->intake_valves[ival].angle_V0), &(this->intake_valves[ival].angle_VC), 
								 &(this->intake_valves[ival].Dv), &(this->intake_valves[ival].Lvmax), 
								 &(this->intake_valves[ival].Cd[0]),&(this->intake_valves[ival].Lv[0]), 
								 &(this->intake_valves[ival].valve_model), &l1, &l2, this->intake_valves[ival].dx_tube,
								 this->intake_valves[ival].Area_tube, this->intake_valves[ival].twall_tube, 
								 this->intake_valves[ival].dAreax_tube, &(this->intake_valves[ival].tube));
	}
	//cout<<"aca"<<endl;
	for(unsigned int ival=0;ival<this->exhaust_valves.size();ival++){
		l1 = this->exhaust_valves[ival].Cd.size();
		l2 = this->exhaust_valves[ival].Lv.size();
		//cout<<"sizes: "<<l1<<", "<<l2<<endl;
		cout<<this->exhaust_valves[ival].Cd[l1-1]<<endl;
		int ivalMod = ival+1;
		initialize_exhaust_valves(&icyl, &ivalMod, &(this->exhaust_valves[ival].Nval), &(this->exhaust_valves[ival].type_dat), 
								  &(this->exhaust_valves[ival].angle_V0), &(this->exhaust_valves[ival].angle_VC), 
								  &(this->exhaust_valves[ival].Dv), &(this->exhaust_valves[ival].Lvmax), 
								  &(this->exhaust_valves[ival].Cd[0]), &(this->exhaust_valves[ival].Lv[0]), 
								  &(this->exhaust_valves[ival].valve_model), &l1, &l2, this->exhaust_valves[ival].dx_tube,
								  this->exhaust_valves[ival].Area_tube, this->exhaust_valves[ival].twall_tube, 
								  this->exhaust_valves[ival].dAreax_tube, &(this->exhaust_valves[ival].tube));
	}
	//cout<<"por here"<<endl;
	l1 = this->prop.size();
	l2 = this->U_crevice.size();
	int l3 = this->data_crevice.size();
	initialize_arrays(&icyl, &(this->prop[0]),&(this->U_crevice[0]), &(this->data_crevice[0]),&l1,&l2,&l3);
	cout<<"ya inicializo arrays fortran"<<endl;
}

/**
	   \brief Makes the struct for send data to Fortran
*/
void Cylinder::makeStruct(dataCylinder &data){
	data.nvi				= this->nvi;
	data.nve				= this->nve;
	data.nnod				= this->nnod;	
	data.nnod_input			= this->nnod_input;	
	data.ndof				= this->ndof;
	data.Bore				= this->Bore;
	data.crank_radius		= this->crank_radius;
	data.Vol_clearance		= this->Vol_clearance;
	data.rod_length			= this->rod_length;
	data.head_chamber_area	= this->head_chamber_area;
	data.piston_area		= this->piston_area;
	data.theta_0			= this->theta_0;
	data.delta_ca			= this->delta_ca;
	data.Twall				= this->Twall[0];
	data.factor_ht			= this->factor_ht;
	data.model_ht			= this->model_ht;
	data.type_ig			= this->type_ig;
	data.full_implicit		= this->full_implicit;
	data.scavenge			= this->scavenge;
	data.nunit 				= this->nunit;
	data.species_model		= this->species_model;
	data.ntemp 				= this->Twall.size();
	if(this->Twall.size() == 1)
		data.nh_temp =  false;
	else
		data.nh_temp =  true;

}

/**
	\brief Actualizes the component data from the struct that received from Fortran
*/
void Cylinder::undoStruct(dataCylinder &data){
	//this->nvi				= data.nvi;
	//this->nve				= data.nve;
	//this->nnod			= data.nnod;	
	//this->nnod_input		= data.nnod_input;	
	//this->ndof			= data.ndof;
	this->Bore				= data.Bore;
	this->crank_radius		= data.crank_radius;
	this->Vol_clearance		= data.Vol_clearance;
	this->rod_length		= data.rod_length;
	this->head_chamber_area	= data.head_chamber_area;
	this->piston_area		= data.piston_area;
	this->theta_0			= data.theta_0;
	this->delta_ca			= data.delta_ca;
	this->Twall[0]			= data.Twall;
	this->factor_ht			= data.factor_ht;
	//this->model_ht		= data.model_ht;
	//this->type_ig			= data.type_ig;
	//this->scavenge		= data.scavenge;
}

/**
	\brief Calcules the cylinder's new state calling Fortran
*/
void Cylinder::calculate(dataSim &globalData){
	dataCylinder myData;
	makeStruct(myData);
	solve_cylinder(&(this->icyl), &myData, &globalData, &(this->state[0]), &(this->new_state[0]), &(this->mass_C[0]), &(this->Twall[0]));
	undoStruct(myData);
}

/**
   \brief Calcules the cylinder's initial state calling Fortran
*/
void Cylinder::calculate_state(double* atm, dataSim &globalData){
	dataCylinder myData;
	makeStruct(myData);
	cout<<"crea struct"<<endl;
	state_initial_cylinder(&(this->icyl), &myData, &atm[0], &globalData, &(this->state_ini[0]), &(this->mass_C[0]), &(this->Twall[0]));
	cout<<"state initital pasa"<<endl;
	undoStruct(myData);
	cout<<"pasa undo"<<endl;
}

/**
   \brief Open Fortran file units for save Extra data
*/
void Cylinder::openFortranUnit(char* file){
	int nfsize = strlen(file);
	open_unit(&(this->nunit),&file[0],&nfsize);
}

/**
   \brief Close Fortran file units
*/
void Cylinder::closeFortranUnit(){
	close_unit(&(this->nunit));
}
