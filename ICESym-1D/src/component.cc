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

#include <sstream>
#include "component.h"

/**
 	\brief Generic Constructor for all Elements
 */
Component::Component(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int implicit, vector<double>& state_ini, vector<int>& histo, char* label) {
	this->nnod			= nnod;
	this->ndof			= ndof;
	this->nnod_input	= nnod_input;
	this->implicit		= implicit;
	this->state_ini		= state_ini;
	this->histo			= histo;
	this->label = new char[strlen(label)];
	strcpy(this->label,label);
	state.resize((nnod+nnod_input)*ndof);
	new_state.resize(nnod*ndof);
}

/**
	\brief Copy Constructor
*/
// Component::Component(const Component& other)
// {
// 	    const Component* c  = &other;
// 		this->nnod			= c->nnod;
// 		this->ndof			= c->ndof;
// 		this->nnod_input	= c->nnod_input;
// 		this->implicit		= c->implicit;
// 		this->state_ini		= c->state_ini;
// 		this->state.resize((nnod+nnod_input)*ndof);
// 		this->new_state.resize(nnod*ndof);
// 		this->histo			= histo;
// }

Component::Component(Component* c){
		this->nnod			= c->nnod;
		this->ndof			= c->ndof;
		this->nnod_input	= c->nnod_input;
		this->implicit		= c->implicit;
		this->state_ini		= c->state_ini;
		this->state.resize((nnod+nnod_input)*ndof);
		this->new_state.resize(nnod*ndof);
		this->histo			= histo;
}

/**
	\brief This function makes the Initial Global State from state_ini's arrays
	\param Xn: Global State in time "n"
	\param Xn1: Global State in time "n+1"
	\param get_state: 0: from state_ini 1: from file, 2: from fortran
 */
void Component::Make(vector<double> &Xn, unsigned int &iXn, int get_state){
	for(unsigned int i=0;i<nnod;i++){
		for(unsigned int j=0;j<ndof;j++){
			if(!(get_state==1))
				Xn.push_back(state_ini[i*ndof+j]);
			i_state.push_back(iXn);
			i_new_state.push_back(iXn);
			iXn++;
		}
	}
}

/**
	\brief This function builds state's for all component from the Global State
	\param Xn1: Global State in time "n+1"
 */
void Component::build(vector<double>&Xn){
	for(unsigned int i=0; i<(nnod+nnod_input); i++){
		for(unsigned int j=0; j<ndof; j++){
				state[i*ndof+j] = Xn[i_state[i*ndof+j]];
		}
	}
}

/**
 	\brief: This function actualizes Global State from the Global New State
	\param Xn: Global State in time "n"
	\param Xn1: Global State in time "n+1"
	\param is_tube: true if the component is tube
 */
void Component::actualize(vector<double>&Xn1, vector<double>&Xn, bool is_tube=false){
	for(unsigned int i=0; i<nnod; i++){
		for(unsigned int j=0;j<ndof;j++){
			Xn1[i_new_state[i*ndof+j]] = new_state[i*ndof+j];
			if(is_tube){
				Xn[i_state[i*ndof+j]] = new_state[i*ndof+j];
			}
		}
	}
	
}

/**
	\brief This function separes solver step for each component in three differents tasks: build, calculate and actualize
	\param globalData: 
	\param Xn: Global State in time "n"
	\param Xn1: Global State in time "n+1"
	\param is_tube: true if the component is tube
 */
void Component::solver(dataSim &globalData, vector<double>& Xn, vector<double>& Xn1, bool is_tube=false){
	build(Xn);
	//cout<<"pasa build"<<endl;
	calculate(globalData);
	//cout<<"pasa calculate"<<endl;
	actualize(Xn1,Xn,is_tube);
	//cout<<"pasa actualize"<<endl;
}

/**
	\brief Copy "s" strchar in "r"
 */
void strcopy(char* &r,char* s){
	int l = strlen(s);
	r = (char *) malloc((l+1)*sizeof(char));
	for(int i=0;i<l;i++){
		r[i]=s[i];
	}
	r[l]='\0';
}

/**
	\brief Compares two strchar (equal or not)
 */
bool strcomp(char* r,char* s){
	int l1 = strlen(r);
	int l2 = strlen(s);
	if(l1!=l2)
		return false;

	for(int i=0;i<l1;i++){
		if(r[i]!=s[i])
			return false;
	}
	
	return true;
}

/**
	\brief Append "s" strchar to "r"
 */
char* strconcat(char* r,char* s){
	int l1 = strlen(r);
	int l2 = strlen(s);
	char* out = (char *) malloc((l1+l2+1)*sizeof(char));

	for(int i=0;i<l1;i++){
		out[i]=r[i];
	}
	for(int i=l1;i<(l1+l2);i++){
		out[i]=s[i-l1];
	}
	out[l1+l2] = '\0'; 
	
	return out;
}

/**
	\brief Convert int "n" in strchar
*/
char* int2char(int n){
		stringstream s;
		s.precision(5);s<<n;
		string aux = s.str();
		const char* strAux = aux.c_str();
		int l = strlen(strAux);
		char* out = (char *) malloc((l+1)*sizeof(char));
		for(int i=0;i<l;i++){
			out[i]=strAux[i];
		}
		out[l] = '\0';
		return out;
} 

/*!
	\brief Saves actual component's state in its file
	\param first: true if you need create the file, false: append file
	\param file: filename
	\param icycle: actual cycle
	\param crank_angle:
	\param time: actual time for this rpm
	\param xnod: position of nodes (only for tubes)
		
*/
void Component::saveHisto(bool first, char* file, int icycle, double crank_angle, double time,vector<double> &xnod){
	ofstream archim;
	if(first)
		archim.open(file,ios::trunc);
	else
		archim.open(file,ios::app);
	if(archim.is_open()){
		for(unsigned int i=0;i<this->histo.size();i++){
			string aux = "";
			stringstream s1,s3,s2,s4,s5;
			s1.precision(0);s3.precision(0);s2.precision(5);s4.precision(5);s5.precision(3);
			s1.width(3);s3.width(3);s2.width(12);s4.width(14);s5.width(10);
			s1<<icycle;s3<<this->histo[i];s2<<crank_angle;s4<<time;
			aux += s3.str();aux += s1.str();aux += s2.str();aux += s4.str();
			for(unsigned int j=0;j<this->ndof;j++){
				double value = this->new_state[this->histo[i]*this->ndof+j];
				stringstream ss;
				ss.precision(10);
				ss.width(20);
				ss<<value;
				aux += ss.str();
			}
			if (xnod.size() != 0){
				s5<<xnod[this->histo[i]];
				aux += s5.str();
			}
			archim<<aux<<"\n";
		}
		archim.close();	
	}
}
