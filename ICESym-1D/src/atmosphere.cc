#include "atmosphere.h"

/**
   	\brief Atmosphere's Constructor
   	\param all each Atmosphere attribute 
*/
Atmosphere::Atmosphere(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int type, vector<double> state_ini, vector<int> histo, char* label):
	Component(nnod,ndof,nnod_input,implicit,state_ini,histo, label){
}

/**
	\brief Atmosphere's Copy Constructor
*/
Atmosphere::Atmosphere(Atmosphere* a):Component(a->nnod,a->ndof,a->nnod_input,a->implicit,a->state_ini, a->histo, a->label){
		
}

/**
	\brief Calcules the atmosphere's new state
*/
void Atmosphere::calculate(dataSim &globalData){
	for(unsigned int i=0; i<(this->nnod*this->ndof);i++)	
		this->new_state[i] = this->state[i];
}

