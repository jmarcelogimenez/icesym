#ifndef _ATMOSPHERE_H_
#define _ATMOSPHERE_H_

#include "component.h"

using namespace std;

/**
	 \class Atmosphere Class
 */
class Atmosphere: public Component 
{
public:
	Atmosphere(unsigned int nnod, unsigned int ndof, unsigned int nnod_input, int type, vector<double> state_ini, vector<int> histo, char* label);
	Atmosphere(){};
	Atmosphere(Atmosphere* a);
protected:

private:
	void calculate(dataSim &globalData);
};

#endif // _ATMOSPHERE_H_

