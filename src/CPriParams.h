#ifndef H_CPRIPARAMS
#define H_CPRIPARAMS

class CPriParams
{
public:

	CPriParams() {base = .95, power=2;}
	double base;
	double power; // p(grow) = base/(1+d)^power
};

#endif
