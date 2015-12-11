#pragma once
#include <lbfgs.h>
#include "Common.h"

class NELbfgs
{
public:
	NELbfgs(void);
	~NELbfgs(void);

	 double local(double *cood, double *force, int N);
};

