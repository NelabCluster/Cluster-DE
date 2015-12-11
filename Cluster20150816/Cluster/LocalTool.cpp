#include "LocalTool.h"
#include <lbfgs.h>
#include "Tool.h"

//µ¥ÀýÄ£Ê½
static LocalTool defaultLocalTool;


LocalTool::LocalTool(void)
{
}

LocalTool::~LocalTool(void)
{
}

LocalTool* LocalTool::shareDefault()
{
	return &defaultLocalTool;
}

void LocalTool::setPE(PotentialEnergyType type)
{
//	LocalTool::PE_E = PE_Tool::GetPEnergy(type);
//	LocalTool::PE_F = PE_Tool::GetPForce(type);
}


void LocalTool::setMethod(LocalOptimizeMethod m)
{
	method = m;
}

static lbfgsfloatval_t evaluate(
	void *instance,
	const lbfgsfloatval_t *x,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
	)
{
	int i;
	lbfgsfloatval_t fx = 0.0;
	int N = n /3;

	for (i = 0; i < n; i++)
		defaultLocalTool._cood[i] = x[i];

	Tool::Distance(defaultLocalTool._cood,defaultLocalTool._R,n/3);
	fx = defaultLocalTool._PE_E(defaultLocalTool._R,n /3);
	defaultLocalTool._PE_F(defaultLocalTool._cood,defaultLocalTool._R,g,n/3);
	for (i = 0 ;i < n; i++)
		g[i] = -g[i];

	return fx;
}

static int progress(
	void *instance,
	const lbfgsfloatval_t *x,
	const lbfgsfloatval_t *g,
	const lbfgsfloatval_t fx,
	const lbfgsfloatval_t xnorm,
	const lbfgsfloatval_t gnorm,
	const lbfgsfloatval_t step,
	int n,
	int k,
	int ls
	)
{
	/*printf("Iteration %d:\n", k);
	printf("  fx = %f\n", fx);
	printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
	printf("\n");*/
	return 0;
}

double LocalTool::localWithCoodAndPE(double *cood,int N, PEnergy PE_E, PForce PE_F)
{
	int ret;
	lbfgsfloatval_t fx;
	lbfgs_parameter_t param;

	_PE_E = PE_E;
	_PE_F = PE_F;

	_cood = new double[N * 3];
	_R = new double[N * N];

	lbfgs_parameter_init(&param);
	ret = lbfgs(N*3, cood, &fx,evaluate, progress, NULL, &param);

	delete(_cood);
	delete(_R);
	return fx;
}