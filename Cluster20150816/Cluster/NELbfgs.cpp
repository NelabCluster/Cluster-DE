#include "NELbfgs.h"
#include "PE_Tool.h"
#include "Tool.h"

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

	double *cood = (double *)malloc(n * sizeof(double));
	for (i = 0; i < n; i++)
		cood[i] = x[i];
	double *R = (double *)calloc(N * N,sizeof(double));
	Tool::Distance(cood,R,n/3);
	fx = PE_Tool::FS_E(R,n / 3);
	//	fx = PE_Tool::Johnson_E(R,n / 3);
	PE_Tool::FS_F(cood,R,g,n/3);
	//	PE_Tool::Johnson_F(cood,R,g,n/3);
	for (i = 0 ;i < n; i++)
		g[i] = -g[i];

	free(cood);
	free(R);
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

double NELbfgs::local(double *cood, double *force, int N)
{
	int ret;
	lbfgsfloatval_t fx;
	lbfgs_parameter_t param;

	lbfgs_parameter_init(&param);
	ret = lbfgs(N*3, cood, &fx,evaluate, progress, NULL, &param);
	printf("L-BFGS optimization terminated with status code = %d\n", ret);
	printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, cood[0], cood[1]);
	return fx;
}


NELbfgs::NELbfgs(void)
{
}


NELbfgs::~NELbfgs(void)
{
}

