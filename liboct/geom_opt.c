#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>

#include "config.h"

struct geom_oct_type_struct
{
	void (* func) (double *x, double *f, double *df);

	int size;
	double *x, f, *df;
};

typedef struct geom_oct_type_struct geom_oct_type;

static void my_work_fdf(const gsl_vector *x, geom_oct_type *p)
{
	int i;

	for(i=0; i<p->size; i++){
		if(gsl_vector_get(x, i) != p->x[i])
			break;
	}

	if(i<p->size){ /* recalculate */
		for(i=0; i<p->size; i++)
			p->x[i] = gsl_vector_get(x, i);
 
		(*(p->func))(p->x, &(p->f), p->df);

		printf("Info: Energy = %14.10lf [H]\n", p->f);
		fflush(stdout);
	}
}

static double my_f(const gsl_vector *x, void *params)
{
  geom_oct_type *p = (geom_oct_type *)params;
  
	my_work_fdf(x, p);
	return p->f;
}

static void my_df(const gsl_vector *x, void *params, gsl_vector *df)
{
	int i;
	geom_oct_type *p = (geom_oct_type *)params;
 
	my_work_fdf(x, p);
	for(i=0; i<p->size; i++)
		gsl_vector_set(df, i, p->df[i]);
}

static void my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) 
{
	int i;
	geom_oct_type *p = (geom_oct_type *)params;
 
	my_work_fdf(x, p);
	*f = p->f;
	for(i=0; i<p->size; i++)
		gsl_vector_set(df, i, p->df[i]);
}

int F90_FUNC_(oct_geom_opt, OCT_GEOM_OPT)
		 (double *x, int *size, int *method, double *tol, int *max_iter,
			void (* cp) (double *, double *, double *))
{
	geom_oct_type c_geom_oct;
	size_t iter = 0;
  int i, status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
	gsl_multimin_function_fdf my_func;
  gsl_vector *vx;
	
	const gsl_multimin_fdfminimizer_type *min_type[] = {
		gsl_multimin_fdfminimizer_steepest_descent,
		gsl_multimin_fdfminimizer_conjugate_pr,
		gsl_multimin_fdfminimizer_conjugate_fr,
		gsl_multimin_fdfminimizer_vector_bfgs
	};

	c_geom_oct.size = *size;
	c_geom_oct.x  = (double *) malloc((*size)*sizeof(double));
	c_geom_oct.df = (double *) malloc((*size)*sizeof(double));
	c_geom_oct.func = *cp;

	/* setup starting point */
  vx    = gsl_vector_alloc (c_geom_oct.size);
	for(i = 0; i<c_geom_oct.size; i++){
		c_geom_oct.x[i] = 0.0;
		gsl_vector_set(vx, i, x[i]);
	}

	/* setup function */
  my_func.f      = &my_f;
  my_func.df     = &my_df;
  my_func.fdf    = &my_fdf;
  my_func.n      = c_geom_oct.size;
  my_func.params = (void *)(&c_geom_oct);
	
	T = min_type[*method - 1];
  s = gsl_multimin_fdfminimizer_alloc (T, c_geom_oct.size);
  gsl_multimin_fdfminimizer_set (s, &my_func, vx, 0.01, *tol);

	printf("Info: Using %s minimiser\n", gsl_multimin_fdfminimizer_name(s));

  do{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);

		if (status)
			break;
		status = gsl_multimin_test_gradient (s->gradient, *tol);

	}while (status == GSL_CONTINUE && iter < *max_iter);

	/* return minimum in x */
	for(i = 0; i<c_geom_oct.size; i++){
		x[i] = gsl_vector_get(s->x, i);
	}

	/* clean up */
	gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(vx);
	free(c_geom_oct.x);
	free(c_geom_oct.df);
	
	return (status == GSL_SUCCESS);
}
