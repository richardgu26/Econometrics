/* @(#)lbfgs_mex.c
 *
 * Copyright 2005 Liam Stewart
 *
 * This file is part of the MATLAB LBFGS wrapper.
 *
 * The MATLAB LBFGS wrapper is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * The MATLAB LBFGS wrapper is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Foobar; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "mex.h"
#include "utils.h"
#include "list.h"

#include <string.h>

#define LBFGS_STRLEN		60

struct lbfgs_opts
{
    int m;
    double factr;
    double pgtol;
    int iprint;
};

typedef struct lbfgs_opts lbfgs_opts_t;

struct lbfgs
{
    int n, m;
    double factr, pgtol;
    int iprint;
    double f, *x, *g;
    double *l, *u;
    int *nbd;
    char task[LBFGS_STRLEN];
    char csave[LBFGS_STRLEN];
    double *wa;
    int *iwa;
    int lsave[4];
    int isave[44];
    double dsave[29];
};

typedef struct lbfgs lbfgs_t;

#define SETULB setulb_
int SETULB(int* n, int* m,
	   double* x, 
	   double* l, double* u, int* nbd,
	   double* f, double* g,
	   double* factr, double* pgtol,
	   double* wa,  int* iwa,
	   char* task,  int* iprint,
	   char* csave, int* lsave, int* isave, double* dsave);

static list_t *list = NULL;
static int mexAtExitIsSet = 0;

lbfgs_t *lbfgs_init(int n, double *x0, double *lbd, double *ubd, double *nbd, lbfgs_opts_t opts);
void lbfgs_destroy(lbfgs_t *l);
void lbfgs_exec(lbfgs_t *l);
lbfgs_t *lbfgs_verify(lid_t id);

void lbfgs_mex_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void lbfgs_mex_step(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void lbfgs_mex_stop(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void lbfgs_mex_destroy(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void exitMex(void);

void
mexFunction(int nlhs, mxArray *plhs[],
	    int nrhs,const mxArray *prhs[])
{
    char *mode = NULL;
    void (*opMode)(int nlhs,
                   mxArray *plhs[],
                   int nrhs,
                   const mxArray *prhs[]) = NULL;

    if ((nrhs < 2) || (nrhs > 7))
        mexErrMsgTxt("Invalid number of input arguments.");

    if (!mexAtExitIsSet)
    {
	list = list_init();
	mexAtExit(exitMex);
	mexAtExitIsSet = 1;
    }

    if (mxIsChar(prhs[0])) 
        mode = mxArrayToString(prhs[0]);
    else
        mexErrMsgTxt("First input to lbfgs_mex must be a string.");

    if (mode == NULL)
        mexErrMsgTxt("Memory allocation failed.");
        
    if (!strcmp(mode,"init"))
        opMode = lbfgs_mex_init;
    else if (!strcmp(mode,"step"))
        opMode = lbfgs_mex_step;
    else if (!strcmp(mode,"stop"))
        opMode = lbfgs_mex_stop;
    else if (!strcmp(mode,"destroy"))
	opMode = lbfgs_mex_destroy;
    else
	mexErrMsgTxt("Unrecognized mode for lbfgs_mex.");

    mxFree(mode);
                
    (*opMode)(nlhs, plhs, nrhs, prhs);
}

void
lbfgs_mex_init(int nlhs, mxArray *plhs[],
	       int nrhs, const mxArray *prhs[])
{
    lid_t id;
    double *nd;
    int n, nf, i;
    double *x0, *lbd, *ubd, *nbd;
    lbfgs_opts_t opts;
    lbfgs_t *l;

    /* Validate inputs */
    if (nrhs != 7)
	mexErrMsgTxt("Invalid number of inputs to lbfgs_init.");

    nd = mxGetPr(prhs[1]);
    n = (int) (*nd);
    x0 = mxGetPr(prhs[2]);

    lbd = mxGetPr(prhs[3]);
    ubd = mxGetPr(prhs[4]);
    nbd = mxGetPr(prhs[5]);

    if (!mxIsStruct(prhs[6]))
	mexErrMsgTxt("Argument 7 must be a structure.");
    nf = mxGetNumberOfFields(prhs[6]);

    /* some defaults */
    opts.m = 5;
    opts.factr = 1e7;
    opts.pgtol = 1e-5;
    opts.iprint = 1;

    if (!mxIsEmpty(prhs[6]))
    {
	for (i = 0; i < nf; i++)
	{
	    const char *fname = mxGetFieldNameByNumber(prhs[6], i);
	    mxArray *field = mxGetFieldByNumber(prhs[6], 0, i);

	    /* ignore unknown fields */
	    if (!strcmp(fname, "m"))
		opts.m = (int)(mxGetScalar(field));
	    else if (!strcmp(fname, "factr"))
		opts.factr = mxGetScalar(field);
	    else if (!strcmp(fname, "pgtol"))
		opts.pgtol = mxGetScalar(field);
	    else if (!strcmp(fname, "iprint"))
		opts.iprint = (int)(mxGetScalar(field));
	}
    }

    l = lbfgs_init(n, x0, lbd, ubd, nbd, opts);

    /* keep track of structures */
    id = list_add(list, (void*) l);

    /* return id */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    *((double *) mxGetPr(plhs[0])) = id;
}

void
lbfgs_mex_step(int nlhs, mxArray *plhs[],
	       int nrhs, const mxArray *prhs[])
{
    lid_t id;
    lbfgs_t *l;
    double f;
    double *g, *x;
    double rv;			/* FIXME: should be int */

    /* verify inputs */
    if (nrhs != 4)
	mexErrMsgTxt("Invalid number of inputs to lbfgs_destroy.");

    id = (lid_t) (*mxGetPr(prhs[1]));
    l = lbfgs_verify(id);

    f = mxGetScalar(prhs[2]);
    g = mxGetPr(prhs[3]);

    /* FIXME: verify size of g */

    l->f = f;
    memcpy(l->g, g, l->n * sizeof(double));

    /* step */
    lbfgs_exec(l);

    /* create outputs */
    if(!strncmp(l->task, "FG", 2))
	rv = 0;
    else if (!strncmp(l->task, "NEW_X", 5))
	rv = 1;
    else if (!strncmp(l->task, "CONV", 4))
	rv = 2;
    else if (!strncmp(l->task, "ABNO", 4))
	rv = 3;
    else if (!strncmp(l->task, "ERROR", 4))
	rv = 4;
    else
	rv = 5;

    plhs[0] = mxCreateDoubleMatrix(l->n, 1, mxREAL);
    x = mxGetPr(plhs[0]);
    memcpy(x, l->x, l->n * sizeof(double));

    plhs[1] = mxCreateScalarDouble(rv);
}

void
lbfgs_mex_stop(int nlhs, mxArray *plhs[],
	       int nrhs, const mxArray *prhs[])
{
    lid_t id;
    lbfgs_t *l;

    if (nrhs != 2)
	mexErrMsgTxt("Invalid number of inputs to lbfgs_destroy.");

    id = (lid_t) *mxGetPr(prhs[1]);
    l = lbfgs_verify(id);

    my_s_copy(l->task, "STOP", LBFGS_STRLEN, 4);

    lbfgs_exec(l);
}

void
lbfgs_mex_destroy(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
    lid_t id;

    if (nrhs != 2)
	mexErrMsgTxt("Invalid number of inputs to lbfgs_destroy.");

    id = (lid_t) *mxGetPr(prhs[1]);
    lbfgs_verify(id);
    list_delete(list, id, (list_item_dtor)lbfgs_destroy);
}

lbfgs_t *
lbfgs_init(int n, double *x0, double *lbd, double *ubd, double *nbd, 
	   lbfgs_opts_t opts)
{
    lbfgs_t *l = XMALLOC_P(lbfgs_t, 1);
    int m = opts.m;
    int i;
    
    l->n = n;
    l->m = m;
    l->factr = opts.factr;
    l->pgtol = opts.pgtol;
    l->iprint = opts.iprint;
    l->f = 0;

    l->x = XMALLOC_P(double, n);
    memcpy(l->x, x0, n*sizeof(double));

    l->g = XCALLOC_P(double, n);

    l->l = XMALLOC_P(double, n);
    memcpy(l->l, lbd, n*sizeof(double));
    
    l->u = XMALLOC_P(double, n);
    memcpy(l->u, ubd, n*sizeof(double));
    
    l->nbd = XMALLOC_P(int, n);
    for (i = 0; i < n; i++)
	l->nbd[i] = (int)nbd[i];

    l->wa = XMALLOC_P(double, 2*m*n+4*n+12*m*m+12*m);
    l->iwa = XMALLOC_P(int, 3*n);

    my_s_copy(l->task, "START", LBFGS_STRLEN, 5);

    return l;
}

void
lbfgs_destroy(lbfgs_t *l)
{
    if (!l)
	return;

    XFREE(l->x);
    XFREE(l->g);
    XFREE(l->l);
    XFREE(l->u);
    XFREE(l->nbd);
    XFREE(l->wa);
    XFREE(l->iwa);
    XFREE(l);
}

void
lbfgs_exec(lbfgs_t *l)
{
    if (l == NULL)
	return;

    SETULB(&(l->n),
	   &(l->m),
	   l->x,
	   l->l, l->u, l->nbd,
	   &(l->f), l->g,
	   &(l->factr), &(l->pgtol),
	   l->wa, l->iwa,
	   l->task,
	   &(l->iprint),
	   l->csave, l->lsave, l->isave, l->dsave);
}

lbfgs_t *
lbfgs_verify(lid_t id)
{
    void *data = list_find(list, id);

    if (data == NULL)
	mexErrMsgTxt("The given identifier is not valid.");

    return (lbfgs_t *)data;
}

void
exitMex(void)
{
    if (list == NULL)
	return;
    if (list->head != NULL)
	mexWarnMsgTxt("Clearing all structures. It is no longer possible to use any of them.");
    list_destroy(list, (list_item_dtor)lbfgs_destroy);
    list = NULL;
}
