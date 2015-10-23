/* @(#)utils.c
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

#include "utils.h"
#include "mex.h"
#include <stdlib.h>

void *
xmalloc_p(size_t size)
{
    void *data = mxMalloc(size);
    if (!data)
	mexErrMsgTxt("Out of memory in xmalloc_p");
    else
	mexMakeMemoryPersistent(data);
    return data;
}

void *
xcalloc_p(size_t nmemb, size_t size)
{
    void *data = mxCalloc(nmemb, size);
    if (!data)
	mexErrMsgTxt("Out of memory in xcalloc_p");
    else
	mexMakeMemoryPersistent(data);
    return data;
}

void
my_s_copy(char *a, char *b, int la, int lb)
{
    char *ae, *be;
    ae = a + la;
    if (la <= lb)
	while (a < ae)
	    *a++ = *b++;
    else
    {
	be = b + lb;
	while (b < be)
	    *a++ = *b++;
	while (a < ae)
	    *a++ = ' ';
    }
}
