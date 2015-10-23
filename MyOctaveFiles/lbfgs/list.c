/* @(#)list.c
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

#include "list.h"
#include "utils.h"

#include <assert.h>

list_t *
list_init()
{
    list_t *l = XMALLOC_P(list_t, 1);

    l->next_id = 0;
    l->head = NULL;

    return l;
}

lid_t
list_add(list_t *l, void *data)
{
    list_item_t *li;

    assert(l != NULL);
    assert(data != NULL);

    li = XMALLOC_P(list_item_t, 1);
    li->data = data;
    li->next = NULL;
    li->id = l->next_id++;

    if (l->head == NULL)
    {
	l->head = li;
    }
    else
    {
	li->next = l->head;
	l->head = li;
    }

    return li->id;
}

void *
list_find(list_t *l, lid_t id)
{
    list_item_t *p;
    
    assert(l != NULL);

    p = l->head;
    while (p != NULL)
    {
	if (p->id == id)
	    return p->data;
	p = p->next;
    }
    return NULL;
}

void *
list_remove(list_t *l, lid_t id)
{
    void *data;
    list_item_t *p, *q;

    assert(l != NULL);

    data = NULL;
    p = l->head, q = NULL;
    while (p != NULL)
    {
	if (p->id == id)
	{
	    if (q != NULL)
		q->next = p->next;
	    else		/* p is head */
		l->head = p->next;
	    p->next = NULL;
	    data = p->data;
	    XFREE(p);
	    break;
	}
	else
	{
	    q = p;
	    p = p->next;
	}
    }

    return data;
}

void
list_delete(list_t *l, lid_t id, void(*dtor)(void *))
{
    void *data;
    
    assert(l != NULL);
    
    data = list_remove(l, id);
    if (data != NULL && dtor != NULL)
	dtor(data);
}

void
list_clean(list_t *l, void(*dtor)(void *))
{
    list_item_t *p;
    
    assert(l != NULL);

    while (l->head != NULL)
    {
	p = l->head;
	l->head = p->next;

	assert(p->data != NULL);
	if (dtor != NULL)
	    dtor(p->data);
	XFREE(p);
    }

    l->next_id = 0;
}

void
list_destroy(list_t *l, void(*dtor)(void *))
{
    list_clean(l, dtor);
    XFREE(l);
}
