/* @(#)list.h
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

#ifndef LBFGS_MEX_H
#define LBFGS_MEX_H 1

typedef unsigned int lid_t;

struct list_item 
{
    void *data;
    lid_t id;
    struct list_item *next;
};

struct list
{
    struct list_item *head;
    lid_t next_id;
};

typedef struct list_item list_item_t;
typedef struct list list_t;

typedef void(*list_item_dtor)(void *);

list_t *list_init();
lid_t list_add(list_t *l, void *data);
void *list_find(list_t *l, lid_t id);
void *list_remove(list_t *l, lid_t id);
void list_delete(list_t *l, lid_t id, list_item_dtor dtor);
void list_clean(list_t *l, list_item_dtor dtor);
void list_destroy(list_t *l, list_item_dtor dtor);

#endif /* _LBFGS_MEX_H */

