// Copyright (C) 2004-2007 Javier Fernández Baldomero, Mancia Anguita López
// Copyright (C) 2009 Riccardo Corradini <riccardocorradini@yahoo.it>
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

/*
 * ----------------------------------------------------
 * Initialize the MPI execution environment
 * info = MPI_Init [ ( 'arg' [, 'arg']... ) ]
 * ----------------------------------------------------
 */

#define NAME MPI_Init

#include "mpi.h"        // mpi.h, oct.h
#include <octave/oct.h>

DEFUN_DLD(NAME, args, nargout,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} @var{INFO} = MPI_Init()\n\
Initialize the MPI execution environment.\n\
\n\
 @example\n\
 @group\n\
    @var{INFO} (int) return code\n\
       0 MPI_SUCCESS    No error\n\
      16 MPI_ERR_OTHER  Attempt was made to call MPI_Init a  second  time\n\
                       MPI_Init may only be called once in a program\n\
                       \n\
SEE ALSO: MPI_Finalize, MPI_Initialized, MPI_Finalized\n\
@end group\n\
@end example\n\
@end deftypefn")
{
  int nargin = args.length();            
  for (int i = 0; i < nargin; i++)
    {
      if (! args(i).is_string ()) 
        {
          error ("MPI_Init: args must be strings");
          return octave_value (MPI_ERR_ARG);    // error returns nothing
        }
    }

    string_vector argvec = args.make_argv ("MPI_Init");
    char **argve= argvec.c_str_vec ();
    char **argv = &argve[1];
    
    int info = MPI_Init (&nargin, &argv);
    free(argve);
    return octave_value (info);
}
 
