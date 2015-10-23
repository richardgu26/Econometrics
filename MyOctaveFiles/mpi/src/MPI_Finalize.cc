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

#define NAME MPI_Finalize

/*
 * ----------------------------------------------------
 * Terminates MPI execution environment
 * info = MPI_Finalize
 * ----------------------------------------------------
 */

#include "mpi.h"       
#include <octave/oct.h>

DEFUN_DLD(NAME, args, nargout,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} @var{INFO} = MPI_Finalize ()\n\
Terminate the MPI execution environment.\n\
\n\
 @example\n\
 @group\n\
    @var{INFO} (int) return code\n\
       0 MPI_SUCCESS    No error\n\
       5 MPI_ERR_COMM   Invalid communicator (NULL?)\n\
      13 MPI_ERR_ARG    Invalid argument (typically a NULL pointer?)\n\
SEE ALSO: MPI_Init\n\
@end group\n\
@end example\n\
@end deftypefn")
{
  int info = MPI_Finalize ();   
  return octave_value (info);
}
