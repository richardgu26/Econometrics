// Copyright (C) 2004-2007 Javier Fernández Baldomero, Mancia Anguita López
// Copyright (C) 2009 Riccardo Corradini <riccardocorradini@yahoo.it>
// Copyright (C) 2009 VZLU Prague
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

#define NAME MPI_Comm_size
/*
 * ----------------------------------------------------
 * Determines the size of the calling process in the communicator
 * [info rank] = MPI_Comm_size (comm)
 * ----------------------------------------------------
 */

#include "simple.h"  
DEFUN_DLD(NAME, args,nargout ,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [@var{RANK} @var{INFO}] = MPI_Comm_size (@var{COMM})\n\
Return the size of a communicator.\n\
The comunicator object @var{COMM} (loaded with MPI_Comm_Load) must be sepcified \n\
or otherwise an error occurs. \n\
 @example\n\
 @group\n\
    @var{exprank} rank of the calling process in group of communicator\n\
    @var{exprinfo} (int) return code\n\
       0 MPI_SUCCESS    No error\n\
       5 MPI_ERR_COMM   Invalid communicator (NULL?)\n\
      13 MPI_ERR_ARG    Invalid argument (typically a NULL pointer?)\n\
SEE ALSO: MPI_Comm_rank\n\
@end group\n\
@end example\n\
@end deftypefn")

{
  octave_value_list results;
  int nargin = args.length ();
  if (nargin != 1)
    print_usage ();
  else
    {
      if (! simple_type_loaded)
        {
          simple::register_type ();
          simple_type_loaded = true;
          mlock ();
        }

      if ((args.length() != 1 )
          || args(0).type_id () != simple::static_type_id ())
        {		
          error ("MPI_Comm_size: Please enter octave comunicator object");
          results(0) = octave_value (-1);
	}
      else
        {
          const octave_base_value& rep = args(0).get_rep ();
          const simple& B = ((const simple &)rep);
          MPI_Comm comm = ((const simple&) B).comunicator_value ();
          
          if (! error_state)
            {
              int my_size;
              int info = MPI_Comm_size (comm, &my_size);
              if (nargout > 1)
                results(1) = info;
              results(0) = my_size;
            }
          else
            print_usage ();
        }
      // comm = NULL;
      /* [size info] = MPI_Comm_size (comm) */   
    }
  return results;
}

