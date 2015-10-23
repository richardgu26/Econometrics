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

#define NAME MPI_Barrier
/*
 * ----------------------------------------------------
 * Blocks until all processes in the communicator have reached this routine
 * [info ] = MPI_Barrier (comm)
 * ----------------------------------------------------
 */

#include "simple.h"    

DEFUN_DLD (NAME, args, ,
"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} @var{INFO} = MPI_Barrier (@var{COMM})\n\
Block processes in the communicator @var{COMM} untill the call to this routine has been reached by all.\n\
The communicator object @var{COMM} (loaded with MPI_Comm_Load) must be specified otherwise \n\
an error occurs. \n\
 @example\n\
 @group\n\
    @var{INFO} (int) return code\n\
       0 MPI_SUCCESS    No error\n\
       5 MPI_ERR_COMM   Invalid communicator (NULL?)\n\
      13 MPI_ERR_ARG    Invalid argument (typically a NULL pointer?)\n\
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

      if((args.length () != 1)
         || args(0).type_id () != simple::static_type_id ())
        {
          error ("MPI_Barrier: Please enter octave comunicator object");
          results(0) = octave_value (-1);
        }
      else
        {	       
          const octave_base_value& rep = args(0).get_rep();
          const simple& B = ((const simple &)rep);
          MPI_Comm comm = ((const simple&) B).comunicator_value ();
          if (! error_state)
            {
              int my_size;
              int info = MPI_Barrier (comm);
              
              results(0) = octave_value (info);
            }
          else
            print_usage ();
        }
    }
   
  return results;
}

