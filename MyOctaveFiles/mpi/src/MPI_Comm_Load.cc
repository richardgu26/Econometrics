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

#include "simple.h"

DEFUN_DLD(MPI_Comm_Load, args, ,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} @var{COMM} = MPI_Comm_Load (@var{DESCRIPTION})\n\
Return @var{COMM} the MPI_Communicator object whose description is  @var{DESCRIPTION}, as a string.\n\
The default value will be MPI_COMM_WORLD. \n\
If @var{DESCRIPTION} is omitted, return anyway an MPI_COMM_WORLD comunicator object \n\
with no decription.\n\
For\n\
example,\n\
\n\
@example\n\
@group\n\
MPI_Init();\n\
X = MPI_Comm_Load(\"description\"); \n\
whos X\n\
MPI_Finalize();\n\
@end group\n\
@end example\n\
@end deftypefn")
{
  octave_value retval;
  if (args.length () != 1 || ! args (0).is_string ())
    print_usage ();
  else
    {
      if (!simple_type_loaded)
        {
          simple::register_type ();
          simple_type_loaded = true;
          mlock ();
        }
      const std::string name = args (0).string_value ();
      retval = new simple (name,MPI_COMM_WORLD);
    }  
  return retval;
}
