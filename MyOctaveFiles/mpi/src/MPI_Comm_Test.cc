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
DEFUN_DLD(MPI_Comm_Test, args, ,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} @var{DESCRIPTION} = MPI_Comm_Test (@var{COMM})\n\
Return @var{DESCRIPTION} string description of the MPI_Communicator  @var{COMM}.\n\
For\n\
example,\n\
\n\
@example\n\
@group\n\
MPI_Init();\n\
X = MPI_Comm_Load(\"description\"); \n\
whos X\n\
MPI_Comm_Test(X) \n\
@result{} \"description\"\n\
MPI_Finalize();\n\
@end group\n\
@end example\n\
@end deftypefn")
{
  octave_value_list results;
  if(args.length() != 1 
     || args(0).type_id () != simple::static_type_id ())
    {
      print_usage ();
      results(0) = octave_value (-1);
    }
  else
    {
      if (! simple_type_loaded)
        {
          simple::register_type ();
          simple_type_loaded = true;
          mlock ();
        }      
      const octave_base_value& rep = args(0).get_rep ();
      const simple& b = ((const simple &)rep);
      //octave_stdout << "MPI_Comm_Test has " << b.name_value()  << " output arguments.\n";
      MPI_Comm res = b.comunicator_value ();
      results(0) = octave_value (b.name_value ());
    }
  return results;
}
