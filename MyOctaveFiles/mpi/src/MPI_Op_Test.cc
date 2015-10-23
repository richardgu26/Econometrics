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

#include "simpleop.h"

DEFUN_DLD(MPI_Op_Test, args, ,"")
{

  if (! simpleop_type_loaded)
    {
      simpleop::register_type ();
      simpleop_type_loaded = true;
      mlock ();
    }

  octave_value retval;
  if (args.length () != 1 
      || args(0).type_id() != simpleop::static_type_id ()
      || error_state)
    {
      print_usage ();
      retval = octave_value (-1);
    }
  else
    {
      const octave_base_value& rep = args(0).get_rep();
      const simpleop& b = ((const simpleop &)rep);
      octave_stdout << "simpleoptest has " << b.name_value()  << " output arguments.\n";
      MPI_Op res = b.operator_value();
    }
  return retval;
}
