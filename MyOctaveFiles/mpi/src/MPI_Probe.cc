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

#define NAME MPI_Probe
/*
 * ----------------------------------------------------
 * Blocking test for a message
 * [info stat] = MPI_Probe (src, tag, comm)
 * ----------------------------------------------------
 */
#include "simple.h"
#include <octave/ov-struct.h>

octave_scalar_map put_MPI_Stat (MPI_Status &stat)
{
  /*---------------------------------------------*/
  octave_scalar_map map;
  octave_value tmp = stat.MPI_SOURCE;
  map.assign ("src", tmp);
  tmp = stat.MPI_TAG;
  map.assign ("tag", tmp );
  tmp = stat.MPI_ERROR;
  map.assign ("err", tmp );
  int itmp;
  MPI_Get_count (&stat, MPI_CHAR, &itmp);
  map.assign ("cnt", itmp);

  MPI_Test_cancelled (&stat, &itmp);
  map.assign ("can", itmp);

  return map;
}

DEFUN_DLD(NAME, args, nargout,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [@var{STAT} @var{INFO}] = MPI_Probe(@var{SRCRANK}, @var{TAG}, @var{COMM})\n \
Blocking test for a message.\n\
 @example\n\
 @group\n\
 \n\
     @var{STAT} struct object\n\
       src (int)       source rank for the accepted message\n\
       tag (int)       message tag for the accepted message\n\
       err(int)        error \n\
       cnt (int)       count\n\
       can (int)       cancel\n\
    @var{INFO} (int) return code\n\
      0 MPI_SUCCESS    No error\n\
     13 MPI_ERR_ARG    Invalid argument\n\
      5 MPI_ERR_COMM   Invalid communicator (null?)\n\
      4 MPI_ERR_TAG    Invalid tag argument (MPI_ANY_TAG, 0..MPI_TAG_UB attr)\n\
      6 MPI_ERR_RANK   Invalid src/dst rank (MPI_ANY_SOURCE, 0..Comm_size-1)\n\
 @end group\n\
 @end example\n\
 \n\
@seealso{MPI_Iprobe, MPI_Recv, and MPI documentation for C examples}\n\
@end deftypefn")
{
  octave_value_list results;
  int nargin = args.length ();
  if (nargin != 3)
    print_usage ();
  else
    {
      if (!simple_type_loaded)
        {
          simple::register_type ();
          simple_type_loaded = true;
          mlock ();
        }

      if (args(2).type_id () == simple::static_type_id ())
        {
          const octave_base_value& rep = args(2).get_rep ();
          const simple& B = ((const simple &)rep);
          MPI_Comm comm = ((const simple&) B).comunicator_value ();
          int src = args(0).int_value ();    
          int tag = args(1).int_value ();    
          
          if (! error_state)
            {
              MPI_Status stat = {0, 0, 0, 0};
              int info = MPI_Probe (src, tag, comm, &stat);
              comm= NULL;
              results(0) = put_MPI_Stat (stat);
              results(1) = info;
            }
        }
      else
        {
          print_usage ();
          results = octave_value (-1);
        }

      return results;
      
      /* [ stat info ] = MPI_Probe (src, tag, comm) */
    }
}


