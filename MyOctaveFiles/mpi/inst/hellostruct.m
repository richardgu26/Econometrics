## Copyright (C) 2009 Riccardo Corradini <riccardocorradini@yahoo.it>
## Copyright (C) 2012 Carlo de Falco
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} = hellostruct ()
## This function demonstrates sending and receiving a struct over MPI.
## Each process will send a a struct to process with rank 0, which will then display it.
## To run this example, set the variables HOSTFILE and NUMBER_OF_MPI_NODES to appropriate values, 
## then type the following command in your shell:
## @example
## mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave --eval 'pkg load mpi; hellostruct ()'
## @end example
## @seealso{hello2dimmat,helloworld,hellocell,hellosparsemat,mc_example,montecarlo,Pi} 
## @end deftypefn

function hellostruct ()

  MPI_Init ();
  ## the string NEWORLD is just a label could be whatever you want 
  CW = MPI_Comm_Load ("NEWORLD");

  my_rank = MPI_Comm_rank (CW);
  p = MPI_Comm_size (CW);

  ## TAG is very important to identify the message
  TAG = 1;

  if (my_rank != 0)
    message = struct ("f1", {1 3; 2 4}, "f2", 25);
    ## Could be a vector containing the list of ranks identifiers; 
    rankvect = 0
    info = MPI_Send (message, rankvect, TAG, CW);
  else
    for source = 1:p-1
      disp ("We are at rank 0 that is master etc..");
      [messagerec, info] = MPI_Recv (source, TAG, CW);
      messagerec
    endfor
  endif

  MPI_Finalize ();

endfunction

%!demo
%! system ("mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave -q --eval 'pkg load mpi; hellostruct ()'");
