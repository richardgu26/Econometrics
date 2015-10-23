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
## @deftypefn {Function File} {} = hello2dimmat ()
## This function demonstrates sending and receiving of a 2-dimensional matrix over MPI.
## Each process in the pool will create a random 90x90 matrix and send it to process with rank 0.
## To run this example, set the variables HOSTFILE and NUMBER_OF_MPI_NODES to appropriate values, 
## then type the following command in your shell:
## @example
## mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave --eval 'pkg load mpi; hello2dimmat ()'
## @end example
## @seealso{hellocell,hellosparsemat,hellostruct,helloworld,mc_example,montecarlo,Pi} 
## @end deftypefn

function hello2dimmat ()

  MPI_SUCCESS = 0;
  MPI_Init ();
  
  ## the string NEWORLD is just a label could be whatever you want    
  CW = MPI_Comm_Load ("NEWORLD");
  my_rank = MPI_Comm_rank (CW);
  p = MPI_Comm_size (CW);
  mytag = 48;
  
  if (my_rank != 0)
    ## Generate a random matrix
    message = rand (90, 90);
    ## load message
    ## rankvect is the vector containing the list of rank  destination process
    rankvect = 0;
    [info] = MPI_Send (message, rankvect, mytag, CW);
  else
    for source = 1:p-1
      disp ("We are at rank 0 that is master etc..");
      [messager, info] = MPI_Recv (source, mytag, CW);
      
      ## You could also save each result and make comparisons if you don't trust MPI
      disp ("Rank 0 is the master receiving ... :");
      if (info == MPI_SUCCESS)
        disp ('OK!');
      endif
    endfor
  endif   

  MPI_Finalize ();

endfunction

%!demo
%! system ("mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave -q --eval 'pkg load mpi; hello2dimmat ()'");

