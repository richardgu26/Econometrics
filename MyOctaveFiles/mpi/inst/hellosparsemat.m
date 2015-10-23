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
## @deftypefn {Function File} {} = hellosparsemat ()
## This function demonstrates sending and receiving a sparse matrix over MPI.
## Each process will send a a sparse matrix to process with rank 0, which will then display it.
## To run this example, set the variables HOSTFILE and NUMBER_OF_MPI_NODES to appropriate values, 
## then type the following command in your shell:
## @example
## mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave --eval 'pkg load mpi; hellosparsemat ()'
## @end example
## @seealso{hello2dimmat,helloworld,hellocell,hellostruct,mc_example,montecarlo,Pi} 
## @end deftypefn

function hellosparsemat ()

  MPI_Init();

  ## the string NEWORLD is just a label could be whatever you want  
  CW = MPI_Comm_Load ("NEWORLD");
  my_rank = MPI_Comm_rank (CW);
  p = MPI_Comm_size (CW);

  ## tag[0] ----> type of octave_value
  ## tag[1] ----> array of three elements 1) num of rows 2) number of columns 3) number of non zero elements
  ## tag[2] ---->  vector of rowindex
  ## tag[3] ---->  vector of columnindex
  ## tag[4] ---->  vector of  non zero elements
  ## These tags will be generated after mytag by the MPI_Send and MPI_Recv (see source code)
  mytag = 48;
  
  ## Fill the sparse matrix
  M = 5;  N = 5;  D = 0.9;

  for one_by_one = p-1:-1:1 ## work one cpu at a time to make the output readable

    if (my_rank != 0)
      
      message = sprand (M, N, D);
      dest = 0;
      info = MPI_Send (message, dest, mytag, CW);
      printf (["on rank %d MPI_Send returned the following ", ...
              "error code (0 = Success): info = %d\n"], my_rank, info)

    endif

    MPI_Barrier (CW);

    if (my_rank == 0)
      
      source = one_by_one;
      [messager, info] = MPI_Recv (source, mytag, CW);
        
      printf (["MPI_Recv returned the following error code (0 = Success)", ...
               " while receving from rank %d : info = %d\n"], source, info);

        
      printf ("This is the matrix received from rank %d: \n", source);
      disp (full (messager))
      
    endif   
    
    MPI_Barrier (CW);
  endfor
  
  MPI_Finalize ();
  
endfunction

%!demo
%! system ("mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave -q --eval 'pkg load mpi; hellosparsemat ()'");
