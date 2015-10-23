## Copyright (C) 2004-2007 Javier Fernández Baldomero, Mancia Anguita López
## Copyright (C) 2009 Riccardo Corradini <riccardocorradini@yahoo.it>
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
## @deftypefn {Function File} {[@var{result}]} = Pi ()
## Classical MPI example that computes @var{PI} by integrating arctan'(x) in [0,1].
##  @var{N} [1e7] #subdivisions of the [0, 1] interval.
##  @var{mod} ['s'] communication modality:  (s)end (r)educe.
##  @var{results} struct contains
##
##  @itemize @minus
##      @item @var{pi}: estimated pi value
##      @item @var{err}: error
##      @item @var{time}: from argument xmit to pi computed
## @end itemize
##
## To run this example, set the variables HOSTFILE and NUMBER_OF_MPI_NODES to appropriate values, 
## then type the following command in your shell:
## @example
## mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave --eval 'pkg load mpi; Pi ()'
## @end example
## @seealso{hello2dimmat,helloworld,hellocell,hellosparsemat,mc_example,montecarlo,hellostruct} 
## @end deftypefn

function results = Pi (N, mod)
        
  ############
  ## ArgChk ##
  ############
  if (nargin < 1)
    N = 1E7;
  end
  if (nargin < 2)
    mod = 's';
  end
  if (nargin > 2)
    print_usage (); # let all ranks complain
  end 
  flag = 0; # code much simpler
  flag = flag || ~isscalar (N) || ~isnumeric (N);
  flag = flag || fix(N) ~= N || N < 1;
  mod = lower (mod); mods = 'sr';
  flag = flag || isempty (findstr (mod, mods)); 
  if (flag)
    print_usage (); # let them all error out
  end

  ####################
  ## Results struct ##
  ####################
  results.pi = 0;
  results.err = 0;
  results.time = 0;

  ####################################################################
  ## PARALLEL: initialization, include MPI_Init time in measurement  ##
  ####################################################################

  T=clock; #
  ############
  MPI_ANY_SOURCE = -1;
  MPI_Init ();  
  MPI_COMM_WORLD = MPI_Comm_Load ("NEWORLD");           
  rnk   = MPI_Comm_rank (MPI_COMM_WORLD); # let it abort if it fails
  siz   = MPI_Comm_size (MPI_COMM_WORLD);

  SLV = logical(rnk);                 # handy shortcuts, master is rank 0
  MST = ~ SLV;                        # slaves are all other

  #####################################################
  ## PARALLEL: computation (depends on rank/size)    ##
  #####################################################

  width=1/N; lsum=0;            # for i=rnk:siz:N-1
  i=rnk:siz:N-1;                #   x=(i+0.5)*width;
  x=(i+0.5)*width;              #   lsum=lsum+4/(1+x^2);
  lsum=sum(4./(1+x.^2));        # end

  #####################################
  ## PARALLEL: reduction and finish  ##
  #####################################

  switch mod
    case 's',                     
      TAG=7;                              # Any tag would do
      if SLV                              # All slaves send result back
        MPI_Send (lsum, 0, TAG, MPI_COMM_WORLD);
      else                                # Here at master
        Sum =lsum;                        # save local result
        for slv=1:siz-1                   # collect in any order
          lsum = MPI_Recv (MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD);
          Sum += lsum;                    # and accumulate
        endfor                            # order: slv or MPI_ANY_SOURCE
      endif
    case 'r',
      disp ("not yet implemented");
  endswitch

  MPI_Finalize ();

  if MST
    Sum      = Sum/N ;                      # better at end: don't loose resolution
    #################################       # stopwatch measurement
    results.time = etime(clock,T);  #       # but only at master after PI computed
    #################################       # all them started T=clock;
    results.err  = Sum - pi;
    results.pi   = Sum # ;
  endif

endfunction

%!demo
%! system ("mpirun --hostfile $HOSTFILE -np $NUMBER_OF_MPI_NODES octave -q --eval 'pkg load mpi; Pi ()'");
