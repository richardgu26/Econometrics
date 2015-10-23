## Copyright (C) 2006, 2009 Michael Creel <michael.creel@uab.es>
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
## @deftypefn {Function File} {[@var{n_received}]} = @
## montecarlo (@var{f}, @var{f_args}, @var{reps}, @var{outfile}, @var{n_pooled}, @var{n_returns}, @var{usempi}, @var{verbose})
## Generate a specified number of replications of a function's
## output and write them to a user-specified output file.
##
## IMPORTANT: @var{f} should return a row vector of output from feval (f, f_args)
##
## For normal evaluation on one core, only the first 4 arguments are required.
##
## @itemize @minus
## @item Arg 1: (required) the function that generates a row vector of output
## @item Arg 2: (required) the arguments of the function, in a cell
## @item Arg 3: (required) the number of replications to generate
## @item Arg 4: (required) the output file name
## @item Arg 5 (optional) number of replications to be pooled together between writes
## @item Arg 6 (optional) verbose: 1 for on, 0 for off
## @end itemize
##
## If using MPI, you should run using ranks equal to number of cores plus 1,
## and should make sure that the core running the frontend is also the one that
## has the second rank. That way the core the frontend is on will also do work.
## @end deftypefn

function n_received = montecarlo (f, f_args, reps, outfile, n_pooled, verbose)

  t0 = clock(); # initialize timing

  ## defaults for optional arguments
  if (nargin < 6) verbose = false; endif
  if (nargin < 5) n_pooled = 1; endif;

  if MPI_Initialized      # check if doing this parallel or serial
    use_mpi = true;
    CW = MPI_Comm_Load("NEWORLD");
    is_node = MPI_Comm_rank(CW);
    nodes = MPI_Comm_size(CW);
    mytag = 48;
  else
    use_mpi = false;
    is_node = 0;
  endif

  if is_node # compute nodes
    more_please = 1;
    while more_please
      for i = 1:n_pooled
        contrib = feval(f, f_args);
        contribs(i,:) = contrib;
      endfor
      MPI_Send(contribs, 0, mytag, CW);
      pause(0.05); # give time for the fronted to send a stop message, if done
                   # check if we're done
      if (MPI_Iprobe(0, is_node, CW)) # check for ping from rank 0
        junk = MPI_Recv(0, is_node, CW);
        break;
      endif
    endwhile
  else # frontend
    received = 0;
    done = false;
    while received < reps
      if use_mpi
        ## retrieve results from compute nodes
        for i = 1:nodes-1
          ## compute nodes have results yet?
          ready = false;
          ready = MPI_Iprobe(i, mytag, CW); # check if message pending
          if ready
            ## get it if it's there
            contribs = MPI_Recv(i, mytag, CW);
            need = reps - received;
            received = received + n_pooled;
            ## truncate?
            if n_pooled  >= need
              contribs = contribs(1:need,:);
              done = true;
            endif
            ## write to output file
            FN = fopen (outfile, "a");
            if (FN < 0) error ("montecarlo: couldn't open output file %s", outfile); endif
            t = etime(clock(), t0);
            for j = 1:rows(contribs)
              fprintf(FN, "%f ", i, t, contribs(j,:));
              fprintf(FN, "\n");
            endfor
            fclose(FN);
            if verbose printf("\nContribution received from node%d.  Received so far: %d\n", i, received); endif
            if done
              ## tell compute nodes to stop loop
              for j = 1:5
                for i = 1:(nodes-1)
                  if (j==1) MPI_Send(" ",i,i,CW); endif # send out message to stop
                  ready = MPI_Iprobe(i, mytag, CW); # get last messages
                  if ready contribs = MPI_Recv(i, mytag, CW); endif
                endfor
              endfor
              break;
            endif   
          endif
        endfor
      else
        for i = 1:n_pooled
          contrib = feval(f, f_args);
          contribs(i,:) = contrib;
        endfor
        need = reps - received;
        received = received + n_pooled;
        ## truncate?
        if n_pooled  >= need
          contribs = contribs(1:need,:);
        endif
        ## write to output file
        FN = fopen (outfile, "a");
        if (FN < 0) error ("montecarlo: couldn't open output file %s", outfile); endif
        t = etime(clock(), t0);
        for j = 1:rows(contribs)
          fprintf(FN, "%f ", 0, t, contribs(j,:));
          fprintf(FN, "\n");
        endfor
        fclose(FN);
        if verbose printf("\nContribution received from node 0.  Received so far: %d\n", received); endif
      endif
    endwhile
  endif
endfunction
