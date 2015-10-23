% Copyright (C) 2003,2004  Michael Creel <michael.creel@uab.es>
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

%% this prints matrices with row and column labels
function prettyprint(mat, rlabels, clabels)

	% left pad the column labels 
	a = size(rlabels,2);
	for i = 1:a
		printf(' ');
	end
	printf('  ');

	% print the column labels
	clabels = ['          ';clabels]; % pad to 8 characters wide
	clabels = strjust(clabels,'right');

	k = columns(mat);
	for i = 1:k
		printf('%s  ',clabels(i+1,:));
	end

	% now print the row labels and rows
	fprintf('\n');
	k = size(mat,1);
	for i = 1:k
		if ischar(rlabels(i,:))
			fprintf(rlabels(i,:));
		else
			fprintf('%i', rlabels(i,:));
		end
		fprintf('  %10.5f', mat(i,:));
		fprintf('\n');
	end
end
