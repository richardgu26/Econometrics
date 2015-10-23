## Copyright (C) 2013 Michael Creel
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## WriteBinaryData
##
## usage: data = WriteBinaryData(data, 'filename');
##
## Writes data in single binary form to an output file.
## If the output file does not exist, the record size is 
## written as the first element, then the data. Otherwise,
## the data is appended to the existing file. There is a 
## check that the new data has the same record size as the
## existing data file, if they differ, the output file is not
## modified. Use ReadBinaryData to load data written by this
## function.

## Author: Michael Creel <michael@pelican>
## Created: 2013-09-04

function ok = WriteBinaryData(data, filename)
	% if the file doesn't exist, write record size as first entry
	if (exist(sprintf("./%s",filename)) !=2)
		recordsize = columns(data);
		FN = fopen (filename, "ab");
		if (FN < 0) error ('WriteBinaryData: could not open output file %s', filename);
		else
			count = fwrite(FN, recordsize, 'double');
			fclose(FN);
			if count==1;
				printf('successfully wrote record size %d to %s\n', recordsize, filename);
			else
				printf('error writing record size to output file %s\n', filename);	
			endif
		endif	
	else
		FN=fopen(filename, 'r');
		recordsize = fread(FN, 1);
		fclose(FN);
	endif	
	
	if recordsize != columns(data)
		printf('existing output file has a different record size than the new data, aborting the write\n');
		ok = 0;	
	else
		FN = fopen (filename, "ab");
		if (FN < 0) error ("WriteBinaryData: couldn't open output file %s", filename);
		else	
			count = fwrite(FN, data);
			fclose(FN);
		endif	
		if (count == recordsize*rows(data));
			printf('successfully wrote %d records of size %d to output file %s\n', rows(data), recordsize, filename);
			ok = 1;
		else	
			error('WriteBinaryData: error writing data to output file %s\n', filename);
			ok = 0;	
		endif
	endif	
endfunction
