// Copyright (C) 2009 Riccardo Corradini <riccardocorradini@yahoo.it>
// Copyright (C) 2009 VZLU Prague
// Copyright (C) 2012, 2013, 2014 Carlo de Falco
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
#include <ov-cell.h>    // avoid errmsg "cell -- incomplete datatype"
#include <oct-map.h>    // avoid errmsg "Oct.map -- invalid use undef type"


// forward declarations
int 
recv_class (MPI_Comm comm, octave_value &ov, int source, int mytag);

int 
recv_cell (MPI_Comm comm, octave_value &ov, int source, int mytag);

int 
recv_struct (MPI_Comm comm, octave_value &ov, int source, int mytag);

int 
recv_string (MPI_Comm comm, octave_value &ov, int source, int mytag);

int
recv_range (MPI_Comm comm, octave_value &ov, int source, int mytag);

template<class AnyElem>
int 
recv_vec (MPI_Comm comm, AnyElem &LBNDA, int nitem, MPI_Datatype TRCV, 
          int source, int mytag);

int 
recv_matrix (bool is_complex, MPI_Datatype TRcv, MPI_Comm comm, 
             octave_value &ov, int source, int mytag);

int
recv_sp_mat (bool is_complex, MPI_Datatype TRcv, MPI_Comm comm, 
             octave_value &ov, int source, int mytag);

template <class Any>
int recv_scalar (MPI_Datatype TRcv, MPI_Comm comm, Any *d, int source, 
                 int mytag);

template <class Any>
int
recv_scalar (MPI_Datatype TRcv, MPI_Comm comm, std::complex<Any> *d, 
             int source, int mytag);


int
recv_range (MPI_Comm comm, octave_value &ov, int source, int mytag)
{

  // Send/Recv as base, limit, incr, nelem
  // just 3 doubles + 1 int   
  // octave_range (double base, double limit, double inc)

  OCTAVE_LOCAL_BUFFER(int, tanktag, 3);
  OCTAVE_LOCAL_BUFFER(double, d, 2);

  MPI_Status stat;
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;


  // first receive
  int info = MPI_Recv (d, 3, MPI_DOUBLE, source, tanktag[1], comm, &stat);
  if (info == MPI_SUCCESS) 
    {
      int nelem = 0;
      info = MPI_Recv 
        (&nelem, 1, MPI_INT, source, tanktag[2], comm, &stat);

      if (info == MPI_SUCCESS) 
        {
          Range r (d[0], d[2], nelem); 
          ov = r;
        }
    }
  return info;
}

// This will get the fortran_vec vector for Any type Octave can handle
template<class AnyElem>
int
recv_vec (MPI_Comm comm, AnyElem &LBNDA, int nitem, MPI_Datatype TRCV, 
          int source, int mytag)
{

  MPI_Status stat;
  int info = MPI_Recv (LBNDA, nitem, TRCV, source, mytag, comm, &stat);
  return (info);
}

// template specialization for complex case
template <class Any>
int 
recv_scalar (MPI_Datatype TRcv, MPI_Comm comm, std::complex<Any> &d, 
             int source, int mytag)
{
  int info;
  MPI_Status stat;
  OCTAVE_LOCAL_BUFFER(int, tanktag, 2);
  tanktag[0] = mytag;
  tanktag[1] = mytag+1;
  OCTAVE_LOCAL_BUFFER(std::complex<Any>, Deco, 2);
  Deco[0] = real (d);
  Deco[1] = imag (d);
		
  info = MPI_Recv (&Deco, 2, TRcv, source, tanktag[1], comm, &stat);

  return info;
}

template <class Any>
int 
recv_scalar (MPI_Datatype TRcv, MPI_Comm comm, Any &d, int source, int mytag)
{
  
  // it's just a value directly MPI_Recv it
 
  OCTAVE_LOCAL_BUFFER(int, tanktag, 2);
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  int info;

  MPI_Status stat;
  info = MPI_Recv (&d, 1, TRcv, source, tanktag[1], comm, &stat);
  return info;
}

int
recv_string (MPI_Comm comm, octave_value &ov, int source, int mytag)
{

  // it's just a  string value directly MPI_Recv it

  std::string cpp_string;
  OCTAVE_LOCAL_BUFFER(int, tanktag, 2);
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;

  int info, nitem;
  MPI_Status stat;
  info = MPI_Recv (&nitem, 1, MPI_INT, source, tanktag[1], comm, &stat);
  OCTAVE_LOCAL_BUFFER(char, mess, nitem + 1);
  if (info == MPI_SUCCESS) 
    {
      info = MPI_Recv (mess, nitem + 1, MPI_CHAR, source, tanktag[2], comm, &stat);
      if (info == MPI_SUCCESS) 
        {
          cpp_string = mess;
          ov = cpp_string;
        }
    }  
  return info; 
}

int 
recv_matrix (bool is_complex, MPI_Datatype TRCV, const MPI_Comm comm, 
             octave_value &ov, int source, int mytag)
{

  OCTAVE_LOCAL_BUFFER(int, tanktag, 6);
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;
  tanktag[3] = mytag + 3;
  tanktag[4] = mytag + 4;
  tanktag[5] = mytag + 5;
  int info;
  int nitem, nd;

  MPI_Status stat;
  dim_vector dv;
 
  // Receive the number of elements
  info = MPI_Recv (&nitem, 1, MPI_INT, source, tanktag[1], comm, &stat);
  if (info != MPI_SUCCESS) return info;
    
  // Receive the number of dimensions
  info = MPI_Recv (&nd, 1, MPI_INT, source, tanktag[2], comm, &stat);
  if (info != MPI_SUCCESS) return info;

  // Create a buffer to store the dim_vector elements and receive data
  dv.resize (nd);
  OCTAVE_LOCAL_BUFFER(int, dimV, nd);

  info = MPI_Recv (dimV, nd, MPI_INT, source, tanktag[3], comm, &stat);
  if (info != MPI_SUCCESS) return info;

  // Now reverse the content of int vector into dim vector
  for (octave_idx_type i = 0; i < nd; i++)
    dv(i) = dimV[i] ;
  
  if (is_complex)
    {
      // FIXME : it should be easy to avoid the extra memory allocation and copying
#define __MAKE_CMPLX_TYPE_BRANCH__(TMPI, T1, A1)                        \
      if (TRCV == TMPI)                                                 \
        {                                                               \
          OCTAVE_LOCAL_BUFFER(T1, LBNDA1, nitem);                       \
          info = recv_vec (comm, LBNDA1, nitem, TRCV,                   \
                           source, tanktag[4]);                         \
          if (info == MPI_SUCCESS)                                      \
            {                                                           \
              A1 myNDA (dv);                                            \
              OCTAVE_LOCAL_BUFFER(T1, LBNDA2, nitem);                   \
              info = recv_vec (comm, LBNDA2, nitem, TRCV,               \
                               source, tanktag[5]);                     \
              if (info == MPI_SUCCESS)                                  \
                {                                                       \
                  for (octave_idx_type i = 0; i < nitem; i++)           \
                    myNDA(i) = real (LBNDA1[i]) + imag (LBNDA2[i]);     \
                                                                        \
                  ov = myNDA;                                           \
                }                                                       \
            }                                                           \
        }                                                               
      
      __MAKE_CMPLX_TYPE_BRANCH__(MPI_DOUBLE, double, ComplexNDArray)
      else 
        __MAKE_CMPLX_TYPE_BRANCH__(MPI_FLOAT, float, FloatComplexNDArray)

#undef __MAKE_CMPLX_TYPE_BRANCH__
    }

  else
    {
      // FIXME : it should be easy to avoid the extra memory allocation and copying
#define __MAKE_TYPE_BRANCH__(TMPI, T1, A1)                      \
      if (TRCV == TMPI)                                         \
        {                                                       \
          A1 myNDA (dv);                                        \
          OCTAVE_LOCAL_BUFFER(T1, LBNDA, nitem);                \
          info = recv_vec                                       \
            (comm, LBNDA, nitem, TRCV, source, tanktag[4]);     \
          if (info == MPI_SUCCESS)                              \
            {                                                   \
              T1 *myNDA_ptr = myNDA.fortran_vec ();             \
              for (octave_idx_type i = 0; i < nitem; i++)       \
                myNDA_ptr[i] = LBNDA[i];                        \
              ov = myNDA;                                       \
            }                                                   \
        }

      __MAKE_TYPE_BRANCH__(MPI_DOUBLE, double, NDArray)
      else 
        __MAKE_TYPE_BRANCH__(MPI_INT, bool, boolNDArray)
      else 
        __MAKE_TYPE_BRANCH__(MPI_FLOAT, float, FloatNDArray)
      else 
        __MAKE_TYPE_BRANCH__(MPI_BYTE, octave_int8, int8NDArray)
      else
        __MAKE_TYPE_BRANCH__(MPI_SHORT, octave_int16, int16NDArray)
      else
        __MAKE_TYPE_BRANCH__(MPI_INT, octave_int32, int32NDArray)
      else 
        __MAKE_TYPE_BRANCH__(MPI_LONG_LONG, octave_int64, int64NDArray)
      else 
        __MAKE_TYPE_BRANCH__(MPI_UNSIGNED_CHAR, octave_uint8, uint8NDArray)
      else
        __MAKE_TYPE_BRANCH__(MPI_UNSIGNED_SHORT, octave_uint16, uint16NDArray)
      else 
        __MAKE_TYPE_BRANCH__(MPI_UNSIGNED, octave_uint32, uint32NDArray)
      else 
        __MAKE_TYPE_BRANCH__(MPI_UNSIGNED_LONG_LONG, octave_uint64, 
                             uint64NDArray)

#undef __MAKE_TYPE_BRANCH__
    }
  return info; 
}


int 
recv_sp_mat (bool is_complex, MPI_Datatype TRcv, MPI_Comm comm, 
             octave_value &ov, int source, int mytag)
{
  int info;   
                
  OCTAVE_LOCAL_BUFFER(int, tanktag,6);
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;
  tanktag[3] = mytag + 3;
  tanktag[4] = mytag + 4;
  tanktag[5] = mytag + 5;

  MPI_Status stat;

  OCTAVE_LOCAL_BUFFER(int, s, 3);  
  
  // receive the shape and capacity of 
  // the matrix in an int vector named s
  info = MPI_Recv (s, 3, MPI_INT, source, tanktag[1], comm, &stat);

  if (info != MPI_SUCCESS) 
    return info;

  // Receive row and column index 
  OCTAVE_LOCAL_BUFFER(int, sridx, s[2]); 
  OCTAVE_LOCAL_BUFFER(int, scidx, s[1] + 1); 

  info = MPI_Recv (sridx, s[2], MPI_INT, source, tanktag[2], comm, &stat);
  if (info != MPI_SUCCESS) 
    return info;

  // receive the vector with column indexes
  info = MPI_Recv (scidx, s[1] + 1, MPI_INT, source, tanktag[3], comm, &stat);

  if (info != MPI_SUCCESS) 
    return info;

  // Now we have a different vector of non zero elements according to datatype
#define __MAKE_TYPE_BRANCH__(TMPI, T1, A1)                              \
  if (TRcv == TMPI)                                                     \
    {                                                                   \
      A1 m (s[0], s[1], s[2]);                                          \
      OCTAVE_LOCAL_BUFFER(T1, LBNDA, s[2]);                             \
      info = recv_vec (comm, LBNDA, s[2], TRcv, source, tanktag[4]);    \
      if (info != MPI_SUCCESS) return info;                             \
                                                                        \
      for (octave_idx_type i = 0; i < s[1]+1; i++)                      \
        m.cidx(i) = scidx[i];                                           \
      for (octave_idx_type i = 0; i < s[2]; i++)                        \
        {                                                               \
          m.ridx(i) = sridx[i];                                         \
          m.data(i) = LBNDA[i];                                         \
        }                                                               \
      ov = m;                                                           \
    }

  if (is_complex)
    {
      if (TRcv == MPI_DOUBLE)
        {  
          TRcv = MPI_DOUBLE;
          SparseComplexMatrix m (s[0], s[1], s[2]);
          OCTAVE_LOCAL_BUFFER(double, LBNDA1, s[2]);
          OCTAVE_LOCAL_BUFFER(double, LBNDA2, s[2]);

          info = recv_vec (comm, LBNDA1, s[2], TRcv, source, tanktag[4]);
          if (info != MPI_SUCCESS) return info;

          info = recv_vec (comm, LBNDA2, s[2], TRcv, source, tanktag[5]);
          if (info != MPI_SUCCESS) return info;		  
          
          for (octave_idx_type i = 0; i < s[1] + 1; i++)
            m.cidx(i) = scidx[i];
          
          for (octave_idx_type i = 0; i < s[2]; i++)
            {
              m.ridx(i) = sridx[i];
              m.data(i) = real (LBNDA1[i]) + imag (LBNDA2[i]);
            }
          
          ov = m;
        }
    }
  else
    {
      __MAKE_TYPE_BRANCH__(MPI_INT, bool, SparseBoolMatrix)
      else 
        __MAKE_TYPE_BRANCH__(MPI_DOUBLE, double, SparseMatrix)
    }

#undef __MAKE_TYPE_BRANCH__

  return info;
}

int 
recv_cell (MPI_Comm comm, octave_value &ov, int source, int mytag)
{
  
  OCTAVE_LOCAL_BUFFER(int, tanktag, 5);
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;
  tanktag[3] = mytag + 3;
  tanktag[4] = mytag + 4;
  
  int info;
  int nitem,nd;
  MPI_Status stat;
  dim_vector dv;
 
  //       nitem is the total number of elements 
  info = MPI_Recv ((&nitem), 1, MPI_INT, source, tanktag[1], comm, &stat);
  if (info != MPI_SUCCESS) return info;

  //      ndims is number of dimensions
  info = MPI_Recv (&nd, 1, MPI_INT, source, tanktag[2], comm, &stat);
  if (info != MPI_SUCCESS) return info;
  
  //  Now create contiguous datatype for dim vector
  dv.resize (nd);
  OCTAVE_LOCAL_BUFFER(int, dimV, nd);
  MPI_Datatype dimvec;
  MPI_Type_contiguous (nd, MPI_INT, &dimvec);
  MPI_Type_commit (&dimvec);
  
  info = MPI_Recv (dimV, 1, dimvec, source, tanktag[3], comm, &stat);
  if (info != MPI_SUCCESS) return info;

  // Now reverse the content of int vector into dim vector
  for (octave_idx_type i=0; i<nd; i++)
    dv(i) = dimV[i] ;

  Cell oc (dv);

  // Now focus on every single octave_value
  int newtag = tanktag[4];
  int ocap;
  for (octave_idx_type i = 0; i < nitem; i++)
    {
      octave_value celem;				
      info = MPI_Recv (&ocap, 1, MPI_INT, source, newtag, comm, &stat);
      if (info != MPI_SUCCESS) return info;

      newtag = newtag + ocap;
      info = recv_class (comm, celem, source, newtag);
      if (info != MPI_SUCCESS) return info;

      oc.Array<octave_value>::elem(i) = celem;
    }
  ov = oc;

  MPI_Type_free (&dimvec);
  return info;
}

int 
recv_struct (MPI_Comm comm, octave_value &ov, int source, int mytag)
{
  octave_scalar_map om;
  int n; // map.fields ();

  OCTAVE_LOCAL_BUFFER(int, tanktag, 2);
  tanktag[0] = mytag;      // t_id
  tanktag[1] = mytag + 1;  // n
  int tagcap = mytag + 2;
  int ntagkey = mytag + 3; // string
  int ctag = mytag + 4;    // cell
  int info;
  MPI_Status stat;
  info = MPI_Recv (&n, 1,MPI_INT, source, tanktag[1], comm, &stat);

  int scap;  
  for (int i = 0; i < n; i++)
    {	
      /* nkeys: foreach, get key */
      octave_value ov_string;
      ntagkey = ntagkey + 3;
      info = recv_class (comm, ov_string, source, ntagkey);
      std::string key = ov_string.string_value ();

      if (info != MPI_SUCCESS) return info;

      /* all elements on this fname */
      octave_value conts;
      
      // Receives capacity
      info = MPI_Recv (&scap, 1, MPI_INT, source, tagcap, comm, &stat);
      tagcap = tagcap + 1;
      ctag = ctag + scap;
      info = recv_class (comm, conts, source, ctag);

      if (! conts.is_cell ()) return MPI_ERR_UNKNOWN;
      
      om.assign (key, conts.cell_value ());
    }
  if (n != om.nfields ())
    {
      error ("MPI_Recv: inconsistent map length");
      return MPI_ERR_UNKNOWN;
    }
  
  ov=om;  
  return MPI_SUCCESS;
}


int
recv_class (MPI_Comm comm, octave_value &ov, int source, int mytag )
{
  /*------------------------------------*/    
  /* varname-strlength 1st, dims[ndim] */
  /* and then appropriate specific info */
  int t_id;
  MPI_Status status;
     
  int info = MPI_Recv (&t_id, 1, MPI_INT, source, mytag, comm, &status);
   
  static string_vector pattern = octave_value_typeinfo::installed_type_names ();
  const std::string tstring = pattern (t_id); 

  if (tstring == "cell")   
    return (recv_cell (comm,  ov, source, mytag));
  
  if (tstring == "struct") 
    return (recv_struct (comm, ov, source, mytag)); 
  
  if (tstring == "scalar")  
    {
      double d = 0; 
      MPI_Datatype TRcv = MPI_DOUBLE;
      info = recv_scalar (TRcv, comm, d, source, mytag);
      ov = d;
      return info;
    }

    if (tstring == "bool")
      {
        bool b; 
        MPI_Datatype TRcv = MPI_INT;
        info = recv_scalar (TRcv,comm, b,source,mytag);   
        ov = b;
        return info;
      }

    if (tstring == "int8 scalar")       
      {
        octave_int8 d; 
        MPI_Datatype TRcv = MPI_BYTE;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov=d ;
        return info;
      }

    if (tstring == "int16 scalar")
      {
        octave_int16 d; 
        MPI_Datatype TRcv = MPI_SHORT;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d;
        return info;
      }

    if (tstring == "int32 scalar")
      {
        octave_int32 d; 
        MPI_Datatype TRcv = MPI_INT;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d;
        return info;
      }

    if (tstring == "int64 scalar")
      {
        octave_int64 d; 
        MPI_Datatype TRcv = MPI_LONG_LONG;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d; 
        return info;
      }

    if (tstring == "uint8 scalar")
      {
        octave_uint8 d; 
        MPI_Datatype TRcv = MPI_UNSIGNED_CHAR;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d; 
        return info;
      }

    if (tstring == "uint16 scalar")
      {
        octave_uint16 d;
        MPI_Datatype TRcv = MPI_UNSIGNED_SHORT;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d; 
        return info;
      }

    if (tstring == "uint32 scalar")
      {
        octave_uint32 d;
        MPI_Datatype TRcv = MPI_UNSIGNED;
        info = recv_scalar (TRcv,comm, d,source,mytag);
        ov = d;
        return info;
      }

    if (tstring == "uint64 scalar")
      {
        octave_uint64 d;
        MPI_Datatype TRcv = MPI_UNSIGNED_LONG_LONG;
        info = recv_scalar (TRcv,comm, d,source,mytag);
        ov = d;
        return info;
      }

    if (tstring == "float scalar")
      {
        float d;
        MPI_Datatype TRcv = MPI_FLOAT;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d;
        return info;
      }

    if (tstring == "complex scalar")
      {
        std::complex<double> d;
        MPI_Datatype TRcv = MPI_DOUBLE;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d;
        return info;
      }

    if (tstring == "float complex scalar")
      { 
        std::complex<float> d; 
        MPI_Datatype TRcv = MPI_FLOAT;
        info = recv_scalar (TRcv, comm, d, source, mytag);
        ov = d;
        return info;
      }

    if (tstring == "string")
      return (recv_string (comm, ov, source, mytag));

    if (tstring == "sq_string")
      return (recv_string (comm, ov, source, mytag));

    if (tstring == "range")
      return (recv_range (comm, ov, source, mytag));
    
    if (tstring == "matrix")
      return (recv_matrix (false, MPI_DOUBLE, comm, ov, source, mytag));
    
    if (tstring == "complex matrix")
      return (recv_matrix (true, MPI_DOUBLE, comm, ov, source, mytag));

    if (tstring == "bool matrix")
      return (recv_matrix (false, MPI_INT, comm, ov, source, mytag));

    if (tstring == "int8 matrix")
      return (recv_matrix (false, MPI_BYTE, comm, ov, source, mytag));
    
    if (tstring == "int16 matrix") 		
        return (recv_matrix (false, MPI_SHORT, comm, ov, source, mytag));

    if (tstring == "int32 matrix")
      return (recv_matrix (false, MPI_INT, comm, ov, source, mytag));

    if (tstring == "int64 matrix") 	
      return (recv_matrix (false, MPI_LONG_LONG, comm, ov, source, mytag));

    if (tstring == "uint8 matrix")
      return (recv_matrix (false, MPI_UNSIGNED_CHAR, comm, ov, source, mytag));

    if (tstring == "uint16 matrix")		
      return (recv_matrix (false, MPI_UNSIGNED_SHORT, comm, ov, source, mytag));

    if (tstring == "uint32 matrix")
      return (recv_matrix (false, MPI_UNSIGNED, comm, ov, source, mytag));

    if (tstring == "uint64 matrix")
      return (recv_matrix (false, MPI_UNSIGNED_LONG_LONG, comm, ov, source, mytag));

    if (tstring == "float matrix")
      return (recv_matrix (false, MPI_FLOAT,comm, ov,source,mytag));
    
    if (tstring == "float complex matrix")
      return (recv_matrix (true, MPI_FLOAT, comm, ov, source, mytag));

    if (tstring == "sparse matrix")
      return (recv_sp_mat (false, MPI_DOUBLE, comm, ov, source, mytag));

    if (tstring == "sparse complex matrix")
      return (recv_sp_mat (true, MPI_DOUBLE, comm, ov, source, mytag));

    if (tstring == "<unknown type>")
      {
        error ("MPI_Recv: unknown class");
        return MPI_ERR_UNKNOWN;
      }
    else
      {
        error ("MPI_Recv: unsupported class %s", ov.type_name ().c_str ());
        return MPI_ERR_UNKNOWN;
      }
}


DEFUN_DLD(MPI_Recv, args, nargout,"-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} [@var{VALUE} @var{INFO}] = \
MPI_Recv(@var{SOURCE},@var{TAG},@var{COMM})\n\
Receive an MPI message containing an Octave variable and extract its value.\n\
The Octave variable being received is returned as @var{VALUE},\n\
while @var{INFO} is an integer indicating success or failure.\n\
 @example\n\
 @group\n\
@var{SOURCE} must be an integer indicating source processes \n\
@var{TAG} must be an integer to identify the message by openmpi \n\
@var{COMM} must be an octave communicator object created by \
MPI_Comm_Load function \n\
@end group\n\
@end example\n\
@seealso{MPI_Comm_Load,MPI_Init,MPI_Finalize,MPI_Send}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();
  if (nargin != 3)
    print_usage ();
  else
    {
      int source = args(0).int_value ();    
      int mytag = args(1).int_value ();
      if (! error_state)
        {
          if (! simple_type_loaded)
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

              octave_value result;
              int info = recv_class (comm, result, source, mytag);

              comm = NULL;
              retval(1) = octave_value (info);
              retval(0) = result;
            }
          else
            {
              error ("Please enter octave comunicator object!");
              retval = octave_value (-1);
            }
        }
    }
  return retval;
}
