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

    
// forward declarations ordered by datatype

int 
send_class (MPI_Comm comm, const octave_value &ov,
            const Array<octave_idx_type> &rankrec, int mytag);

int
send_string (int t_id, MPI_Comm comm, std::string  oi8, 
             const Array<octave_idx_type> &rankrec, int mytag);

int
send_cell (int t_id, MPI_Comm comm, Cell cell, 
           const Array<octave_idx_type> &rankrec, int mytag);

int
send_struct (int t_id, MPI_Comm comm, octave_map map, 
             const Array<octave_idx_type> &rankrec, int mytag);

template <class Any>
int
send_scalar (int t_id, MPI_Datatype TSnd, MPI_Comm comm, std::complex<Any> d, 
             const Array<octave_idx_type> &rankrec, int mytag);

template <class Any>
int
send_scalar (int t_id, MPI_Datatype TSnd, MPI_Comm comm, Any d, 
             const Array<octave_idx_type> &rankrec, int mytag);

int
send_range (int t_id, MPI_Comm comm, Range range, 
            const Array<octave_idx_type> &rankrec, int mytag);

int
send_matrix (int t_id, MPI_Datatype TSnd, MPI_Comm comm, 
             const octave_value &myOv, const Array<octave_idx_type> &rankrec, 
             int mytag);

int
send_sp_mat (int t_id, MPI_Datatype TSnd, MPI_Comm comm, 
             const octave_value &MyOv, const Array<octave_idx_type> &rankrec, 
             int mytag);

// template specialization for complex case
template <class Any>
int
send_scalar (int t_id, MPI_Datatype TSnd, MPI_Comm comm, std::complex<Any> d, 
             const Array<octave_idx_type> &rankrec, int mytag)
{
  OCTAVE_LOCAL_BUFFER(int, tanktag, 2);
  OCTAVE_LOCAL_BUFFER(std::complex<Any>, Deco, 2);

  int info;
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  Deco[0] = real (d);
  Deco[1] = imag (d);

  const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();
  for (octave_idx_type i = 0; i < rankrec.nelem (); i++)
    {
      info = MPI_Send 
        (&t_id, 1, MPI_INT, rankrec_ptr[i], tanktag[0], comm);
      if (info != MPI_SUCCESS) return info;

      info = MPI_Send
        (&Deco, 2, TSnd, rankrec_ptr[i], tanktag[1], comm);
      if (info != MPI_SUCCESS) return info;
    }

  return (info);
}


template <class Any>
int
send_scalar (int t_id, MPI_Datatype TSnd, MPI_Comm comm, Any d, 
             const Array<octave_idx_type> &rankrec, int mytag)
{

  OCTAVE_LOCAL_BUFFER(int,tanktag,2);

  int info;
  tanktag[0] = mytag;
  tanktag[1] = mytag+1;
  
  const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();
  for (octave_idx_type  i = 0; i< rankrec.nelem(); i++)
    {
      info = MPI_Send 
        (&t_id, 1, MPI_INT, rankrec_ptr[i], tanktag[0], comm);
      if (info != MPI_SUCCESS) return info;
      
      info = MPI_Send 
        (&d, 1, TSnd, rankrec_ptr[i], tanktag[1], comm);
      if (info != MPI_SUCCESS) return info;
    }
  
  return (info);
}

int 
send_range (int t_id, MPI_Comm comm, Range range, 
            const Array<octave_idx_type> &rankrec, int mytag)
{

  OCTAVE_LOCAL_BUFFER(int, tanktag, 3);
  OCTAVE_LOCAL_BUFFER(double, d, 3);

  // send as: base, limit, incr, nelem 
  // just 3 doubles + 1 int    
  
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;

  //   Range (double  b, double  l, double  i)
  d[0]= range.base ();
  d[1]= range.limit ();
  d[2]= range.inc ();

  int nele = range.nelem ();
  int info = MPI_SUCCESS;

  const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();
  for (octave_idx_type i = 0; i< rankrec.nelem (); i++)
    {
      info = MPI_Send 
        (&t_id, 1, MPI_INT, rankrec_ptr[i], tanktag[0], comm);
      if (info != MPI_SUCCESS) return info;

      info = MPI_Send 
        (d, 2, MPI_DOUBLE, rankrec_ptr[i], tanktag[1], comm);
      if (info != MPI_SUCCESS) return info;

      info = MPI_Send
        (&nele, 1, MPI_INT, rankrec_ptr[i], tanktag[2], comm);
      if (info != MPI_SUCCESS) return info;
    }
   
  return (info);
}


int
send_matrix (int t_id, MPI_Datatype TSnd, MPI_Comm comm, 
             const octave_value &myOv, const Array<octave_idx_type> &rankrec, 
             int mytag)
{

  OCTAVE_LOCAL_BUFFER(int, tanktag, 6);
  // real branch:
  // OCTAVE_LOCAL_BUFFER(int, dimV, nd);                             
  // complex branch:
  // OCTAVE_LOCAL_BUFFER(T3, LBNDA1, nitem);
  // OCTAVE_LOCAL_BUFFER(T3, LBNDA2, nitem);

  int info;
  int nitem;
  dim_vector dv;
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;
  tanktag[3] = mytag + 3;
  tanktag[4] = mytag + 4;
  tanktag[5] = mytag + 5;
  int nd;

#define __MAKE_TYPE_BRANCH__(TMPI, T1, T2, T3, A1, A2)                  \
  if (TSnd == TMPI && myOv.T1  && myOv.T2 )                             \
    {                                                                   \
      A1 myNDA = myOv.A2;                                               \
      nitem = myNDA.nelem ();                                           \
      dv = myNDA.dims ();                                               \
      nd = myNDA.ndims ();                                              \
      OCTAVE_LOCAL_BUFFER(int, dimV, nd);                               \
                                                                        \
      for (octave_idx_type i = 0; i < nd; i++)                          \
        dimV[i] = dv(i) ;                                               \
                                                                        \
      T3 * LBNDA = myNDA.fortran_vec ();                                \
                                                                        \
      const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();      \
      for (octave_idx_type  i = 0; i < rankrec.nelem (); i++)           \
        {                                                               \
          info = MPI_Send                                               \
            (&t_id, 1, MPI_INT, rankrec_ptr[i], tanktag[0], comm);      \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (&nitem, 1, MPI_INT, rankrec_ptr[i], tanktag[1], comm);     \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (&nd, 1, MPI_INT, rankrec_ptr[i], tanktag[2], comm);        \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (dimV, nd, MPI_INT, rankrec_ptr[i], tanktag[3], comm);      \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info =  MPI_Send                                              \
            (LBNDA, nitem, TSnd, rankrec_ptr[i], tanktag[4], comm);     \
          if (info != MPI_SUCCESS) return info;                         \
        }                                                               \
    }                                                                 
  
#define __MAKE_CMPLX_TYPE_BRANCH__(TMPI,T1,T2,T3,A1,A2)                 \
  if (TSnd == TMPI && myOv.T1 && myOv.T2)                               \
    {                                                                   \
      A1 myNDA = myOv.A2;                                               \
      nitem = myNDA.nelem ();                                           \
      OCTAVE_LOCAL_BUFFER(T3, LBNDA1, nitem);                           \
      OCTAVE_LOCAL_BUFFER(T3, LBNDA2, nitem);                           \
                                                                        \
      dv = myNDA.dims ();                                               \
      nd = myNDA.ndims ();                                              \
      OCTAVE_LOCAL_BUFFER(int, dimV, nd);                               \
                                                                        \
      for (octave_idx_type i = 0; i < nd; ++i)                          \
        dimV[i] = dv(i);                                                \
                                                                        \
      for (octave_idx_type i = 0; i < nitem; ++i)                       \
        {                                                               \
          LBNDA1[i] = real (myNDA(i));                                  \
          LBNDA2[i] = imag (myNDA(i));                                  \
        }                                                               \
                                                                        \
                                                                        \
      const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();      \
      for (octave_idx_type  i = 0; i< rankrec.nelem (); i++)            \
        {                                                               \
          info = MPI_Send                                               \
            (&t_id, 1, MPI_INT, rankrec_ptr[i], tanktag[0], comm);      \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (&nitem, 1, MPI_INT, rankrec_ptr[i], tanktag[1], comm);     \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (&nd, 1, MPI_INT, rankrec_ptr[i], tanktag[2], comm);        \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (dimV, nd, MPI_INT, rankrec_ptr[i], tanktag[3], comm);      \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info =  MPI_Send                                              \
            (LBNDA1, nitem, TSnd, rankrec_ptr[i], tanktag[4], comm);    \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info =  MPI_Send                                              \
            (LBNDA2, nitem, TSnd, rankrec_ptr[i], tanktag[5], comm);    \
          if (info != MPI_SUCCESS) return info;                         \
        }                                                               \
    }

  __MAKE_CMPLX_TYPE_BRANCH__(MPI_DOUBLE, is_complex_type (), is_double_type (), 
                             double,ComplexNDArray,complex_array_value())
  else 
    __MAKE_CMPLX_TYPE_BRANCH__(MPI_FLOAT, is_complex_type (), is_single_type (),
                               double, FloatComplexNDArray, 
                               float_complex_array_value ())
  else 
    __MAKE_TYPE_BRANCH__(MPI_DOUBLE, is_real_type (), is_double_type (),
                         double, NDArray, array_value ())
  else
    __MAKE_TYPE_BRANCH__(MPI_INT, is_bool_type (), is_real_type (), 
                         bool, boolNDArray, bool_array_value ())
  else 
    __MAKE_TYPE_BRANCH__(MPI_FLOAT, is_single_type (), is_real_type (),
                         float, FloatNDArray, float_array_value ())
  else 
    __MAKE_TYPE_BRANCH__(MPI_BYTE, is_int8_type (), is_real_type (), 
                         octave_int8, int8NDArray, int8_array_value ())
  else 
    __MAKE_TYPE_BRANCH__(MPI_SHORT, is_int16_type (), is_real_type (),
                         octave_int16, int16NDArray, int16_array_value ())
  else
    __MAKE_TYPE_BRANCH__(MPI_INT, is_int32_type (), is_real_type (),
                         octave_int32, int32NDArray, int32_array_value ())
  else
    __MAKE_TYPE_BRANCH__(MPI_LONG_LONG, is_int64_type (), is_real_type (),
                         octave_int64, int64NDArray, int64_array_value ())
  else 
    __MAKE_TYPE_BRANCH__(MPI_UNSIGNED_CHAR, is_uint8_type (), is_real_type (),
                         octave_uint8, uint8NDArray, uint8_array_value ())
  else 
    __MAKE_TYPE_BRANCH__(MPI_UNSIGNED_SHORT, is_uint16_type (), is_real_type (), 
                         octave_uint16, uint16NDArray, uint16_array_value ())
  else
    __MAKE_TYPE_BRANCH__(MPI_UNSIGNED, is_uint32_type (), is_real_type (),
                         octave_uint32, uint32NDArray, uint32_array_value ())
  else
    __MAKE_TYPE_BRANCH__(MPI_UNSIGNED_LONG_LONG, is_uint64_type (), 
                         is_real_type (), octave_uint64, uint64NDArray, 
                         uint64_array_value ())
  else 
    error ("MPI_Send: matrix with unexpected data type");

  return (info);

#undef __MAKE_TYPE_BRANCH__
#undef __MAKE_CMPLX_TYPE_BRANCH__
}


int
send_sp_mat (int t_id, MPI_Datatype TSnd, MPI_Comm comm, 
             const octave_value &MyOv, const Array<octave_idx_type> &rankrec, 
             int mytag)
{

  OCTAVE_LOCAL_BUFFER(int,tanktag,6);

  int info;
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;
  tanktag[3] = mytag + 3;
  tanktag[4] = mytag + 4;
  tanktag[5] = mytag + 5;
  
#define __MAKE_TYPE_BRANCH__(TMPI,T0,T1,T2,A1)                          \
  if (TSnd == TMPI and MyOv.T1)                                         \
    {                                                                   \
      OCTAVE_LOCAL_BUFFER(int, s, 3);                                   \
      T2 m = MyOv.A1;                                                   \
      OCTAVE_LOCAL_BUFFER(int, sridx, m.capacity ());                   \
      OCTAVE_LOCAL_BUFFER(int, scidx, m.cols () + 1);                   \
      s[0]= m.rows ();                                                  \
      s[1]= m.cols ();                                                  \
      s[2]= m.capacity ();                                              \
                                                                        \
      for (octave_idx_type ix = 0; ix < m.cols () + 1; ix++)            \
        scidx[ix]= m.cidx (ix);                                         \
                                                                        \
      OCTAVE_LOCAL_BUFFER(T0, sdata, m.capacity ());                    \
                                                                        \
      for (octave_idx_type ix = 0; ix < m.capacity (); ix++)            \
        {                                                               \
          sdata[ix]= m.data (ix);                                       \
          sridx[ix]= m.ridx (ix);                                       \
        }                                                               \
                                                                        \
      const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();      \
      for (octave_idx_type i = 0; i < rankrec.nelem (); i++)            \
        {                                                               \
          info = MPI_Send                                               \
            (&t_id, 1, MPI_INT, rankrec_ptr[i], tanktag[0], comm);      \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (s, 3, MPI_INT, rankrec_ptr[i], tanktag[1], comm);          \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (sridx, m.capacity (), MPI_INT, rankrec_ptr[i],             \
             tanktag[2], comm);                                         \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info = MPI_Send                                               \
            (scidx, m.cols () + 1, MPI_INT, rankrec_ptr[i],             \
             tanktag[3], comm);                                         \
          if (info != MPI_SUCCESS) return info;                         \
                                                                        \
          info =  MPI_Send                                              \
            (sdata, m.capacity (), TSnd, rankrec_ptr[i],                \
             tanktag[4], comm);                                         \
          if (info != MPI_SUCCESS) return info;                         \
        }                                                               \
    }                                                                   

  __MAKE_TYPE_BRANCH__(MPI_INT, bool, is_bool_type (), 
                       SparseBoolMatrix, sparse_bool_matrix_value ())
  else
    __MAKE_TYPE_BRANCH__(MPI_DOUBLE, double, is_real_type (), 
                         SparseMatrix, sparse_matrix_value ())
  else if (TSnd == MPI_DOUBLE and MyOv.is_complex_type ())
    { 
      SparseComplexMatrix m = MyOv.sparse_complex_matrix_value ();
      OCTAVE_LOCAL_BUFFER(int,s,3);  
      s[0]= m.rows ();
      s[1]= m.cols ();
      s[2]= m.capacity ();

      OCTAVE_LOCAL_BUFFER(int,sridx,m.capacity());
      OCTAVE_LOCAL_BUFFER(int,scidx,m.cols()+1);
                  
      for (octave_idx_type ix = 0; ix < m.cols () + 1; ix++)
          scidx[ix]= m.cidx(ix);   

      OCTAVE_LOCAL_BUFFER(double, sdata1, m.capacity());
      OCTAVE_LOCAL_BUFFER(double, sdata2, m.capacity());

      // Fill them with their respective value
      for (octave_idx_type ix = 0; ix < m.capacity (); ix++)
        {
          sdata1[ix] = real (m.data(ix));
          sdata2[ix] = imag (m.data(ix));
          sridx[ix]  = m.ridx(ix);
        }

      const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();
      for (octave_idx_type i = 0; i < rankrec.nelem (); i++)
        {
          info = MPI_Send 
            (&t_id, 1, MPI_INT, rankrec_ptr[i], tanktag[0], comm);
          if (info != MPI_SUCCESS) return info;

          info = MPI_Send
            (s, 3, MPI_INT, rankrec_ptr[i], tanktag[1], comm);
          if (info != MPI_SUCCESS) return info;

          info =  MPI_Send
            (sridx, m.capacity (), MPI_INT, rankrec_ptr[i], tanktag[2], comm);
          if (info != MPI_SUCCESS) return info;

          info =  MPI_Send 
            (scidx, m.cols () + 1, MPI_INT, rankrec_ptr[i], tanktag[3], comm);
          if (info != MPI_SUCCESS) return info;

          info =  MPI_Send 
            (sdata1, m.capacity (), TSnd, rankrec_ptr[i], tanktag[4], comm);
          if (info != MPI_SUCCESS) return info;

          info =  MPI_Send
            (sdata2, m.capacity (), TSnd, rankrec_ptr[i], tanktag[5], comm);
          if (info != MPI_SUCCESS) return info;
        }
    }
  return(info);
}

int
send_string (int t_id, MPI_Comm comm, std::string oi8, 
             const Array<octave_idx_type> &rankrec, int mytag)
{
  int info;
  int nitem = oi8.length ();
  int tanktag[3];
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;

  char i8[nitem+1];
  strcpy (i8, oi8.c_str ());

  const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();
  for (octave_idx_type i = 0; i < rankrec.nelem (); i++)
    {
      info = MPI_Send 
        (&t_id, 1, MPI_INT, rankrec_ptr[i], mytag, comm);
      if (info != MPI_SUCCESS) return info;

      info = MPI_Send
        (&nitem, 1, MPI_INT, rankrec_ptr[i], tanktag[1], comm);
      if (info != MPI_SUCCESS) return info;

      info =  MPI_Send 
        (&i8, nitem + 1, MPI_CHAR, rankrec_ptr[i], tanktag[2], comm);
      if (info != MPI_SUCCESS) return info;
    }

  return(info);
}

int
send_cell (int t_id, MPI_Comm comm, Cell cell, 
           const Array<octave_idx_type> &rankrec, int mytag)
{    
  /*-------------------------------------*/ 
  /* we first store nelems and then      */
  /* recursively the elements themselves */

  // Lists of items to send
  // type_id to identify octave_value
  // n for the cell capacity
  // nd for number of dimensions
  // dimvec derived datatype
  // item of cell
  int n = cell.capacity ();
  int info;
  int tanktag[5];
  tanktag[0] = mytag;
  tanktag[1] = mytag + 1;
  tanktag[2] = mytag + 2;
  tanktag[3] = mytag + 3;
  tanktag[4] = mytag + 4;
  int newtag = tanktag[4];
  dim_vector vdim = cell.dims ();
  int nd = cell.ndims ();

  // Declare here the octave_local_buffers
  OCTAVE_LOCAL_BUFFER(int,dimV,nd);
  for (octave_idx_type i = 0; i < nd; i++)
    dimV[i] = vdim(i) ;

  // Now start the big loop
  const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();
  for (octave_idx_type i = 0; i < rankrec.nelem (); i++)
    {
      info = MPI_Send (&t_id, 1, MPI_INT, rankrec_ptr[i],
                       tanktag[0], comm);
      if (info != MPI_SUCCESS) 
        return info;
      info = MPI_Send (&n, 1, MPI_INT, rankrec_ptr[i], 
                       tanktag[1], comm);
      if (info != MPI_SUCCESS) 
        return info;
      info = MPI_Send (&nd, 1, MPI_INT, rankrec_ptr[i], 
                       tanktag[2], comm);
      if (info != MPI_SUCCESS) 
        return info;
      // send the dim vector
      info = MPI_Send (dimV, nd, MPI_INT, rankrec_ptr[i], 
                       tanktag[3], comm);
      if (info != MPI_SUCCESS) 
        return info;
      
      int cap;
      // Now focus on every single octave_value
      for (octave_idx_type j = 0; j < n; j++)
        {
          octave_value ov = cell.data ()[j];
          cap = ov.capacity ();
          info = MPI_Send (&cap, 1, MPI_INT, rankrec_ptr[i], 
                           newtag, comm);
          if (info != MPI_SUCCESS) 
            return info;
          newtag = newtag + ov.capacity ();
          info = send_class (comm, ov, rankrec, newtag);
          if (info != MPI_SUCCESS) 
            return info;
        }                                       
    }

  return (info); 
}

int
send_struct (int t_id, MPI_Comm comm, octave_map map, 
             const Array<octave_idx_type> &rankrec, int mytag)
{        /* we store nkeys, */
  int n = map.nfields (); 
  int info;
  OCTAVE_LOCAL_BUFFER(int,tanktag,2);
  tanktag[0] = mytag; //t-id
  tanktag[1] = mytag + 1; // n
  int tagcap = mytag + 2;
  int ntagkey = mytag + 3; // string

  // struct array dimensions (ND)
  // each key stores ND field-values
  dim_vector struc_dims = map.dims ();    
  dim_vector conts_dims;         

  // Now we start the big loop
  const octave_idx_type *rankrec_ptr = rankrec.fortran_vec ();
  for (octave_idx_type i = 0; i < rankrec.nelem (); i++)
    {
      info = MPI_Send (&t_id, 1, MPI_INT, rankrec_ptr[i], 
                       tanktag[0], comm);
      if (info != MPI_SUCCESS) 
        return info;
      info = MPI_Send (&n, 1, MPI_INT, rankrec_ptr[i], 
                       tanktag[1], comm);
      if (info != MPI_SUCCESS) 
        return info;

      // This is to avoid confusion between tags of strings and tags of Cells
      int   ntagCell = ntagkey + 1;

      // iterate through keys(fnames)
      int scap;
      for (octave_map::const_iterator p = map.begin (); p != map.end (); p++)
        {
          // field name
          std::string key = map.key (p);

          // Cell w/ND contents
          Cell        conts = map.contents (p);    

          // each elemt should have same ND
          conts_dims = conts.dims ();              

          if (struc_dims != conts_dims)
            {
              error ("MPI_Send: inconsistent map dims\n"); 
              return (MPI_ERR_UNKNOWN);
            }

          // Sending capacity of octave_cell
          scap = conts.capacity (); 
          info = MPI_Send 
            (&scap, 1, MPI_INT, rankrec_ptr[i], tagcap, comm);
          if (info != MPI_SUCCESS) return info;
          
          tagcap = tagcap + 1;
          ntagkey = ntagkey + 3;
          info = send_class (comm, key, rankrec, ntagkey);
          if (info != MPI_SUCCESS) return info;
          
          // Sending Cell
          ntagCell = ntagCell + conts.capacity();
          info = send_class (comm, conts, rankrec, ntagCell);
          if (info != MPI_SUCCESS) return info;
        }

      if (n != map.nfields ()) 
        {
          error ("MPI_Send: inconsistent map length\n"); 
          return (MPI_ERR_UNKNOWN);
        }
    }
  return (info);
}

int
send_class (MPI_Comm comm, const octave_value &ov, 
            const Array<octave_idx_type> &rankrec, 
            int mytag)
{
  /*----------------------------------  */
  /* varname-strlength 1st, dims[ndim]  */
  /* and then appropriate specific info */

  int t_id = ov.type_id ();
  const std::string mystring = ov.type_name ();
  MPI_Datatype TSnd;

  // range
  if (mystring == "range")
    return (send_range (t_id, comm, ov.range_value (), rankrec, mytag));

  // scalar
  else if (mystring == "scalar")  
    return (send_scalar (t_id, MPI_DOUBLE, comm, ov.scalar_value (), 
                         rankrec, mytag)); 

  else if (mystring == "int8 scalar") 
      return (send_scalar (t_id, MPI_BYTE, comm, ov.int8_scalar_value (), 
                           rankrec, mytag));

  else if (mystring == "int16 scalar")  
    return (send_scalar (t_id, MPI_SHORT, comm, ov.int16_scalar_value (), 
                         rankrec, mytag));

  else if (mystring == "int32 scalar")  
      return (send_scalar (t_id, MPI_INT, comm, ov.int32_scalar_value (), 
                           rankrec, mytag));

  else if (mystring == "int64 scalar")  
    return (send_scalar (t_id, MPI_LONG_LONG, comm, ov.int64_scalar_value (), 
                         rankrec, mytag));

  else if (mystring == "uint8 scalar")
    return (send_scalar (t_id, MPI_UNSIGNED_CHAR, comm, ov.uint8_scalar_value (), 
                         rankrec, mytag));

  else if (mystring == "uint16 scalar")
    return (send_scalar (t_id, MPI_UNSIGNED_SHORT, comm, 
                         ov.uint16_scalar_value (), rankrec, mytag));
  
  else if (mystring == "uint32 scalar")
    return (send_scalar (t_id, MPI_UNSIGNED, comm, ov.uint32_scalar_value (), 
                         rankrec, mytag));

  else if (mystring == "uint64 scalar")
    return (send_scalar (t_id, MPI_UNSIGNED_LONG_LONG, comm, 
                         ov.uint64_scalar_value (), rankrec, mytag));

  else if (mystring == "bool")
    return (send_scalar (t_id, MPI_INT, comm, ov.int_value (), rankrec, mytag));

  else if (mystring == "float scalar")
    return (send_scalar (t_id, MPI_FLOAT, comm, ov.float_value (), 
                         rankrec, mytag));

  else if (mystring == "complex scalar")
    return (send_scalar (t_id, MPI_DOUBLE, comm, ov.complex_value (), rankrec, 
                         mytag));

  else if (mystring == "float complex scalar")
    return (send_scalar (t_id, MPI_FLOAT, comm, ov.float_complex_value (), 
                         rankrec, mytag));

  // matrix
  else if (mystring ==  "matrix") 
    return (send_matrix (t_id, MPI_DOUBLE, comm, ov, rankrec, mytag));

  else if (mystring ==  "bool matrix") 
    return (send_matrix (t_id, MPI_INT, comm, ov, rankrec, mytag));

  else if (mystring ==  "int8 matrix")
    return (send_matrix (t_id, MPI_BYTE, comm, ov, rankrec, mytag));

  else if (mystring ==  "int16 matrix")
    return (send_matrix (t_id, MPI_SHORT, comm, ov, rankrec, mytag));

  else if (mystring == "int32 matrix")
    return (send_matrix (t_id, MPI_INT, comm, ov, rankrec, mytag));

  else if (mystring == "int64 matrix")
    return (send_matrix (t_id, MPI_LONG_LONG, comm, ov, rankrec, mytag));

  else if (mystring == "uint8 matrix")
    return (send_matrix (t_id, MPI_UNSIGNED_CHAR, comm, ov, rankrec, mytag));

  else if (mystring == "uint16 matrix")
    return (send_matrix (t_id, MPI_UNSIGNED_SHORT, comm, ov, rankrec, mytag));

  else if (mystring == "uint32 matrix")
    return (send_matrix (t_id, MPI_UNSIGNED, comm, ov, rankrec, mytag));

  else if (mystring == "uint64 matrix")
    return (send_matrix (t_id, MPI_UNSIGNED_LONG_LONG, comm, ov, rankrec, 
                         mytag));

  //       complex matrix
  else if (mystring ==  "complex matrix")
    return (send_matrix (t_id, MPI_DOUBLE, comm, ov, rankrec, mytag));

  else if (mystring ==  "float complex matrix")
    return (send_matrix (t_id, MPI_FLOAT, comm, ov, rankrec, mytag));

  //       sparse matrix
  else if (mystring ==  "sparse bool matrix")
    return (send_sp_mat (t_id, MPI_INT, comm, ov, rankrec, mytag));

  else if (mystring ==  "sparse matrix")
    return (send_sp_mat (t_id, MPI_DOUBLE, comm, ov, rankrec, mytag));

  else if (mystring ==  "sparse complex matrix")
    return (send_sp_mat (t_id, MPI_DOUBLE, comm, ov, rankrec, mytag));
      
  else if (mystring == "string")
    return (send_string (t_id, comm, ov.string_value (), rankrec, mytag));

  else if (mystring == "sq_string")
    return (send_string (t_id, comm, ov.string_value (), rankrec, mytag));
      
  else if (mystring == "struct")
    return (send_struct (t_id, comm, ov.map_value (), rankrec, mytag));

  else if (mystring == "cell")
    return (send_cell (t_id, comm, ov.cell_value (), rankrec, mytag));
      
  else if (mystring ==  "<unknown type>")  
    {
      error ("MPI_Send: unknown class\n");
      return (MPI_ERR_UNKNOWN);
    } 
  else
    {
      error ("MPI_Send: unsupported class %s\n", ov.type_name ().c_str ());
      return (MPI_ERR_UNKNOWN);
    }
}


DEFUN_DLD(MPI_Send, args, nargout, 
          "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{INFO} =} MPI_Send(@var{VALUE},@var{RANKS},@var{TAG},@var{COMM})\n\
Transmit an Octave variable as a set of MPI message.\n\
Return an integer @var{INFO} to indicate success or failure.\n\
@example\n\
@group\n\
@var{VALUE} must be an octave variable \n\
@var{RANKS} must be a vector containing the list of rank destination processes \n\
@var{TAG} must be an integer to identify the message by openmpi \n\
@var{COMM} must be an octave communicator object created by MPI_Comm_Load function \n\
@end group\n\
@end example\n\
\n\
@seealso{MPI_Comm_Load,MPI_Init,MPI_Finalize,MPI_Recv}\n\
@end deftypefn")
{
  octave_value retval;

  int nargin = args.length ();
  if (nargin != 4 )
    print_usage ();
  else
    {
      Array<octave_idx_type> tankrank = 
        args(1).octave_idx_type_vector_value ();
      int mytag = args(2).int_value ();
      if (! error_state)
        {
          if (! simple_type_loaded)
            {
              simple::register_type ();
              simple_type_loaded = true;
              mlock ();
            }

          if (args(3).type_id () == simple::static_type_id ())
            {
              const octave_base_value& rep = args(3).get_rep ();
              const simple& B = ((const simple &)rep);
              MPI_Comm comm = ((const simple&) B).comunicator_value ();
              int info = send_class (comm, args(0), tankrank, mytag);
              comm = NULL;
              retval = info;
            }
          else 
            error ("MPI_Send: Please enter octave comunicator object!");
        }
    }
  return retval;
}

