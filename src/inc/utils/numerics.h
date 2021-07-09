// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/numerics.h
///
/// \author Walter Dehnen
/// \date   1994-2008
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1994-2008 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// CONTENTS                                                                     
//                                                                              
// sorting                                              v0.2 & v0.3             
// find position in ordered table                       v0.0                    
// polynomial interpolation on grids                    v0.0                    
// root finding                                         v0.0                    
// bracket a minimum                                    v0.0                    
// function minimization                                v0.0                    
// Burlisch-Stoer integration of 1D real integrals      v0.0                    
// Runge-Kutte 4th order integrator                     v0.0                    
// Legendre polynomials and their derivatives           v0.1                    
// cubic spline (#included from walter/spln.h)          v0.1                    
// Gauss-Legendre integration: points & weights         v0.1                    
// eigenvalues and vectors of symmetric matrices        v0.4                    
//                                                                              
////////////////////////////////////////////////////////////////////////////////
#ifndef WDutils_included_numerics_h
#define WDutils_included_numerics_h

#ifndef WDutils_included_iostream
#  include <iostream>
#  define WDutils_included_iostream
#endif
#ifndef WDutils_included_cstdlib
#  include <cstdlib>
#  define WDutils_included_cstdlib
#endif
#ifndef WDutils_included_algorithm
#  include <algorithm>
#  define WDutils_included_algorithm
#endif
#ifndef WDutils_included_inline_h
#  include <inline.h>
#endif
#ifndef WDutils_included_memory_h
#  include <memory.h>
#endif

////////////////////////////////////////////////////////////////////////////////
#ifdef __INTEL_COMPILER
#pragma warning (disable:1418) /* intel: "this warning can safely be ignored" */
#endif
////////////////////////////////////////////////////////////////////////////////
namespace WDutils {
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name finding position in ordered table                                   
  //@{                                                                          
  //----------------------------------------------------------------------------
  /// find position in ordered table (see NR)                                   
  ///                                                                           
  /// \note template parameter @c T must support < and > operators
  /// \param[in] xarr  array of T, must be ordered (ascending or descending)
  /// \param[in] n     size of array xarr
  /// \param[in] x     value to find position fo
  /// \param[in] j     initial guess for position
  /// \return position jlo such that xarr[jlo] <= x < xarr[jlo+1]
  ///
  /// The ordered table xarr is hunted for jlo such that
  ///    xarr[jlo] <= x < xarr[jlo+1].
  /// For an ascendingly ordered array, we return -1 if x < xarr[0], n-1 if
  /// x == xarr[n-1], and n if x > xarr[n--1].
  template<typename T>
  int hunt(const T*xarr, int n, T x, int j) {
    int  jlo=j,jhi,l=n-1;
    bool ascnd=(xarr[l]>xarr[0]);
    if( !ascnd && xarr[l]>=xarr[0] ) return -1;	    // x_0 = x_l
    if( (ascnd && x<xarr[0]) || (!ascnd && x>xarr[0]) ) return -1;
    if( (ascnd && x>xarr[l]) || (!ascnd && x<xarr[l]) ) return  n;

    if(jlo<0 || jlo>l) {                            // input guess not useful,
      jlo = -1;                                     //   go to bisection below
      jhi = n;
    } else {
      int inc = 1;
      if((x>=xarr[jlo]) == ascnd) {                 // hunt upward
	if(jlo == l) return (x<=xarr[l] && x>=xarr[l])? l : n;
	jhi = jlo+1;
	while((x>=xarr[jhi]) == ascnd) {            // not done hunting
	  jlo =jhi;
	  inc+=inc;                                 // so double the increment
	  jhi =jlo+inc;
	  if(jhi>l) {                               // off end of table
	    jhi=n;
	    break;
	  }
	}
      } else {                                      // hunt downward
	if(jlo == 0) return ascnd? -1 : 0;
	jhi = jlo;
	jlo-= 1;
	while((x<xarr[jlo]) == ascnd) {             // not done hunting
	  jhi = jlo;
	  inc+= inc;                                // so double the increment
	  jlo = jhi-inc;
	  if(jlo < 0) {                             // off end of table 
	    jlo = 0;
	    break;
	  }
	}
      }
    }
    while (jhi-jlo != 1) {                          // bisection phase
      int jm=(jhi+jlo) >> 1;
      if((x>=xarr[jm]) == ascnd) jlo=jm;
      else jhi=jm;
    }
    return jlo;
  }
  //----------------------------------------------------------------------------
  /// find position in ordered table (see NR)
  ///
  /// \note template param T must support < and > operators
  /// \param[in,out] k  jlo such that xarr[jlo] <= x < xarr[jlo+1]
  /// \param[in]     x  array of T, must be ordered (ascending or descending) 
  /// \param[in]     n  size of array xarr
  /// \param[in]     xi value to find position for
  ///
  /// If the original value for k already gives the position, we return.
  /// Otherwise, we guess k from linear interpolation and then invoke hunt().
  /// If x is not in range, we throw an error.
  template<typename T>
  inline void find(int&k, int n, const T*x, T xi)
  {
    if(k<0 || k>=n-1 || x[k]>xi || x[k+1]<xi) {
      k = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
      k = hunt(x,n,xi,k);
      if(k<0 || k>=n) 
	WDutils_THROW("find(): x=%f out of range [%f,%f]\n", xi,x[0],x[n-1]);
    }
  }
  //@}
  // ///////////////////////////////////////////////////////////////////////////
  //                                                                            
  /// \name sorting and related                                                 
  //{@                                                                          
  //----------------------------------------------------------------------------
  /// the numbers 0 to n-1 are ordered in ascending order of A[i]
  /// using the heap-sort algorithm; based on a routine given in NR.
  /// \note a routine HeapSort() can be found in file heap.h.
  /// \note type @c sortit must have @c sortable sortit::operator[](int)
  /// \param A (input) array or values to be sorted
  /// \param n (input) number of elements
  /// \param indx (output) index table so that A[indx[i]] is sorted
  template<typename sortable, class sortit>
  void HeapIndex(const sortit&A, int n, int*indx)
  {
    if(n<=0) return;
    if(n==1) { indx[0]=0; return; }
    for(int j=0; j!=n; ++j) indx[j] = j;
    int l = n>>1;
    int ir= n-1;
    sortable q;
    for(;;) {
      int indxt;
      if(l>0)
	q = A[indxt=indx[--l]];
      else {
	q = A[indxt=indx[ir]];
	indx[ir] = indx[0];
	if(--ir == 0) {
	  indx[0] = indxt;
	  return;
	}
      }
      int i = l;
      int j = (l<<1) + 1;
      while(j<=ir) {
	if(j < ir && A[indx[j]] < A[indx[j+1]] ) j++;
	if(q < A[indx[j]] ) {
	  indx[i] = indx[j];
	  j+= 1+(i=j);
	} else
	  j = ir+1;
      }
      indx[i] = indxt;
    }
  }
  //----------------------------------------------------------------------------
  /// the numbers 0 to n-1 are ordered in ascending order of A[i]
  /// using the heap-sort algorithm; based on a routine given in NR.
  /// \note a routine HeapSort() can be found in file heap.h.
  /// \param A (input) array or values to be sorted
  /// \param n (input) number of elements
  /// \param indx (output) index table so that A[indx[i]] is sorted
  template<typename sortable, typename integer>
  void HeapIndex(const sortable*A, size_t n, integer*indx)
    // based on a routine given in NR                                           
    // the numbers 0 to n-1 are ordered in ascending order of A[i]              
  {
#if __cplusplus >= 201103L
    WDutilsStaticAssert(std::is_integral<integer>::value);
#endif
    WDutilsAssert      (n <= size_t(std::numeric_limits<integer>::max()));
    if(n<=0) return;
    if(n==1) { indx[0]=0; return; }
    for(size_t j=0; j!=n; ++j) indx[j] = integer(j);
    size_t l = n>>1;
    size_t ir= n-1;
    for(;;) {
      const sortable*q;
      integer indxt;
      if(l>0)
	q = A+(indxt=indx[--l]);
      else {
	q = A+(indxt=indx[ir]);
	indx[ir] = indx[0];
	if(--ir == 0) {
	  indx[0] = indxt;
	  return;
	}
      }
      size_t i = l;
      size_t j = (l<<1) + 1;
      while(j<=ir) {
	if(j  < ir && A[indx[j]] < A[indx[j+1]] ) j++;
	if(*q < A[indx[j]] ) {
	  indx[i] = indx[j];
	  j+= 1+(i=j);
	} else
	  j = ir+1;
      }
      indx[i] = indxt;
    }
  }
  //----------------------------------------------------------------------------
  /// the numbers 0 to n-1 are ordered in ascending order of A[i]
  /// using the heap-sort algorithm; based on a routine given in NR
  /// \param A (input) Array or values to be sorted
  /// \param I (output) index table so that A[I[i]] is sorted
  template<typename sortable, typename integer>
  void HeapIndex(Array<sortable,1>const&A, Array<integer,1>&I)
  {
    if(A.size() != I.size()) WDutils_THROW("size mismatch in HeapIndex()\n");
    HeapIndex(A.array(),A.size(),I.array());
  }
  //----------------------------------------------------------------------------
  /// given an array of values, produce a table of their ranks
  template<typename sortable>
  void HeapRank(const sortable*A, int n, int*rank)
  {
    if(n<=0) return;
    int*indx = WDutils_NEW(int,size_t(n));
    HeapIndex(A,n,indx);
    for(int r=0; r!=n; ++r) rank[indx[r]] = r;
    WDutils_DEL_A(indx);
  }
  //----------------------------------------------------------------------------
  /// given an array of values, produce a table of their ranks
  template<typename sortable>
  void HeapRank(Array<sortable,1>const&A, Array<int,1>&I)
  {
    if(A.size() != I.size()) WDutils_THROW("size mismatch in HeapRank()\n");
    HeapRank(A.array(),A.size(),I.array());
  }
  //----------------------------------------------------------------------------
  /// \brief
  /// Find index or position of a 1D distribution given rank or percentile.
  /// \details
  /// Given a 1D arrays of positions and weights (optional), we facilitate the
  /// search of the element with given rank w.r.t. the position and return its
  /// original index or position. If weights are provided, we also facilitate
  /// the search of the element just below a given cumulative weight and
  /// return its original index or position.
  /// \note The position need not to be ordered originally (that's the point).
  /// \note Instantinations for float and double
  template<typename scalar> class FindPercentile {
    void *DATA;
    /// copy constructor disabled; use references instead of copies
    FindPercentile(const FindPercentile&);
    void setup(const scalar*, unsigned, const scalar*, unsigned);
    void setup(unsigned, scalar(*)(unsigned), unsigned);
    void setup(unsigned, void(*)(unsigned, scalar&, scalar&), unsigned);
  public:
    /// ctor: setup sort tree
    /// \param[in] X  array with positions
    /// \param[in] N  number of elements
    /// \param[in] W  (optional) array with positive weights
    /// \param[in] K  (optional) expected number of calls to functional members
    /// \note If no weights are provided, each element is assumed to have
    ///       weight 1. This implies that the cumulative weight and the rank
    ///       for each element are identical.
    FindPercentile(const scalar*X, unsigned N, const scalar*W=0, unsigned K=0)
      : DATA(0)
    {
      setup(X,N,W,K);
    }
    /// ctor: setup sort tree
    /// \param[in] X  Array with positions
    /// \param[in] K  (optional) expected number of calls to Index() or Value()
    /// \note Each element is assumed to have weight 1. This implies that the
    ///       cumulative weight and the rank for each element are identical.
    FindPercentile(Array<scalar,1>const&X, unsigned K=0)
      : DATA(0)
    {
      setup(X.array(),X.size(),0,K);
    }
    /// ctor: setup sort tree
    /// \param[in] X  Array with positions
    /// \param[in] W  Array with weights
    /// \param[in] K  (optional) expected number of calls to Index() or Value()
    FindPercentile(Array<scalar,1>const&X, Array<scalar,1>const&W, unsigned K=0)
      : DATA(0)
    {
      if(X.size() != W.size())
	WDutils_THROW("FindPercentile: positions vs weight number mismatch "
		      "(%d vs %d)\n",X.size(),W.size());
      setup(X.array(),X.size(),W.array(),K);
    }
    /// ctor: setup sort tree
    /// \param[in] N  number of elements
    /// \param[in] X  function retuning position given index i in [0,N-1]
    /// \param[in] K  (optional) expected number of calls to Index() or Value()
    /// \note Each element is assumed to have weight 1. This implies that the
    ///       cumulative weight and the rank for each element are identical.
    FindPercentile(unsigned N, scalar(*X)(unsigned), unsigned K=0)
      : DATA(0)
    {
      setup(N,X,K);
    }
    /// ctor: setup sort tree
    /// \param[in] N  number of elements
    /// \param[in] F  function setting position & weight given i in[0,N-1]
    /// \param[in] K  (optional) expected number of calls to Index() or Value()
    FindPercentile(unsigned N, void(*F)(unsigned,scalar&,scalar&), unsigned K=0)
      : DATA(0)
    {
      setup(N,F,K);
    }
    /// dtor: de-allocate sorttree
    ~FindPercentile();
    /// total weight of all points: beyond maximum for cumulative weight
    scalar   TotalWeight() const;
    /// total number of points: beyond maximum for ranks
    unsigned TotalNumber() const;
    /// handle for search result
    typedef const void* handle;
    /// find element of given rank
    /// \param[in] r  positional rank in [0,N-1]
    /// \return    handle to element of given rank
    handle FindRank(unsigned r) const;
    /// find element of given cumulative weight
    /// \param[in] w  cumulative weight in [0, TotalWeight()[
    /// \return    handle to element at which cumulative weight equals w
    /// \note If no weights were given, this is equivalent to FindRank(r=w).
    handle FindCumulativeWeight(scalar w) const;
    /// given a handle, find the next element up
    /// \param[in] h element handle
    /// \return    handle to next element in ascending position
    /// \note If h refers to the last element, NULL is returned
    /// \note Looping through all elments in ascending order using this
    ///       method, whilst in principle possible, is not efficient. If
    ///       the rank of all elements is required, use a sorting algorithm.
    handle Next(handle h) const
    {
      unsigned r = Rank(h) + 1;
      return r == TotalNumber()? 0 : FindRank(r);
    }
    /// given a handle, find the next element down
    /// \param[in] h element handle
    /// \return    handle to previous element in ascending position
    /// \note If h refers to the first element, NULL is returned
    /// \note Looping through all elments in descending order using this
    ///       method, whilst in principle possible, is not efficient. If
    ///       the rank of all elements is required, use a sorting algorithm.
    handle Previous(handle h) const
    {
      unsigned r = Rank(h);
      return r? FindRank(r-1) : 0;
    }
    /// index given handle
    /// \param[in] h handle as returned from FindRank() or FindWeight()
    /// \return index of element referred to by handle
    unsigned Index(handle h, bool=true) const;
    /// rank given handle
    /// \param[in] h element handle as returned from FindRank() or FindWeight()
    /// \return rank of element referred to by handle
    unsigned Rank(handle h, bool=true) const; 
    /// position given handle
    /// \param[in] h element handle as returned from FindRank() or FindWeight()
    /// \return position of element referred to by handle
    scalar Position(handle h, bool=true) const;
    /// cumulative weight given handle
    /// \param[in] h handle as returned from FindRank() or FindWeight()
    /// \return cumulative weight at element referred to by handle
    scalar CumulativeWeight(handle h, bool=true) const;
    /// local weight given handle
    /// \param[in] h handle as returned from FindRank() or FindWeight()
    /// \return weight of element referred to by handle
    scalar Weight(handle h, bool=true) const;
    /// index given rank
    /// \param[in] r  positional rank in [0,N-1]
    /// \return index of element with that rank
    unsigned IndexOfRank(unsigned r) const
    {
      return Index(FindRank(r),0);
    }
    /// position given rank
    /// \param[in] r  positional rank in [0,N-1]
    /// \return position of element with that rank
    scalar PositionOfRank(unsigned r) const
    {
      return Position(FindRank(r),0);
    }
    /// index given cumulative weight
    /// \param[in] w  cumulative weight in [0, TotalWeight()]
    /// \return index of element with that cumulative weight
    unsigned IndexOfCumulativeWeight(scalar w) const
    {
      return Index(FindCumulativeWeight(w),0);
    }
    /// position given cumulative weight
    /// \param[in] w  cumulative weight in [0, TotalWeight()]
    /// \return position of element with that cumulative weight
    scalar PositionOfCumulativeWeight(scalar w) const
    {
      return Position(FindCumulativeWeight(w),0);
    }
  };
  //@}
  //----------------------------------------------------------------------------
  /// class to support reporting file and line number on error
  class FileLineFind {
  protected:
    const char* file;
    const int   line;
    /// find the index of the first point to be used in interpolation
    /// \param[in] n size of table == total number of points
    /// \param[in] m number of points used for interpolation
    /// \param[in] x table of points
    /// \param[in] xi point to be interpolated at
    /// \param[in,out] j index of first point to be used in interpolation
    /// \return number of points required in interpolation (1 or min{m,n})
    template<typename X>
    static int find(int&j, int n, int m, const X*x, X xi)
    {
      int M=m<n? m:n;
      j = int( (xi-x[0]) / (x[n-1]-x[0]) * (n-1) );
      j = hunt(x,n,xi,j) - (M+1)/2 + 1;
      if(j>=0 && j<n && x[j]<=xi && x[j]>=xi) M = 1;
      else if(j<0)		              j = 0;
      else if(j>n-M)		              j = n-M;
      return M;
    }
    FileLineFind() : file(0), line(0) {}
    /// constructor: take file and line
    FileLineFind(const char*f, int l) : file(f), line(l) {}
  };
  //----------------------------------------------------------------------------
  /// \name polynomial interpolation in 1D
  //@{
  /// supporting macro Polev.
  class PolynomialEvaluation : private FileLineFind {
    /// polynomial interpolation using n values; adapted from NR
    /// \param[in] x  array of points
    /// \param[in] y  array of values
    /// \param[in] n  number of points
    /// \param[in] P  auxialiary array of size n
    /// \param[in] xi position to be interpolated at
    /// \return y(xi) as interpolated
    template<typename X, typename Y>
    Y polint(int n, const X*x, const Y*y, Y*P, X xi)
    {
      for(int i=0;i!=n;++i)
	P[i]=y[i];
      for(int m=1;m!=n;++m)
	for(int i=0;i<n-m;++i) {
	  if(x[i]<=x[i+m] && x[i]>=x[i+m]) {
	    if(file)
	      WDutils_THROWN("[%s:%d]: x's not distinct in Polev(): "
			     "x[%d]=%g=x[%d]=%g (xi=%g, x=%p)\n",
			     file,line,i,x[i],i+m,x[i+m],xi,x);
	    else
	      WDutils_THROW ("x's not distinct in polev(): "
			     "x[%d]=%g=x[%d]=%g (xi=%g, x=%p)\n",
			     i,x[i],i+m,x[i+m],xi,x);
	  }
	  P[i]= ( (xi-x[i+m])*P[i] + (x[i]-xi)*P[i+1] ) / (x[i] - x[i+m]);
	}
      return P[0];
    }
    //..........................................................................
    template<int M, typename X, typename Y> inline
    Y interpol(X xi, const X*x, const Y*y, int n)
    { 
      Y P[M];
      int j;
      return find(j,n,M,x,xi)==1? y[j] : polint(M,x+j,y+j,P,xi);
    }
  public:
    /// default constructor
    PolynomialEvaluation() : FileLineFind() {}
    /// constructor: take file and line
    PolynomialEvaluation(const char*f, int l) : FileLineFind(f,l) {}
    /// polynomial interpolation using m of n values
    /// \param[in] n  total size of arrays
    /// \param[in] xi position to find function value at
    /// \param[in] x  array of points
    /// \param[in] y  array of values
    /// \param[in] m  number of points to use in interpolation
    /// \note Together with the macro Polev this implements the function
    ///       polev()
    template<typename X, typename Y> inline
    Y operator()(X xi, const X*x, const Y*y, int n, int m=4)
    {
      switch(m) {
      case 2: return interpol<2,X,Y>(xi,x,y,n);
      case 3: return interpol<3,X,Y>(xi,x,y,n);
      case 4: return interpol<4,X,Y>(xi,x,y,n);
      case 5: return interpol<5,X,Y>(xi,x,y,n);
      case 6: return interpol<6,X,Y>(xi,x,y,n);
      default: 
	if(file) 
	  WDutils_THROW("[%s:%d]: m=%d not supported in polev\n",file,line,m);
	else
	  WDutils_THROW("m=%d not supported in polev\n",m);
      }
    }
    /// polynomial interpolation using m of n values, taking Array<T> arguments
    /// \param[in] xi position to find function value at
    /// \param[in] x  array of points
    /// \param[in] y  array of values
    /// \param[in] m  number of points to use in interpolation
    /// \note Together with the macro polev this implements the function
    ///       polev() Array<> arguments
    template<typename X, typename Y>
    Y operator()(X xi, const Array<X,1>&x, const Array<Y,1>&y, int m=4)
    {
      if(x.size() != y.size()) {
	if(file)
	  WDutils_THROWN("[%s:%d]: Array size mismatch in Polev()",file,line);
	else
	  WDutils_THROW ("Array size mismatch in polev()");
      }
      return operator()(xi,x.array(),y.array(),x.size(),m);
    }
  };
  /// macro Polev: implements functions Polev() like polev().
  /// The idea is to implement the "functions" Polev() via a macro such that
  /// on error the file and line of the call to Polev() can be reported. \n
  /// The trick is simple: the macro expands code like
  /// \code Polev(xi,x,y,n); \endcode into
  /// \code PolynomialEvaluation(__FILE__,__LINE__)(xi,x,y,n); \endcode
  /// The first argument list invokes the constructor and the second the
  /// operator() members of class PolynomialEvaluation.\n
  /// \note We provide ordinary functions polev() within the namespace WDutils
  ///       with the same semantics as Polev(). The only difference being (i)
  ///       the error reporting (with file and line number in Polev()) and (ii)
  ///       the fact that Polev as a macro is in the global namespace
#define Polev WDutils::PolynomialEvaluation(__FILE__,__LINE__)
  /// polynomial interpolation using @a m of @a n values
  /// \param[in] xi position to find function value at
  /// \param[in] x  array of points
  /// \param[in] y  array of values
  /// \param[in] n  total size of arrays
  /// \param[in] m  number of points to use in interpolation
  template<typename X, typename Y>
  inline Y polev(X xi, const X*x, const Y*y, int n, int m=4)
  { return PolynomialEvaluation(0,0)(xi,x,y,n,m); }
  /// polynomial interpolation using @a m of @a n values,
  /// taking Array<T> arguments
  /// \param[in] xi position to find function value at
  /// \param[in] x  array of points
  /// \param[in] y  array of values
  /// \param[in] m  number of points to use in interpolation
  template<typename X, typename Y>
  inline Y polev(X xi, const Array<X,1>&x, const Array<Y,1>&y, int m=4)
  { return PolynomialEvaluation()(xi,x,y,m); }
  //----------------------------------------------------------------------------
  /// like polev(), but no extrapolation; gives boundary values instead
  template<typename X, typename Y>
  inline Y ipolev(X xi, const X*x, const Y*y, int n, int m=4)
  {
    if(xi < x[0]   && x[0]   < x[n-1]) return y[0];
    if(xi > x[0]   && x[0]   > x[n-1]) return y[0];
    if(xi > x[n-1] && x[n-1] > x[0]  ) return y[n-1];
    if(xi < x[n-1] && x[n-1] < x[0]  ) return y[n-1];
    return polev(xi,x,y,n,m);
  }
  //@}
  //----------------------------------------------------------------------------
  /// \name root finding
  //@{
  /// encodes NR rtsafe
  class RootSafe: private FileLineFind {
  public:
    RootSafe() : FileLineFind() {}
    RootSafe(const char*f, int l) : FileLineFind(f,l) {}
    template <typename X>
    X operator() (void(*func)(X,X&,X&), X x1, X x2, X xacc)
    {
      const int maxit=100;
      X xh,xl,dx,dxo,f,df,fh,fl,rts,temp;
      func(x1,fl,df);
      func(x2,fh,df);
      if(fl*fh >= 0.) {
	if(file)
	  throw exception("[%s.%d]: root must be bracketed in Rtsafe()",
			  file,line);
	else
	  throw exception("root must be bracketed rtsafe()");
      }
      if(fl<0.) {
	xl  = x1;
	xh  = x2;
      } else {
	xh  = x1;
	xl  = x2;
	temp= fl;
	fl  = fh;
	fh  = temp;
      }
      rts = 0.5*(x1+x2);
      dxo = abs(x2-x1);
      dx  = dxo;
      func(rts,f,df);
      for(int j=0; j!=maxit; ++j) {
	if((((rts-xh)*df-f)*((rts-xl)*df-f)>= 0.) || (abs(2.*f)>abs(dxo*df))) {
	  dxo = dx;
	  dx  = 0.5*(xh-xl);
	  rts = xl+dx;
	  if(xl==rts) return rts;
	} else {
	  dxo = dx;
	  dx  = f/df;
	  temp=rts;
	  rts-= dx;
	  if(temp==rts) return rts;
	}
	if(abs(dx)<xacc) return rts;
	func(rts,f,df);
	if(f<0.) {
	  xl  = rts;
	  fl  = f;
	} else {
	  xh  = rts;
	  fh  = f;
	}
      }
      if(file)
	throw exception("[%s.%d]: "
			"maximum number of %d iterations exceeded in Rtsafe()",
			file,line,maxit);
      else
	throw exception("maximum number of %d iterations exceeded in rtsafe()",
			maxit);
      return rts;
    }
  };
  /// like NR rtsafe
  /// \return root: x value at which f(x)=0
  /// \param[in] func function returning f & df/dx (args 2&3) given x (1st arg)
  /// \param[in] x1   left boundary of interval
  /// \param[in] x2   right boundary of interval
  /// \param[in] xacc desired accuracy for root
  /// \note f(x1)*f(x2) \b must not be positive on input
  template <typename X>
  X rtsafe(void(*func)(X,X&,X&), X x1, X x2, X xacc)
  {
    return RootSafe()(func,x1,x2,xacc);
  }
  /// macro with same functionality as function rtsafe() above, except that on
  /// error file and line number of the call are reported.
#define Rtsafe WDutils::RootSafe(__FILE__,__LINE__)
  //@}
  //----------------------------------------------------------------------------
  // bracketing a minimum                                                       
  //----------------------------------------------------------------------------
  template<typename scalar_type>
  void minbrak(scalar_type& ax,                    // I/O ax                    
	       scalar_type& bx,                    // I/O bx                    
	       scalar_type& cx,                    // O:  cx                    
	       scalar_type(*f)(const scalar_type)) // I: f(x)                   
  {
    const scalar_type GLIMIT=100.0, GOLD=1.618034, TINY=1.e-20;
    scalar_type ulim,u,r,q,fu,fa=f(ax),fb=f(bx),fc;
    if(fb<fa) { swap(ax,bx); swap(fa,fb); }
    cx = bx+GOLD*(bx-ax);
    fc = f(cx);
    while(fb > fc) {
      r    = (bx-ax)*(fb-fc);
      q    = (bx-cx)*(fb-fa);
      u    = bx-((bx-cx)*q-(bx-ax)*r)/(2*sign(max(abs(q-r),TINY),q-r));
      ulim = bx+GLIMIT*(cx-bx);
      if((bx-u)*(u-cx) > 0.0) {                    // try parabolic u in [b,c]  
	fu = f(u);
	if(fu < fc) {
	  ax = bx;
	  bx = u;
	  return;
	} else if(fu>fb) {
	  cx = u;
	  return;
	}
	u = cx + GOLD * (cx-bx);
	fu= f(u);
      } else if((cx-u)*(u-ulim) > 0.0) {
	fu = f(u);
	if(fu < fc) {
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
	  SHFT(bx,cx,u,cx+GOLD*(cx-bx))
	  SHFT(fb,fc,fu,f(u))
	}
      } else if((u-ulim)*(ulim-cx) >= 0.0) {
	u  = ulim;
	fu = f(u);
      } else {
	u  = cx + GOLD * (cx-bx);
	fu = f(u);
      }
      SHFT(ax,bx,cx,u)
      SHFT(fa,fb,fc,fu)
#undef  SHFT
    }
  }
  //----------------------------------------------------------------------------
  // function minimization                                                      
  //----------------------------------------------------------------------------
  template<typename scalar_type> scalar_type
  brent(                                           // R: f_min = f(x_min)       
	scalar_type ax,                            // I: ax   these bracket the 
	scalar_type bx,                            // I: bx   minimum, e.g. out-
	scalar_type cx,                            // I: cx   put of minbrak()  
	scalar_type(*f)(scalar_type),              // I: function f(x)          
	scalar_type tol,                           // I: accuracy wanted        
	scalar_type& xmin)                         // O: x_min                  
  {
    const int   itmax=100;
    const scalar_type cgold=0.381966, zeps=1.e-10;
    int         iter;
    scalar_type a,b,d=0.,e=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
    a = (ax<cx) ? ax:cx;
    b = (ax>cx) ? ax:cx;
    x = w = v = bx;
    fw= fv= fx= f(x);
    for(iter=0; iter!=itmax; ++iter) {
      xm   = 0.5*(a+b);
      tol1 = tol*abs(x)+zeps;
      tol2 = 2.0* tol1;
      if(abs(x-xm) <= (tol2-0.5*(b-a))) {
	xmin = x;
	return fx;
      }
      if(abs(e) > tol1) {
	r = (x-w) * (fx-fv);
	q = (x-v) * (fx-fw);
	p = (x-v) * q - (x-w) * r;
	q = 2 * (q-r);
	if(q>0.) p =-p;
	q     = abs(q);
	etemp = e;
	e     = d;
	if(abs(p)>=abs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	  d = cgold * (e = (x>=xm) ? a-x : b-x);
	else {
	  d = p / q;
	  u = x + d;
	  if ((u-a)<tol2 || (b-u) < tol2) d=sign(tol1,xm-x);
	}
      } else
	d = cgold * (e= (x>=xm) ? a-x : b-x);
      u  = (abs(d)>=tol1) ? x+d : x+sign(tol1,d);
      fu = f(u);
      if(fu<=fx) {
	if(u>=x) a=x;
	else     b=x;
#define SHFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
	SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
#undef SHFT
      } else {
	if(u<x) a=u;
	else    b=u;
	if(fu<=fw || w==x) {
	  v = w;
	  w = u;
	  fv= fw;
	  fw= fu;
	} else if (fu<=fv || v==x || v==w) {
	  v = u;
	  fv= fu;
	}
      }
    }
    WDutils_Error("in brent(): exceeding iterations");
    xmin   =x;
    return fx;
  }
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// Burlisch-Stoer integration of 1D real integrals
  //
  /// \return approximated value for integral
  /// \param[in]  func  pointer to function to be integrated
  /// \param[in]  a     lower boundary of integration interval
  /// \param[in]  b     upper boundary of integration interval
  /// \param[in]  eps   desired relative accuracy
  /// \param[out] err   (optional) actual relative error of the return value
  /// \param[in]  abort (optional) abort if exceeding maximum of iterations?
  /// \param[in]  miter (optional) maximum number of iterations
  ///                                                                           
  /// Quadrature program using the Bulirsch sequence and rational extrapolation.
  /// The algorithm is puplished in Bulirsch & Stoer, Num. Math. 9, 271-278
  /// (1967), where a routine in ALGOL is given. This is a straightforward
  /// translation into C++.
  ///
  /// \warning Do not use this routine for integrating low order polynomials (up
  /// to fourth order) or periodic functions with period equal to the interval
  /// of integration or linear combinations of both.
  //
  // ///////////////////////////////////////////////////////////////////////////
  double qbulir(double(*func)(double),
		double  a,
		double  b,
		double  eps,
		double* err  =0,
                bool    abort=true,
		int     miter=25);
  //----------------------------------------------------------------------------
  // Runge-Kutta 4th order integrator for ODEs                                  
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename vector_type>
  inline
  vector_type rk4(vector_type const&y,
		  vector_type const&dy0,
		  scalar_type x,
		  scalar_type h,
		  vector_type(*derivs)(scalar_type, vector_type const&))
  {
    const scalar_type hh=0.5*h, h6=h/6., xh=x+hh;
    vector_type yt,dym,dyt;
    dyt = derivs(xh,y+hh*dy0);
    dym = derivs(xh,y+hh*dyt);
    yt  = y+h*dym;
    dym+= dyt;
    dyt = derivs(x+h,yt);
    return y+h6*(dy0+dyt+dym+dym);
  }
  //----------------------------------------------------------------------------
  template<typename scalar_type, typename vector_type>
  inline
  vector_type rk4(vector_type const&y,
		  scalar_type x,
		  scalar_type h,
		  vector_type(*derivs)(scalar_type, vector_type const&))
  {
    return rk4(y,derivs(x,y),x,h,derivs);
  }
  //----------------------------------------------------------------------------
  // Legendre polynomials and their derivatives                                 
  //----------------------------------------------------------------------------
  template<typename S, int N>
  void LegendrePeven(S*p, double x) 
  {
    // based on a routine from J.J. Binney                                      
    // evaluates even Legendre Polys up to l=2*(N-1) at x                       
    double x2=x*x;
    p[0] = 1.;
    p[1] = 1.5*x2-0.5;
    for(int n=2; n<N; n++) {
      int l  = 2*(n-1);
      int l2 = l+l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
  }
  //----------------------------------------------------------------------------
  template<typename S, int N>
  void dLegendrePeven(S*p, S*d, double x)
    // based on a routine from J.J. Binney                                      
    // evaluates even Legendre Polys and its derivs up to l=2*(N-1) at x        
  {
    double x2=x*x;
    p[0] = 1.;
    d[0] = 0.;
    p[1] = 1.5*x2-0.5;
    d[1] = 1.5;
    for(int n=2; n<N; ++n) {
      int l  = 2*(n-1);
      int l2 = l+l;
      p[n] = - p[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + p[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) );
      p[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
      d[n] = - d[n-2] * l*(l-1)/double((l2+1)*(l2-1))
	     + d[n-1] * (x2 - (l2*l+l2-1)/double((l2-1)*(l2+3)) )
	     + p[n-1];
      d[n]*= (l2+1)*(l2+3) / double((l+1)*(l+2));
    }
    x2 = x+x;
    for(int n=0; n<N; n++)
      d[n] *= x2;
  }
  //----------------------------------------------------------------------------
#ifdef WDutils_included_tupel_h
  template<typename S, int N> inline
  void LegendrePeven(tupel<N,S>& p, double x) 
  { return LegendrePeven<S,N>(p,x); }
  template<typename S, int N> inline
  void dLegendrePeven(tupel<N,S>& p, tupel<N,S>& d, double x)
  { return dLegendrePeven<S,N>(p,d,x); }
#endif
#ifdef WDutils_included_vector_h
  template<typename S, int N> inline
  void LegendrePeven(vector<N,S>& p, double x) 
  { return LegendrePeven<S,N>(p,x); }
  template<typename S, int N> inline
  void dLegendrePeven(vector<N,S>& p, vector<N,S>& d, double x)
  { return dLegendrePeven<S,N>(p,d,x); }
#endif
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Gauss-Legendre integration: points & weights                               
  //----------------------------------------------------------------------------
  void GaussLegendre(double*, double*, unsigned);
  // ///////////////////////////////////////////////////////////////////////////
  //
  /// \name eigensystem of symmetric matrix using Jacobi transformation
  /// \note not fully tested
  //@{
  // ---------------------------------------------------------------------------
  /// Eigen values and vectors for symmetric matrix using Jacobi transformation
  /// \note    @a N  size of matrix
  /// \note    @c X  scalar type (either float or double)
  /// \param[in]  M  symmetric matrix
  /// \param[out] V  matrix with Eigenvectors of M
  /// \param[out] D  vector with Eigenvalues of M
  /// \param[out] R  number of rotations required
  ///
  /// The eigen values and vectors of a symmetric matrix are computed using
  /// Jacobi transformations (see NR section 11.1). This is not efficient for
  /// large N, hence we have coded N as template parameter.
  /// --------------------------------------------------------------------------
  template<int N, typename X>
  void EigenSymJacobi(const X M[N][N], X V[N][N], X D[N], int&R)
  {
    const X   zero    = X(0);
    const X   half    = X(0.5);
    const X   one     = X(1);
    const int MaxIter = sizeof(X)*13;
    const X   eps     = std::numeric_limits<X>::epsilon();
    // copy M to A, set V to unity, copy diagonal of A to B and D
    X A[N][N], B[N], Z[N];
    for(int ip=0; ip!=N; ++ip) {
      for(int iq=0; iq!=N; ++iq) {
	A[ip][iq] = M[ip][iq];
	V[ip][iq] = zero;
      }
      V[ip][ip] = one;
      B[ip]     = A[ip][ip];
      D[ip]     = A[ip][ip];
      Z[ip]     = zero;
    }
    // perform iteration
    R = 0;
    for(int iter=0; iter!=MaxIter; ++iter) {
      X sm(zero);
      for(int ip=0; ip!=N-1; ++ip)
	for(int iq=ip+1; iq!=N; ++iq)
	  sm += abs(A[ip][iq]);
      if(sm <= zero) return;
      X tresh = iter<3? sm/X(5*N*N) : zero;
      for(int ip=0; ip!=N-1; ++ip) {
	for(int iq=ip+1; iq!=N; ++iq) {
	  X a = A[ip][iq];
	  X g = 100 * abs(a);
	  if(iter>3 && g <= eps * abs(D[ip]) && g <= eps * abs(D[iq]))
	    A[ip][iq] = zero;
	  else if(abs(a) > tresh) {
	    X h = D[iq]-D[ip], t;
	    if(g <= eps*abs(h))
	      t = A[ip][iq]/h;
	    else {
	      X theta = half*h/a;
	      t = one/(abs(theta)+sqrt(one+theta*theta));
	      if(theta < zero) t = -t;
	    }
	    X c   = one/sqrt(one+t*t);
	    X s   = t*c;
	    X tau = s/(one+c);
	    h     = t*a;
	    Z[ip] -= h;
	    Z[iq] += h;
	    D[ip] -= h;
	    D[iq] += h;
	    A[ip][iq] = zero;
#define Rotate(M,i,j,k,l)			\
  {						\
    g = M[i][j];				\
    h = M[k][l];				\
    M[i][j] = g-s*(h+g*tau);			\
    M[k][l] = h+s*(g-h*tau);			\
  }
	    for(int j=0;    j<ip; ++j) Rotate(A,j,ip,j,iq);
	    for(int j=ip+1; j<iq; ++j) Rotate(A,ip,j,j,iq);
	    for(int j=iq+1; j<N;  ++j) Rotate(A,ip,j,iq,j);
	    for(int j=0;    j<N;  ++j) Rotate(V,j,ip,j,iq);
#undef Rotate
	    ++R;
	  }
	}
      }
      for(int ip=0; ip!=N; ++ip) {
	B[ip] += Z[ip];
	D[ip]  = B[ip];
	Z[ip]  = zero;
      }
    }
    WDutils_THROW("EigenSymJacobi(): number iteration exceeds %d\n",MaxIter);
  }
  // ---------------------------------------------------------------------------
  /// function template sorting the eigenvalues & vectors by straight insertion
  ///
  /// \note @a       N  size of matrix
  /// \note @c       X  scalar type (float or double)
  /// \param[in,out] V  matrix with Eigenvectors
  /// \param[in,out] D  vector with Eigenvalues
  // ---------------------------------------------------------------------------
  template<int N, typename X>
  void EigenSort(X V[N][N], X D[N])
  {
    for(int k,i=0; i!=N-1; ++i) {
      X p = D[k=i];
      for(int j=i+1; j!=N; ++j)
	if(D[j] >= p) p=D[k=j];
      if(k!=i) {
	D[k] = D[i];
	D[i] = p;
	for(int j=0; j!=N; ++j) {
	  p       = V[j][i];
	  V[j][i] = V[j][k];
	  V[j][k] = p;
	}
      }
    }
  }
  // ---------------------------------------------------------------------------
  /// Sorted eigen values and vectors of symmetric matrix with Jacobi transform
  /// \note    @a N  size of matrix
  /// \note    @c X  scalar type (float or double)
  /// \param[in]  M  symmetric matrix
  /// \param[out] V  matrix with Eigenvectors, sorted
  /// \param[out] D  vector with Eigenvalues = columns, sorted
  /// \param[out] R  number of rotations required
  /// This routine simply combines EigenSymJacobi() and EigenSort().
  template<int N, typename X>
  void EigenSymJacobiSorted(const X M[N][N], X V[N][N], X D[N], int&R)
  {
    EigenSymJacobi(M,V,D,R);
    EigenSort(V,D);
  }
  //@}                                                                          
  // ---------------------------------------------------------------------------
  /// replaces matrix M with its transposed
  template<int N, typename X>
  void Transpose(X M[N][N])
  {
    for(int i=1; i!=N; ++i)
      for(int j=0; j!=i; ++j)
	std::swap(M[i][j], M[j][i]);
  }
  // ---------------------------------------------------------------------------
  /// \name eigensystem of symmetric matrix using Householder transformation
  /// \note not fully tested
  //@{
  // ---------------------------------------------------------------------------
  /// reduce real symmetric matrix to tridiagonal form
  /// \note       @c T prepare for eigenvector extraction?
  /// \note       @c X scalar type (float or double)
  /// \param[in]     N size of matrix
  /// \param[in,out] A on input: real symmetric matrix,
  ///                  on putput: input required by \a EigenSystemTridiagonal()
  /// \param[out]    D diagonal elements of tridiagonal form
  /// \param[out]    E off-diagonal elements of tridiagonal form
  ///
  /// Householder reduction of a real symmetric matrix
  /// \a A[0..\a N -1][0..\a N -1]. On output, \a A is replaced by an orthogonal
  /// matrix effecting the transformation. \a D[0..\a N -1] returns the diagonal
  /// elements of the tridiagonal matrix, and \a E[0..\a N -1] the off-diagonal
  /// elements, with \a E[0]=0.
  /// If \a T is set to false, \a A has no sensible meaning on output. See NR
  /// section 11.2 for details.
  /// \warning not tested, may be buggy
  template<bool T, typename X> void HouseholderReduction(int N, X**A, X*D, X*E);
  // ---------------------------------------------------------------------------
  /// compute eigensystem of tridiagonal symmetric matrix
  ///
  /// \note       @c X  only X=float and X=double are implemented
  /// \param[in]     N  size of matrix
  /// \param[in,out] D  on input: diagonal elements; on output: eigenvalues
  /// \param[in,out] E  on input: off-diagonal elements; on output: destroyed
  /// \param[in,out] Z  on input: see below; on output: EVs corresponding to D
  ///
  /// QL algorithm with implicit shifts to determine the eigenvalues and eigen-
  /// vectors of a real symmetric tridiagonal matrix, or of a real symmetric
  /// matrix previously reduced by \a HouseholderReduction(). In the first case,
  /// \a Z on input must be the unit matrix. In the second case, \a Z on input
  /// must be the matrix returned by \a HouseholderReduction(). For details, see
  /// NR section 11.3.
  /// \warning not tested, may be buggy
  template<typename X> void EigenSystemTridiagonal(int N, X*D, X*E, X**Z);
  // ---------------------------------------------------------------------------
  /// compute eigenvalues of tridiagonal symmetric matrix
  ///
  /// \note     @c   X  only X=float and X=double are implemented
  /// \param[in]     N  size of matrix
  /// \param[in,out] D  on input: diagonal elements; on output: eigenvalues
  /// \param[in,out] E  on input: off-diagonal elements; on output: destroyed
  ///
  /// QL algorithm with implicit shifts to determine the eigenvalues of a real
  /// symmetric tridiagonal matrix, or of a real symmetric matrix previously
  /// reduced by \a HouseholderReduction(). For details, see NR section 11.3.
  /// \warning not tested, may be buggy
  template<typename X> void EigenValuesTridiagonal(int N, X*D, X*E);
  // ---------------------------------------------------------------------------
  /// eigensystem of symmetric matrix, replaces original matrix
  ///
  /// \note       @a EIGENVECTORS want eigenvectors?
  /// \note       @c X  scalar type (either float or double)
  /// \param[in]     N  size of matrix
  /// \param[in,out] A  on input: symmetric matrix; on output: eigenvectors
  /// \param[out]    D  vector with eigenvalues 
  ///
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() or
  /// EigenValuesTridiagonal() for \a N known at run time. If \a N is known at
  /// compile time, use \a EigenSymmetricFixed() below.
  /// \warning not tested, may be buggy
  template<bool EIGENVECTORS, typename X>
  void EigenSymmetricReplace(int N, X**A, X*D)
  {
    X*E = WDutils_NEW(X,N);
    HouseholderReduction<EIGENVECTORS>(N,A,D,E);
    if(EIGENVECTORS) EigenSystemTridiagonal(N,D,E,A);
    else             EigenValuesTridiagonal(N,D,E);
    WDutils_DEL_A(E);
  }
  // ---------------------------------------------------------------------------
  /// eigensystem of symmetric matrix, keeps original matrix
  ///
  /// \note    @c X  scalar type (either float or double)
  /// \param[in]  N  size of matrix
  /// \param[in]  M  symmetric matrix
  /// \param[out] D  vector with eigenvalues
  /// \param[out] V  (optional) vector with eigenvectors
  ///
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() or
  /// EigenValuesTridiagonal() for \a N known at run time. If \a N is known at
  /// compile time, use \a EigenSymmetricFixed() below.
  /// \warning not tested, may be buggy
  template<typename X>
  void EigenSymmetricKeep(int N, const X**M, X*D, X**V=0)
  {
    X*E = WDutils_NEW(X,N);
    if(V) {
      for(int i=0; i!=N; ++i)
	for(int j=0; j!=N; ++j)
	  V[i][j] = M[i][j];
      HouseholderReduction<1>(N,V,D,E);
      EigenSystemTridiagonal(N,D,E,V);
    } else {
      V = WDutils_NEW(X*,N);
      X* m = WDutils_NEW(X,N*N);
      for(int i=0; i!=N; ++i, m+=N) {
	V[i] = m;
	for(int j=0; j!=N; ++j)
	  V[i][j] = M[i][j];
      }
      HouseholderReduction<0>(N,V,D,E);
      EigenValuesTridiagonal(N,D,E);
      WDutils_DEL_A(V[0]);
      WDutils_DEL_A(V);
    }
    WDutils_DEL_A(E);
  }
  // ---------------------------------------------------------------------------
  /// eigensystem of symmetric matrix, keeps original matrix, N template param
  ///
  /// \note    @a N  size of matrix
  /// \note    @c X  scalar type (float or double)
  /// \param[in]  M  matrix
  /// \param[out] D  vector with Eigenvalues
  /// \param[out] V  (optional) matrix with Eigenvectors
  ///
  /// This combines HouseholderReduction() and EigenSystemTridiagonal() for
  /// fixed \a N (known as template parameter at compile time).
  /// \warning not tested, may be buggy
  template<int N, typename X>
  void EigenSymmetricFixed(const X M[N][N], X D[N], X**V=0)
  {
    X E[N];
    if(V) {
      for(int i=0; i!=N; ++i)
	for(int j=0; j!=N; ++j)
	  V[i][j] = M[i][j];
      HouseholderReduction<1>(N,V,D,E);
      EigenSystemTridiagonal (N,D,E,V);

    } else {
      X Z[N][N];
      for(int i=0; i!=N; ++i)
	for(int j=0; j!=N; ++j)
	  Z[i][j] = M[i][j];
      HouseholderReduction<0>(N,Z,D,E);
      EigenValuesTridiagonal (N,D,E);
    }
  }
  //@}
  //////////////////////////////////////////////////////////////////////////////
} // namespace WDutils {
////////////////////////////////////////////////////////////////////////////////
#endif // WDutils_included_numerics_h
