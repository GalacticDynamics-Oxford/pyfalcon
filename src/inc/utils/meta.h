// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   utils/inc/meta.h
///
/// \brief  support for (mostly numerical) metaprogramming
///
/// \author Walter Dehnen
///                                                                             
/// \date   2008-2013
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008-2013 Walter Dehnen
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
#ifndef WDutils_included_meta_h
#define WDutils_included_meta_h

#ifndef WDutils_included_exception_h
# include <exception.h>
#endif

#ifndef WD_HOT
#  if defined(__GNUC__) && (__GNUC__ > 4 && __GNUC_MINOR__ > 2)
#    define WD_HOT __attribute__((hot))
#  else
#    define WD_HOT
#  endif
#endif

#if __cplusplus >= 201103L
# ifndef WDutils_included_type_traits
#  include <type_traits>
#  define WDutils_included_type_traits
# endif
#endif

#if __cplusplus < 201103L
#  define noexcept
#endif

namespace WDutils {

#if __cplusplus >= 201103L
  using std::enable_if;
  using std::is_same;
  using std::is_floating_point;
  using std::remove_volatile;
  using std::remove_const;
  using std::remove_cv;
  ///
  /// is_callable<Func,Args>::value is true if F::operator(Args) exists
  ///
  template<typename Func, typename... Args>
  class is_callable
  {
    template <typename U>
    static decltype(std::declval<U>()(std::declval<Args>()...), void(), 0)
      test(int);
    template <typename>
    static void test(...);
  public:
    static const bool value = !std::is_void<decltype(test<Func>(0))>::value;
  };
#else
  template<bool B, class T = void>
  struct enable_if {};
  template<class T>
  struct enable_if<true, T> { typedef T type; };

  template<class A, class B = void>
  struct is_same { static const bool value = false; };
  template<class A>
  struct is_same<A,A> { static const bool value = true; };

  template<class T> struct remove_const          { typedef T type; };
  template<class T> struct remove_const<const T> { typedef T type; };
 
  template<class T> struct remove_volatile             { typedef T type; };
  template<class T> struct remove_volatile<volatile T> { typedef T type; };

  template<class T>
  struct remove_cv {
    typedef typename 
    remove_volatile<typename remove_const<T>::type>::type type;
  };
  template<class T>
  struct is_floating_point {
    static const bool value = 
      is_same<float, typename remove_cv<T>::type>::value  ||
      is_same<double, typename remove_cv<T>::type>::value  ||
      is_same<long double, typename remove_cv<T>::type>::value;
  };
#endif

  ///
  /// support for (mostly numerical) metaprogramming
  ///
  namespace meta {

    /// to be used instead of bool variables true
    struct True {};
    /// to be used instead of bool variables false
    struct False {};
    ////////////////////////////////////////////////////////////////////////////
    // operational dependence on even or odd integer
    template<int> struct _m_EvenOdd;
    template<> struct _m_EvenOdd<0> {
      template<typename R> static void neg (R& ) noexcept {}
      template<typename R> static R    equ (R x) noexcept { return x; }
      template<typename R> static R    pow (R  ) noexcept { return R(1); }
      template<typename R> static R    powp(R ,  R y) noexcept { return y; }
      template<typename R> static void ass (R&x, R y) noexcept { x =y; }
      template<typename R> static void add (R&x, R y) noexcept { x+=y; }
      template<typename R> static void sub (R&x, R y) noexcept { x-=y; }
      template<typename R> static void add2(R&x, R y) noexcept { x+=y; x+=y; }
      template<typename R> static void sub2(R&x, R y) noexcept { x-=y; x-=y; }
      template<typename R> static R    sum (R x, R y) noexcept { return x+y; }
      template<typename R> static R    diff(R x, R y) noexcept { return x-y; }
    };
    template<> struct _m_EvenOdd<1> {
      template<typename R> static void neg (R&x) noexcept { x=-x; }
      template<typename R> static R    equ (R x) noexcept { return -x; }
      template<typename R> static R    pow (R x) noexcept { return x; }
      template<typename R> static R    powp(R x, R y) noexcept { return x*y; }
      template<typename R> static void ass (R&x, R y) noexcept { x =-y; }
      template<typename R> static void add (R&x, R y) noexcept { x-=y; }
      template<typename R> static void sub (R&x, R y) noexcept { x+=y; }
      template<typename R> static void add2(R&x, R y) noexcept { x-=y; x-=y; }
      template<typename R> static void sub2(R&x, R y) noexcept { x+=y; x+=y; }
      template<typename R> static R    sum (R x, R y) noexcept { return x-y; }
      template<typename R> static R    diff(R x, R y) noexcept { return x+y; }
    };
    /// returns +/- x, depending on whether N is even or odd
    /// \param[in] x  rvalue
    /// \return    x  if N is even, -x if n is odd
    template<int N, typename R> inline
    R Equate(R x) noexcept { return _m_EvenOdd<N&1>::equ(x); }
    /// assigns positive or negative, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by y if N is even, -y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename R> inline
    void Assign(R&x, R y) noexcept { _m_EvenOdd<N&1>::ass(x,y); }
    /// switches sign of x if N is odd
    /// \param[in,out] x lvalue, replaced by -x if N is odd
    template<int N, typename R> inline
    void Negate(R&x) noexcept { _m_EvenOdd<N&1>::neg(x); }
    /// switches sign of x if N is odd
    /// \param[in,out] x lvalue, replaced by -x if N is odd
    template<int N, typename R> inline
    void negate_if_odd(R&x) noexcept { _m_EvenOdd<N&1>::neg(x); }
    /// adds or subtracts, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x+y if N is even, x-y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename R> inline
    void Add(R&x, R y) noexcept { _m_EvenOdd<N&1>::add(x,y); }
    /// adds or subtracts twice, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x+2y if N is even, x-2y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename R> inline
    void AddTwice(R&x, R y) noexcept { _m_EvenOdd<N&1>::add2(x,y); }
    /// subtracts or adds, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x-y if N is even, x+y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename R> inline
    void Sub(R&x, R y) noexcept { _m_EvenOdd<N&1>::sub(x,y); }
    /// subtracts or adds twice, depending on whether N is even or odd
    /// \param[in,out] x lvalue, replaced by x-2y if N is even, x+2y if N is odd
    /// \param[in]     y rvalue
    template<int N, typename R> inline
    void SubTwice(R&x, R y) noexcept { _m_EvenOdd<N&1>::sub2(x,y); }
    /// sum or difference, depending on whether N is even or odd
    /// \param[in] x rvalue
    /// \param[in] y rvalue
    /// \return x+y if N is even, x-y if N is odd
    template<int N, typename R> inline
    R Sum(R x, R y) noexcept { return _m_EvenOdd<N&1>::sum(x,y); }
    /// difference or sum, depending on whether N is even or odd
    /// \param[in] x rvalue
    /// \param[in] y rvalue
    /// \return x-y if N is even, x+y if N is odd
    template<int N, typename R> inline
    R Diff(R x, R y) noexcept { return _m_EvenOdd<N&1>::diff(x,y); }
    ////////////////////////////////////////////////////////////////////////////
    /// sign of an integer
    template<int N> struct Sign
    { const static int S = N<0? -1 : N>0? 1 : 0; };
    ///  x=-x
    inline float &negate(float &x) noexcept { return x=-x; }
    inline double&negate(double&x) noexcept { return x=-x; }
    /// Inverse of an integer: replace integer divisions of real numbers by
    /// multiplications to generate faster code.
    /// \note instantinated for all integers up to 100
    template<int N> struct IntegerInverse
#ifndef WDutilsMetaNoDefaultIntegerInverse
    {
      WDutilsStaticAssert((N>1));
      /// inverse of integer
      template<typename R> static R Inverse() noexcept
      { return R(1)/R(N); }
      /// divide scalar
      /// \param[in,out] X scalar, replaced by X/N on output
      template<typename R>
      static void Divide(R&X) noexcept
      { X *= Inverse<R>(); }
      /// ratio of scalar and N
      /// \param[in] X scalar
      /// \return    X/N
      template<typename R>
      static R Ratio(R X) noexcept
      { return X*Inverse<R>(); }
    }
#endif
    ;
    ///
    template<> struct IntegerInverse<-1> {
      template<typename R> static R Inverse() noexcept { return R(-1); }
      template<typename R> static void Divide(R&X) noexcept { negate(X); }
      template<typename R> static R Ratio(R X) noexcept { return -X; }
    };
    template<> struct IntegerInverse<0> {};
    template<> struct IntegerInverse<1> {
      template<typename R> static R Inverse() noexcept { return R(1); }
      template<typename R> static void Divide(R&) noexcept {}
      template<typename R> static R Ratio(R X) noexcept { return X; }
    };
#define DEFINVERSE(NUM,INVERSE)					\
    template<> struct IntegerInverse<NUM> {			\
      template<typename R> static R Inverse() noexcept		\
      { return R(INVERSE); }					\
      template<typename R> static void Divide(R&X) noexcept	\
      { X *= R(INVERSE); }					\
      template<typename R> static R Ratio(R X) noexcept		\
      { return X * R(INVERSE); }				\
    };								\
    template<> struct IntegerInverse<- NUM> {			\
      template<typename R> static R Inverse() noexcept		\
      { return R(- INVERSE); }					\
      template<typename R> static void Divide(R&X) noexcept	\
      { X *= R(- INVERSE); }					\
      template<typename R> static R Ratio(R X) noexcept		\
      { return X * R(- INVERSE); }				\
    }
    // all numbers up to 100
    DEFINVERSE(   2, 0.5);
    DEFINVERSE(   3, 0.333333333333333333333);
    DEFINVERSE(   4, 0.25);
    DEFINVERSE(   5, 0.2);
    DEFINVERSE(   6, 0.166666666666666666667);
    DEFINVERSE(   7, 0.142857142857142857143);
    DEFINVERSE(   8, 0.125);
    DEFINVERSE(   9, 0.111111111111111111111);
    DEFINVERSE(  10, 0.1);
    DEFINVERSE(  11, 0.090909090909090909091);
    DEFINVERSE(  12, 0.083333333333333333333);
    DEFINVERSE(  13, 0.076923076923076923077);
    DEFINVERSE(  14, 0.071428571428571428571);
    DEFINVERSE(  15, 0.066666666666666666667);
    DEFINVERSE(  16, 0.0625);
    DEFINVERSE(  17, 0.058823529411764705882);
    DEFINVERSE(  18, 0.055555555555555555556);
    DEFINVERSE(  19, 0.052631578947368421053);
    DEFINVERSE(  20, 0.05);
    DEFINVERSE(  21, 0.047619047619047619048);
    DEFINVERSE(  22, 0.045454545454545454545);
    DEFINVERSE(  23, 0.043478260869565217391);
    DEFINVERSE(  24, 0.041666666666666666667);
    DEFINVERSE(  25, 0.04);
    DEFINVERSE(  26, 0.038461538461538461538);
    DEFINVERSE(  27, 0.037037037037037037037);
    DEFINVERSE(  28, 0.035714285714285714286);
    DEFINVERSE(  29, 0.034482758620689655172);
    DEFINVERSE(  30, 0.033333333333333333333);
    DEFINVERSE(  31, 0.032258064516129032258);
    DEFINVERSE(  32, 0.03125);
    DEFINVERSE(  33, 0.030303030303030303030);
    DEFINVERSE(  34, 0.029411764705882352941);
    DEFINVERSE(  35, 0.028571428571428571429);
    DEFINVERSE(  36, 0.027777777777777777778);
    DEFINVERSE(  37, 0.027027027027027027027);
    DEFINVERSE(  38, 0.026315789473684210526);
    DEFINVERSE(  39, 0.025641025641025641026);
    DEFINVERSE(  40, 0.025);
    DEFINVERSE(  41, 0.024390243902439024390);
    DEFINVERSE(  42, 0.023809523809523809524);
    DEFINVERSE(  44, 0.022727272727272727273);
    DEFINVERSE(  45, 0.022222222222222222222);
    DEFINVERSE(  46, 0.021739130434782608696);
    DEFINVERSE(  47, 0.021276595744680851064);
    DEFINVERSE(  48, 0.020833333333333333333);
    DEFINVERSE(  49, 0.020408163265306122449);
    DEFINVERSE(  50, 0.02);
    DEFINVERSE(  51, 0.019607843137254901961);
    DEFINVERSE(  52, 0.019230769230769230769);
    DEFINVERSE(  53, 0.018867924528301886792);
    DEFINVERSE(  54, 0.018518518518518518519);
    DEFINVERSE(  55, 0.018181818181818181818);
    DEFINVERSE(  56, 0.017857142857142857143);
    DEFINVERSE(  57, 0.017543859649122807018);
    DEFINVERSE(  58, 0.017241379310344827586);
    DEFINVERSE(  59, 0.016949152542372881356);
    DEFINVERSE(  60, 0.016666666666666666667);
    DEFINVERSE(  61, 0.016393442622950819672);
    DEFINVERSE(  62, 0.016129032258064516129);
    DEFINVERSE(  63, 0.015873015873015873016);
    DEFINVERSE(  64, 0.015625);
    DEFINVERSE(  65, 0.015384615384615384615);
    DEFINVERSE(  66, 0.015151515151515151515);
    DEFINVERSE(  67, 0.014925373134328358209);
    DEFINVERSE(  68, 0.014705882352941176471);
    DEFINVERSE(  69, 0.014492753623188405797);
    DEFINVERSE(  70, 0.014285714285714285714);
    DEFINVERSE(  71, 0.014084507042253521127);
    DEFINVERSE(  72, 0.013888888888888888889);
    DEFINVERSE(  73, 0.013698630136986301370);
    DEFINVERSE(  74, 0.013513513513513513514);
    DEFINVERSE(  75, 0.013333333333333333333);
    DEFINVERSE(  76, 0.013157894736842105263);
    DEFINVERSE(  77, 0.012987012987012987013);
    DEFINVERSE(  78, 0.012820512820512820513);
    DEFINVERSE(  79, 0.012658227848101265823);
    DEFINVERSE(  80, 0.0125);
    DEFINVERSE(  81, 0.012345679012345679012);
    DEFINVERSE(  82, 0.012195121951219512195);
    DEFINVERSE(  83, 0.012048192771084337349);
    DEFINVERSE(  84, 0.011904761904761904762);
    DEFINVERSE(  85, 0.011764705882352941176);
    DEFINVERSE(  86, 0.011627906976744186047);
    DEFINVERSE(  87, 0.011494252873563218391);
    DEFINVERSE(  88, 0.011363636363636363636);
    DEFINVERSE(  89, 0.011235955056179775281);
    DEFINVERSE(  90, 0.011111111111111111111);
    DEFINVERSE(  91, 0.010989010989010989011);
    DEFINVERSE(  92, 0.010869565217391304348);
    DEFINVERSE(  93, 0.010752688172043010753);
    DEFINVERSE(  94, 0.010638297872340425532);
    DEFINVERSE(  95, 0.010526315789473684211);
    DEFINVERSE(  96, 0.010416666666666666667);
    DEFINVERSE(  97, 0.010309278350515463918);
    DEFINVERSE(  98, 0.010204081632653061224);
    DEFINVERSE(  99, 0.010101010101010101010);
    DEFINVERSE( 100, 0.01);
    // all products of two numbers up to 16 and some more
    DEFINVERSE( 104, 0.0096153846153846153846);
    DEFINVERSE( 105, 0.0095238095238095238095);
    DEFINVERSE( 108, 0.0092592592592592592593);
    DEFINVERSE( 110, 0.0090909090909090909091);
    DEFINVERSE( 112, 0.0089285714285714285714);
    DEFINVERSE( 115, 0.0086956521739130434783);
    DEFINVERSE( 117, 0.0085470085470085470085);
    DEFINVERSE( 119, 0.0084033613445378151260);
    DEFINVERSE( 120, 0.0083333333333333333333);
    DEFINVERSE( 121, 0.0082644628099173553719);
    DEFINVERSE( 125, 0.008);
    DEFINVERSE( 126, 0.0079365079365079365079);
    DEFINVERSE( 128, 0.0078125);
    DEFINVERSE( 130, 0.0076923076923076923077);
    DEFINVERSE( 131, 0.0076335877862595419847);
    DEFINVERSE( 132, 0.0075757575757575757576);
    DEFINVERSE( 133, 0.0075187969924812030075);
    DEFINVERSE( 135, 0.0074074074074074074074);
    DEFINVERSE( 140, 0.0071428571428571428571);
    DEFINVERSE( 143, 0.0069930069930069930070);
    DEFINVERSE( 144, 0.0069444444444444444444);
    DEFINVERSE( 147, 0.0068027210884353741497);
    DEFINVERSE( 150, 0.0066666666666666666667);
    DEFINVERSE( 153, 0.0065359477124183006536);
    DEFINVERSE( 154, 0.0064935064935064935065);
    DEFINVERSE( 156, 0.0064102564102564102564);
    DEFINVERSE( 160, 0.00625);
    DEFINVERSE( 161, 0.0062111801242236024845);
    DEFINVERSE( 165, 0.0060606060606060606061);
    DEFINVERSE( 168, 0.0059523809523809523810);
    DEFINVERSE( 169, 0.0059171597633136094675);
    DEFINVERSE( 171, 0.0058479532163742690058);
    DEFINVERSE( 175, 0.0057142857142857142857);
    DEFINVERSE( 176, 0.0056818181818181818182);
    DEFINVERSE( 180, 0.0055555555555555555556);
    DEFINVERSE( 182, 0.0054945054945054945055);
    DEFINVERSE( 187, 0.0053475935828877005348);
    DEFINVERSE( 189, 0.0052910052910052910053);
    DEFINVERSE( 192, 0.0052083333333333333333);
    DEFINVERSE( 195, 0.0051282051282051282051);
    DEFINVERSE( 196, 0.0051020408163265306122);
    DEFINVERSE( 200, 0.005);
    DEFINVERSE( 207, 0.0048309178743961352657);
    DEFINVERSE( 208, 0.0048076923076923076923);
    DEFINVERSE( 209, 0.0047846889952153110048);
    DEFINVERSE( 210, 0.0047619047619047619048);
    DEFINVERSE( 216, 0.0046296296296296296296);
    DEFINVERSE( 224, 0.0044642857142857142857);
    DEFINVERSE( 220, 0.0045454545454545454545);
    DEFINVERSE( 221, 0.0045248868778280542986);
    DEFINVERSE( 225, 0.0044444444444444444444);
    DEFINVERSE( 231, 0.0043290043290043290043);
    DEFINVERSE( 240, 0.0041666666666666666667);
    DEFINVERSE( 247, 0.0040485829959514170040);
    DEFINVERSE( 250, 0.004);
    DEFINVERSE( 252, 0.0039682539682539682540);
    DEFINVERSE( 255, 0.0039215686274509803922);
    DEFINVERSE( 256, 0.00390625);
    // more numbers needed in fmm
    DEFINVERSE( 286, 0.0034965034965034965035);
    DEFINVERSE( 330, 0.0030303030303030303030);
    DEFINVERSE( 364, 0.0027472527472527472527);
    DEFINVERSE( 455, 0.0021978021978021978022);
    DEFINVERSE( 462, 0.0021645021645021645022);
    DEFINVERSE( 495, 0.0020202020202020202020);
    DEFINVERSE( 560, 0.0017857142857142857143);
    DEFINVERSE( 715, 0.0013986013986013986014);
    DEFINVERSE( 792, 0.0012626262626262626263);
    DEFINVERSE( 924, 0.0010822510822510822511);
    DEFINVERSE(1001, 0.00099900099900099900100);
    DEFINVERSE(1287, 0.00077700077700077700078);
    DEFINVERSE(1365, 0.00073260073260073260073);
    DEFINVERSE(1716, 0.00058275058275058275058);
    DEFINVERSE(1820, 0.00054945054945054945055);
    DEFINVERSE(2002, 0.00049950049950049950050);
    DEFINVERSE(3003, 0.00033300033300033300033);
    DEFINVERSE(3432, 0.00029137529137529137529);
    DEFINVERSE(4368, 0.00022893772893772893773);
    DEFINVERSE(5005, 0.00019980019980019980020);
    DEFINVERSE(6435, 0.00015540015540015540016);
    DEFINVERSE(8008, 0.00012487512487512487512);
    DEFINVERSE(11440,0.000087412587412587412587);
    DEFINVERSE(12870,0.000077700077700077700078);
    // further powers of two
    DEFINVERSE( 512, 0.001953125);
    DEFINVERSE(1024, 0.0009765625);
    DEFINVERSE(2048, 0.00048828125);
    DEFINVERSE(4096, 0.000244140625);
    DEFINVERSE(8192, 0.0001220703125);
    // some low factorials
    DEFINVERSE( 720, 0.0013888888888888888889);
    DEFINVERSE(5040, 0.00019841269841269841270);
    DEFINVERSE(40320,0.000024801587301587301587);
#undef DEFINVERSE
    /// meta programming arithmetic with small non-negative integers
    /// \note specialisations for N=0,1,2,3
    template<int N> struct Integer : public IntegerInverse<N> {
      using IntegerInverse<N>::Inverse;
      using IntegerInverse<N>::Divide;
      using IntegerInverse<N>::Ratio;
      /// modulus with 2
      static const int Odd = N & 2;
      /// half rounded down to integer
      static const int Half = N>>1;
      /// largest power of two not greater than N
      static const int LargestPow2 = 2*Integer<Half>::LargestPow2;
      /// logarithm of two rounded down to integer
      static const int Log2 = 1+Integer<Half>::Log2;
      /// multiply scalar
      /// \param[in,out] x scalar, replaced by N*x on output
      template<typename R>
      static void Multiply(R&x) noexcept { x *= R(N); }
      /// scalar to power N
      /// \param[in] x scalar
      /// \return x^N
      template<typename R>
      static R Power(R x) noexcept
      {
	R t = Integer<(N>>1)>::Power(x);
	return _m_EvenOdd<N&1>::powp(x, t*t);
      }
      /// replace x with x^N
      template<typename R> static void Exponentiate(R&x) noexcept
      { x = Power(x); }
      /// product of scalar with N
      /// \param[in] x scalar
      /// \return    x*N
      template<typename R>
      static R Product(R x) noexcept
      { return x*R(N); }
    };
    //
    template<> struct Integer<0> {
      static const int Odd = 0;
      static const int Half = 0;
      static const int Factorial = 1;
      template<typename R> static void Multiply    (R&x) noexcept
      { x=0; }
      template<typename R> static R    Product     (R  ) noexcept
      { return 0; }
      template<typename R> static void Exponentiate(R&x) noexcept
      { x=1; }
      template<typename R> static R    Power       (R  ) noexcept
      { return 1; }
      template<typename R> static R    PowerProduct(R, R Y) noexcept
      { return Y; }
    };
    //
    template<> struct Integer<1> : IntegerInverse<1> {
      static const int Odd = 1;
      static const int Half = 0;
      static const int LargestPow2 = 1;
      static const int Log2 = 0;
      static const int Factorial = 1;
      using IntegerInverse<1>::Inverse;
      using IntegerInverse<1>::Divide;
      using IntegerInverse<1>::Ratio;
      template<typename R>
      static void Multiply    (R& ) noexcept {}
      template<typename R>
      static R Product        (R x) noexcept { return x; }
      template<typename R>
      static void Exponentiate(R& ) noexcept {}
      template<typename R>
      static R Power          (R x) noexcept { return x; }
      template<typename R>
      static R PowerProduct(R x, R Y) noexcept { return x*Y; }
    };
    //
    template<> struct Integer<2> : IntegerInverse<2> {
      static const int Odd = 0;
      static const int Half = 1;
      static const int LargestPow2 = 2;
      static const int Log2 = 1;
      static const int Factorial = 2;
      using IntegerInverse<2>::Inverse;
      using IntegerInverse<2>::Divide;
      using IntegerInverse<2>::Ratio;
      template<typename R>
      static void Multiply    (R&x) noexcept { x+=x; }
      template<typename R>
      static R    Product     (R x) noexcept { return x+x; }
      template<typename R>
      static void Exponentiate(R&x) noexcept { x*=x; }
      template<typename R>
      static R    Power       (R x) noexcept { return x*x; }
    };
    //
    template<> struct Integer<3> : IntegerInverse<3> {
      static const int Odd = 1;
      static const int Half = 1;
      static const int LargestPow2 = 2;
      static const int Log2 = 1;
      static const int Factorial = 6;
      using IntegerInverse<3>::Inverse;
      using IntegerInverse<3>::Divide;
      using IntegerInverse<3>::Ratio;
      template<typename R>
      static void Multiply    (R&x) noexcept { x+=x+x; }
      template<typename R>
      static R    Product     (R x) noexcept { return x+x+x; }
      template<typename R>
      static void Exponentiate(R&x) noexcept { x*=x*x; }
      template<typename R>
      static R    Power       (R x) noexcept { return x*x*x; }
    };
    /// Ratio with integer: convert into real-valued multiplication
    template<int N, typename R> inline
    R Over(R x) noexcept { return Integer<N>::Ratio(x); }
    /// Product with integer: convert to sum for N=0,1,2,3
    template<int N, typename R> inline
    R Times(R x) noexcept { return Integer<N>::Product(x); }
#if(0)
    ////////////////////////////////////////////////////////////////////////////
    /// \name simple functors for assign-type operations
    //@{
    /// functor base class
    template<typename _lVal, typename _rVal>
    struct assign_function {
      typedef _lVal result_type;
      typedef _rVal argument_type;
    };
    /// =
    template<typename _Tp>
    struct assign : assign_function<_Tp&,_Tp>
    {
      _Tp& operator()(_Tp& _m_x, _Tp const&_m_y) const noexcept
      { return _m_x = _m_y; }
      static _Tp& operate(_Tp& _m_x, _Tp const&_m_y) noexcept
      { return _m_x = _m_y; }
    };
    /// +=
    template<typename _Tp>
    struct add_assign : assign_function<_Tp&,_Tp>
    {
      _Tp& operator()(_Tp& _m_x, _Tp const&_m_y) const noexcept
      { return _m_x += _m_y; }
      static _Tp& operate(_Tp& _m_x, _Tp const&_m_y) noexcept
      { return _m_x += _m_y; }
    };
    /// -=
    template<typename _Tp>
    struct subtract_assign : assign_function<_Tp&,_Tp>
    {
      _Tp& operator()(_Tp& _m_x, _Tp const&_m_y) const noexcept
      { return _m_x -= _m_y; }
      static _Tp& operate(_Tp& _m_x, _Tp const&_m_y) noexcept
      { return _m_x -= _m_y; }
    };
    /// *=
    template<typename _Tp>
    struct multiply_assign : assign_function<_Tp&,_Tp>
    {
      _Tp& operator()(_Tp& _m_x, _Tp const&_m_y) const noexcept
      { return _m_x *= _m_y; }
      static _Tp& operate(_Tp& _m_x, _Tp const&_m_y) noexcept
      { return _m_x *= _m_y; }
    };
    /// /=
    template<typename _Tp>
    struct divide_assign : assign_function<_Tp&,_Tp>
    {
      _Tp& operator()(_Tp& _m_x, _Tp const&_m_y) const noexcept
      { return _m_x /= _m_y; }
      static _Tp& operate(_Tp& _m_x, _Tp const&_m_y) noexcept
      { return _m_x /= _m_y; }
    };
    //@}
#endif
    ////////////////////////////////////////////////////////////////////////////
    //
    // struct ONE<N>
    //
    // F: factorial N!
    // G: N!!
    // H: (2*N-1)!!
    // P: 2^N
    //
    ////////////////////////////////////////////////////////////////////////////
    template<int N> struct ONE {
      enum {
	F  = N*ONE<N-1>::F,
	G  = N*ONE<N-2>::G,
	H  = (N+N-1)*ONE<N-1>::H,
	P  = 4*ONE<N-2>::P
      }; };
    template<> struct ONE<1> { enum { F=1, G=1, H=1, P=2 }; };
    template<> struct ONE<0> { enum { F=1, G=1, H=1, P=1 }; };
    ////////////////////////////////////////////////////////////////////////////
    //
    // struct TWO<L,M>
    //
    // B: binomial (L,M)
    // I: index of (L,M)
    //
    ////////////////////////////////////////////////////////////////////////////
    template<int L, int M> struct TWO {
      enum { B = TWO<L-1,M-1>::B + TWO<L-1,M>::B,
	     I = (L*(L+1))/2+M
      }; };
    template<int M> struct TWO<0,M> { enum { B=1, I=0 }; };
    template<int L> struct TWO<L,0> { enum { B=1, I=(L*(L+1))/2 }; };
    template<int L> struct TWO<L,L> { enum { B=1, I=(L*(L+3))/2 }; };
    template<>      struct TWO<0,0> { enum { B=1, I=0 }; };
    //
    /// array operations to unroll at compile time
    template<int I> struct ArrayLoop {
      /// f(a[i],x)
      template<typename X, typename AssignFunctor>
      static void assign(X*a, X x, AssignFunctor f)
      { ArrayLoop<I-1>::assign(a,x,f); f(a[I],x); }
      /// f(a[i],b[i])
      template<typename X, typename AssignFunctor>
      static void assign(X*a, const X*b, AssignFunctor f)
      { ArrayLoop<I-1>::assign(a,b,f); f(a[I],b[I]); }
    };
    template<> struct ArrayLoop<0> {
      template<typename X, typename AssignFunctor>
      static void assign(X*a, X x, AssignFunctor f) { f(a[0],x); }
      template<typename X, typename AssignFunctor>
      static void assign(X*a, const X*b, AssignFunctor f) { f(a[0],b[0]); }
    };
    //
    template<bool> struct _m_Bool {
      static bool OR (bool  ) { return 1; }
      static bool AND(bool x) { return x; }
    };
    template<> struct _m_Bool<0> {
      static bool OR (bool x) { return x; }
      static bool AND(bool  ) { return 0; }
    };
    /// boolean OR between compile-time and run-time argument
    /// \param[in] Y  boolean expression which may be ignored
    /// \return    @a X || @a Y
    template<bool X> bool Or(bool Y) { return _m_Bool<X>::OR(Y); }
    /// boolean AND between compile-time and run-time argument
    /// \param[in] Y  boolean expression which may be ignored
    /// \return    @a X && @a Y
    template<bool X> bool And(bool Y) { return _m_Bool<X>::AND(Y); }
    /// \name functors useful as template arguments
    //@{
    /// x=y
    struct assign {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x=y; }
    };
    /// x+=y
    struct add {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x+=y; }
    };
    /// x-=y
    struct subtract {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x-=y; }
    };
    /// x*=y
    struct multiply {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x*=y; }
    };
    /// x/=y
    struct divide {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x/=y; }
    };
    /// x <-> y
    struct swap {
      template<typename X> static void operate(X&x, X&y) noexcept
      { X tmp(x); x=y; y=tmp; }
    };
    /// x=sqrt(y)
    struct sqrt {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x=sqrt(y); }
    };
    /// x=y*y
    struct square {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x=y*y; }
    };
    /// x=1/y
    struct reciprocal {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x=X(1)/y; }
    };
    /// x=max(x,y)
    struct maximum {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x=max(x,y); }
    };
    /// x=min(x,y)
    struct minimum {
      template<typename X> static X& operate(X&x, X y) noexcept
      { return x=min(x,y); }
    };
#if(0)
    /// functor base class
    template<typename _arg1, typename _arg2, typename _result>
    struct assign_function {
      typedef _arg1 arg1_type;
      typedef _arg2 arg2_type;
      typedef _result result_type;
    };
    /// x=y
    template<typename T>
    struct assign : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x = y; }
    };
    /// x+=y
    template<typename T>
    struct add : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x += y; }
    };
    /// x-=y
    template<typename T>
    struct subtract : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x -= y; }
    };
    /// x*=y
    template<typename T>
    struct multiply : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x *= y; }
    };
    /// x/=y
    template<typename T>
    struct divide : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x /= y; }
    };
    /// x <-> y
    template<typename T>
    struct swap : assign_function<T&,T&,void>
    {
      void operator()(T&x, T&y) const noexcept 
      { T tmp(x); x=y; y=tmp; }
    };
    /// x=sqrt(y)
    template<typename T>
    struct sqrt : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x=sqrt(y); }
    };
    /// x=y*y
    template<typename T>
    struct square : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x=y*y; }
    };
    /// x=1/y
    template<typename T>
    struct reciprocal : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x=T(1)/y; }
    };
    /// x=max(x,y)
    template<typename T>
    struct maximum : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x=max(x,y); }
    };
    /// x=min(x,y)
    template<typename T>
    struct minimum : assign_function<T&,T,T&>
    {
      T& operator()(T& x, T const&y) const noexcept { return x=min(x,y); }
    };
#endif
    //@}
  } // namespace WDutils::meta
  ///
  /// \name methods supporting SFINAE
  ///
  //@{
#if __cplusplus >= 201103L
  /// is_equal<l,r>::value = l==r
  template<int Left, int Right>
  struct is_equal : std::integral_constant<bool, (Left==Right) > {};
  /// is_not_equal<l,r>::value = l!=r
  template<int Left, int Right>
  struct is_not_equal : std::integral_constant<bool, (Left!=Right) > {};
  /// is_less<l,r>::value =  l< r
  template<int Left, int Right>
  struct is_less : std::integral_constant<bool, (Left<Right) > {};
  /// is_greater<l,r>::value =  l> r
  template<int Left, int Right>
  struct is_greater : std::integral_constant<bool, (Left>Right) > {};
  /// is_less_or_equal<l,r>::value =  l<=r
  template<int Left, int Right>
  struct is_less_or_equal : std::integral_constant<bool, (Left<=Right) > {};
  /// is_greater_or_equal<l,r>::value =  l>=r
  template<int Left, int Right>
  struct is_greater_or_equal : std::integral_constant<bool, (Left>=Right) > {};
  /// is_in_range<x,l,r>::value =  l<=x && x<=r
  template<int X, int Low, int High>
  struct is_in_range : std::integral_constant<bool, (Low<=X && X<=High) > {};
  /// is_not_in_range<x,l,r>::value =  l<=x && x<=r
  template<int X, int Low, int High>
  struct is_not_in_range : std::integral_constant<bool, (High<X || X<Low) > {};
  //
  namespace detail {
    enum class enabler {};
  }
  /// nifty template magic due to R.Martinho Fernandes
  /// http://rmartinho.github.com/2012/05/29/type-traits-galore.html


  /// Then or Else, depending on If::value
  template <typename If, typename Then, typename Else>
  using Conditional = typename std::conditional<If::value,Then,Else>::type;
  /// dependent boolean type
  template <bool B, typename...>
  struct dependent_bool_type : std::integral_constant<bool, B> {};
// #ifndef __INTEL_COMPILER
  //  icpc 13 doesn't have variadic templates
  /// and an alias
  template <bool B, typename... T>
  using Bool = typename dependent_bool_type<B, T...>::type;
  /// meta-logical negation
  template <typename T>
  using Not = Bool<!T::value>;
  /// meta-logical disjunction
  template <typename... T>
  struct Any : Bool<false> {};
  template <typename Head, typename... Tail>
  struct Any<Head, Tail...> : Conditional<Head, Bool<true>, Any<Tail...>> {};
  /// meta-logical conjunction
  template <typename... T>
  struct All : Bool<true> {};
  template <typename Head, typename... Tail>
  struct All<Head, Tail...> : Conditional<Head, All<Tail...>, Bool<false>> {};
// #endif
  ///
  template <bool If, typename Then, typename Else>
  using Condition = typename std::conditional<If,Then,Else>::type;
  template <bool Condition, typename T = void>
  using EnableIf  = typename std::enable_if<Condition, T>::type;
  template <bool Condition, typename T = void>
  using DisableIf = typename std::enable_if<!Condition, T>::type;
#endif// __cplusplus >= 201103L
  //@} 
} // namespace WDutils
//
#endif //WDutils_included_meta_h
