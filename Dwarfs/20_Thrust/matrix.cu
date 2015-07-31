#include <iostream>
#include <vector>
#include "cuda_runtime.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>


using namespace std;


template<size_t... I>
struct indexSeq {};


template<size_t N, size_t... S>
struct makeIndexSeqImpl : makeIndexSeqImpl<N-1,N-1,S...> {};
template<size_t... S>
struct makeIndexSeqImpl<0,S...>
{
    using type = indexSeq<S...>;
};
template<size_t N>
using makeIndexSeq = typename makeIndexSeqImpl<N>::type;




template <int r, int c>
class mymatrix
{
  public:
    __host__ __device__
      mymatrix(){};
    __host__ __device__
      ~mymatrix(){};


    template<typename  Tuple  ,int which, size_t... I>
      __host__ __device__ 
      auto unpack(Tuple t,indexSeq<I...>) -> std::initializer_list<double>
      {
        // NO IDEA HOW TO GET AROUND THE HARD CODING HERE
        auto foo = {data[I]=thrust::get<I>(thrust::get<which>(t))...};
        return foo; /* removes compiler warning */
      }


    template<typename Tuple ,int which, size_t... I>
      __host__ __device__ 
      auto pack(Tuple t,indexSeq<I...>) -> std::initializer_list<double>
      {
        // NO IDEA HOW TO GET AROUND THE HARD CODING HERE
        auto foo = {thrust::get<I>(thrust::get<which>(t))=data[I]...};
        return foo; /* removes compiler warning */
      }
    
    double data[r*c];
};


template<int rowsleft,int colsleft,int rowsright,int colsright>
  __host__ __device__ 
void multiply_mat(mymatrix<rowsleft,colsleft>&m1,mymatrix<rowsleft,colsleft>&m2,mymatrix<rowsleft,colsright>&m3)
{
  for (int row=0;row<rowsleft;row++)
    for (int col=0;col<colsright;col++)
    {
      m3.data[col+row*colsright]=0.;
      for (int sum_i=0;sum_i<colsleft;sum_i++)
        m3.data[col+row*colsright]+=m1.data[sum_i+row*colsleft]*m2.data[col+sum_i*colsright]; //THIS IS NICE (PROBABLY SOME BLAS STUFF WOULD BE EVEN NICER, THATS NOT THE POINT THO)
    }
};


template<int r1,int c1,int r2,int c2>
struct my_matrix_functor
{
  mymatrix<r1,c1> m1;
  mymatrix<r2,c2> m2;
  mymatrix<r1,c2> m3;
  my_matrix_functor(){};
  template<typename Tuple>
    __host__ __device__
    void operator()(Tuple t)
    {
      static_assert(r1*c2 == r2*c2, "Matrices are not of the same size");
      m1.template unpack<Tuple,0>(t,makeIndexSeq<r1*c1>());
      m2.template unpack<Tuple,1>(t,makeIndexSeq<r1*c1>());
      multiply_mat<r1,c1,r2,c2>(m1,m2,m3);
      m3.template pack<Tuple,2>(t,makeIndexSeq<r1*c1>());  
    }
};


template<size_t... I, typename T, typename Func>
void my_foreach(const T& matrices1_SoA, const T& matrices2_SoA, T& matrices3_SoA, Func func, indexSeq<I...>)
{
  thrust::for_each(
      thrust::make_zip_iterator(
        thrust::make_tuple(
          thrust::make_zip_iterator(thrust::make_tuple(matrices1_SoA[I].begin()...)),
          thrust::make_zip_iterator(thrust::make_tuple(matrices2_SoA[I].begin()...)),
          thrust::make_zip_iterator(thrust::make_tuple(matrices3_SoA[I].begin()...))
          )
        ),
      thrust::make_zip_iterator(
        thrust::make_tuple(
          thrust::make_zip_iterator(thrust::make_tuple(matrices1_SoA[I].end()...)),
          thrust::make_zip_iterator(thrust::make_tuple(matrices2_SoA[I].end()...)),
          thrust::make_zip_iterator(thrust::make_tuple(matrices3_SoA[I].end()...))
          )
        )
      , func);
}


int  main()
{
  //do  A*B -> C for 10000 matrices
  // size(A) = M*M
  // size(B) = M*M
  // size(C) = M*M
  
  constexpr auto M = size_t{2};


  vector<thrust::device_vector<double> > matrices1_SoA(M*M, thrust::device_vector<double>(10000,1.)); //left matrices
  vector<thrust::device_vector<double> > matrices2_SoA(M*M, thrust::device_vector<double>(10000,2.)); //right matrices
  vector<thrust::device_vector<double> > matrices3_SoA(M*M, thrust::device_vector<double>(10000,0.)); //result matrices


  //NO IDEA HOW THE CALL COULD BE PARAMETRIZED
  // use variadic templates...
  my_foreach(matrices1_SoA, matrices2_SoA, matrices3_SoA,  my_matrix_functor<M,M,M,M>(), makeIndexSeq<M*M>());


  cout<<"OUTPUT MATRIX:"<<endl;
  cout<<matrices3_SoA[0][0]<<" "<<matrices3_SoA[1][0]<<endl;
  cout<<matrices3_SoA[2][0]<<" "<<matrices3_SoA[3][0]<<endl;


  return 1;
}
