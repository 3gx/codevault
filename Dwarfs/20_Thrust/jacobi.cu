#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include "timer.h"


void jacobi_gold(
    thrust::host_vector<float> &output,
    const thrust::host_vector<float> &input, 
    const thrust::host_vector<float> &source, 
    float constant)
{
  const auto N = static_cast<int>(std::sqrt(input.size()));
  assert((int)input.size()== N*N);

  auto u = [&](int i, int j) -> float& { return output[j*N+i]; };
  auto g = [&](int i, int j) { return input[j*N+i]; };
  auto f = [&](int i, int j) { return source[j*N+i]; };


  //u(i,j) = 0.25*(g(i+1, j) + g(i-1, j) + g(i, j+1) + g(i, j-1) + constant * f(i,j));
  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
    {
      auto im = i-1;
      auto ip = i+1;
      auto jm = j-1;
      auto jp = j+1;
      if (im <  0) im += N;
      if (ip >= N) ip -= N;
      if (jm <  0) jm += N;
      if (jp >= N) jp -= N;
      u(i,j) = 0.25f*(g(ip,j)+g(im,j)+g(i,jp)+g(i,jm) + constant*f(i,j));
    }
}

bool diff(const thrust::host_vector<float> &x, const thrust::host_vector<float> &y)
{
  using namespace thrust;
  const auto n = x.size();
  for (size_t i = 0; i < n; i++)
  {
    const auto z = x[i] - y[i];
    if (abs(z) > 1.0e-7)
      return false;
  }
  return true;
}

struct jacobi_functor
{
  int N;
  float constant;
  thrust::device_ptr<float> input,source;
  thrust::device_ptr<float> output;
  __host__ __device__
  void operator()(int idx)
  {
    using namespace thrust;
    auto u = [&](int i, int j) { return output[j*N+i]; };
    auto g = [&](int i, int j) { return  input[j*N+i]; };
    auto f = [&](int i, int j) { return source[j*N+i]; };
    
    auto j = idx/N;
    auto i = idx - j*N;
    auto im = i-1;
    auto ip = i+1;
    auto jm = j-1;
    auto jp = j+1;
    if (im <  0) im += N;
    if (ip >= N) ip -= N;
    if (jm <  0) jm += N;
    if (jp >= N) jp -= N;
    u(i,j) = 0.25f*(g(ip,j)+g(im,j)+g(i,jp)+g(i,jm) + constant*f(i,j));
  }
};



int main(void)
{
  using namespace thrust;

  const auto N = 2048;
  host_vector<float> input(N*N), source(N*N), result(N*N);

  for (int j = 0; j < N; j++)
    for (int i = 0; i < N; i++)
    {
       input[j*N+i] = (i+j)&2;
      source[j*N+i] = (i+j)&4;
    }
  const auto constant = 0.25f;


  jacobi_gold(result, input,source,constant);

  int nrep = 10;

  Timer timer_gold  ("gold  ");
  Timer timer_thrust("thrust");
  for (int rep = 0; rep < nrep; rep++)
  {
    timer_gold.tbeg();
    jacobi_gold(result, input,source,rep/8.0f+constant);
    timer_gold.tend();
  }


  device_vector<float> d_input  = input;
  device_vector<float> d_source = source;
  device_vector<float> d_result = input;

  using CountingIterator = typename thrust::counting_iterator<int>;

  for_each(CountingIterator(0), CountingIterator(N*N),
      jacobi_functor{N,constant,
      d_input.data(),
      d_source.data(),
      d_result.data()});
      cudaThreadSynchronize();
  for (int rep = 0; rep < nrep; rep++)
  {
    timer_thrust.tbeg();
    for_each(CountingIterator(0), CountingIterator(N*N),
        jacobi_functor{N,rep/8.0f + constant,
        d_input.data(),
        d_source.data(),
        d_result.data()});
    cudaThreadSynchronize();
    timer_thrust.tend();
  }

  
  assert(diff(result, d_result));

  std::cout << "OK\n";

  timer_gold.finalize();
  timer_thrust.finalize();
  auto bw = [&](double dt, double ddt) { std::cout << " - BW= " << N*N*sizeof(float)*3/dt/1e9 << " GB/s "; };
  timer_gold  .print(bw);
  timer_thrust.print(bw);
  std::cout << " speedup: " << timer_gold.dtmean()/timer_thrust.dtmean() << "x \n";

  return 0;
}
