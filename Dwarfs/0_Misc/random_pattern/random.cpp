#include <iostream>
#include <parse_arguments.h>
#include <vector>
#include <chrono>
#include <algorithm>

template<size_t N_, typename T>
struct Data
{
  static constexpr auto N = N_;
  std::array<T,N> data;
  Data() {}
  Data(const T value)
  {
    for (auto& d : data)
      d = value;
  }
};

template<typename T>
struct numa_allocator
{
  using value_type = T;
  T* allocate(size_t size)
  {
    std::cout << "numa_allocator allocate(" << size << ") : " << size*sizeof(T)/1e6 << "MB \n";
    auto ptr = malloc(size * sizeof(T));
    auto p = static_cast<char*>(ptr);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < sizeof(T)*size; i++)
      p[i] = 0;
    return static_cast<T*>(ptr);
  }
  void deallocate(T* ptr, size_t size)
  {
    free(ptr);
  }
};
template<typename T>
using numa_vector = std::vector<T,numa_allocator<T>>;

template<size_t N, typename T>
void bw_test(int ndata, bool seq)
{
  using namespace std;
  using data_type = Data<N,T>;
  numa_vector<data_type> x(ndata), y(ndata);
  numa_vector<int> idx(ndata);

  std::iota(idx.begin(), idx.end(),0);
  if (!seq)
    std::random_shuffle(idx.begin(), idx.end());

#pragma omp parallel for schedule(static)
  for (int i = 0; i < ndata; i++)
  {
    x[i] = (ndata*0.5-i)*2.0/ndata;
    y[i] = 0.0;
  }

  auto start = std::chrono::steady_clock::now();

  const size_t nloop = 10;
  for (size_t iloop = 0; iloop < nloop; iloop++)
  {
    cout << "iloop= " << iloop << " / " << nloop << endl;
#pragma omp parallel for schedule(static)
    for (int j = 0; j < ndata; j++)
    {
      auto i = idx[j];
      auto im1 = max(i-1,0);
      auto ip1 = min(i+1,ndata-1);
      for (size_t k = 0; k < N; k++)
        y[j].data[k] = x[ip1].data[k] - 2*x[i].data[k] + x[im1].data[k];
    }
#pragma omp parallel for schedule(static)
    for (int i = 0; i < ndata; i++)
      x[i] = y[i];
  }

  
  const auto end = std::chrono::steady_clock::now();
  const auto dt_in_sec= 1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()/nloop;

  const auto bw = ndata*sizeof(data_type)*4/dt_in_sec/1e9;
  cout << "dt_per_loop= " << dt_in_sec << "sec\n";
  cout << "BW= " << bw << "GB/s\n";

}


template<typename T>
void bw_test(int ndata, int nelems, bool seq)
{
  switch(nelems)
  {
    case  1:  bw_test< 1,T>(ndata,seq); break;
    case  2:  bw_test< 2,T>(ndata,seq); break;
    case  4:  bw_test< 4,T>(ndata,seq); break;
    case  8:  bw_test< 8,T>(ndata,seq); break;
    case 12:  bw_test<12,T>(ndata,seq); break;
    case 16:  bw_test<16,T>(ndata,seq); break;
    case 20:  bw_test<20,T>(ndata,seq); break;
    case 24:  bw_test<24,T>(ndata,seq); break;
    case 28:  bw_test<28,T>(ndata,seq); break;
    case 32:  bw_test<32,T>(ndata,seq); break;
    default:
          std::cerr << "Unknown number of elements\n";
        exit(-1);
  }

}

int main(int argc, char * argv[])
{
#ifndef _OPENMP
  static_assert(false, "OpenMP is not defined. Please rebuild with OpenMP");
#endif

  using namespace std;
  using namespace parse_arguments;

  auto nelems     = 4;
  auto ndata      = 1;
  auto sequential = false;
  auto use_double = false;

  auto params = pack(argc, argv, 
      param("data size in M", ndata, "n", "ndata"),
      param("number of elements", nelems, "e", "nelems"),
      param("use sequential", sequential, "seq",""),
      param("use double", use_double, "d","double")
      );


  cerr << params.parse_all();

  ndata *= 1024*1024;

  if (use_double)
    bw_test<double>(ndata, nelems,sequential);
  else
    bw_test<float>(ndata, nelems,sequential);

  return 0;
}

