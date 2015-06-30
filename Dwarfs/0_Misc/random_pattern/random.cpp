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
    void* ptr;
    posix_memalign(&ptr, 64, size*sizeof(T));

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

  auto dt_in_sec= 0.0;

  const size_t nloop = 10;
  for (size_t iloop = 0; iloop < nloop; iloop++)
  {
    cout << "iloop[" << N << "]= " << iloop << " / " << nloop;
    const auto start = std::chrono::steady_clock::now();
    asm("#loopBeg");
#pragma omp parallel for schedule(static)
    for (int j = 0; j < ndata; j++)
    {
      if (!seq)
        for (int k = 0; k < 2; k++)
        {
          auto inew = idx[min(j+k,ndata-1)];
          auto hint = _MM_HINT_T0;
          auto ptr = (char*)&x[inew-1];
          for (int kk = 0; kk < sizeof(data_type)*3; kk += 64)
            _mm_prefetch(ptr+kk, hint);
          for (int kk = 0; kk < sizeof(data_type); kk += 64)
            _mm_prefetch(((char*)&y[j+k])+kk,_MM_HINT_T0);
        }
      auto i = idx[j];
      auto im1 = max(i-1,0);
      auto ip1 = min(i+1,ndata-1);
      auto const & xm1 = x[im1];
      auto const & xi  = x[i];
      auto const & xp1 = x[ip1];
#pragma simd
      for (size_t k = 0; k < N; k++)
        y[j].data[k] += xp1.data[k] - 2*xi.data[k] + xm1.data[k];
    }
    asm("#loopEnd");
    const auto end = std::chrono::steady_clock::now();
    auto dt = 1.0e-6*std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    dt_in_sec += dt;
    cout << " [ " << dt << " sec]\n";
#if 0
    if (!seq)
      std::random_shuffle(idx.begin(), idx.end());
#endif
  }

  

  dt_in_sec *= 1.0/nloop;
  const auto bw = ndata*(sizeof(data_type)*2+sizeof(int))/dt_in_sec/1e9;
  cout << "sizeof(data_type)= " << sizeof(data_type) << " B\n";
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

  auto nelems     = 8;
  auto ndata      = 4;
  auto sequential = false;
  auto use_single = false;

  auto params = pack(argc, argv, 
      param("data size in M", ndata, "n", "ndata"),
      param("number of elements", nelems, "e", "nelems"),
      param("use sequential", sequential, "seq",""),
      param("use fp32", use_single, "","fp32")
      );


  cerr << params.parse_all();

  ndata *= 1024*1024;

  if (!use_single)
    bw_test<double>(ndata, nelems,sequential);
  else
    bw_test<float>(ndata, nelems,sequential);

  return 0;
}

