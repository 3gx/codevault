#include <thrust/iterator/constant_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "timer.h"

typedef double value_type;

class varset {
public: 
	value_type v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;
	__host__ __device__ varset() : v0(0.0), v1(1.0), v2(1.5), v3(1.75), v4(1.875), v5(1.9375), v6(1.9688), v7(1.9922), v8(1.9961), v9(1.9980), v10(1.9990) { }
};

//conditionally compile for the GPU or the CPU
#ifdef GPU

typedef thrust::device_vector< value_type > state_type;
typedef thrust::device_vector< varset > var_type;

#else

typedef thrust::host_vector< value_type > state_type;
typedef thrust::host_vector< varset > var_type;

#endif

#if 0
#include <trove/block.h>

template<typename T, int s>
__global__ void test_block_copy(const T* x, T* r, int l) {
    typedef trove::array<T, s> s_ary;
    int global_index = threadIdx.x + blockIdx.x * blockDim.x;

    for(int index = global_index; index < l; index += gridDim.x * blockDim.x) {

        //The block memory accesses only function
        //correctly if the warp is converged. Here we check.
        if (trove::warp_converged()) {
            //Warp converged, indices are contiguous, call the fast
            //load and store
            s_ary d = trove::load_array_warp_contiguous<s>(x, index);
            trove::store_array_warp_contiguous(r, index, d);
        } else {
            //Warp not converged, call the slow load and store
            s_ary d = trove::load_array<s>(x, index);
            trove::store_array(r, index, d);
        }
    }
}
#endif


struct updatevars_functor {

   	template< class Tuple >
	__host__ __device__
	void operator()( Tuple tuple_in ) {
		// dynamics at dv_n/dt = v_{n-1} - v_n
		varset var      = thrust::get<0>(tuple_in);
		value_type dt   = thrust::get<1>(tuple_in);
    	    	varset varnew;

		varnew.v0  = var.v0  + dt*( var.v10 - var.v0  );
		varnew.v1  = var.v1  + dt*( var.v0  - var.v1  );
		varnew.v2  = var.v2  + dt*( var.v1  - var.v2  );
		varnew.v3  = var.v3  + dt*( var.v2  - var.v3  );
		varnew.v4  = var.v4  + dt*( var.v3  - var.v4  );
		varnew.v5  = var.v5  + dt*( var.v4  - var.v5  );
		varnew.v6  = var.v6  + dt*( var.v5  - var.v6  );
		varnew.v7  = var.v7  + dt*( var.v6  - var.v7  );
		varnew.v8  = var.v8  + dt*( var.v7  - var.v8  );
		varnew.v9  = var.v9  + dt*( var.v8  - var.v9  );
		varnew.v10 = var.v10 + dt*( var.v9  - var.v10 );
		
		thrust::get<0>(tuple_in) = varnew;
	}
};

int main( int arc , char* argv[] )
{	
	if (arc!=2) exit(-1);
	int N = atoi(argv[1]); // number of copies to simulation

	// initial condition
	var_type var(N);
	
	// numerical parameters
	value_type t  = 0;
	value_type dt = 0.001;
	value_type tf = 100;

	//setup
	thrust::constant_iterator<value_type> dtval(dt);

	// solve the ode
  {
    Timer ts_("simulation", Timer::verbose_destructor{});
    ts_.tbeg();
    Timer t_("step",Timer::verbose_destructor{});
    while ( t < tf ) {
      // update the variables
      t_.tbeg();
      thrust::for_each_n(
          thrust::make_zip_iterator(thrust::make_tuple(var.begin(),dtval)),
          N,
          updatevars_functor()
          );
      cudaDeviceSynchronize();
      t_.tend();
      t += dt;
    }
    ts_.tend();
  }

	
	return 0;
}
