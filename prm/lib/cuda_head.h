#undef qualifier
#ifdef CUDA_IMPLEMENTATION
  #include <cuda.h>
  #define qualifier __host__ __device__
  #define cudaonly(x) x
  #define hostonly(x)
#else
  #define qualifier inline
  #define cudaonly(x)
  #define hostonly(x) x
#endif

#ifdef CUDA_IMPLEMENTATION
  #define sin_ sinf
  #define cos_ cosf
  #define sqrt_ sqrtf
#else
  #include <math.h>
  #define sin_ sin
  #define cos_ cos
  #define sqrt_ sqrt
#endif
