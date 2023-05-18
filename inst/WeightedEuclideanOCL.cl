__kernel void weighted_euclidean(
  __global numeric* output,
  const unsigned int count,
  __global numeric* input,
  __global numeric* weights,
  const int N,
  const int DIM){
  
  size_t i = get_global_id(0);
  size_t j = get_global_id(1);
  
  if(i < N && i == j){
    output[i*N+i] = 0;
  }
  
  if(i < N && j < i){
    double tmpRes = 0;
    for(int k = 0; k < DIM; ++k){
      tmpRes = tmpRes + weights[k] * (input[i+k*N] - input[j+k*N]) * (input[i+k*N] - input[j+k*N]);
    }
    tmpRes = sqrt(tmpRes);
    
    output[j*N+i] = tmpRes;
    output[i*N+j] = tmpRes;
  }
}
