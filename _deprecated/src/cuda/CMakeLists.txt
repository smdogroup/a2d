add_library(a2dcuda a2dcuda.cu)
set_target_properties(a2dcuda PROPERTIES CUDA_ARCHITECTURES all)