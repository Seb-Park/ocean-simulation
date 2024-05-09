#include "halftone.h"
#include <cmath>
#include <vector>

halftone::halftone()
{

}

//void halftone::init(int n, int m, float density){
//    N = n;
//    M = m;
//   // std::vector<std::vector<int>> m_halftone;
//    // 1. intiialize pattern array with density
//    for (int i=0; i<N; i++){
//        for (int j=0; j<M; j++){
//            // determine whether to place a 1 or 0 here based on density

//            // generate random number 0-1
//            float r = ((float) rand() / (float) RAND_MAX);
//            if (r < density){
//                m_halftone[i][j] = 1;
//            } else {
//                m_halftone[i][j] = 0;
//            }
//        }
//    }

//}

//void halftone::apply_gaussian(){
//    for (int i=0; i<N; i++){
//        for (int j=0; j<M; j++){
//           m_halftone[i][j] = radial_gaussian(i, j);
//        }
//    }

//}

//float halftone::radial_gaussian(int i, int j){
//    float r = sqrt((i-.5*N)*(i-.5*N) + (j-.5*M)*(j-.5*M));
//    float result = exp(-(r*r) / (2.f*m_sigma*m_sigma));

//    return result;

//}

