#include "vect.h"

double len2(vect inV) {
   return dot(inV,inV);
}

double len(vect inV) {
   return sqrt(len2(inV));
}
 
//double len2(vectInt inV) {
//   return dot(inV,inV);
//}

//double len(vectInt inV) {
//   return sqrt(len2(inV));
//}
