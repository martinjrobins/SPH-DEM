#include <blitz/tinyvec-et.h>
#include <iostream>
using namespace blitz;

int main() {
   TinyVector<double,2> test1(1);
   TinyVector<double,2> test2(2);

   cout << all(test1<test2);
   //if (test1<test2) cout << "first one worked";
   //test2[1] = 0.5;
   //if (test1<test2) cout << "second one failed";
   return 1;
}
    
   
