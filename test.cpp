#include <iostream>
int main() {
   int i = 1;
   int &j = i;
   int k = 2;
   j = k;
   std::cout << j;
   return 1;
}
