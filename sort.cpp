#include <vector>

template <typename T>
void insertionSort(std::vector<T> *v,std::vector<T> *keys)
{
  const int size = v.size();
  for(int i=1;i<size;++i)
  {
    for(int j=0;j<i;++j)
    {
      if (keys[j]>keys[i])
      {
        const int tempv=v[j];
	const int tempkey=keys[j];
        v[j]=v[i];
	keys[j]=keys[i];
        for(int k=i;k>j;--k) { v[k]=v[k-1];keys[k]=keys[k-1]; }
        v[j+1]=tempv;
	keys[j+1]=tempkey;
      }
    }
  }
}

template <typename T>
void insertionSort(std::list<T> *v) {
   const int size = v.size();
   for (std::list<T>::iterator i=++v.begin();i!=v.end();i++) {
      for (std::list<T>::iterator j=v.begin();j!=i;j++) {
         if (j->key>i->key) {
            v.splice(j,v,i);
         }
      }
   }
}
