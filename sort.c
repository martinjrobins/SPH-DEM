#include <vector>

template <typename T>
void InsertionSort(std::vector<T> *v,std::vector<T> *keys)
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
