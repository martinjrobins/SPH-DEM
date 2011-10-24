#ifndef SORT_H
#define SORT_H

#include <vector>

template <typename T>
void insertionSort(std::list<T> &v) {
   list<CpInfo>::iterator second = v.begin();
   second++;
   for (list<CpInfo>::iterator i = second;i!=v.end();i++) {
      for (list<CpInfo>::iterator j=v.begin();j!=i;j++) {
         if (j->key > i->key) {
            v.splice(i,v,j);
         }
      }
   }
}

template <typename T1,typename T2>
void insertionSort(std::vector<T1> &v,std::vector<T2> &keys)
{
  const int size = v.size();
  for(int i=1;i<size;++i)
  {
    for(int j=0;j<i;++j)
    {
      if (keys[j]>keys[i])
      {
        T1 tempv=v[j];
	     T2 tempkey=keys[j];
        v[j]=v[i];
	keys[j]=keys[i];
        for(int k=i;k>j;--k) { v[k]=v[k-1];keys[k]=keys[k-1]; }
        v[j+1]=tempv;
	keys[j+1]=tempkey;
      }
    }
  }
}

template <typename T1,typename T2,typename T3>
void insertionSort(std::vector<T1> &v,std::vector<T2> &v2,std::vector<T3> &keys)
{
  const int size = v.size();
  for(int i=1;i<size;++i)
  {
    for(int j=0;j<i;++j)
    {
      if (keys[j]>keys[i])
      {
        T1 tempv=v[j];
	T2 tmpv2=v2[j];
	     T3 tempkey=keys[j];
	
        v[j]=v[i];
	v2[j]=v2[i];
	keys[j]=keys[i];
        for(int k=i;k>j;--k) { v[k]=v[k-1];v2[k]=v2[k-1];keys[k]=keys[k-1]; }
        v[j+1]=tempv;
	v2[j+1]=tmpv2;
	keys[j+1]=tempkey;
      }
    }
  }
}
 

#endif
