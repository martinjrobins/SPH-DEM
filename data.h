
#ifndef DATA_H
#define DATA_H

#include <vector>
#include <stack>
#include <ext/hash_map>
#include <functional>
namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later


#include "particle.h"
#include "sort.h"
#include "customConstants.h"

#ifdef _3D_
#define KEYSIZE 64
typedef unsigned long long keytype;
#define NO_HASH_BITS 16
#endif
#ifdef _2D_
#define KEYSIZE 32
typedef unsigned long keytype;
#define NO_HASH_BITS 8
#endif
#ifdef _1D_
#define KEYSIZE 32
typedef unsigned long keytype;
#define NO_HASH_BITS 8
#endif

using namespace std;

class Ctree_node {
    public:
        keytype key;
        vect COM;
        vect centre;
        double cellr;
        double cell2h;
        double mass;
        Cparticle *p;
        vect rmin;
        vect rmax;
        bool leaf;
        bool noInit;

       Ctree_node() {
          COM = 0;
          cellr = 0;
          cell2h = 0;
          mass = 0;
          p = NULL;
          rmin = 100000;
          rmax = -100000;
          leaf = true;
          noInit = true;
       }

       inline void add_particle(Cparticle &_p) {
          p = &_p;
          leaf = false;
       }
       void remove_particle();
       inline bool is_in(Cparticle &_p) {
          return all(_p.r > rmin) && all(_p.r < rmax);
       } 
       inline bool has_neighbrs(Cparticle &_p) {
         //cout << "Checking node "<<key<<" has neighbours. cellr="<<cellr<<"cell2h="<<cell2h<<endl;
         double _r2 = len2(_p.r-centre);
         return (_r2 <= pow(max(2*_p.h,cell2h),2));
       }
       inline void initNode(keytype _key, vect _centre, vect _rmin, vect _rmax) {
          key = _key;
          centre = _centre;
          rmin = _rmin;
          rmax = _rmax;
          noInit = false;
       }
       inline void updateNodeWithParticle(Cparticle &_p) {
          double _thisCellr = len(_p.r-centre);
          if (_thisCellr > cellr) { cellr = _thisCellr; }
          double _thisCell2h = _thisCellr+2*_p.h;
          if (_thisCell2h > cell2h) { cell2h = _thisCell2h; }
          mass += _p.mass;
       }
};

struct hashFunc {
   size_t operator()(keytype _key) const {
      return size_t(_key & ((2<<NO_HASH_BITS)-1));
   }
};
typedef Sgi::hash_map<keytype,Ctree_node,hashFunc,equal_to<keytype> > hashTableType;

class Cdata {
public:
      Cdata(vector<Cparticle> &in_ps);
      void newTimeStep();
      template <void Thefunct(Cparticle &)>
      void traverse();
      template <void Thefunct(Cparticle &,double)>
      void traverse(double);
      template <void Thefunct(Cparticle &), bool Iffunct(Cparticle &)>
      void traverse();
      template <void Thefunct(Cparticle &,Cparticle &)>
      void neighbours();
      template <void Thefunct(Cparticle &,Cparticle &), bool Iffunct(Cparticle &)>
      void neighbours();
      void printMe();

private:
      void printMe(Ctree_node *_ptn);
      void calc_neighbours(int i);
      Ctree_node *parent(Ctree_node *_tree_node);

      inline Ctree_node *daughter(Ctree_node *_ptree_node, Cparticle &_p) {
  
         unsigned int _mask = 0x01<<(KEYSIZE-1);
         if ((_ptree_node->key & _mask) != 0) {
            cerr << "Ctree_node::daughter. Reached bottom of possible tree! key="<<_ptree_node->key<<endl;;
            exit(-1);
         }
         keytype _newKey = _ptree_node->key;
         vect _newCentre = _ptree_node->centre;
         vect _newRmin = _ptree_node->rmin;
         vect _newRmax = _ptree_node->rmax;
   
         for (int i=0;i<NDIM;i++) {
            _newKey <<= 1;
            if (_p.r[i] > _ptree_node->centre[i]) {
               _newKey |= 0x01;
               _newCentre[i] = 0.5*(_ptree_node->centre[i]+_ptree_node->rmax[i]);
               _newRmin[i] = _ptree_node->centre[i];
            }
            else {
               _newCentre[i] = 0.5*(_ptree_node->centre[i]+_ptree_node->rmin[i]);
               _newRmax[i] = _ptree_node->centre[i];
            }
         }
         Ctree_node &_d = hashTable[_newKey];

         if (_d.noInit) {
            _d.initNode(_newKey,_newCentre,_newRmin,_newRmax);
         }
         //cout << "Finding daughter for node with key="<<_ptree_node->key<<". Found daughter with key="<<_d.key<<endl;
         return &_d;
      }
      
      inline vector<Ctree_node *> &find_daughters(Ctree_node *_ptree_node) {
         vector<Ctree_node *> *_pdaughters = new vector<Ctree_node *>();
         keytype _newkey = _ptree_node->key;
         _newkey <<= NDIM;
         for (unsigned int i=0;i<(1<<NDIM);i++) {
            hashTableType::iterator _pdaughter_node = hashTable.find(_newkey|i);
            if (_pdaughter_node != hashTable.end()) {
	            _pdaughters->push_back(&(_pdaughter_node->second));
            } 
         }
         return *_pdaughters;
      }

      void calcDomainLimits();
      keytype calc_key(vect r);
      Ctree_node *get_root_node();
      bool insert(Cparticle &_p); 
      hashTableType hashTable;
      vector<Cparticle> &ps;
      vect rmax;
      vect rmin;
      vector<keytype> keys;
      vector<vector<Cparticle *> > neighbrs;
      int n;
};

template <void T_thefunct(Cparticle &)>
void Cdata::traverse() {
    for (int i=0;i<n;i++) {
       T_thefunct(ps[i]);
    }
}

template <void T_thefunct(Cparticle &, double)>
void Cdata::traverse(double param) {
    for (int i=0;i<n;i++) {
       T_thefunct(ps[i],param);
    }
}

template <void Thefunct(Cparticle &,Cparticle &)>
void Cdata::neighbours() {
   for (int i=0;i<n;i++) {
      if (neighbrs[i].empty()) {
         calc_neighbours(i);
         //cout << "Found "<<neighbrs[i].size()<<" neighbours for particle "<<i<<endl;
      }
      //for each neighbours run function
      for (int j=0;j<neighbrs[i].size();j++) {
         Cparticle *_pthep = neighbrs[i][j];
         Thefunct(ps[i],*_pthep);
      }
   }
}

template <void Thefunct(Cparticle &),bool Iffunct(Cparticle &)>
void Cdata::traverse() {
    for (int i=0;i<n;i++) {
       if (Iffunct(ps[i])) {
          Thefunct(ps[i]);
       }
    }
}

template <void Thefunct(Cparticle &,Cparticle &),bool Iffunct(Cparticle &)>
void Cdata::neighbours() {
   for (int i=0;i<n;i++) {
      if (Iffunct(ps[i])) {
         if (neighbrs[i].empty()) {
            calc_neighbours(i);
            //cout << "Found "<<neighbrs[i].size()<<" neighbours for particle "<<i<<endl;
         }
         //for each neighbours run function
         for (int j=0;j<neighbrs[i].size();j++) {
            Cparticle *_pthep = neighbrs[i][j];
            Thefunct(ps[i],*_pthep);
         }
      }
   }
}
#endif
