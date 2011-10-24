Cdata::Cdata() {
   //make sure key type is correct
   if (size_of(keytype)<KEYSIZE) then
      printf("Size of keytype not big enough for KEYSIZE (%d)",KEYSIZE);
      exit();
   }
}


Cdata::construct(particle *in_ps,int in_n) {
    //point to particles
    ps = in_ps;
    n = in_n;
    //allocate memory for keys and hash tree
    for (int i=0;i<n;i++) {
        keys[i] = calc_key(ps[i].x);
    }
    insertionSort(ps,keys);
    //create root node of tree 
    Ctree_node *ptree_node = get_root_node();
    for (i=0;i<n;i++) {
        ptree_node = octinsert(&ps[i],ptree_node);
    }
}

Cdata::reconstruct() {
    hashTable.clear(); 
    insertionSort(ps,keys);
    Ctree_node *ptree_node = get_root_node();
    for (i=0;i<n;i++) {
        ptree_node = ptree_node.pinsert(&ps[i]);
    }
}

Ctree_node *Ctree_node::insert(particle *pp) {
    
   while (!ptree_node.is_in(pp)) {
       ptree_node = ptree_node.parent();
   }
   while (1) {
      if (ptree_node.is_leaf()) {
         ptree_node.add_particle(pp);
         return(ptree_node);
      }
	   else if ((pexist = ptree_node.get_particle())!= NULL) {
         ptree_node.remove_particle();
         dexist = ptree_node.daughter(pexist);
	      dp = ptree_node.daughter(pp);
	      while (dexist==dp) { 
            dexist = dexist.daughter(pexist);
	         dp = dp.daughter(pp);
	      }
	      dexist.add_particle(pexist);
	      dp.add_particle(pp);
	      return(dp);
      }
      ptree_node = ptree_node.daughter(pp);
   }
}

keytype Cdata::calc_key(vect x) {
   unsigned int xint[NDIM];
   keytype key = 0;
   for (i=0;i<NDIM;i++) {
      xint[i] = static_cast<unsigned int>((x[i]-xmin[i])*UINT_MAX/(xmax[i]-xmin[i]));
   }
   unsigned int mask = 0x80;
   for (i=0;i<(KEYSIZE-1)/NDIM;i++) {
      for (j=0;j<NDIM;j++) {
         key |= (xint[j] & mask) >> (sizeof(unsigned int)-1-i);
         key <<= 1;
         key |= xint[j] & 0x01;
         key >>= 1;
         xint[j] <<= 1;
      }
      mask >>= 1;
   } 
   return key; 
}

Ctree_node *Ctree_node::parent() {
   return hashTable[key>>NDIM];
}

Ctree_node *Ctree_node::daughter(particle *pp) {
   keytype newkey = key;
   for (int i=0;i<NDIM;i++) {
      newkey <<= 1;
      if (pp->x[i] > centre[i]) {
         newkey |= 0x01;
      }
   }
   return hashTable[newkey];
}

Ctree_node *Cdata::get_root_node() {
   keytype rkey = 0x1;
   return hashTable[rkey];
}

particle *Ctree_node::get_particle() {
   return p;
}

bool Ctree_node::is_leaf() {
   return leaf;
}

bool Ctree_node::is_in(particle *pp) {
   return pp->x > rmin && pp->x < rmax
} 

void Ctree_node::add_particle(particle *pp) {
   p = pp;
   leaf = FALSE;
}

template <void thefunct(&particle)>
Cdata::traverse() {
    for (int i=0;i<n;i++) {
       thefunct(&ps[i])
    }
}

template <void thefunct(&particle,&particle)>
Cdata::neighbours() {
   list<Ctree_ddnode *> nodes_to_search;  
   particle *p;
   for (int i=0;i<n;i++) {
      neighbrs[i].empty();
      //find neighbours
      nodes_to_search.add(get_root_node());
      while (!nodes_to_search.empty) {
         ptree_node = nodes_to_search.pop;
         //check if its a particle node
         if ((p = ptree_node.get_particle) != NULL) {
            neighbrs[i].add(p);
         }
         //check to see if it contains potential neighbours
         else if (ptree_node.has_neighbrs(&ps[i])) {
            Ctree_node *all_daughters[2^NDIM];
            ptree_node.find_daughters(all_daughters);
            for (int j=0;j<2^NDIM;j++) {
               nodes_to_search.add(all_daughters[j]);
            }
         }
      }
      //for each neighbours run function
      for (int j=0;j<neighbrs[i].length;j++) {
         thefunct(&ps[i],neighbrs[i][j]);
      }
   }
}

bool Ctree_node::has_neighbrs(particle *_p) {
   real _r2 = lensqr(_p->x-centre);
   return (_r2 <= (cellr + max(2*_p->h,2*cellh))^2);
}

void Ctree_node::find_daughters(Ctree_node daughters[2^NDIM]) {
   keytype newkey = key;
   newkey <<= NDIM;
   for (unsigned int i=0;i<2^NDIM;i++) {
      daughters[i] = hashTable[newkey|i];
   }
}
