


Cdata::Cdata(vector<Cparticle> &in_ps) : ps(in_ps) {
    //make sure key type is correct
    if (sizeof(keytype)*8<KEYSIZE) {
       cout << "Size of keytype not big enough for KEYSIZE (" << KEYSIZE << "). sizeof(keytype)*8="<<sizeof(keytype)*8<<endl;
       exit(-1);
    }
    newTimeStep();
}

void Cdata::newTimeStep() {
   //point to particles
    hashTable.clear();
    n = ps.size();
    neighbrs.resize(n);
    cout << "n = "<<n<<endl;
    cout << "neighbrs size is "<<neighbrs.size()<<endl;
    calcDomainLimits();
    //keys.resize(n);
    //for (int i=0;i<n;i++) {
    //    keys[i] = calc_key(ps[i].r);
    //}
    //insertionSort<Cparticle,keytype>(ps,keys);
    //create root node of tree 
    for (int i=0;i<n;i++) {
        neighbrs[i].clear();
        insert(ps[i]);
    }
}    

void Cdata::calcDomainLimits() {
   rmin = 100000;
   rmax = -100000;
   for (int i=0;i<n;i++) {
      for (int j=0;j<NDIM;j++) {
         if (ps[i].r[j] < rmin[j]) { rmin[j] = ps[i].r[j]; }
         if (ps[i].r[j] > rmax[j]) { rmax[j] = ps[i].r[j]; }
      }
   }
   cout << "domain extents: rmin="<<rmin<<". rmax = "<<rmax<<endl;
}

void Cdata::calc_neighbours(int i) {
   stack<Ctree_node *> nodes_to_search;  
   Cparticle *_pp;
   cout << "finding neighbours for particle, r="<<ps[i].r<<" and tag="<<ps[i].tag<<endl;

   //find neighbours
   nodes_to_search.push(get_root_node());
   while (!nodes_to_search.empty()) {
      Ctree_node *_ptree_node = nodes_to_search.top();
      cout << "Checking for neighbours in node with key="<<_ptree_node->key<<"and centre="<<_ptree_node->centre<<endl;
      nodes_to_search.pop();
      //check if its a particle node
      if ((_pp = _ptree_node->p) != NULL) {
         cout << "its got a particle, r ="<<_pp->r<<"and tag="<<_pp->tag<<endl;
         if (ps[i].is_neighbr(*_pp)) {
            cout << "particle is a neighbour, adding to neighbrs list" << endl;
            neighbrs[i].push_back(_pp);
         }
         //check to see if it contains potential neighbours
      } else {
         vector<Ctree_node *> &_all_daughters = find_daughters(_ptree_node);
         for (int j=0;j<_all_daughters.size();j++) {
            if (_all_daughters[j]->has_neighbrs(ps[i])) {
               nodes_to_search.push(_all_daughters[j]);
            }
         }
      }
      cout << "all done for this one."<<endl<<endl;
   }
}

keytype Cdata::calc_key(vect r) {
   unsigned int rint[NDIM];
   keytype key = 0;
   for (int i=0;i<NDIM;i++) {
      rint[i] = static_cast<unsigned int>((r[i]-rmin[i])*UINT_MAX/(rmax[i]-rmin[i]));
      cout << "integer from dim "<<i<<" is "<<rint[i]<<endl;
   }
   unsigned int mask = 0x01<<(KEYSIZE-1);
   for (int i=0;i<(KEYSIZE-1)/NDIM;i++) {
      for (int j=0;j<NDIM;j++) {
         key <<= 1;
         key |= (rint[j] & mask) >> (sizeof(unsigned int)*8-1-i);
      }
      mask >>= 1;
   } 
   cout << "calculating key from r="<<r<<". Result is key="<<key<<endl;
   return key; 
}



Ctree_node *Cdata::get_root_node() {
   keytype _rkey = 0x1;
   Ctree_node &_rootNode = hashTable[_rkey];
   if (_rootNode.noInit) {
      _rootNode.initNode(_rkey,0.5*(rmax+rmin),rmin,rmax);
   }
   cout << "getting root node. This node has key="<<_rootNode.key<<" and centre="<<_rootNode.centre<<endl;
   return &_rootNode;
}



vector<Ctree_node *> &Cdata::find_daughters(Ctree_node *_ptree_node) {
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


Ctree_node *Cdata::parent(Ctree_node *_ptree_node) {
   return &(hashTable[_ptree_node->key>>NDIM]);
}

Ctree_node *Cdata::daughter(Ctree_node *_ptree_node, Cparticle &_p) {
  
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
   cout << "Finding daughter for node with key="<<_ptree_node->key<<". Found daughter with key="<<_d.key<<endl;
   return &_d;
}


void Ctree_node::updateNodeWithParticle(Cparticle &_p) {
      double _thisCellr = len(_p.r-centre);
      if (_thisCellr > cellr) { cellr = _thisCellr; }
      double _thisCell2h = _thisCellr+2*_p.h;
      if (_thisCell2h > cell2h) { cell2h = _thisCell2h; }
      mass += _p.mass;
}

bool Cdata::insert(Cparticle &_p) {
    
   Ctree_node *_ptree_node = get_root_node();
   cout << "inserting particle with r="<<_p.r<<" and tag="<<_p.tag<<endl;
   cout << "starting at root node with key="<<_ptree_node->key<<endl;
   while (1) {
      _ptree_node->updateNodeWithParticle(_p);

      Cparticle *_ppexist = _ptree_node->p;
      if (_ppexist != NULL) {
         Cparticle &_pexist = *_ppexist;
         cout << "this node has a particle. removing particle with r="<<_pexist.r<<"and tag="<<_pexist.tag<<endl;
         _ptree_node->p = NULL;
         Ctree_node *_pdexist = daughter(_ptree_node,_pexist);
         _pdexist->updateNodeWithParticle(_pexist);
         Ctree_node *_pdp = daughter(_ptree_node,_p);
         _pdp->updateNodeWithParticle(_p);
	 while (_pdexist==_pdp) { 
            cout << "both daughters are the same. going deeper...(key="<<_pdp->key<<")"<<endl;
            _pdexist->leaf = false;
            _pdexist = daughter(_pdexist,_pexist);
            _pdexist->updateNodeWithParticle(_pexist);       
            _pdp = daughter(_pdp,_p);
            _pdp->updateNodeWithParticle(_p);
         }
         cout << "found two separate daughters. orig key="<<_pdp->key<<"orig centre="<<_pdp->centre<<". exist key="<<_pdexist->key<<" exist centre="<<_pdexist->centre<<endl;
         _pdexist->add_particle(_pexist);
         _pdp->add_particle(_p);
         return true;
      } else if (_ptree_node->leaf) {
         cout << "this node is a leaf. adding particle..."<<endl;
         _ptree_node->add_particle(_p);
         return true;
      }
      _ptree_node = daughter(_ptree_node,_p);
      cout << "moving down tree... new key = "<<_ptree_node->key<<" new centre="<<_ptree_node->centre<<endl;
   }
}





bool Ctree_node::has_neighbrs(Cparticle &_p) {
   cout << "Checking node "<<key<<" has neighbours. cellr="<<cellr<<"cell2h="<<cell2h<<endl;
   double _r2 = len2(_p.r-centre);
   return (_r2 <= pow(max(2*_p.h,cell2h),2));
}



Ctree_node::Ctree_node() {
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

void Ctree_node::initNode(keytype _key, vect _centre, vect _rmin, vect _rmax) {
   key = _key;
   centre = _centre;
   rmin = _rmin;
   rmax = _rmax;
   noInit = false;
}


bool Ctree_node::is_in(Cparticle &_p) {
   return all(_p.r > rmin) && all(_p.r < rmax);
} 

void Ctree_node::add_particle(Cparticle &_p) {
   p = &_p;
   leaf = false;
}

