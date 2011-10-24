#ifndef CUSTOMSIMBASE_H
#define CUSTOMSIMBASE_H

#include "customConstants.h"
#include "dataLL.h"
#include "vect.h"

class CcustomSimBase {
   public:
      CcustomSimBase(CdataLL *_data,double _time): data(_data),time(_time) {}
      double getTime() { return time; }
      virtual void beforeStart(double newTime) { time = newTime; }
      virtual void beforeMiddle(double newTime) { time = newTime; }
      virtual void beforeEnd(double newTime) { time = newTime; }
      virtual void afterEnd(double newTime) { time = newTime; }
      virtual vect vFilter(const vect vin, const vect rin) { return vin; } 
      virtual vect rFilter(const vect rin) { return rin; }
      
   protected:
      double time;
      CdataLL *data;
};


#endif
