#ifndef DBINDEXITEM_H
#define DBINDEXITEM_H

/*! \class  StEEmcDbIndexItem1
   Helper class to gather local  Db info about
   any ADC channel (PMT + MAPMT)
*/

#define StEEmcNameLen 16  // to avoid dependency on "cstructs/eemcConstDB.hh"

class StEEmcDbIndexItem1 {

 public:
  char name[StEEmcNameLen]; ///< ASCII name of the channel, see Readme 
  int crate, chan; ///< hardware channel
  float gain, hv; 
  float ped, sig;

  void clear();
  void print();
  void setName(char *text);

};

#endif

