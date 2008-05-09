#ifndef EMCdsm2Tree_h
#define EMCdsm2Tree_h
/**************************************************************
 * $Id: EMCdsm2Tree.h,v 1.1 2007/08/17 01:15:38 balewski Exp $
 * Emulates functionality of  Endcap DSM1-tree
 **************************************************************/
#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
class BEdsm2; 

class EMCdsm2Tree  {  // DSM0 tree emulators
 public:
  enum { Nbe2=4, Nbe2Cha=2, mxTxt=16}; // # of DSM2 for Endcap & Barrel
 private:
  BEdsm2 *be2; //  Level-2 Endcap[0], Barrel[1,2,3] boards
  char name[mxTxt];

  int mYear; // DSM algo changes w/ years
  int BEsumthr8bit, EEsumthr6bit, JPSIthrSelc2bit, BarreSide2bit, EtotThr8bit;
  int OutEndcapJP2bit;
  int OutEndcapHT2bit;
  int OutEndcapSum1bit; 
  int OutEndcapHTTP1bit;
  int OutEndcapTP1bit;
  int OutBarreJP2bit;
  int OutBarreHT2bit;
  int OutBarreSum1bit;
  int OutBarreJPSi1bit;
  int OutBarreHTTP1bit;
  int OutBarreTP1bit;
  int OutEtot1bit;

  int intBarreSum;
  int intEndcapSum;
  int intEtot;


 public:
  
  EMCdsm2Tree(char *);
  void setYear(int x, int BEsumthr, int EEsumthr, int JPSIthrSelc, int BarreSide, int EtotThr);
  ~EMCdsm2Tree(); 
  void  print(int k=0);
  void  clear();
  void setInput16bit(int ibrd, int ch,  ushort val); // words
  void compute();
  //  int getNboards() { return Nee1;}
  //int getInpTPsum(int ch /*0...11*/) ;// halfPatches, for both DSM1 boards

  //...... input..........
  int getInpHT2bit(int ibr, int ch); // ibr=0=Endcap
  int getInpHT2bit_2(int ibr, int ch); // ibr=0=Endcap
  int getInpTP1bit(int ibr, int ch); // ibr=0=Endcap
  int getInpHTTP1bit(int ibr, int ch); // ibr=0=Endcap
  int getInpEsum5bit(int ibr, int ch); // ibr=0=Endcap
  int getInpJP2bit(int ibr, int ch); // ibr=0=Endcap

  //........output........
  int getOutEndcapJP2bit(){return OutEndcapJP2bit;}
  int getOutEndcapHT2bit(){return OutEndcapHT2bit;}
  int getOutEndcapSum1bit(){return OutEndcapSum1bit;}
  int getOutEndcapHTTP1bit(){return OutEndcapHTTP1bit;}
  int getOutEndcapTP1bit(){return OutEndcapTP1bit;}
  int getIntEndcapSum(){return intEndcapSum;}

  int getOutBarreJP2bit(){return OutBarreJP2bit;}
  int getOutBarreHT2bit(){return OutBarreHT2bit;}
  int getOutBarreSum1bit(){return OutBarreSum1bit;}
  int getOutBarreJPSi1bit(){return OutBarreJPSi1bit;}
  int getOutBarreHTTP1bit(){return OutBarreHTTP1bit;}
  int getOutBarreTP1bit(){return OutBarreTP1bit;}
  int getIntBarreSum(){return intBarreSum;}

  int getOutEtot1bit(){return OutEtot1bit;}
  int getIntEtot(){return intEtot;}
};

#endif

/*
 * $Log: EMCdsm2Tree.h,v $
 * Revision 1.1  2007/08/17 01:15:38  balewski
 * full blown Endcap trigger simu, by Xin
 *
 *
 **************************************************************/

