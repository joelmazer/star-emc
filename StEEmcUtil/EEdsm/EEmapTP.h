
#ifdef EEmapTP_USE
// hardocded map of TP for EEMC DSM-0

#define  EEnTPeta 3  // numbers TP patches in eta  (from center)
#define  EEnTPphi 30  // numbers TP patches in phi (clock wise)

struct  EEmapTP { 
  int JPid, TPid; // Jet Patch ID [1-6] & Trig Patch ID [1-15] according to Steve's numbering scheme 
  int brdIn,chIn; //DSM level-0 input board [1-9] & channel number 0-9]
  int chOut; // DSM level-0 output channel [0-11] 
}; 
//

static EEmapTP eeMapTP[EEnTPphi][EEnTPeta] ={
 { {1, 1, 7, 0, 8}, {1, 2, 7, 1, 8}, {1, 3, 8, 0, 9} }, 
 { {1, 4, 7, 2, 8}, {1, 5, 7, 3, 8}, {1, 6, 8, 1, 9} }, 
 { {1, 7, 7, 4, 8}, {1, 8, 7, 5, 8}, {1, 9, 8, 2, 9} }, 
 { {1,10, 7, 6, 8}, {1,11, 7, 7, 8}, {1,12, 8, 3, 9} }, 
 { {1,13, 7, 8, 8}, {1,14, 7, 9, 8}, {1,15, 8, 4, 9} }, 
 { {2, 1, 9, 0,11}, {2, 2, 9, 1,11}, {2, 3, 8, 5,10} }, 
 { {2, 4, 9, 2,11}, {2, 5, 9, 3,11}, {2, 6, 8, 6,10} }, 
 { {2, 7, 9, 4,11}, {2, 8, 9, 5,11}, {2, 9, 8, 7,10} }, 
 { {2,10, 9, 6,11}, {2,11, 9, 7,11}, {2,12, 8, 8,10} }, 
 { {2,13, 9, 8,11}, {2,14, 9, 9,11}, {2,15, 8, 9,10} }, 
 { {3, 1, 1, 0, 0}, {3, 2, 1, 1, 0}, {3, 3, 2, 0, 1} }, 
 { {3, 4, 1, 2, 0}, {3, 5, 1, 3, 0}, {3, 6, 2, 1, 1} }, 
 { {3, 7, 1, 4, 0}, {3, 8, 1, 5, 0}, {3, 9, 2, 2, 1} }, 
 { {3,10, 1, 6, 0}, {3,11, 1, 7, 0}, {3,12, 2, 3, 1} }, 
 { {3,13, 1, 8, 0}, {3,14, 1, 9, 0}, {3,15, 2, 4, 1} }, 
 { {4, 1, 3, 0, 3}, {4, 2, 3, 1, 3}, {4, 3, 2, 5, 2} }, 
 { {4, 4, 3, 2, 3}, {4, 5, 3, 3, 3}, {4, 6, 2, 6, 2} }, 
 { {4, 7, 3, 4, 3}, {4, 8, 3, 5, 3}, {4, 9, 2, 7, 2} }, 
 { {4,10, 3, 6, 3}, {4,11, 3, 7, 3}, {4,12, 2, 8, 2} }, 
 { {4,13, 3, 8, 3}, {4,14, 3, 9, 3}, {4,15, 2, 9, 2} }, 
 { {5, 1, 4, 0, 4}, {5, 2, 4, 1, 4}, {5, 3, 5, 0, 5} }, 
 { {5, 4, 4, 2, 4}, {5, 5, 4, 3, 4}, {5, 6, 5, 1, 5} }, 
 { {5, 7, 4, 4, 4}, {5, 8, 4, 5, 4}, {5, 9, 5, 2, 5} }, 
 { {5,10, 4, 6, 4}, {5,11, 4, 7, 4}, {5,12, 5, 3, 5} }, 
 { {5,13, 4, 8, 4}, {5,14, 4, 9, 4}, {5,15, 5, 4, 5} }, 
 { {6, 1, 6, 0, 7}, {6, 2, 6, 1, 7}, {6, 3, 5, 5, 6} }, 
 { {6, 4, 6, 2, 7}, {6, 5, 6, 3, 7}, {6, 6, 5, 6, 6} }, 
 { {6, 7, 6, 4, 7}, {6, 8, 6, 5, 7}, {6, 9, 5, 7, 6} }, 
 { {6,10, 6, 6, 7}, {6,11, 6, 7, 7}, {6,12, 5, 8, 6} }, 
 { {6,13, 6, 8, 7}, {6,14, 6, 9, 7}, {6,15, 5, 9, 6} }
};


#endif