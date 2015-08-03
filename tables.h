/*

Header file for S-boxes and MDS-matrices for the reference implementation of the Kupyna hash function (DSTU 7564:2014)

Authors: Ruslan Kiianchuk, Ruslan Mordvinov, Roman Oliynykov

*/


#ifndef KUPYNA_TABLES_H_
#define KUPYNA_TABLES_H_


#include "kupyna.h"

extern uint8_t mds_matrix[8][8];

extern uint8_t sboxes[4][256];

#endif  /* KUPYNA_TABLES_H_ */

