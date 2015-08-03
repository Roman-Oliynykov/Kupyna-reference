/*

Header file for the reference implementation of the Kupyna hash function (DSTU 7564:2014), all hash length variants

Authors: Ruslan Kiianchuk, Ruslan Mordvinov, Roman Oliynykov

*/


#ifndef SRC_KUPYNA_H_
#define SRC_KUPYNA_H_

#include <stdlib.h>
#include <limits.h>


#define ROWS 8
#define NB_512 8  ///< Number of 8-byte words in state for <=256-bit hash code.
#define NB_1024 16  ///< Number of 8-byte words in state for <=512-bit hash code. 
#define STATE_BYTE_SIZE_512 (ROWS * NB_512)
#define STATE_BYTE_SIZE_1024 (ROWS * NB_1024)
#define NR_512 10  ///< Number of rounds for 512-bit state.
#define NR_1024 14  ///< Number of rounds for 1024-bit state.
#define REDUCTION_POLYNOMIAL 0x011d  /* x^8 + x^4 + x^3 + x^2 + 1 */

#if (ULLONG_MAX != 0xFFFFFFFFFFFFFFFFULL)
#error "Architecture not supported. Required type to fit 64 bits."
#endif

#define BITS_IN_WORD 64

#if (CHAR_BIT != 8)
#error "Architecture not supported. Required type to fit 8 bits."
#endif

#define BITS_IN_BYTE 8


typedef unsigned char uint8_t;
typedef unsigned long long uint64_t;
typedef struct {
    uint8_t state[NB_1024][ROWS];  ///< Hash function internal state (of maximum possible size to fit for all modes of operation).
    size_t nbytes;  ///< Number of bytes currently located in state.
    size_t data_nbytes;  ///< Number of bytes in input data sequence.
    uint8_t padding[STATE_BYTE_SIZE_1024 * 2];  ///< Space for extra bytes and padding.
    size_t pad_nbytes;  ///< Number of bytes currently located in padding buffer.
    size_t hash_nbits;  ///< Hash code bit length.
    int columns;  ///< Number of columns (8-byte vectors) located in internal state.
    int rounds;  ///< Number of rounds for current mode of operation.
} kupyna_t;


/*!
 * Initialize Kupyna hash function context.
 * Choose appropriate state size.
 * @param hash_hbits Bit size of the hash code to generate.
 * @param ctx Context of the Kupyna hash function.
 * @return Non-zero value in case of error, 0 in case of success.
 */
int KupynaInit(size_t hash_nbits, kupyna_t* ctx);


/*!
 * Generate hash code for input message of given bit length.
 * @param ctx Context of the Kupyna hash function.
 * @param data Pointer to input message bytes array to be digested.
 * @param msg_nbits Input message bit size.
 * @param hash_code Pointer to memory, allocated for storing generated hash
 * code.
 * It is user's responsibility to allocate enough memory for the hash code.
 * @return value
 */
void KupynaHash(kupyna_t* ctx, uint8_t* data, size_t msg_nbits, uint8_t* hash_code);

/*!
 * Generate message authentication code with key usage for input message of given bit length.
 * @param ctx Context of the Kupyna hash function.
 * @param key Pointer to input key bytes array.
 * @param digist_nbits Input key bit size.
 * @param data Pointer to input message bytes array to be digested.
 * @param msg_nbits Input message bit size.
 * @param digest Pointer to memory, allocated for storing generated message authentication code.
 * It is user's responsibility to allocate enough memory for the message authentication code.
 * @return return Non-zero value in case of error, 0 in case of success.
 */
int KupynaKmac(kupyna_t* ctx, uint8_t* key, size_t digest_nbits, uint8_t* data, size_t msg_nbits, uint8_t* digest);

#endif  /* SRC_KUPYNA_H_ */
