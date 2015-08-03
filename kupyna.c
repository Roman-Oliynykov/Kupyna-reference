/*

Reference implementation of the Kupyna hash function (DSTU 7564:2014), all hash code length variants

Authors: Ruslan Kiianchuk, Ruslan Mordvinov, Roman Oliynykov

*/


#include <memory.h>
#include <stdio.h>
#include "kupyna.h"
#include "tables.h"


int KupynaInit(size_t hash_nbits, kupyna_t* ctx) {
    if ((hash_nbits % 8 != 0) || (hash_nbits > 512)) {
        return -1;
    }
    if (hash_nbits <= 256) {
        ctx->rounds = NR_512;
        ctx->columns = NB_512;
        ctx->nbytes = STATE_BYTE_SIZE_512;
    } else {
        ctx->rounds = NR_1024;
        ctx->columns = NB_1024;
        ctx->nbytes = STATE_BYTE_SIZE_1024;
    }
    ctx->hash_nbits = hash_nbits;
    memset(ctx->state, 0, ctx->nbytes);
    // Set init value according to the specification.
    ctx->state[0][0] = ctx->nbytes;
    return 0;
}


void SubBytes(uint8_t state[NB_1024][ROWS], int columns) {
    int i, j;
    uint8_t temp[NB_1024];
    for (i = 0; i < ROWS; ++i) {
        for (j = 0; j < columns; ++j) {
            state[j][i] = sboxes[i % 4][state[j][i]];
        }
    }
}

void ShiftBytes(uint8_t state[NB_1024][ROWS], int columns) {
    int i, j;
    uint8_t temp[NB_1024];
    int shift = -1;
    for (i = 0; i < ROWS; ++i) {
        if ((i == ROWS - 1) && (columns == NB_1024)) {
            shift = 11;
        } else {
            ++shift;
        }
        for (j = 0; j < columns; ++j) {
            temp[(j + shift) % columns] = state[j][i];
        }
        for (j = 0; j < columns; ++j) {
            state[j][i] = temp[j];
        }
    }
}


uint8_t MultiplyGF(uint8_t x, uint8_t y) {
    int i;
    uint8_t r = 0;
    uint8_t hbit = 0;
    for (i = 0; i < BITS_IN_BYTE; ++i) {
        if ((y & 0x1) == 1)
            r ^= x;
        hbit = x & 0x80;
        x <<= 1;
        if (hbit == 0x80)
            x ^= REDUCTION_POLYNOMIAL;
        y >>= 1;
    }
    return r;
}

void MixColumns(uint8_t state[NB_1024][ROWS], int columns) {
    int i, row, col, b;
    uint8_t product;
    uint8_t result[ROWS];
    for (col = 0; col < columns; ++col) {
        memset(result, ROWS, 0);
        for (row = ROWS - 1; row >= 0; --row) {
            product = 0;
            for (b = ROWS - 1; b >= 0; --b) {
                product ^= MultiplyGF(state[col][b], mds_matrix[row][b]);
            }
            result[row] = product;
        }    
        for (i = 0; i < ROWS; ++i) {
            state[col][i] = result[i];
        }
    }
}


void AddRoundConstantP(uint8_t state[NB_1024][ROWS], int columns, int round) {
    int i;
    for (i = 0; i < columns; ++i) {
        state[i][0] ^= (i * 0x10) ^ round;
    }
}

void AddRoundConstantQ(uint8_t state[NB_1024][ROWS], int columns, int round) {
    int j;
    uint64_t* s = (uint64_t*)state;
    for (j = 0; j < columns; ++j) {
        s[j] = s[j] + (0x00F0F0F0F0F0F0F3ULL ^ 
                ((((columns - j - 1) * 0x10ULL) ^ round) << (7 * 8)));
    }
}


void P(kupyna_t* ctx, uint8_t state[NB_1024][ROWS]) {
    int i;
    for (i = 0; i < ctx->rounds; ++i) {
        AddRoundConstantP(state, ctx->columns, i);
        SubBytes(state, ctx->columns);
        ShiftBytes(state, ctx->columns);
        MixColumns(state, ctx->columns);
    }
}

void Q(kupyna_t* ctx, uint8_t state[NB_1024][ROWS]) {
    int i;
    for (i = 0; i < ctx->rounds; ++i) {
        AddRoundConstantQ(state, ctx->columns, i);
        SubBytes(state, ctx->columns);
        ShiftBytes(state, ctx->columns);
        MixColumns(state, ctx->columns);
    }
}


int Pad(kupyna_t* ctx, uint8_t* data, size_t msg_nbits) {
    int i;
    int mask;
    int pad_bit;
    int extra_bits;
    int zero_nbytes;
    size_t msg_nbytes = msg_nbits / BITS_IN_BYTE;
    size_t nblocks = msg_nbytes / ctx->nbytes;
    ctx->pad_nbytes = msg_nbytes - (nblocks * ctx->nbytes);
    ctx->data_nbytes = msg_nbytes - ctx->pad_nbytes;
    uint8_t* pad_start = data + ctx->data_nbytes;
    extra_bits = msg_nbits % BITS_IN_BYTE;
    if (extra_bits) {
        ctx->pad_nbytes += 1;
    }
    memcpy(ctx->padding, pad_start, ctx->pad_nbytes);
    extra_bits = msg_nbits % BITS_IN_BYTE;
    if (extra_bits) {
        mask = ~(0xFF >> (extra_bits));
        pad_bit = 1 << (7 - extra_bits);
        ctx->padding[ctx->pad_nbytes - 1] = (ctx->padding[ctx->pad_nbytes - 1] & mask) | pad_bit;
    } else {
        ctx->padding[ctx->pad_nbytes] = 0x80;
        ctx->pad_nbytes += 1;
    }
    zero_nbytes = ((-msg_nbits - 97) % (ctx->nbytes * BITS_IN_BYTE)) / BITS_IN_BYTE;
    memset(ctx->padding + ctx->pad_nbytes, 0, zero_nbytes);
    ctx->pad_nbytes += zero_nbytes;
    for (i = 0; i < (96 / 8); ++i, ++ctx->pad_nbytes) {
        if (i < sizeof(size_t)) {
            ctx->padding[ctx->pad_nbytes] = (msg_nbits >> (i * 8)) & 0xFF;
        } else {
            ctx->padding[ctx->pad_nbytes] = 0;
        }
    }
    return 0;
}


void Digest(kupyna_t* ctx, uint8_t* data) {
    int b, i, j;
    uint8_t temp1[NB_1024][ROWS];
    uint8_t temp2[NB_1024][ROWS];
    for (b = 0; b < ctx->data_nbytes; b += ctx->nbytes) {
        for (i = 0; i < ROWS; ++i) {
            for (j = 0; j < ctx->columns; ++j) {
                temp1[j][i] = ctx->state[j][i] ^ data[b + j * ROWS + i];
                temp2[j][i] = data[b + j * ROWS + i];
            }
        }
        P(ctx, temp1);
        Q(ctx, temp2);
        for (i = 0; i < ROWS; ++i) {
            for (j = 0; j < ctx->columns; ++j) {
                ctx->state[j][i] ^= temp1[j][i] ^ temp2[j][i];
            }
        }
    }
    /* Process extra bytes in padding. */
    for (b = 0; b < ctx->pad_nbytes; b += ctx->nbytes) {
        for (i = 0; i < ROWS; ++i) {
            for (j = 0; j < ctx->columns; ++j) {
                temp1[j][i] = ctx->state[j][i] ^ ctx->padding[b + j * ROWS + i];
                temp2[j][i] = ctx->padding[b + j * ROWS + i];
            }
        }
        P(ctx, temp1);
        Q(ctx, temp2);
        for (i = 0; i < ROWS; ++i) {
            for (j = 0; j < ctx->columns; ++j) {
                ctx->state[j][i] ^= temp1[j][i] ^ temp2[j][i];
            }
        }
    }
}


void Trunc(kupyna_t* ctx, uint8_t* hash_code) {
    int i;
    size_t hash_nbytes = ctx->hash_nbits / BITS_IN_BYTE;    
    memcpy(hash_code, (uint8_t*)ctx->state + ctx->nbytes - hash_nbytes, hash_nbytes);
}


void OutputTransformation(kupyna_t* ctx, uint8_t* hash_code) {
    int i, j;
    uint8_t temp[NB_1024][ROWS];
    memcpy(temp, ctx->state, ROWS * NB_1024);
    P(ctx, temp);
    for (i = 0; i < ROWS; ++i) {
        for (j = 0; j < ctx->columns; ++j) {
            ctx->state[j][i] ^= temp[j][i];
        }
    }
    Trunc(ctx, hash_code);
}


void KupynaHash(kupyna_t* ctx, uint8_t* data, size_t msg_bit_len, uint8_t* hash_code) {
    /* Reinitialize internal state. */
    memset(ctx->state, 0, ctx->nbytes);
    ctx->state[0][0] = ctx->nbytes;

    Pad(ctx, data, msg_bit_len);
    Digest(ctx, data);
    OutputTransformation(ctx, hash_code);
}

int KupynaKmac(kupyna_t* ctx, uint8_t* key, size_t digest_nbits, uint8_t* data, size_t msg_nbits, uint8_t* mac) {
    size_t total_nbytes;
    uint8_t* input;
    size_t i = 0;
    kupyna_t kpad;
    kupyna_t mpad;
    /* Reinitialize internal state. */
    memset(ctx->state, 0, ctx->nbytes);
    ctx->state[0][0] = ctx->nbytes;

    if (digest_nbits != 256 && digest_nbits != 384 && digest_nbits != 512) {
        /* Invalid key and digest size. */
        printf("Error: MAC has invalid bit size.");
        return -1;
    }

    KupynaInit(digest_nbits, &kpad);
    KupynaInit(digest_nbits, &mpad);

    Pad(&kpad, key, digest_nbits);
    Pad(&mpad, data, msg_nbits);

    total_nbytes = kpad.pad_nbytes + mpad.pad_nbytes + (digest_nbits / 8);
    if (kpad.data_nbytes > 0) {
        total_nbytes += kpad.data_nbytes;
    }
    if (mpad.data_nbytes > 0) {
        total_nbytes += kpad.data_nbytes;
    }
    input = calloc(total_nbytes, sizeof(uint8_t));

    if (kpad.data_nbytes > 0) {
        memcpy(input, key, digest_nbits / 8);
        i += digest_nbits / 8;
    }
    memcpy(&input[i], kpad.padding, kpad.pad_nbytes);
    i += kpad.pad_nbytes;
    if (mpad.data_nbytes > 0) {
        memcpy(&input[i], data, mpad.data_nbytes);
        i += mpad.data_nbytes;
    }
    memcpy(&input[i], mpad.padding, mpad.pad_nbytes);
    i += mpad.pad_nbytes;
    memcpy(&input[i], key, digest_nbits / 8);
    /* Invert key. */
    while(i < total_nbytes) {
        input[i] = input[i] ^ 0xFF;
        ++i;
    }

    KupynaHash(ctx, input, total_nbytes * 8, mac);
    free(input);
    return 0;
}
