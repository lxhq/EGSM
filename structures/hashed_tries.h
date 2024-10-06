#ifndef STRUCTURES_HASHED_TRIES_H
#define STRUCTURES_HASHED_TRIES_H

#include <cstdint>


struct HashedTries
{
    // buffer_ is used for num_candidates_,
    // compacted_vs_sizes_, compacted_vs_offs_,
    // and hash_keys_offs_.
    uint32_t *buffer_;

    // tries._num_candidates_,               Number of candidated data vertices for each query vertex      length = NUM_VQ
    // tries._compacted_vs_sizes_,           Number of candidates keys for each query edges                length = NUM_EQ * 2
    // tries._num_buckets_,                  Number of buckets for each query edge: _num_candidates_ / 8   length = NUM_EQ * 2
    // tries._hash_table_offs_,              Offset of each query edge's buckets                           length = NUM_EQ * 2, _hash_table_offs_[i + 1] = _hash_table_offs_[i] + _num_buckets_[i] * BUCKET_DIM;
    // compacted_vs_temp_,                   Candidate data vertices for each query vertex:                length =  NUM_VQ * C_MAX_L_FREQ, all candidates are stored in the column (C_MAX_L_FREQ)

    uint32_t *num_candidates_;
    uint32_t *compacted_vs_sizes_;
    uint32_t *num_buckets_;
    uint32_t *hash_table_offs_;

    // hash tables
    uint32_t *compacted_vs_;
    uint32_t *keys_[2];
    uint32_t *values_[2];
    uint32_t *neighbors_[2];

    HashedTries();

    __device__ uint32_t BIdx(
        uint32_t table_index, 
        uint32_t element_index) const;

    __host__ __device__ uint32_t CIdx(
        uint32_t trie_index, 
        uint32_t table_index, 
        uint32_t const_index) const;

    __device__ uint2 HashSearch(
        uint32_t table_index, 
        const uint32_t key) const;
};

extern __constant__ HashedTries tries;

#endif //STRUCTURES_HASHED_TRIES_H
