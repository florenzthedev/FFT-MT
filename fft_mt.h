// Copyright (c) 2023 Zachary Todd Edwards
// MIT License

#ifndef FFT_MT_H_INCLUDED
#define FFT_MT_H_INCLUDED

#include <complex.h>
#include <stdbool.h>

/**
 * @brief Performs the Fourier transform. This interface is intended to be used
 * with dynamic loading. This function handles putting the input into
 * bit-reversal-permutation order and spawning threads.
 *
 * @param X Input dataset, will be overwritten by results.
 * @param N Size of input dataset, must be a power of two greater than or equal
 * to 2.
 * @param aux the number of threads to spawn, if 0 two threads will be spawned,
 * if negative the reverse Fourier transform will be calculated.
 *
 * @return  0 on success, 1 on failure (usually malloc failure).
 */
int fourier_transform(double complex* X, long N, int aux);

// BELOW ARE INTERNAL FUNCTIONS

/**
 * @brief Struct containing the partition table information. Sets that are a
 * power of two in size can always be partitioned into two groups.
 *
 */
struct partition_s {
  long N;
  long portion_a;
  long portion_b;
  int count_a;
  int count_b;
};

/**
 * @brief Creates a radix-2 partition table.
 *
 * @param N The size of the input data.
 * @param threads The number of threads to spread the data across.
 * @param parts The struct to store the partition table in.
 */
void partition_pow2(long N, int threads, struct partition_s* parts);

/**
 * @brief Reorders the input dataset into bit-reversal-permutation for use with
 * the FFT.
 *
 * @param x Input dataset.
 * @param N Size of input dataset.
 */
void bit_reversal_permutation(double complex* x, long N);

/**
 * @brief The fast Fourier transform. This implementation uses a lookup table,
 * by changing the lookup table this can be used for the forward or inverse FFT.
 *
 * @param omegas The lookup table for first roots of unity from 2 to at least
 * log2(n).
 * @param X The input dataset, will be overwritten by results.
 * @param N The size of the dataset, must be a power of two in size.
 */
void fft(double complex omegas[], double complex X[], long N);

/**
 * @brief A single FFT butterfly pass, use for consolidating smaller parts other
 * nodes worked on.
 *
 * @param omegas The lookup table for first roots of unity from 2 to at least
 * log2(n).
 * @param X The input dataset, will be overwritten by results.
 * @param n The size of the dataset, must be a power of two in size.
 */
void fft_butterfly(double complex omegas[], double complex X[], long n);

#endif  // FFT_MT_H_INCLUDED