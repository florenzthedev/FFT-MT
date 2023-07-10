#include "fft_mt.h"

#include <assert.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdlib.h>

/**
 * @brief A struct containing all of the information needed by all threads. All
 * of these are unchanging after they are set so no locks are needed to read
 * from these. This does not necessarily apply to the variables they point-to
 * however. This way a bunch of pointers and structs don't need to be copied for
 * every thread.
 *
 */
struct global_args_s {
  sem_t* semaphores;
  double complex* X;
  double complex* omegas;
  struct partition_s parts;
};

/**
 * @brief Struct actually passed to each thread as its argument. The ID here is
 * independent from the OS or pthread ID, and is used to determine the domain of
 * each thread and which locks are needed when.
 *
 */
struct thread_args_s {
  int id;
  struct global_args_s* global;
};

#ifdef __GNUC__
#include <limits.h>
#define bit_length(x) (sizeof(long) * CHAR_BIT - __builtin_clzl(x))
#else
unsigned int bit_length(unsigned long x) {
  unsigned int bits = 0;
  while (x) {
    bits++;
    x >>= 1;
  }
  return bits;
}
#endif  // __GNUC__

/**
 * @brief Struct passed to the omega thread as its argument.
 *
 */
struct omega_args_s {
  bool inverse;
  long N;
};

/**
 * @brief Generates the first-nth-roots-of-unity lookup table. During this the
 * base fourier_transform function handles the bit-reversal-permutation.
 *
 * Originally this thread also handled building the communication tree for the
 * threads, but a more elegant solution was found and that is no longer needed.
 * This might not take long enough to warrant its own thread, more testing is
 * needed.
 *
 * @param p_args The omega_args_s struct with the arguments needed by
 * omega_thread.
 * @return void* -> double complex* The omegas lookup table.
 */
void* omega_thread(void* p_args) {
  struct omega_args_s* args = p_args;
  int Nth = (bit_length(args->N) - 1);  // Log2(N)
  double complex* omegas = malloc(sizeof(double complex) * Nth);
  if (omegas == NULL) pthread_exit(NULL);
  omegas[0] = CMPLX(-1, 0);
  omegas[1] = CMPLX(0, -1);
  for (int j = 2; j < Nth; j++) omegas[j] = csqrt(omegas[j - 1]);
  if (args->inverse)
    for (int j = 1; j < Nth; j++) omegas[j] = conj(omegas[j]);
  pthread_exit(omegas);
}

/**
 * @brief Threading function for the fast Fourier transform. Each thread based
 * on its ID and the partition table determines what size of data it is
 * initially responsible for and where. After performing the FFT on that data it
 * will double the size its responsible for and check to see if, based on its
 * ID, it is still needed. If it is not it will unlock its ID mutex. If it is it
 * will lock the ID for the first part of the next chunk of data it will try to
 * operate on, which will be the last to unlock for that data chunk size. Then
 * it will perform a FFT butterfly pass on that data and repeat. Eventually only
 * thread ID 0 will be left performing the final butterfly operation on the
 * whole dataset.
 *
 * @param p_args
 * @return void*
 */
void* fast_fourier_thread(void* p_args) {
  struct thread_args_s* args = p_args;
  long myPortion, myStart;

  // get our starting size and index
  if (args->id < args->global->parts.count_a) {
    myPortion = args->global->parts.portion_a;
    myStart = args->id * myPortion;
  } else {
    myPortion = args->global->parts.portion_b;
    myStart = (args->global->parts.count_a * args->global->parts.portion_a) +
              ((args->id - args->global->parts.count_a) * myPortion);
  }
  if (myPortion == 0) {
    pthread_exit(NULL);
  }

  fft(args->global->omegas, &args->global->X[myStart], myPortion);

  int id_offset = 1;
  while (((myStart + myPortion) % (myPortion * 2))) {
    myPortion *= 2;
    sem_wait(&args->global->semaphores[args->id + id_offset]);
    id_offset *= 2;
    fft_butterfly(args->global->omegas, &args->global->X[myStart], myPortion);
    if (myPortion == args->global->parts.N) break;
  }

  sem_post(&args->global->semaphores[args->id]);
  pthread_exit(NULL);
}

int fourier_transform(double complex* X, long N, int aux) {
  assert((N & (N - 1)) == 0);  // Must be a power of two
  assert(N > 1);               // Must be at least 2
  if (aux == 0) aux = 2;
  bool inverse = false;
  if (aux < 0) {
    aux = -aux;
    inverse = true;
  }

  // Generate lookup table and perform bit-reversal-permutation concurrently.
  struct omega_args_s omega_args;
  omega_args.N = N;
  omega_args.inverse = inverse;
  pthread_t omega_thread_id;
  pthread_create(&omega_thread_id, NULL, &omega_thread, &omega_args);
  bit_reversal_permutation(X, N);
  double complex* omegas;
  pthread_join(omega_thread_id, (void*)&omegas);
  if (omegas == NULL) return 1;

  // Get function arguments in order
  struct global_args_s global;
  partition_pow2(N, aux, &global.parts);
  if (global.parts.portion_b == 0)
    aux = global.parts.count_a;  // do not spawn threads with nothing to do!
  global.X = X;
  global.omegas = omegas;
  global.semaphores = malloc(sizeof(sem_t) * aux);
  if (global.semaphores == NULL) return 1;
  for (int j = 0; j < aux; j++) sem_init(&global.semaphores[j], 0, 0);
  pthread_t* threads = malloc(sizeof(pthread_t) * aux);
  struct thread_args_s* args = malloc(sizeof(struct thread_args_s) * aux);
  if (threads == NULL) return 1;
  if (args == NULL) return 1;

  // Spawn threads
  for (int j = 0; j < aux; j++) {
    args[j].id = j;
    args[j].global = &global;

    pthread_create(&threads[j], NULL, &fast_fourier_thread, &args[j]);
  }

  // Wait for threads to finish
  for (int j = aux; j > 0; j--) {
    pthread_join(threads[j - 1], NULL);
  }

  // For the inverse FFT one final division is needed
  if (inverse)
    for (int j = 0; j < N; j++) X[j] /= N;

  // Free data and exit
  for (int j = 0; j < aux; j++) sem_destroy(&global.semaphores[j]);
  free(global.semaphores);
  free(omegas);
  free(args);
  free(threads);
  return 0;
}

void fft_butterfly(double complex omegas[], double complex X[], long n) {
  double complex omega = omegas[bit_length(n) - 2];
  double complex root_of_unity = 1;
  for (long j = 0; j < n / 2; j++) {
    double complex product = root_of_unity * X[j + n / 2];
    X[j + n / 2] = X[j] - product;
    X[j] = X[j] + product;
    root_of_unity *= omega;
  }
}

void fft(double complex omegas[], double complex X[], long N) {
  assert((N & (N - 1)) == 0);  // Must be a power of two
  for (long n = 2; n <= N; n *= 2) {
    double complex omega = omegas[bit_length(n) - 2];
    double complex root_of_unity = 1;
    for (long j = 0; j < n / 2; j++) {
      for (long k = 0; k < N; k += n) {
        double complex product = root_of_unity * X[k + j + n / 2];
        X[k + j + n / 2] = X[k + j] - product;
        X[k + j] = X[k + j] + product;
      }
      root_of_unity *= omega;
    }
  }
}

void partition_pow2(long N, int threads, struct partition_s* parts) {
  assert((N & (N - 1)) == 0);                    // Must be a power of two
  int balance = 1 << (bit_length(threads) - 1);  // pow(2,floor(log2(threads)))
  long halfN = N >> 1;
  parts->N = N;
  if (balance >= halfN) {  // Everyone gets 2 or 0
    parts->count_a = halfN;
    parts->portion_a = 2;
    parts->count_b = threads - parts->count_a;
    parts->portion_b = 0;
    return;
  }
  parts->portion_a = N / balance;           // Only ever need two values
  parts->portion_b = parts->portion_a / 2;  // a & b
  int difference = threads - balance;       // num above a power of two
  parts->count_a = balance - difference;    // num of nodes with higher value
  parts->count_b = 2 * difference;          // num of nodes with lower value
}

#ifdef __clang__
#define bit_reverse(x, n) (__builtin_bitreverse64(x) >> (64 - n))
#else
unsigned long bit_reverse(unsigned long x, unsigned int n) {
  unsigned long r = 0;
  for (unsigned int i = 0; i < n; ++i) {
    if (x & (1 << i)) r |= 1 << ((n - 1) - i);
  }
  return r;
}
#endif  // __clang__

void bit_reversal_permutation(double complex* x, long N) {
  assert((N & (N - 1)) == 0);  // Must be a power of two

  // Don't forget bit_length is one indexed!
  long bl = bit_length(N) - 1;

  // We can skip the first and last index, they never need to be moved
  for (long i = 1; i < N - 1; i++) {
    long ri = bit_reverse(i, bl);
    if (i < ri) {
      double complex temp = x[i];
      x[i] = x[ri];
      x[ri] = temp;
    }
  }
}