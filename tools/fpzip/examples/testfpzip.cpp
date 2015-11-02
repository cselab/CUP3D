#include <climits>
#include <cstdio>
#include <cstdlib>
#include <sys/time.h>
#include "fpzip.h"

#ifndef TYPE
  #error "must define TYPE as float or double"
#endif

#define string_(s) #s
#define string(s) string_(s)

static double
now()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  return tv.tv_sec + 1e-6 * tv.tv_usec;
}

int main(int argc, char* argv[])
{
  char* infile = 0;
  char* outfile = 0;
  unsigned dim[4] = { 1, 1, 1, 1 };
  switch (argc) {
    case 7:
      outfile = argv[6];
    case 6:
      if (sscanf(argv[5], "%u", &dim[3]) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 5:
      if (sscanf(argv[4], "%u", &dim[2]) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 4:
      if (sscanf(argv[3], "%u", &dim[1]) != 1)
        goto usage;
      /*FALLTHROUGH*/
    case 3:
      if (sscanf(argv[2], "%u", &dim[0]) != 1)
        goto usage;
      infile = argv[1];
      break;
    default:
    usage:
      fprintf(stderr, "Usage: %s <infile> <nx> [ny [nz [nf [outfile]]]]\n", *argv);
      return EXIT_FAILURE;
  }
  unsigned size = dim[0] * dim[1] * dim[2] * dim[3];
  unsigned inbytes = size * sizeof(TYPE);
  double t;
  int dp = (sizeof(TYPE) == sizeof(double));
  int* prec = new int[dim[3]];
  for (unsigned i = 0; i < dim[3]; i++)
    prec[i] = CHAR_BIT * sizeof(TYPE);

  fprintf(stderr, "testing %s data type\n", string(TYPE));

  // read raw data
  fprintf(stderr, "reading\n");
  FILE* file = fopen(infile, "rb");
  if (!file) {
    fprintf(stderr, "open failed\n");
    return EXIT_FAILURE;
  }
  TYPE* data = new TYPE[size];
  TYPE* copy = new TYPE[size];
  if (fread(data, sizeof(*data), size, file) != size) {
    fprintf(stderr, "read failed\n");
    return EXIT_FAILURE;
  }
  fclose(file);

  if (outfile) {
    // compress to file
    fprintf(stderr, "compressing to file\n");
    t = now();
    file = fopen(outfile, "wb");
    if (!file) {
      fprintf(stderr, "open failed\n");
      return EXIT_FAILURE;
    }
    unsigned outbytes = fpzip_file_write(file, data, prec, dp, dim[0], dim[1], dim[2], dim[3]);
    fclose(file);
    if (!outbytes) {
      fprintf(stderr, "compression failed\n");
      return EXIT_FAILURE;
    }
    t = now() - t;
    fprintf(stderr, "in=%u out=%u ratio=%.2f seconds=%.3f\n", inbytes, outbytes, (double)inbytes / outbytes, t);

    // decompress from file
    fprintf(stderr, "decompressing from file\n");
    t = now();
    file = fopen(outfile, "rb");
    if (!file) {
      fprintf(stderr, "open failed\n");
      return EXIT_FAILURE;
    }
    if (!fpzip_file_read(file, copy, prec, dp, dim[0], dim[1], dim[2], dim[3])) {
      fprintf(stderr, "decompression failed\n");
      return EXIT_FAILURE;
    }
    fclose(file);
    t = now() - t;
    fprintf(stderr, "seconds=%.3f\n", t);

    // validate data
    fprintf(stderr, "validating\n");
    for (unsigned i = 0; i < size; i++)
      if (data[i] != copy[i]) {
        fprintf(stderr, "validation failed at %u\n", i);
        return EXIT_FAILURE;
      }
  }

  // compress to memory
  fprintf(stderr, "compressing to memory\n");
  t = now();
  unsigned char* buffer = new unsigned char[size * sizeof(*data)];
  unsigned outbytes = fpzip_memory_write(buffer, size * sizeof(*data), data, prec, dp, dim[0], dim[1], dim[2], dim[3]);
  if (!outbytes) {
    fprintf(stderr, "compression failed\n");
    return EXIT_FAILURE;
  }
  t = now() - t;
  fprintf(stderr, "in=%u out=%u ratio=%.2f seconds=%.3f\n", inbytes, outbytes, (double)inbytes / outbytes, t);

  // decompress from memory
  fprintf(stderr, "decompressing from memory\n");
  t = now();
  if (!fpzip_memory_read(buffer, copy, prec, dp, dim[0], dim[1], dim[2], dim[3])) {
    fprintf(stderr, "decompression failed\n");
    return EXIT_FAILURE;
  }
  t = now() - t;
  fprintf(stderr, "seconds=%.3f\n", t);

  // validate data
  fprintf(stderr, "validating\n");
  for (unsigned i = 0; i < size; i++)
    if (data[i] != copy[i]) {
      fprintf(stderr, "validation failed\n");
      return EXIT_FAILURE;
    }

  // clean up
  delete[] prec;
  delete[] data;
  delete[] copy;
  delete[] buffer;

  fprintf(stderr, "OK\n");

  return 0;
}
