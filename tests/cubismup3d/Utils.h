#ifndef CUBISMUP3D_TESTS_UTILS_H
#define CUBISMUP3D_TESTS_UTILS_H

#include <cstdio>
#include <string>

namespace cubismup3d {

void initMPI(int *argc, char ***argv);
void finalizeMPI(void);
std::string computeNumBlocksArg(int cells, int blockSize);

#define CUP_RUN_TEST(test) do { \
      if(!(test)()) { \
        fprintf(stderr, "Failed on test \"%s\"!\n", #test); \
        exit(1); \
      } \
    } while (0)
#define CUP_CHECK(condition, ...) do { \
      if (!(condition)) {  \
        fprintf(stderr, __VA_ARGS__); \
        exit(1); \
      } \
    } while (0);

}  // cubismup3d

#endif
