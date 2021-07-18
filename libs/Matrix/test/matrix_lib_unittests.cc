#include "gtest/gtest.h"
#include <Matrix/matrix_base_class.h>

namespace {

TEST(MatrixTests, DefaultConstructor) {
  MatrixBase<float> test_matrix;

  EXPECT_EQ(true,  test_matrix.is_empty());

}

}