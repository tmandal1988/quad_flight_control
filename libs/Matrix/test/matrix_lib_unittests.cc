#include "gtest/gtest.h"
#include <Matrix/matrix_base_class.h>

namespace {

TEST(MatrixTests, TestConstructors) {

  // Create empty matrix	
  MatrixBase<float> test_matrix1;

  EXPECT_EQ(true,  test_matrix1.is_empty());
  EXPECT_EQ(0, test_matrix1.get_ncols());
  EXPECT_EQ(0, test_matrix1.get_nrows());

  // Create a matrix filled with 0s
  MatrixBase<float> test_matrix2(3, 5);

  EXPECT_EQ(5, test_matrix2.get_ncols());
  EXPECT_EQ(3, test_matrix2.get_nrows());
  for (size_t idx_r = 0; idx_r < test_matrix2.get_nrows(); idx_r++){
  	for(size_t idx_c = 0;  idx_c < test_matrix2.get_ncols(); idx_c++)
  		EXPECT_FLOAT_EQ(0.0f, test_matrix2.get_element(idx_r, idx_c));
  }

  // Create an identity matrix
  MatrixBase<float> test_matrix3(3, 3, "eye");
  EXPECT_EQ(3, test_matrix3.get_ncols());
  EXPECT_EQ(3, test_matrix3.get_nrows());
  for (size_t idx_r = 0; idx_r < test_matrix3.get_nrows(); idx_r++){
  	for(size_t idx_c = 0;  idx_c < test_matrix3.get_ncols(); idx_c++){
  		if (idx_r == idx_c)
  			EXPECT_FLOAT_EQ(1.0f, test_matrix3.get_element(idx_r, idx_c));
  		else
  			EXPECT_FLOAT_EQ(0.0f, test_matrix3.get_element(idx_r, idx_c));
  	}
  }

  // Copy constructor
  MatrixBase<float> test_matrix4 = test_matrix3;
  EXPECT_EQ(3, test_matrix4.get_ncols());
  EXPECT_EQ(3, test_matrix4.get_nrows());
  for (size_t idx_r = 0; idx_r < test_matrix4.get_nrows(); idx_r++){
  	for(size_t idx_c = 0;  idx_c < test_matrix4.get_ncols(); idx_c++){
  		if (idx_r == idx_c)
  			EXPECT_FLOAT_EQ(1.0f, test_matrix4.get_element(idx_r, idx_c));
  		else
  			EXPECT_FLOAT_EQ(0.0f, test_matrix4.get_element(idx_r, idx_c));
  	}
  }
}

}