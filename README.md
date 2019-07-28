# Class encapsulating simple real matrices and fundamental matrix operations
This is the basic implementation of the fundamental matrix operation. 

# Overview of available operations (with examples)
## Creating of the new matrix
For a zero matrix (each element equals to zero):
```
Matrix zeroMatrixInst = Matrix.generateZeroMatrix(COLUMN, ROWS);
```
For a matrix with random real numbers in interval [lower bound, upper bound):
```
Matrix randomMatrixInst = Matrix.generateRandomMatrix(COLUMN, ROWS, LOWER_BOUND, UPPER_BOUND);
```
For a matrix with random integers in interval [lower bound, upper bound):
```
Matrix randomIntMatrixInst = Matrix.generateRandomIntegerMatrix(COLUMN, ROWS, LOWER_BOUND, UPPER_BOUND);
```
For generating one matrix (all elements in diagonal equal to 1, all elements everywhere else 0):
```
Matrix oneMatrix = Matrix.generateOneDiagonalMatrix(COLUMN, ROWS);
```
For generating of the matrix that has each element equal to desired value (constant):
```
Matrix oneMatrix = Matrix.generateConstantMatrix(COLUMN, ROWS, CONSTANT);
```
For cloning of the existing matrix (deep cloning):
```
Matrix newMatrix = oldMatrix.deepCopy();
```
## Operations on the row space
