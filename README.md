# Class encapsulating simple real matrices and fundamental matrix operations
This is the basic implementation of the fundamental matrix operations.

The motivation for this repository is to have a small replacement for Python's NumPy library in Java or C#.

Available in Java in OpenJDK compatible and in C# in Mono compatible format.

# Overview of available operations (with examples)
The overall logic is that indices of rows or columns are *indexed from zero* (0). 

*In-place* operation means that no new matrix is created and return as a result of 
the operation.

## Creating of the new matrix
### Zero matrix
For a zero matrix (each element equals to zero):
```
Matrix zeroMatrixInst = Matrix.generateZeroMatrix(int COLUMN, int ROWS);
```
### Random real matrix
For a matrix with random real numbers in interval [lower bound, upper bound):
```
Matrix randomMatrixInst = Matrix.generateRandomMatrix(int COLUMN, int ROWS, double LOWER_BOUND, double UPPER_BOUND);
```
### Random integer matrix
For a matrix with random integers in interval [lower bound, upper bound):
```
Matrix randomIntMatrixInst = Matrix.generateRandomIntegerMatrix(int COLUMN, int ROWS, int LOWER_BOUND, int UPPER_BOUND);
```
### One diagonal matrix (identity matrix)
For generating one matrix (all elements in diagonal equal to 1, all elements everywhere else 0):
```
Matrix oneMatrix = Matrix.generateOneDiagonalMatrix(int COLUMN, int ROWS);
```
### Constant matrix
For generating of the matrix that has each element equal to desired value (constant):
```
Matrix oneMatrix = Matrix.generateConstantMatrix(int COLUMN, int ROWS, double CONSTANT);
```
### Clonning of the matrix
For cloning of the existing matrix (deep cloning):
```
Matrix newMatrix = oldMatrix.deepCopy();
```
### Parsing of the Matlab or Wolfram Alpha strings
Convert the Matlab (or rather Python) or WolframAlpha string to the Matrix instance:
```
Matrix newMatrix = Matrix.parseMatlabString(String input);
Matrix newMatrix = Matrix.parseWolframString(String input);
```
## Operations on the row space
### Get rows count 
Return the number of rows in the matrix:
```
int rowCount = matrixInstance.getRowCount();
```

### Set the row value using vector
To set (in-place) the values in row at given position with new values given by input vector (position of the row is indexed from zero):
```
matrixInstance.setRow(int rowPosition, double [] vector);
matrixInstance.setRow(int rowPosition, Matrix vector);
```
### Swap two rows
To swap two rows at given position (in-place):
```
matrixInstance.swapRowsInplace(int firstRowIndex, int anotherRowIndex);
```
To swap two rows at given position with creating a new Matrix as a result:
```
Matrix swappedRows = matrixInstance.swapRows(int firstRowIndex, int anotherRowIndex);
```
### Row weight
The row weight is the special value that returns the possition of the first non-zero
element (indexed from left, starting with index zero):
```
int rowWeight = matrixInstance.getRowWeight(int rowIndex);
```
To receive of weights of all rows:
```
int [] rowWeights = matrixInstance.getRowsWeight();
```
### Selecting of the row
To selecting of the row (returns the new Matrix instance):
```
Matrix row = matrixInstance.chooseRow(int rowIndex);
```
To select of the multiple rows (defined by list of the row indices):
```
Matrix rows = matrixInstance.chooseRows(int [] rowIndices);
```
### Adding the multiplication of the chosen row to selected row
Standard step in Gaussian Elimination (skipped if row's positions are equals). Available as in-place operation or not in-place operation.
```
matrixInstance.addMultiplicationOfRowToAnotherRowInplace(int sourceRow, int targetRow, double multiplier);
Matrix afterOperation = matrixInstance.addMultiplicationOfRowToAnotherRow(int sourceRow, int targetRow, double multiplier);
```
### Multiply row by value
Multiply all elements in row by given value (essential for Gaussian Elimination):
```
matrixInstance.multiplyRowByValueInplace(int rowIndex, double multiplier);
Matrix afterOperation = matrixInstance.multiplyRowByValue(int rowIndex, double multiplier);
```
### Mean of the values in rows
To compute the mean of the values in row/rows use following (works for the one row, many rows, all matrix):
```
double [] meanOfSelectedRows = matrixInstance.meanOnRowSpace(int [] rowIndices);
double [] meanOfAllRows = matrixInstance.meanOnRowSpace();
double [] meanOfOneRow = matrixInstance.meanOnRowSpace(int rowIndex);
```

## Operations on the column space
### Get columns count 
Return the number of rows in the matrix:
```
int columnCount = matrixInstance.getColumnCount();
```

### Set the row value using vector
To set (in-place) the values in row at given position with new values given by input vector (position of the row is indexed from zero):
```
matrixInstance.setColumn(int columnPosition, double [] vector);
matrixInstance.setColumn(int columnPosition, Matrix vector);
```
### Swap two columns
To swap two columns at given position (in-place):
```
matrixInstance.swapColumnsInplace(int firstColumnIndex, int anotherColumnIndex);
```
To swap two columns at given position with creating a new Matrix as a result:
```
Matrix swappedColumns = matrixInstance.swapColumns(int firstColumnIndex, int anotherColumnIndex);
```

### Selecting of the column
To selecting of the column (returns the new Matrix instance):
```
Matrix column = matrixInstance.chooseColumn(int columnIndex);
```
To select of the multiple columns (defined by list of the column indices):
```
Matrix columns = matrixInstance.chooseColumns(int [] columnIndices);
```

### Mean of the values in columns
To compute the mean of the values in column/columns use following (works for the one column, many columns, all matrix):
```
double [] meanOfSelectedColumns = matrixInstance.meanOnColumnSpace(int [] columnIndices);
double [] meanOfAllColumns = matrixInstance.meanOnColumnSpace();
double [] meanOfOneColumn = matrixInstance.meanOnColumnSpace(int columnIndex);
```

## Matrix operations
Following operations are available for the whole matrix.

### Get the dimension (shape) of the matrix (first index is rows count, second columns count):
```
int [] dimension = matrixInstance.getDimension();
int [] dimension = matrixInstance.shape();
```
_shape is the same as dimension (it is there only because of compatibility with numpy)_

### Row echelon form (pivot form) of the matrix
To compute row echelon form (without swapping rows). The same opperations that are done to the source
matrix are reflected to mirror matrix in argument (if wanted).
```
matrixInstance.pivotMatrixInplace([Matrix mirror]);
Matrix reducedForm = matrixInstance.pivotMatrix([Matrix mirror]);
```
### Transpose matrix
Compute the matrix transposition ():
```
matrixInstance.transposeInplace();
Matrix transpose = matrixInstance.transpose();
```
### Matrix inversion
Get the inverse matrix:
```
Matrix inverse = matrixInstance.getInverseMatrix();
```
### Matrix multiplication (standard matrix multiplication)
Multiply two matrices (not available in-place):
```
Matrix result = matrixInstance.multiply(Matrix multiplier);
```
### Add matrix to another matrix
Add the values of another matrix (same dimension) to current matrix (in-place or not in-place):
```
matrixInstance.addInplace(Matrix adder);
Matrix result = matrixInstance.add(Matrix adder);
```
### Subract matrix to another matrix
Subract the values of another matrix (same dimension) to current matrix (in-place or not in-place):
```
matrixInstance.subtractInplace(Matrix subtracter);
Matrix result = matrixInstance.subtract(Matrix subtracter);
```
### Multiply matrix by constant
Multiply all elements in the matrix by given constant (there are two methods doing the same
kept just for naming compatibility).
```
matrixInstance.multiplyByConstantInplace(double value);
matrixInstance.multiplyScalarToEachElementInplace(double value);
Matrix result = matrixInstance.multiplyByConstant(double value);
Matrix result = matrixInstance.multiplyScalarToEachElement(double value);
```

### Divide matrix by constant
Divide each element in the matrix by constant value:
```
matrixInstance.divideScalarToEachElementInplace(double value);
Matrix result = matrixInstance.divideScalarToEachElement(double value);
```

### Hadamard product
Compute the Hadamard product (also known as the Schur product or the entrywise product) of two matrices. 
```
matrixInstance.hadamardProductInpace(Matrix mulitiplier);
Matrix result = matrixInstance.hadamardProduct(Matrix mulitiplier);
```

### Entrywise division
Divide each element in source matrix by corresponding element (at the same position) of the argument.
```
matrixInstance.entrywiseDivisionInplace(Matrix divider);
Matrix result = matrixInstance.entrywiseDivisionInplace(Matrix divider);
```

### Power to the scalar
Compute the power of each element in matrix to given scalar value
```
matrixInstance.powerToScalarInplace(double scalar);
Matrix result = matrixInstance.powerToScalarInplace(double scalar);
```

### Add constant to each element in the matrix
Add some constant to each element in the matrix.
```
matrixInstance.addScalarToEachElementInplace(double value);
Matrix result = matrixInstance.addScalarToEachElementInplace(double value);
```
### Subtract constant from each element in the matrix
Subtract some constant from each element in the matrix.
```
matrixInstance.subtractScalarToEachElementInplace(double value);
Matrix result = matrixInstance.subtractScalarToEachElement(double value);
```


### Get the upper triangular matrix
Get the upper triangular matrix of current matrix:
```
matrixInstance.getUpperTriangularInplace();
Matrix result = matrixInstance.getUpperTriangular();
```
### Get the trace of the matrix
Trace is the sum of all element on the diagonal of the matrix
```
double trace = matrixInstance.trace();
```
### Compute the determinant of the matrix
```
double determinant = matrixInstance.determinant();
```
### Compute the rank of the matrix
Rank is the number of linearly independent rows.
```
int rank = matrixInstance.rank();
```

### Compute the null space (kernel)
Find the kernel (null space) of the matrix:
```
Matrix nullspace = matrixInstance.nullSpace();
```


## Format conversions

### Get the 2D double array representation
Create the 2D (asymetric) array representation of the matrix:
```
double [][] matrix = matrixInstance.get2DArrayRepresentation();
```
### Convert the Matrix to the String representation
Convert the matrix to the string representation (compatible with Matlab and Python):
```
String result = matrixInstance.toString();
```
Convert to the string compatible with WolframAlpha (with rounding of values to given number of decimal places, default 4):
```
String result = matrixInstance.toWolframAlphaString([int decimalPrecision]);
```


## Operations with elements
### Get the element
To get element on given position:
```
double value = matrixInstance.getElement(int COLUMN, int ROW);
```
### Set the element
To set (in-place) element at given position:
```
matrixInstance.setElement(int COLUMN, int ROW, double VALUE);
```


# License

MIT License

Copyright (c) 2019 David Salac

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
