/*
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
 */
using System;
namespace Matrix
{
    public class Matrix
    {
        // Representation of the matrix as 1D array (indexing logic: get(row, col) = matrix[colNumber * row + col])
        private double[] matrix;
        // The number of columns in the matrix
        private int columnNumber;
        // The number of rows in the matrix
        private int rowNumber;
        // The delat for comparison of integer for the logic: compare(a,b) = absolute_value(a - b) < DELTA
        private double DELTA_FOR_COMPARE;

        /**
        * Convert the double value to string with given number of decimal places.
        * @param number input number.
        * @param precision the precision (defined number of decimal places).
        * @return The string with rounded double value.
        */
        private static string doubleToString(double number, int precision)
        {
            string start = number.ToString();
            string result = "";
            for (int i = 0; i < start.Length; i++)
            {
                if (start[i] == '.' || start[i] == ',')
                {
                    result += start[i];
                    for (int j = i + 1; j < precision + i + 1; j++)
                    {
                        if (j < start.Length)
                        {
                            result += start[j];
                        }
                        else
                        {
                            result += '0';
                        }
                    }
                    break;
                }
                result += start[i];
            }
            return result;
        }

        /**
        * Sorting an array and swapping the mirror array in a same manner as an source array.
        * @param arr the input array to be sorted.
        * @param mirror the mirror array which elements are swapped same as an input array.
        */
        private static void sortWithMirror(int [] arr, int [] mirror)
        {
            int n = arr.Length;

            // Build heap (rearrange array) 
            for (int i = n / 2 - 1; i >= 0; i--)
                heapifyWithMirror(arr, mirror, n, i);

            // One by one extract an element from heap 
            for (int i = n - 1; i >= 0; i--)
            {
                // Move current root to end 
                int temp = arr[0];
                arr[0] = arr[i];
                arr[i] = temp;

                // Swap mirror
                int tempM = mirror[0];
                mirror[0] = mirror[i];
                mirror[i] = tempM;

                // call max heapify on the reduced heap 
                heapifyWithMirror(arr, mirror, i, 0);
            }
        }

        /**
        * Sub-method for sorting an array using heap sort
        * @param arr input array
        * @param mirror mirror array
        * @param n number of elements
        * @param i start
        */
        private static void heapifyWithMirror(int [] arr, int [] mirror, int n, int i)
        {
            int largest = i; // Initialize largest as root 
            int l = 2 * i + 1; // left = 2*i + 1 
            int r = 2 * i + 2; // right = 2*i + 2 

            // If left child is larger than root 
            if (l < n && arr[l] > arr[largest])
                largest = l;

            // If right child is larger than largest so far 
            if (r < n && arr[r] > arr[largest])
                largest = r;

            // If largest is not root 
            if (largest != i)
            {
                int swap = arr[i];
                arr[i] = arr[largest];
                arr[largest] = swap;
                // Swap mirror
                int swapM = mirror[i];
                mirror[i] = mirror[largest];
                mirror[largest] = swapM;

                // Recursively heapify the affected sub-tree 
                heapifyWithMirror(arr, mirror, n, largest);
            }
        }

        /**
         * Create an instance with defined delta for comparison of doubles.
         */
        private Matrix()
        {
            // Compare doubles in given delta;
            DELTA_FOR_COMPARE = 0.0001;
        }

        /**
         * Create an instance with zeros everywhere
         * @param columnNumber The required number of columns in the matrix
         * @param rowNumber  The required number of rows in the matrix
         */
        private Matrix(int columnNumber, int rowNumber) : this()
        {
            this.matrix = new double[columnNumber * rowNumber];
            this.columnNumber = columnNumber;
            this.rowNumber = rowNumber;
        }

        /**
         * Get the element in matrices on given position
         * @param column the column index
         * @param row the row index
         * @return element at given position
         */
        public double getElement(int column, int row) 
        {
            if (row < 0 || row > rowNumber)
            {
                throw new IndexOutOfRangeException("Index row out of bound!");
            }
            if (column < 0 || column > columnNumber)
            {
                throw new IndexOutOfRangeException("Index column out of bound!");
            }
            return matrix[row * columnNumber + column];
        }


        /**
        * Set the element at given position.
        * @param column the column position of the element
        * @param row the row position of the element
        * @param value the new value for the element
        */
        public void setElement(int column, int row, double value)
        {
            if (row < 0 || row > rowNumber)
            {
                throw new IndexOutOfRangeException("Index row out of bound!");
            }
            if (column < 0 || column > columnNumber)
            {
                throw new IndexOutOfRangeException("Index column out of bound!");
            }
            matrix[row * columnNumber + column] = value;
        }

        /**
         * Set the row for the new values given by array of double values.
         * @param row The row where an array is inserted
         * @param values The values for the row
         * @throws IndexOutOfRangeException in the case that row index is out of bound
         * @throws ArgumentException  in the case that dimension of the array does not fit to matrix
         */
        public void setRow(int row, double [] values) {
            if (row < 0 || row > rowNumber) {
                throw new IndexOutOfRangeException("Index row out of bound!");
            }
            if (values.Length != columnNumber) {
                throw new ArgumentException("The number of columns in the matrix has to match the length of the input vector!");
            }
            for (int col = 0; col < columnNumber; col++) {
                setElement(col, row, values[col]);
            }
        }

        /**
        * Set the row using row vector represented as Matrix class instance.
        * @param row the row where the new values are inserted.
        * @param values the vector of values.
        * @throws IndexOutOfRangeException in the case that row index is out of bound
        * @throws ArgumentException  in the case that dimension of the array does not fit to matrix
        */
        public void setRow(int row, Matrix values) {
            if (row< 0 || row> rowNumber) {
                throw new IndexOutOfRangeException("Index row out of bound!");
            }
            if (values.getRowCount() != 1 || values.getColumnCount() != columnNumber) {
                throw new ArgumentException("The number of columns in the matrix has to match the length of the input vector or the number of rows in input vector has to be just 1!");
            }
            for (int col = 0; col<columnNumber; col++) {
                setElement(col, row, values.getElement(col, 0));
            }
        }

        /**
         * Set the values of given column with input double values.
         * @param column The position of the column which will be rewritten.
         * @param values The new values that are inserted to the matrix
         * @throws IndexOutOfRangeException in the case that column index is out of bound
         * @throws ArgumentException in the case that dimension of the array does not fit to matrix
         */
        public void setColumn(int column, double[] values) {
            if (column< 0 || column> columnNumber) {
                throw new IndexOutOfRangeException("Index column out of bound!");
            }
            if (values.Length != rowNumber) {
                throw new ArgumentException("The number of rows in the matrix has to match the length of the input vector!");
            }
            for (int row = 0; row<rowNumber; row++) {
                setElement(column, row, values[row]);
            }
        }
        
        /**
         * Set the values of given column with input double values represented by instance of the Matrix class.
         * @param column The position of the column which will be rewritten
         * @param values The new values
         * @throws IndexOutOfRangeException in the case that column index is out of bound
         * @throws ArgumentException in the case that dimension of the array does not fit to matrix
         */
        public void setColumn(int column, Matrix values) {
            if (column< 0 || column> columnNumber) {
                throw new IndexOutOfRangeException("Index column out of bound!");
            }
            if (values.getColumnCount() != 1 || values.getRowCount()!= rowNumber) {
                throw new ArgumentException("The number of rows in the matrix has to match the length of the input vector or the number of columns in input vector has to be just 1!");
            }
            for (int row = 0; row<rowNumber; row++) {
                setElement(column, row, values.getElement(0, row));
            }
        }
        
        /**
         * Generate new matrix with random double values in given range
         * @param columnNumber Number of columns in the new matrix
         * @param rowNumber Number of rows in the new matrix
         * @param valuesFrom The lower bound for values
         * @param valuesTo The upper bound for values
         * @return The new matrix with given values in defined bound.
         */
        public static Matrix generateRandomMatrix(int columnNumber, int rowNumber, double valuesFrom, double valuesTo)
        {
            Matrix tempMatrix = generateZeroMatrix(columnNumber, rowNumber);
            for (int row = 0; row < rowNumber; row++)
            {
                for (int col = 0; col < columnNumber; col++)
                {
                    Random rnd = new Random();
                    tempMatrix.setElement(col, row, rnd.NextDouble() * (valuesTo - valuesFrom) + valuesFrom);
                }
            }
            return tempMatrix;
        }

        /**
     * Generate a new matrix with integer values in given range.
     * @param columnNumber Number of columns in the new matrix
     * @param rowNumber Number of rows in the new matrix
     * @param valuesFrom The lower bound for values
     * @param valuesTo The upper bound for values
     * @return The new matrix with given values in defined bound.
     */
        public static Matrix generateRandomIntegerMatrix(int columnNumber, int rowNumber, int valuesFrom, int valuesTo)
        {
            Matrix tempMatrix = generateZeroMatrix(columnNumber, rowNumber);
            for (int row = 0; row < rowNumber; row++)
            {
                for (int col = 0; col < columnNumber; col++)
                {
                    Random rnd = new Random();
                    tempMatrix.setElement(col, row, (int)(rnd.NextDouble() * (valuesTo - valuesFrom) + valuesFrom));
                }
            }
            return tempMatrix;
        }

        /**
         * Generate new matrix that has only zeros as elements
         * @param columnNumber The number of columns in a matrix.
         * @param rowNumber The number of rows in a matrix.
         * @return New matrix with given size.
         */
        public static Matrix generateZeroMatrix(int columnNumber, int rowNumber)
        {
            if (columnNumber < 1 || rowNumber < 1)
            {
                throw new ArgumentException("The number of columns and rows has to be greater than 1!");
            }
            return new Matrix(columnNumber, rowNumber);
        }

        /**
         * Generate the matrix that has element equal to 1 on diagonal and 0 everywhere else
         * @param columnNumber The number of columns in a matrix.
         * @param rowNumber The number of rows in a matrix.
         * @return New matrix with given size.
         */
        public static Matrix generateOneDiagonalMatrix(int columnNumber, int rowNumber)
        {
            Matrix tempMatrix = generateZeroMatrix(columnNumber, rowNumber);
            for (int i = 0; i < Math.Min(rowNumber, columnNumber); ++i)
            {
                tempMatrix.setElement(i, i, 1.0);
            }
            return tempMatrix;
        }

        /**
         * Generate the matrix that has element equal to given constant everywhere.
         * @param columnNumber The number of columns in a matrix.
         * @param rowNumber The number of rows in a matrix.
         * @param scalar The number that will be everywhere in the matrix.
         * @return New matrix with given size.
         */
        public static Matrix generateConstantMatrix(int columnNumber, int rowNumber, double scalar)
        {
            Matrix tempMatrix = generateZeroMatrix(columnNumber, rowNumber);
            for (int i = 0; i < tempMatrix.matrix.Length; ++i)
            {
                tempMatrix.matrix[i] = scalar;
            }
            return tempMatrix;
        }

        /**
         * Swap the rows in the matrix without creating a new Matrix.
         * @param rowA The position of the first row.
         * @param rowB The position of another row.
         */
        public void swapRowsInplace(int rowA, int rowB)
        {
            if (rowA < 0 || rowA >= rowNumber)
            {
                throw new IndexOutOfRangeException("Index rowA out of bound!");
            }
            if (rowB < 0 || rowB >= rowNumber)
            {
                throw new IndexOutOfRangeException("Index rowB out of bound!");
            }
            for (int col = 0; col < columnNumber; col++)
            {
                double temp = getElement(col, rowA);
                setElement(col, rowA, getElement(col, rowB));
                setElement(col, rowB, temp);
            }
        }

        /**
         * Swap the rows in the matrix together with creating a new Matrix for a result.
         * @param rowA The position of the first row.
         * @param rowB The position of another row.
         * @return The new matrix with swapped rows.
         */
        public Matrix swapRows(int rowA, int rowB)
        {
            Matrix temporaryMatrix = this.deepCopy();
            temporaryMatrix.swapRowsInplace(rowA, rowB);
            return temporaryMatrix;
        }

        /**
         * Swap the columns in the matrix without creating a new Matrix.
         * @param columnA The position of the first column.
         * @param columnB The position of another column.
         */
        public void swapColumnsInplace(int columnA, int columnB)
        {
            if (columnA < 0 || columnA >= columnNumber)
            {
                throw new IndexOutOfRangeException("Index columnA out of bound!");
            }
            if (columnB < 0 || columnB >= columnNumber)
            {
                throw new IndexOutOfRangeException("Index columnB out of bound!");
            }
            for (int row = 0; row < rowNumber; row++)
            {
                double temp = getElement(columnA, row);
                setElement(columnA, row, getElement(columnB, row));
                setElement(columnB, row, temp);
            }
        }

        /**
         * Swap the columns in the matrix together with creating a new Matrix for a result.
         * @param columnA The position of the first column.
         * @param columnB The position of another column.
         * @return The new matrix with swapped column.
         */
        public Matrix swapColumns(int columnA, int columnB)
        {
            Matrix temporaryMatrix = this.deepCopy();
            temporaryMatrix.swapColumnsInplace(columnA, columnB);
            return temporaryMatrix;
        }



        /**
         * Compare two double values with respect to given delta values.
         * @param A First double value
         * @param B Another double value
         * @param delta The delta parameters for comparison
         * @return abs(A-B) < delta
         */
        public static bool deltaDoubleCompare(double A, double B, double delta)
        {
            // Usable because of weakness of IEEE double implementation
            return Math.Abs(A - B) < delta;
        }

        /**
         * Get the weight of the row (the column position of the first non-zero element in row)
         * @param row Probed row
         * @return The position (from left) of the first non-zero element in the row or the number of rows if none such element exists.
         */
        public int getRowWeight(int row)
        {
            // Find the position (from left) of the first non-zero element in the matrix, return 0 if it does not exists.
            // Usable for implementation of Gauss Elimination
            for (int col = 0; col < columnNumber; col++)
            {
                if (!deltaDoubleCompare(0, getElement(col, row), DELTA_FOR_COMPARE))
                {
                    return col;
                }
            }
            return columnNumber;
        }

        /**
         * Get the weight of all rows in the matrix.
         * @return Weights of all rows in the matrix (weight represents the position of the first non-zero element in the row or the number of rows if none such element exists.)
         */
        public int[] getRowsWeight()
        {
            int[] weights = new int[rowNumber];
            for (int row = 0; row < rowNumber; row++)
            {
                weights[row] = getRowWeight(row);
            }
            return weights;
        }

        /**
         * Choose the column of the matrix at given position
         * @param column Index of desired column
         * @return The column on desired position represented as an Matrix type
         */
        public Matrix chooseColumn(int column)
        {
            Matrix tempCol = generateZeroMatrix(1, rowNumber);
            for (int row = 0; row < rowNumber; row++)
            {
                tempCol.setElement(0, row, getElement(column, row));
            }
            return tempCol;
        }

        /**
         * Select specified columns
         * @param columns An array with indices of required columns.
         * @return The columns on desired positions represented as an Matrix type
         */
        public Matrix chooseColumns(int[] columns)
        {
            Matrix temp = generateZeroMatrix(columns.Length, rowNumber);
            for (int col = 0; col < columns.Length; col++)
            {
                temp.setColumn(col, this.chooseColumn(columns[col]));
            }
            return temp;
        }

        /**
         * Choose the row of the matrix at given position
         * @param row Index of desired row
         * @return The row on desired position represented as an Matrix type
         */
        public Matrix chooseRow(int row)
        {
            Matrix tempRow = generateZeroMatrix(columnNumber, 1);
            for (int col = 0; col < columnNumber; col++)
            {
                tempRow.setElement(col, 0, getElement(col, row));
            }
            return tempRow;
        }

        /**
         * Select specified rows
         * @param rows An array with indices of required rows.
         * @return The rows on desired positions represented as an Matrix type
         */
        public Matrix chooseRows(int[] rows)
        {
            Matrix temp = generateZeroMatrix(columnNumber, rows.Length);
            for (int row = 0; row < rows.Length; row++)
            {
                temp.setColumn(row, this.chooseColumn(rows[row]));
            }
            return temp;
        }

        /**
         * Add the multiplication of the chosen row to selected row.
         * @param rowSource The index of the source row.
         * @param rowTarget The index of the target row.
         * @param multiplier The multiplication constant.
         */
        public void addMultiplicationOfRowToAnotherRowInplace(int rowSource, int rowTarget, double multiplier)
        {
            if (rowSource < 0 || rowSource >= rowNumber)
            {
                throw new IndexOutOfRangeException("Index rowSource out of bound!");
            }
            if (rowTarget < 0 || rowTarget >= rowNumber)
            {
                throw new IndexOutOfRangeException("Index rowTarget out of bound!");
            }
            if (rowSource == rowTarget)
            {
                return;
            }
            for (int col = 0; col < columnNumber; col++)
            {
                setElement(col, rowTarget, getElement(col, rowTarget) + multiplier * getElement(col, rowSource));
            }
        }

        /**
         * Add the multiplication of the chosen row to selected row and return a new object.
         * @param rowSource The index of the source row.
         * @param rowTarget The index of the target row.
         * @param multiplier The multiplication constant.
         * @return The matrix after adding the multiplication of the chosen row to selected row 
         */
        public Matrix addMultiplicationOfRowToAnotherRow(int rowSource, int rowTarget, double multiplier)
        {
            Matrix temp = deepCopy();
            temp.addMultiplicationOfRowToAnotherRowInplace(rowSource, rowTarget, multiplier);
            return temp;
        }

        /**
         * Compute the pivotal representation of the matrix without swapping rows
         * @param mirror The matrix where exactly the same operations are done (or null)
         */
        public void pivotMatrixInplace(Matrix mirror)
        {
            for (int i = 0; i < rowNumber; i++)
            {
                int[] rowsWeights = getRowsWeight();
                int[] indicesOfRows = new int[rowNumber];
                for (int idx = 0; idx < rowNumber; idx++)
                {
                    indicesOfRows[idx] = idx;
                }
                sortWithMirror(rowsWeights, indicesOfRows);
                /*for (int kk = 0; kk < rowsWeights.length; kk++) {
                    System.out.print(rowsWeights[kk] + ", ");
                }
                System.out.println(" ");*/
                int row = indicesOfRows[i];
                int rowWeight = rowsWeights[i];

                if (rowWeight >= columnNumber)
                {
                    break;
                }
                double valueOnPivot = getElement(rowWeight, row);
                double inverse = 1.0 / valueOnPivot;
                multiplyRowByValueInplace(row, inverse);
                if (mirror != null)
                {
                    mirror.multiplyRowByValueInplace(row, inverse);
                }
                //setElement(rowWeight, row, 1.0);

                // Subtract row value from each other row
                for (int j = 0; j < rowNumber; j++)
                {
                    double mulFactor = -1.0 * getElement(rowWeight, j);
                    addMultiplicationOfRowToAnotherRowInplace(row, j, mulFactor);
                    // For reducing of rounding error set directly to zero:
                    if (j != row)
                    {
                        setElement(rowWeight, j, 0.0);
                    }
                    else
                    {
                        setElement(rowWeight, j, 1.0);
                    }
                    if (mirror != null)
                    {
                        mirror.addMultiplicationOfRowToAnotherRowInplace(row, j, mulFactor);
                    }
                }
            }

        }

        /**
         * Compute the pivotal representation of the matrix without swapping rows.
         */
        public void pivotMatrixInplace()
        {
            pivotMatrixInplace(null);
        }

        /**
         * Compute the pivotal representation of the matrix without swapping rows and return modified matrix as a new object.
         * @param mirror The matrix where exactly the same operations are done (or null)
         * @return Matrix in pivotal form.
         */
        public Matrix pivotMatrix(Matrix mirror)
        {
            Matrix temp = deepCopy();
            temp.pivotMatrixInplace(mirror);
            return temp;
        }

        /**
         * Pivot the matrix and return it as a new object
         * @return Matrix in pivotal form.
         */
        public Matrix pivotMatrix()
        {
            return pivotMatrix(null);
        }

        /**
         * Transposition of the matrix
         */
        public void transposeInplace()
        {
            for (int row = 0; row < rowNumber; row++)
            {
                for (int col = 0; col < columnNumber; col++)
                {
                    setElement(row, col, getElement(col, row));
                }
            }
        }

        /**
         * Transposition of the matrix
         * @return Transposed matrix
         */
        public Matrix transpose()
        {
            Matrix temp = deepCopy();
            temp.transposeInplace();
            return temp;
        }

        /**
         * Compute the inversion of the matrix if exists
         * @return Inverse matrix to current matrix.
         */
        public Matrix getInverseMatrix()
        {
            if (rowNumber != columnNumber)
            {
                throw new ArgumentException("The matrix must be rectangular!");
            }
            Matrix copy = deepCopy();
            Matrix inverse = generateOneDiagonalMatrix(columnNumber, columnNumber);
            copy.pivotMatrixInplace(inverse);
            int[] rowWeights = copy.getRowsWeight();
            int[] rowIndices = new int[columnNumber];
            for (int i = 0; i < rowWeights.Length; i++)
            {
                rowIndices[i] = i;
            }
            sortWithMirror(rowWeights, rowIndices);
            Matrix inverseFind = generateZeroMatrix(columnNumber, rowNumber);
            for (int row = 0; row < rowNumber; row++)
            {
                /*System.out.println();*/
                inverseFind.setRow(row, inverse.chooseRow(rowIndices[row]));
            }
            return inverseFind;
        }

        /**
         * Compute the multiplication of the current matrix with given matrix in the argument.
         * @param multiplier Another matrix
         * @return current matrix multiplied by given matrix in the argument.
         */
        public Matrix multiply(Matrix multiplier)
        {
            if (columnNumber != multiplier.rowNumber)
            {
                throw new ArgumentException("The dimension of input matrices does not match!");
            }
            Matrix result = generateZeroMatrix(multiplier.columnNumber, rowNumber);

            for (int resRow = 0; resRow < result.rowNumber; resRow++)
            {
                for (int resCol = 0; resCol < result.columnNumber; resCol++)
                {
                    double mulRes = 0.0;
                    for (int i = 0; i < columnNumber; i++)
                    {
                        mulRes += this.getElement(i, resRow) * multiplier.getElement(resCol, i);
                        if (deltaDoubleCompare(mulRes, 0, DELTA_FOR_COMPARE))
                        {
                            mulRes = 0;
                        }
                        if (deltaDoubleCompare(mulRes, 1, DELTA_FOR_COMPARE))
                        {
                            mulRes = 1;
                        }
                    }
                    result.setElement(resCol, resRow, mulRes);
                }
            }
            return result;
        }

        /**
         * Multiply all element in given row by value.
         * @param row Index of the row
         * @param value The multiplication constant
         */
        public void multiplyRowByValueInplace(int row, double value)
        {
            if (row < 0 || row >= rowNumber)
            {
                throw new IndexOutOfRangeException("Index row out of bound!");
            }
            for (int col = 0; col < columnNumber; col++)
            {
                setElement(col, row, value * getElement(col, row));
            }
        }
        /**
         * Multiply all element in given row by value and return result as a new matrix
         * @param row Index of the row
         * @param value The multiplication constant
         * @return New matrix after operation
         */
        public Matrix multiplyRowByValue(int row, double value)
        {
            Matrix temp = deepCopy();
            temp.multiplyRowByValueInplace(row, value);
            return temp;
        }

        /**
         * Add the values of another matrix to this matrix.
         * @param adder The matrix to be added
         */
        public void addInplace(Matrix adder)
        {
            if (columnNumber != adder.getColumnCount() || rowNumber != adder.getRowCount())
            {
                throw new ArgumentException("The dimension of input matrices does not match!");
            }
            for (int row = 0; row < rowNumber; row++)
            {
                for (int col = 0; col < columnNumber; col++)
                {
                    setElement(col, row, this.getElement(col, row) + adder.getElement(col, row));
                }
            }
        }


        /**
         * Add the values of another matrix to this matrix.
         * @param adder The matrix to be added
         * @return New matrix after operation.
         */
        public Matrix add(Matrix adder)
        {
            Matrix temp = deepCopy();
            temp.addInplace(adder);
            return temp;
        }

        /**
         * Subtract the matrix values from current matrix
         * @param subtracter The matrix to be subtracted.
         */
        public void subtractInplace(Matrix subtracter)
        {
            if (columnNumber != subtracter.getColumnCount() || rowNumber != subtracter.getRowCount())
            {
                throw new ArgumentException("The dimension of input matrices does not match!");
            }

            for (int row = 0; row < rowNumber; row++)
            {
                for (int col = 0; col < columnNumber; col++)
                {
                    setElement(col, row, this.getElement(col, row) - subtracter.getElement(col, row));
                }
            }
        }

        /**
         * Subtract the matrix values from current matrix
         * @param subtracter The matrix to be subtracted.
         * @return The given matrix subtracted by argument.
         */
        public Matrix subtract(Matrix subtracter)
        {
            Matrix temp = deepCopy();
            temp.subtractInplace(subtracter);
            return temp;
        }

        /**
         * Multiply all elements in the matrix by given constant.
         * @param constant The constant.
         */
        public void multiplyByConstantInplace(double constant)
        {
            for (int row = 0; row < rowNumber; row++)
            {
                for (int col = 0; col < columnNumber; col++)
                {
                    setElement(col, row, constant * getElement(col, row));
                }
            }
        }

        /**
         * Multiply all elements in the matrix by given constant.
         * @param constant The constant.
         * @return The new matrix with values multiplied by constant.
         */
        public Matrix multiplyByConstant(double constant)
        {
            Matrix temp = deepCopy();
            temp.multiplyByConstantInplace(constant);
            return temp;
        }


        /**
         * Return the matrix representation in MATLAB(R) compatible format.
         * @return The string matrix representation in MATLAB(R) compatible format
         */
        public string toString()
        {
            String temp = "[";
            for (int row = 0; row < rowNumber; row++)
            {
                temp += "[";
                for (int col = 0; col < columnNumber; col++)
                {
                    temp += (getElement(col, row)).ToString();
                    if (col != columnNumber - 1)
                    {
                        temp += ",\t ";
                    }
                }
                temp += "]";
                if (row != (rowNumber - 1))
                {
                    temp += ",\n";
                }
            }
            temp += "]";
            return temp;
        }

        /**
         * C# Version
         */ 
        public override string ToString()
        {
            return toString();
        }


        /**
         * Convert the matrix to the WolframAlpha(R) compatible format
         * @param decimalPlacePrecision The rounding of decimal place to given position
         * @return WolframAlpha string representation of the matrix.
         */
        public string toWolframAlphaString(int decimalPlacePrecision)
        {
            string temp = "{";
            for (int row = 0; row < rowNumber; row++)
            {
                temp += "{";
                for (int col = 0; col < columnNumber; col++)
                {
                    temp += doubleToString(getElement(col, row), decimalPlacePrecision);
                    if (col != columnNumber - 1)
                    {
                        temp += ", ";
                    }
                }
                temp += "}";
                if (row != (rowNumber - 1))
                {
                    temp += ", ";
                }
            }
            temp += "}";
            return temp;
        }

        /**
         * Convert the matrix to the WolframAlpha compatible format
         * @return WolframAlpha string representation of the matrix.
         */
        public string toWolframAlphaString()
        {
            return toWolframAlphaString(4);
        }

        /**
         * Create a deep copy of the element (copy each element to the newly created matrix).
         * @return Newly created copy of the matrix.
         */
        public Matrix deepCopy()
        {
            Matrix temp = generateZeroMatrix(columnNumber, rowNumber);
            for (int col = 0; col < columnNumber; col++)
            {
                for (int row = 0; row < rowNumber; row++)
                {
                    temp.setElement(col, row, getElement(col, row));
                }
            }
            return temp;
        }

        /**
         * Get the dimension of the matrix
         * @return Dimension of the matrix in logic [rows, columns]
         */
        public int[] getDimension()
        {
            return new int[] { getRowCount(), getColumnCount() };
        }

        /**
         * Return the shape of the matrix in form (row count, column count)
         * @return the shape of the matrix in form (row count, column count)
         */
        public int[] shape()
        {
            return getDimension();
        }

        /**
         * Get the number of rows in the matrix
         * @return Number of rows in the matrix
         */
        public int getRowCount()
        {
            return rowNumber;
        }

        /**
         * Get the number of columns in the matrix.
         * @return The number of columns in the matrix.
         */
        public int getColumnCount()
        {
            return columnNumber;
        }

        /**
         * Compute 2D double array representation of the matrix
         * @return The 2D representation of the matrix as double[,] value (first index is row, second is column)
         */
        public double[,] get2DArrayRepresentation()
        {
            double[,] values = new double[rowNumber,columnNumber];
            for (int row = 0; row < rowNumber; row++)
            {
                for (int col = 0; col < columnNumber; col++)
                {
                    values[row,col] = getElement(col, row);
                }
            }
            return values;
        }

        /**
         * Find the upper triangular form of the matrix (usable not only for determinant computation).
         */
        public void getUpperTriangularInplace()
        {
            // Follow the logic of the reduced Gaussian Elimination with swapping

            for (int i = 0; i < rowNumber; i++)
            {
                int[] rowsWeights = getRowsWeight();
                int[] indicesOfRows = new int[rowNumber];
                for (int idx = 0; idx < rowNumber; idx++)
                {
                    indicesOfRows[idx] = idx;
                }
                sortWithMirror(rowsWeights, indicesOfRows);
                swapRows(i, indicesOfRows[i]);
                /*for (int kk = 0; kk < rowsWeights.length; kk++) {
                    System.out.print(rowsWeights[kk] + ", ");
                }
                System.out.println(" ");*/
                int row = i;
                int rowWeight = rowsWeights[i];

                if (rowWeight >= columnNumber)
                {
                    break;
                }
                double valueOnPivot = getElement(rowWeight, row);
                double inverse = 1.0 / valueOnPivot;
                //multiplyRowByValueInplace(row, inverse);
                //setElement(rowWeight, row, 1.0);

                // Subtract row value from each other row
                for (int j = rowWeight + 1; j < rowNumber; j++)
                {
                    double mulFactor = -1.0 * getElement(rowWeight, j) * inverse;
                    addMultiplicationOfRowToAnotherRowInplace(row, j, mulFactor);
                    // For reducing of rounding error set directly to zero:
                    if (j != row)
                    {
                        setElement(rowWeight, j, 0.0);
                    }
                    else
                    {
                        setElement(rowWeight, j, 1.0);
                    }
                }
            }

        }

        /**
         * Find the upper triangular form of the matrix (usable not only for determinant computation)
         * @return The upper triangular form of the matrix
         */
        public Matrix getUpperTriangular()
        {
            Matrix temp = deepCopy();
            temp.getUpperTriangularInplace();
            return temp;
        }

        /**
         * Compute the trace of the matrix (sum of all elements on diagonal)
         * @return The trace of the matrix
         */
        public double trace()
        {
            double trace = 0.0;
            for (int i = 0; i < Math.Min(getColumnCount(), getRowCount()); i++)
            {
                trace += getElement(i, i);
            }
            return trace;
        }

        /**
         * Compute the determinant of the matrix
         * @return The value of the determinant.
         */
        public double determinant()
        {
            if (rowNumber != columnNumber)
            {
                return 0.0;
            }
            Matrix detTemp = getUpperTriangular();
            double determinant = 1.0;
            for (int i = 0; i < columnNumber; i++)
            {
                determinant *= detTemp.getElement(i, i);
            }
            return determinant;
        }

        /**
         * Compute the rank of the matrix (number of linearly independent rows)
         * @return the rank of current matrix
         */
        public int rank()
        {
            Matrix temp = pivotMatrix();

            int[] weights = temp.getRowsWeight();
            int[] temp1 = new int[weights.Length];
            sortWithMirror(weights, temp1);
            bool nonZeroExists = false;
            for (int i = 0; i < weights.Length; i++)
            {
                if (weights[i] >= columnNumber)
                {
                    return i;
                }
                if (weights[i] != 0)
                {
                    nonZeroExists = true;
                }
            }
            if (!nonZeroExists)
            {
                return 0;
            }
            if (rowNumber < columnNumber)
            {
                return weights[weights.Length - 1] + 1;
            }
            return columnNumber;
        }

        /**
         * Get the null space (kernel of matrix linear projection) of the matrix
         * @return Null space of the matrix
         */
        public Matrix nullSpace()
        {
            int matrixRank = rank();
            Matrix nullspace = generateZeroMatrix(columnNumber - matrixRank, columnNumber);

            Matrix extendedMatrix = generateZeroMatrix(columnNumber, columnNumber);
            Matrix reduced = deepCopy();
            reduced.pivotMatrixInplace();
            int[] weights = reduced.getRowsWeight();
            int[] indices = new int[weights.Length];
            for (int i = 0; i < rowNumber; i++)
            {
                indices[i] = i;
            }

            sortWithMirror(weights, indices);
            for (int i = 0; i < weights.Length; i++)
            {
                extendedMatrix.setRow(weights[i], reduced.chooseRow(indices[i]));
            }
            extendedMatrix.multiplyByConstantInplace(-1.0);
            int nullIdx = 0;
            for (int diag = 0; diag < extendedMatrix.rowNumber; diag++)
            {
                if (deltaDoubleCompare(extendedMatrix.getElement(diag, diag), 0, DELTA_FOR_COMPARE))
                {
                    nullspace.setColumn(nullIdx, extendedMatrix.chooseColumn(diag));
                    nullspace.setElement(nullIdx, diag, 1);
                    nullIdx++;
                }
            }

            return nullspace;
        }

        /**
         * Add given constant to each element in matrix.
         * @param scalar Given constant
         */
        public void addScalarToEachElementInplace(double scalar)
        {
            for (int i = 0; i < matrix.Length; ++i)
            {
                matrix[i] = matrix[i] + scalar;
            }
        }

        /**
         * Add given constant to each element in matrix.
         * @param scalar Given constant
         * @return New matrix after operation.
         */
        public Matrix addScalarToEachElement(double scalar)
        {
            Matrix temp = deepCopy();
            temp.addScalarToEachElementInplace(scalar);
            return temp;
        }

        /**
         * Subtract given constant to each element in matrix.
         * @param scalar Given constant
         */
        public void subtractScalarToEachElementInplace(double scalar)
        {
            for (int i = 0; i < matrix.Length; ++i)
            {
                matrix[i] = matrix[i] - scalar;
            }
        }

        /**
         * Subtract given constant to each element in matrix.
         * @param scalar Given constant
         * @return New matrix after operation.
         */
        public Matrix subtractScalarToEachElement(double scalar)
        {
            Matrix temp = deepCopy();
            temp.subtractScalarToEachElementInplace(scalar);
            return temp;
        }

        /**
         * Multiply each element in matrix by given scalar.
         * @param scalar Given constant
         */
        public void multiplyScalarToEachElementInplace(double scalar)
        {
            multiplyByConstantInplace(scalar);
        }

        /**
         * Multiply each element in matrix by given scalar.
         * @param scalar Given constant
         * @return New matrix after operation.
         */
        public Matrix multiplyScalarToEachElement(double scalar)
        {
            return multiplyByConstant(scalar);
        }

        /**
         * Divide each element in matrix by given scalar.
         * @param scalar Given constant
         */
        public void divideScalarToEachElementInplace(double scalar) 
        {
            if (deltaDoubleCompare(scalar, 0.0, DELTA_FOR_COMPARE)) {
                throw new ArgumentException("Given scalar constant must not be equal to zero!");
            }
            for (int i = 0; i<matrix.Length; ++i) {
                matrix[i] = matrix[i] * scalar;
            }
        }

        /**
         * Divide each element in matrix by given scalar.
         * @param scalar Given constant
         * @return New matrix after operation.
         */
        public Matrix divideScalarToEachElement(double scalar)
        {
            Matrix temp = deepCopy();
            temp.multiplyScalarToEachElementInplace(scalar);
            return temp;
        }

        /**
         * Compute the Hadamard product (also known as the Schur product or the entrywise product) - multiply each element in matrix by corresponding element in given matrix.
         * @param m Matrix of the same size use in Hadamard product computation.
         */
        public void hadamardProductInpace(Matrix m)
        {
            if (m.getRowCount() != getRowCount() || m.getColumnCount() != getColumnCount()) {
                throw new ArgumentException("The shape (number of rows and columns) of both matrices has to be same.");
            }
            for (int i = 0; i < matrix.Length; i++) {
                this.matrix[i] = this.matrix[i] * m.matrix[i];
            }
        }


        /**
         * Compute the Hadamard product (also known as the Schur product or the entrywise product) - multiply each element in matrix by corresponding element in given matrix.
         * @param m Matrix of the same size use in Hadamard product computation.
         * @return Hadamard product of two matrices.
         */
        public Matrix hadamardProduct(Matrix m) 
        {
            Matrix temp = deepCopy();
            temp.hadamardProductInpace(m);
            return temp;
        }

        /**
         * Divide each element in matrix by corresponding element in input matrix.
         * @param m Input matrix
         * @throws ArithmeticException if division by zero occurs
         */
        public void entrywiseDivisionInplace(Matrix m) {
            if (m.getRowCount() != getRowCount() || m.getColumnCount() != getColumnCount()) {
                throw new ArgumentException("The shape (number of rows and columns) of both matrices has to be same.");
            }
            for (int i = 0; i<matrix.Length; i++) {
                this.matrix[i] = this.matrix[i] / m.matrix[i];
            }
        }
    
        /**
         * Divide each element in matrix by corresponding element in input matrix.
         * @param m Input matrix
         * @throws ArithmeticException if division by zero occurs
         * @return The new matrix after computation.
         */
        public Matrix entrywiseDivision(Matrix m) {
            Matrix temp = deepCopy();
            temp.entrywiseDivisionInplace(m);
            return temp;
        }

        /**
     * Compute the power of each element to given scalar
     * @param scalar Given exponent
     */
        public void powerToScalarInplace(double scalar)
        {
            for (int i = 0; i < matrix.Length; i++)
            {
                this.matrix[i] = Math.Pow(this.matrix[i], scalar);
            }
        }

        /**
         * Compute the power of each element to given scalar
         * @param scalar Given exponent
         * @return New matrix after computation.
         */
        public Matrix powerToScalar(double scalar)
        {
            Matrix temp = deepCopy();
            temp.powerToScalarInplace(scalar);
            return temp;
        }

        /**
         * Return the mean of selected columns as a double array.
         * @param columnIndices The list of columns for computation of mean.
         * @return The array of mean values of elements in chosen columns.
         */
        public double[] meanOnColumnSpace(int[] columnIndices)
        {
            double [] average = new double[columnIndices.Length];
            for (int i = 0; i < columnIndices.Length; i++)
            {
                average[i] = meanOnColumnSpace(columnIndices[i]);
            }
            return average;
        }

        /**
         * Return the mean of selected columns as a double array.
         * @param columnIndex The column for computation of mean.
         * @return The mean values of elements in chosen column.
         */
        public double meanOnColumnSpace(int columnIndex)
        {
            double sum = 0.0;
            for (int row = 0; row < getRowCount(); ++row)
            {
                sum += getElement(columnIndex, row);
            }
            return sum / ((double)getRowCount());
        }

        /**
         * Return the mean of all columns as a double array.
         * @return The array of mean values of elements in all columns.
         */
        public double[] meanOnColumnSpace()
        {
            double [] average = new double[getColumnCount()];
            for (int i = 0; i < average.Length; i++)
            {
                average[i] = meanOnColumnSpace(i);
            }
            return average;
        }

        /**
         * Return the mean of selected rows as a double array.
         * @param rowIndices The list of rows for computation of mean.
         * @return The array of mean values of elements in chosen rows.
         */
        public double[] meanOnRowSpace(int[] rowIndices)
        {
            double [] average = new double[rowIndices.Length];
            for (int i = 0; i < rowIndices.Length; i++)
            {
                average[i] = meanOnRowSpace(rowIndices[i]);
            }
            return average;
        }

        /**
         * Return the mean of selected rows as a double array.
         * @param rowIndex The row for computation of mean.
         * @return The mean values of elements in chosen row.
         */
        public double meanOnRowSpace(int rowIndex)
        {
            double sum = 0.0;
            for (int column = 0; column < getRowCount(); ++column)
            {
                sum += getElement(column, rowIndex);
            }
            return sum / ((double)getColumnCount());
        }

        /**
         * Return the mean of all rows as a double array.
         * @return The array of mean values of elements in all rows.
         */
        public double[] meanOnRowSpace()
        {
            double [] average= new double[getRowCount()];
            for (int i = 0; i < average.Length; i++)
            {
                average[i] = meanOnRowSpace(i);
            }
            return average;
        }

        /**
         * Parse the input string to the matrix values.
         * @param input Input string.
         * @param opening The opening character (typically [ or {).
         * @param closing The closing character (typically ] or }).
         * @param separator The separator between two values (typically ,).
         * @return The parsed matrix.
         * @throws ArgumentException If the string is impaired.
         */
        private static Matrix parsingString(string input, char opening, char closing, char separator)
        {
            int rowsCount = 0;
            int columnCounts = 0;
            int elementCounts = 0;
            bool openedLevel0 = false;
            bool openedLevel1 = false;

            for (int i = 0; i < input.Length; i++)
            {
                if (input[i] == opening)
                {
                    if (!openedLevel0)
                    {
                        openedLevel0 = true;
                    }
                    else if (!openedLevel1)
                    {
                        rowsCount++;
                        openedLevel1 = true;
                    }
                    else
                    {
                        throw new ArgumentException("Wrong input string!");
                    }
                }
                if (openedLevel1 && input[i] == separator)
                {
                    elementCounts++;
                }
                if (input[i] == closing)
                {
                    if (openedLevel1)
                    {
                        openedLevel1 = false;
                        elementCounts++;
                    }
                    else if (openedLevel0)
                    {
                        openedLevel0 = false;
                    }
                    else
                    {
                        throw new ArgumentException("Wrong input string!");
                    }
                }
            }
            columnCounts = elementCounts / rowsCount;

            Matrix result = generateZeroMatrix(columnCounts, rowsCount);
            int elementPosition = 0;
            string value = "";
            for (int i = 0; i < input.Length; i++)
            {
                if (input[i] == opening)
                {
                    if (!openedLevel0)
                    {
                        openedLevel0 = true;
                    }
                    else if (!openedLevel1)
                    {
                        openedLevel1 = true;
                    }
                }
                if (openedLevel1 && input[i] == separator)
                {
                    // PARSE
                    double valueInt = Double.Parse(value);
                    result.setElement(elementPosition % columnCounts, elementPosition / columnCounts, valueInt);
                    elementPosition++;
                    value = "";
                }
                else if (openedLevel1 && (Char.IsDigit(input[i])) || input[i] == '.')
                {
                    value += input[i];
                }
                if (input[i] == closing)
                {
                    if (openedLevel1)
                    {
                        openedLevel1 = false;
                        // PARSE
                        double valueInt = Double.Parse(value);
                        result.setElement(elementPosition % columnCounts, elementPosition / columnCounts, valueInt);
                        elementPosition++;
                        value = "";
                    }
                    else if (openedLevel0)
                    {
                        openedLevel0 = false;
                    }
                }
            }

            return result;
        }

        /**
         * Parse the Matlab (Python) compatible string
         * @param input The input string
         * @return Parsed Matrix
         * @throws ArgumentException If there is an error in the string.
         */
            public static Matrix parseMatlabString(string input)
            {
                return parsingString(input, '[', ']', ',');
            }

        /**
         * Parse the WolframAlpha compatible string
         * @param input The input string
         * @return Parsed Matrix
         * @throws ArgumentException If there is an error in the string.
         */
        public static Matrix parseWolframString(string input)
        {
            return parsingString(input, '{', '}', ',');
        }

    }
}
