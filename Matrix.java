/*
Copyright (c) 2019 David Salac.

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

package eu.davidsalac.matrix;

import java.text.DecimalFormat;

/**
 * The class for dealing with fundamental matrix operation for matrices in R(n,m)
 * @author David Salac
 */
public class Matrix {
    private double [] matrix;
    private int columnNumber;
    private int rowNumber;
    private final double DELTA_FOR_COMPARE;
    private int rank;
    
    /**
     * Sorting an array and swapping the mirror array in a same manner as an source array.
     * @param arr the input array to be sorted.
     * @param mirror the mirror array which elements are swapped same as an input array.
     */
    private static void sortWithMirror(int arr[], int mirror[]) 
    { 
        int n = arr.length; 
  
        // Build heap (rearrange array) 
        for (int i = n / 2 - 1; i >= 0; i--) 
            heapifyWithMirror(arr, mirror, n, i); 
  
        // One by one extract an element from heap 
        for (int i=n-1; i>=0; i--) 
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
    private static void heapifyWithMirror(int arr[], int mirror[], int n, int i) 
    { 
        int largest = i; // Initialize largest as root 
        int l = 2*i + 1; // left = 2*i + 1 
        int r = 2*i + 2; // right = 2*i + 2 
  
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
     * Create an instance with defined delta for comparison of doubles
     */
    private Matrix () {
        // Compare doubles in given delta;
        DELTA_FOR_COMPARE = 0.0001;
    }
    
    /**
     * Create an instance with zeros everywhere
     * @param columnNumber The required number of columns in the matrix
     * @param rowNumber  The required number of rows in the matrix
     */
    private Matrix (int columnNumber, int rowNumber) {
        this.matrix = new double[columnNumber * rowNumber];
        this.columnNumber = columnNumber;
        this.rowNumber = rowNumber;
        // For comparison of doubles in given delta
        DELTA_FOR_COMPARE = 0.0001;
    }
    
    /**
     * Get the element in matrices on given position
     * @param column the column index
     * @param row the row index
     * @return element at given position
     */
    public double getElement(int column, int row) {
        if (row < 0 || row > rowNumber) {
            throw new IndexOutOfBoundsException("Index row out of bound!");
        }
        if (column < 0 || column > columnNumber) {
            throw new IndexOutOfBoundsException("Index column out of bound!");
        }
        return matrix[row * columnNumber + column];
    }
    
    /**
     * Set the element value
     * @param column the column position of the element
     * @param row the row position of the element
     * @param value the new value for the element
     */
    public void setElement(int column, int row, double value) {
        if (row < 0 || row > rowNumber) {
            throw new IndexOutOfBoundsException("Index row out of bound!");
        }
        if (column < 0 || column > columnNumber) {
            throw new IndexOutOfBoundsException("Index column out of bound!");
        }
        matrix[row * columnNumber + column] = value;
    }
    
    /**
     * Set the row for the new values given by array of delta values
     * @param row The row where an array is inserted
     * @param values The values for the row
     * @throws IndexOutOfBoundsException in the case that row index is out of bound
     * @throws IllegalArgumentException  in the case that dimension of the array does not fit to matrix
     */
    public void setRow(int row, double [] values) throws IndexOutOfBoundsException, IllegalArgumentException {
        if (row < 0 || row > rowNumber) {
            throw new IndexOutOfBoundsException("Index row out of bound!");
        }
        if (values.length != columnNumber) {
            throw new IllegalArgumentException("The number of columns in the matrix has to match the length of the input vector!");
        }
        for (int col = 0; col < columnNumber; col++) {
            setElement(col, row, values[col]);
        }
    }
    
    /**
     * Set the row using row vector represented as Matrix class instance
     * @param row the row where the new values are inserted.
     * @param values the vector of values.
     * @throws IndexOutOfBoundsException in the case that row index is out of bound
     * @throws IllegalArgumentException  in the case that dimension of the array does not fit to matrix
     */
    public void setRow(int row, Matrix values) throws IndexOutOfBoundsException, IllegalArgumentException {
        if (row < 0 || row > rowNumber) {
            throw new IndexOutOfBoundsException("Index row out of bound!");
        }
        if (values.getRowCount() != 1 || values.getColumnCount() != columnNumber) {
            throw new IllegalArgumentException("The number of columns in the matrix has to match the length of the input vector or the number of rows in input vector has to be just 1!");
        }
        for (int col = 0; col < columnNumber; col++) {
            setElement(col, row, values.getElement(col, 0));
        }
    }
    
    /**
     * Set the column of values to the matrix
     * @param column The position of the column which will be rewritten.
     * @param values The new values that are inserted to the matrix
     * @throws IndexOutOfBoundsException in the case that column index is out of bound
     * @throws IllegalArgumentException in the case that dimension of the array does not fit to matrix
     */
    public void setColumn(int column, double [] values) throws IndexOutOfBoundsException, IllegalArgumentException {
        if (column < 0 || column > columnNumber) {
            throw new IndexOutOfBoundsException("Index column out of bound!");
        }
        if (values.length != rowNumber) {
            throw new IllegalArgumentException("The number of rows in the matrix has to match the length of the input vector!");
        }
        for (int row = 0; row < rowNumber; row++) {
            setElement(column, row, values[row]);
        }
    }
    
    /**
     * Set the values in given column to the new input values represented by Matrix vector
     * @param column The position of the column which will be rewritten
     * @param values The new values
     * @throws IndexOutOfBoundsException in the case that column index is out of bound
     * @throws IllegalArgumentException in the case that dimension of the array does not fit to matrix
     */
    public void setColumn(int column, Matrix values) throws IndexOutOfBoundsException, IllegalArgumentException {
        if (column < 0 || column > columnNumber) {
            throw new IndexOutOfBoundsException("Index column out of bound!");
        }
        if (values.getColumnCount() != 1 || values.getRowCount()!= rowNumber) {
            throw new IllegalArgumentException("The number of rows in the matrix has to match the length of the input vector or the number of columns in input vector has to be just 1!");
        }
        for (int row = 0; row < rowNumber; row++) {
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
    public static Matrix generateRandomMatrix(int columnNumber, int rowNumber, double valuesFrom, double valuesTo) {
        Matrix tempMatrix = generateZeroMatrix(columnNumber, rowNumber);
        for (int row = 0; row < rowNumber; row++) {
            for (int col = 0; col < columnNumber; col++) {
                tempMatrix.setElement(col, row, Math.random() * (valuesTo - valuesFrom)+valuesFrom);
            }
        }
        return tempMatrix;
    }
    
    /**
     * Generate a new matrix with integer values.
     * @param columnNumber Number of columns in the new matrix
     * @param rowNumber Number of rows in the new matrix
     * @param valuesFrom The lower bound for values
     * @param valuesTo The upper bound for values
     * @return The new matrix with given values in defined bound.
     */
    public static Matrix generateRandomIntegerMatrix(int columnNumber, int rowNumber, int valuesFrom, int valuesTo) {
        Matrix tempMatrix = generateZeroMatrix(columnNumber, rowNumber);
        for (int row = 0; row < rowNumber; row++) {
            for (int col = 0; col < columnNumber; col++) {
                tempMatrix.setElement(col, row, (int)(Math.random() * (valuesTo - valuesFrom)+valuesFrom));
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
    public static Matrix generateZeroMatrix(int columnNumber, int rowNumber) {
        if (columnNumber < 1 || rowNumber < 1) {
            throw new IllegalArgumentException("The number of columns and rows has to be greater than 1!");
        }
        return new Matrix(columnNumber, rowNumber);
    }
    
    /**
     * Generate the matrix that has element equal to 1 on diagonal and 0 everywhere else
     * @param columnNumber The number of columns in a matrix.
     * @param rowNumber The number of rows in a matrix.
     * @return New matrix with given size.
     */
    public static Matrix generateOneDiagonalMatrix(int columnNumber, int rowNumber) {
        Matrix tempMatrix = generateZeroMatrix(columnNumber, rowNumber);
        for (int i = 0; i < Math.min(rowNumber, columnNumber); ++i) {
            tempMatrix.setElement(i, i, 1.0);
        }
        return tempMatrix;
    }
    
    /**
     * Swap the rows in the matrix without creating a new one
     * @param rowA The position of the first row.
     * @param rowB The position of another row.
     */
    public void swapRowsInplace(int rowA, int rowB) {
        if (rowA < 0 || rowA >= rowNumber) {
            throw new IndexOutOfBoundsException("Index rowA out of bound!");
        }
        if (rowB < 0 || rowB >= rowNumber) {
            throw new IndexOutOfBoundsException("Index rowB out of bound!");
        }
        for (int col = 0; col < columnNumber; col ++) {
            double temp = getElement(col, rowA);
            setElement(col, rowA, getElement(col, rowB));
            setElement(col, rowB, temp);
        }
    }
    
    /**
     * Swap the rows in the matrix together with creating a new one
     * @param rowA The position of the first row.
     * @param rowB The position of another row.
     * @return The new matrix with swapped rows.
     */
    public Matrix swapRows(int rowA, int rowB) {
        Matrix temporaryMatrix = this.deepCopy();
        temporaryMatrix.swapRowsInplace(rowA, rowB);
        return temporaryMatrix;
    }
    
    /**
     * Swap the columns in the matrix without creating a new one
     * @param columnA The position of the first column.
     * @param columnB The position of another column.
     */
    public void swapColumnsInplace(int columnA, int columnB) {
        if (columnA < 0 || columnA >= columnNumber) {
            throw new IndexOutOfBoundsException("Index columnA out of bound!");
        }
        if (columnB < 0 || columnB >= columnNumber) {
            throw new IndexOutOfBoundsException("Index columnB out of bound!");
        }
        for (int row = 0; row < rowNumber; row ++) {
            double temp = getElement(columnA, row);
            setElement(columnA, row, getElement(columnB, row));
            setElement(columnB, row, temp);
        }
    }
    
    /**
     * Swap the columns in the matrix together with creating a new one
     * @param columnA The position of the first column.
     * @param columnB The position of another column.
     * @return The new matrix with swapped column.
     */
    public Matrix swapColumns(int columnA, int columnB) {
        Matrix temporaryMatrix = this.deepCopy();
        temporaryMatrix.swapColumnsInplace(columnA, columnB);
        return temporaryMatrix;
    }
    
    /**
     * Compare two double values with respect to given delta values
     * @param A First double value
     * @param B Another double value
     * @param delta The delta parameters for comparison
     * @return abs(A-B) < delta
     */
    public static boolean deltaDoubleCompare(double A, double B, double delta) {
        return Math.abs(A - B) < delta;
    }
    
    /**
     * Get the weight of the row (the column position of the first non-zero element in row)
     * @param row Probed row
     * @return The position of the first non-zero element in the row or the number of rows if none such element exists.
     */
    public int getRowWeight(int row) {
        // Find the position (from left) of the first non-zero element in the matrix, return 0 if it does not exists.
        for (int col = 0; col < columnNumber; col ++) {
            if (!deltaDoubleCompare(0, getElement(col, row), DELTA_FOR_COMPARE)) {
                return col;
            }
        }
        return columnNumber;
    }
    
    /**
     * Get the weight of all rows in the matrix.
     * @return Weights of all rows in the matrix (weight represents the position of the first non-zero element in the row or the number of rows if none such element exists.)
     */
    public int [] getRowsWeight() {
        int [] weights = new int[rowNumber];
        for(int row = 0; row < rowNumber; row++) {
            weights[row] = getRowWeight(row);
        }
        return weights;
    }
    
    /**
     * Choose the column of the matrix at given position
     * @param column Index of desired column
     * @return The column on desired position represented as an Matrix type
     */
    public Matrix chooseColumn(int column) {
        Matrix tempCol = generateZeroMatrix(1, rowNumber);
        for (int row = 0; row< rowNumber; row++) {
            tempCol.setElement(0, row, getElement(column, row));
        }
        return tempCol;
    }
    
    /**
     * Select specified columns
     * @param columns An array with indices of required columns.
     * @return The columns on desired positions represented as an Matrix type
     */
    public Matrix chooseColumns(int [] columns) {
        Matrix temp = generateZeroMatrix(columns.length, rowNumber);
        for (int col = 0; col < columns.length; col++) {
            temp.setColumn(col, this.chooseColumn(columns[col]));
        }
        return temp;
    }
    
    /**
     * Choose the row of the matrix at given position
     * @param row Index of desired row
     * @return The row on desired position represented as an Matrix type
     */
    public Matrix chooseRow(int row) {
        Matrix tempRow = generateZeroMatrix(columnNumber, 1);
        for (int col = 0; col< columnNumber; col++) {
            tempRow.setElement(col, 0, getElement(col, row));
        }
        return tempRow;
    }
    
    /**
     * Select specified rows
     * @param rows An array with indices of required rows.
     * @return The rows on desired positions represented as an Matrix type
     */
    public Matrix chooseRows(int [] rows) {
        Matrix temp = generateZeroMatrix(columnNumber, rows.length);
        for (int row = 0; row < rows.length; row++) {
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
    public void addMultiplicationOfRowToAnotherRowInplace(int rowSource, int rowTarget, double multiplier) {
        if (rowSource < 0 || rowSource >= rowNumber) {
            throw new IndexOutOfBoundsException("Index rowSource out of bound!");
        }
        if (rowTarget < 0 || rowTarget >= rowNumber) {
            throw new IndexOutOfBoundsException("Index rowTarget out of bound!");
        }
        if(rowSource == rowTarget) {
            return;
        }
        for (int col = 0; col < columnNumber; col++) {
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
    public Matrix addMultiplicationOfRowToAnotherRow(int rowSource, int rowTarget, double multiplier) {
        Matrix temp = deepCopy();
        temp.addMultiplicationOfRowToAnotherRowInplace(rowSource, rowTarget, multiplier);
        return temp;
    }
    
    /**
     * Compute the pivotal representation of the matrix without swapping rows
     * @param mirror The matrix where exactly the same operations are done (or null)
     */
    public void pivotMatrixInplace(Matrix mirror) {
        for (int i = 0; i < rowNumber; i++) {
            int [] rowsWeights = getRowsWeight();
            int [] indicesOfRows = new int[rowNumber];
            for (int idx = 0; idx < rowNumber; idx++) {
                indicesOfRows[idx] = idx;
            }
            sortWithMirror(rowsWeights, indicesOfRows);
            /*for (int kk = 0; kk < rowsWeights.length; kk++) {
                System.out.print(rowsWeights[kk] + ", ");
            }
            System.out.println(" ");*/
            int row = indicesOfRows[i];
            int rowWeight = rowsWeights[i];
            
            if (rowWeight >= columnNumber) {
                break;
            }
            double valueOnPivot = getElement(rowWeight, row);
            double inverse = 1.0/valueOnPivot;
            multiplyRowByValueInplace(row, inverse);
            if (mirror != null) {
                mirror.multiplyRowByValueInplace(row, inverse);
            }
            //setElement(rowWeight, row, 1.0);
            
            // Subtract row value from each other row
            for (int j = 0; j < rowNumber; j++) {
                double mulFactor = -1.0 * getElement(rowWeight, j);
                addMultiplicationOfRowToAnotherRowInplace(row, j, mulFactor);
                // For reducing of rounding error set directly to zero:
                if (j != row) {
                    setElement(rowWeight, j, 0.0);
                } else {
                    setElement(rowWeight, j, 1.0);
                }
                if (mirror != null) {
                    mirror.addMultiplicationOfRowToAnotherRowInplace(row, j, mulFactor);
                }
            }
        }
        
    }
    
    /**
     * Compute the pivotal representation of the matrix without swapping rows
     */
    public void pivotMatrixInplace() {
        pivotMatrixInplace(null);
    }
    
    /**
     * Compute the pivotal representation of the matrix without swapping rows and return modified matrix as a new object.
     * @param mirror The matrix where exactly the same operations are done (or null)
     * @return Matrix in pivotal form.
     */
    public Matrix pivotMatrix(Matrix mirror) {
        Matrix temp = deepCopy();
        temp.pivotMatrixInplace(mirror);
        return temp;
    }
    
    /**
     * Pivot the matrix and return it as a new object
     * @return Matrix in pivotal form.
     */
    public Matrix pivotMatrix() {
        return pivotMatrix(null);
    }
    
    /**
     * Transposition of the matrix
     * @return Transposed matrix
     */
    public Matrix transpose() {
        Matrix temp = generateZeroMatrix(rowNumber, columnNumber);
        for (int row = 0; row < rowNumber; row++) {
            for (int col = 0; col < columnNumber; col++) {
                temp.setElement(row, col, getElement(col, row));
            }
        }
        return temp;
    }
    
    /**
     * Compute the inversion of the matrix if exists
     * @return Inverse matrix to current matrix
     */
    public Matrix getInverseMatrix() {
        if (rowNumber != columnNumber) {
            throw new IllegalArgumentException("The matrix must be rectangular!");
        }
        Matrix copy = deepCopy();
        Matrix inverse = generateOneDiagonalMatrix(columnNumber, columnNumber);
        copy.pivotMatrixInplace(inverse);
        int [] rowWeights = copy.getRowsWeight();
        int [] rowIndices = new int[columnNumber];
        for (int i = 0; i < rowWeights.length; i++) {
            rowIndices[i] = i;
        }
        sortWithMirror(rowWeights, rowIndices);
        Matrix inverseFind = generateZeroMatrix(columnNumber, rowNumber);
        for (int row = 0; row < rowNumber; row++) {
            /*System.out.println();*/
            inverseFind.setRow(row, inverse.chooseRow(rowIndices[row]));
        }
        return inverseFind;
    }
    
    /**
     * Compute the multiplication of the current matrix with inserted matrix
     * @param multiplier Another matrix
     * @return current matrix multiplied by given matrix in the argument.
     */
    public Matrix multiply(Matrix multiplier) {
        if (columnNumber != multiplier.rowNumber) {
            throw new IllegalArgumentException("The dimension of input matrices does not match!");
        }
        Matrix result = generateZeroMatrix(multiplier.columnNumber, rowNumber);
        
        for (int resRow = 0; resRow < result.rowNumber; resRow++) {
            for (int resCol = 0; resCol < result.columnNumber; resCol++) {
                double mulRes = 0.0;
                for (int i = 0; i < columnNumber; i++) {
                    mulRes += this.getElement(i, resRow) * multiplier.getElement(resCol, i);
                    if (deltaDoubleCompare(mulRes, 0, DELTA_FOR_COMPARE)) {
                        mulRes = 0;
                    }
                    if (deltaDoubleCompare(mulRes, 1, DELTA_FOR_COMPARE)) {
                        mulRes = 1;
                    }
                }
                result.setElement(resCol, resRow, mulRes);
            }
        }
        return result;
    }
    
    /**
     * Multiply all element in given row by value
     * @param row Index of the row
     * @param value The multiplication constant
     */
    public void multiplyRowByValueInplace(int row, double value) {
        if (row < 0 || row >= rowNumber) {
            throw new IndexOutOfBoundsException("Index row out of bound!");
        }
        for (int col = 0; col < columnNumber; col++) {
            setElement(col, row, value * getElement(col, row));
        }
    }
    /**
     * Multiply all element in given row by value and return result as a new matrix
     * @param row Index of the row
     * @param value The multiplication constant
     * @return New matrix after operation
     */
    public Matrix multiplyRowByValue(int row, double value) {
        Matrix temp = deepCopy();
        temp.multiplyRowByValueInplace(row, value);
        return temp;
    }
    
    /**
     * Add the values of another matrix to this matrix
     * @param adder The matrix to be added
     * @return Result after adding two matrices
     */
    public Matrix add(Matrix adder) {
        if (columnNumber != adder.getColumnCount() || rowNumber != adder.getRowCount()) {
            throw new IllegalArgumentException("The dimension of input matrices does not match!");
        }
        Matrix temp = this.deepCopy();
        for (int row = 0; row < rowNumber; row++) {
            for (int col = 0; col < columnNumber; col++) {
                temp.setElement(col, row, this.getElement(col, row) + adder.getElement(col, row));
            }
        }
        return temp;
    }
    
    /**
     * Subtract the matrix values from current matrix
     * @param subtracter The matrix to be subtracted.
     * @return The given matrix subtracted by argument.
     */
    public Matrix subtract(Matrix subtracter) {
        if (columnNumber != subtracter.getColumnCount() || rowNumber != subtracter.getRowCount()) {
            throw new IllegalArgumentException("The dimension of input matrices does not match!");
        }
        Matrix temp = this.deepCopy();
        for (int row = 0; row < rowNumber; row++) {
            for (int col = 0; col < columnNumber; col++) {
                temp.setElement(col, row, this.getElement(col, row) - subtracter.getElement(col, row));
            }
        }
        return temp;
    }
    
    /**
     * Multiply all elements in the matrix by given constant.
     * @param constant The constant.
     */
    public void multiplyByConstantInplace(double constant) {
        for (int row = 0; row < rowNumber; row++) {
            for (int col = 0; col < columnNumber; col++) {
                setElement(col, row, constant*getElement(col, row));
            }
        }
    }
    
    /**
     * Multiply all elements in the matrix by given constant.
     * @param constant The constant.
     * @return The new matrix with values multiplied by constant.
     */
    public Matrix multiplyByConstant(double constant) {
        Matrix temp = deepCopy();
        temp.multiplyByConstantInplace(constant);
        return temp;
    }
    
    /**
     * Return the matrix representation in MATLAB(R) compatible format.
     * @return The string matrix representation in MATLAB(R) compatible format
     */
    @Override
    public String toString() {
        String temp = "[";
        for (int row = 0; row < rowNumber; row ++) {
            temp += "[";
            for (int col = 0; col < columnNumber; col++) {
                temp += Double.toString(getElement(col, row));
                if (col != columnNumber - 1) {
                    temp += ",\t ";
                }
            }
            temp += "]";
            if (row != (rowNumber - 1)) {
                temp += ",\n";
            }
        }
        temp += "]";
        return temp;
    }
    
    /**
     * Convert the matrix to the WolframAlpha compatible format
     * @return WolframAlpha string representation of the matrix.
     */
    public String toWolframAlphaString() {
        String temp = "{";
        for (int row = 0; row < rowNumber; row ++) {
            temp += "{";
            for (int col = 0; col < columnNumber; col++) {
                DecimalFormat formatTo2DecDigit = new DecimalFormat(".####"); //for being able to format your output to 2 decimal places
                temp += (formatTo2DecDigit.format(getElement(col, row)));
                if (col != columnNumber - 1) {
                    temp += ", ";
                }
            }
            temp += "}";
            if (row != (rowNumber - 1)) {
                temp += ", ";
            }
        }
        temp += "}";
        return temp;
    }
    
    /**
     * Create a deep copy of the element (copy each element to the newly created matrix).
     * @return Newly created copy of the matrix.
     */
    public Matrix deepCopy() {
        Matrix temp = generateZeroMatrix(columnNumber, rowNumber);
        for (int col = 0; col < columnNumber; col++) {
            for (int row = 0; row < rowNumber; row ++) {
                temp.setElement(col, row, getElement(col, row));
            }
        }
        return temp;
    }
    
    /**
     * Get the dimension of the matrix
     * @return Dimension of the matrix in logic [rows, columns]
     */
    public int [] getDimension() {
        return new int [] { getRowCount(), getColumnCount()};
    }
    
    /**
     * Get the number of rows in the matrix
     * @return Number of rows in the matrix
     */
    public int getRowCount() {
        return rowNumber;
    }
    
    /**
     * Get the number of columns in the matrix.
     * @return The number of columns in the matrix.
     */
    public int getColumnCount() {
        return columnNumber;
    }
    
    /**
     * Compute 2D double array representation of the matrix
     * @return The 2D representation of the matrix as double[][] value (first index is row, second is column)
     */
    public double [][] get2DArrayRepresentation() {
        double [][] values = new double[rowNumber][columnNumber];
        for (int row = 0; row < rowNumber; row++) {
            for (int col = 0; col < columnNumber; col++) {
                values[row][col] = getElement(col, row);
            }
        }
        return values;
    }
    
    /**
     * Find the upper triangular form of the matrix (usable not only for determinant computation)
     */
    public void getUpperTriangularInplace() {
        // Follow the logic of the reduced Gaussian Elimination with swapping
        
        for (int i = 0; i < rowNumber; i++) {
            int [] rowsWeights = getRowsWeight();
            int [] indicesOfRows = new int[rowNumber];
            for (int idx = 0; idx < rowNumber; idx++) {
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
            
            if (rowWeight >= columnNumber) {
                break;
            }
            double valueOnPivot = getElement(rowWeight, row);
            double inverse = 1.0/valueOnPivot;
            //multiplyRowByValueInplace(row, inverse);
            //setElement(rowWeight, row, 1.0);
            
            // Subtract row value from each other row
            for (int j = rowWeight + 1; j < rowNumber; j++) {
                double mulFactor = -1.0 * getElement(rowWeight, j) * inverse;
                addMultiplicationOfRowToAnotherRowInplace(row, j, mulFactor);
                // For reducing of rounding error set directly to zero:
                if (j != row) {
                    setElement(rowWeight, j, 0.0);
                } else {
                    setElement(rowWeight, j, 1.0);
                }
            }
        }
        
    }
    
    /**
     * Find the upper triangular form of the matrix (usable not only for determinant computation)
     * @return The upper triangular form of the matrix
     */
    public Matrix getUpperTriangular() {
        Matrix temp = deepCopy();
        temp.getUpperTriangularInplace();
        return temp;
    }
    
    /**
     * Compute the trace of the matrix (sum of all elements on diagonal)
     * @return The trace of the matrix
     */
    public double trace() {
        if (rowNumber != columnNumber) {
            throw new IllegalArgumentException("The matrix must be rectangular!");
        }
        double trace = 0;
        for (int i = 0; i < rowNumber; i++) {
            trace += getElement(i, i);
        }
        return trace;
    }
    
    /**
     * Compute the determinant of the matrix
     * @return The value of the determinant.
     */
    public double determinant() {
        if (rowNumber != columnNumber) {
            return 0;
        }
        Matrix detTemp = getUpperTriangular();
        double determinant = 1;
        for (int i =0; i < columnNumber; i ++) {
            determinant *= detTemp.getElement(i, i);
        }
        return determinant;
    }
    
    /**
     * Compute the rank of the matrix (number of linearly independent rows)
     * @return the rank of current matrix
     */
    public int rank() {
        Matrix temp = pivotMatrix();
        
        int [] weights = temp.getRowsWeight();
        int [] temp1 = new int[weights.length];
        sortWithMirror(weights, temp1);
        boolean nonZeroExists = false;
        for (int i = 0; i < weights.length; i ++) {
            if (weights[i] >= columnNumber) {
                return i;
            }
            if (weights[i] != 0) {
                nonZeroExists = true;
            }
        }
        if (!nonZeroExists) {
            return 0;
        }
        return columnNumber;
    }
    
    public Matrix parseMatlabString(String input) {
        throw new java.lang.UnsupportedOperationException("Not supported yet.");
    }
    public Matrix parseWolframString(String input) {
        throw new java.lang.UnsupportedOperationException("Not supported yet.");
    }
    
    public Matrix nullSpace() {
        throw new java.lang.UnsupportedOperationException("Not supported yet.");
    }
    public double [] eigenvalues() {
        throw new java.lang.UnsupportedOperationException("Not supported yet.");
    }
    public Matrix [] eigenvectors() {
        throw new java.lang.UnsupportedOperationException("Not supported yet.");
    }
    public void rotateInplace(double angleInRadians) {
        throw new java.lang.UnsupportedOperationException("Not supported yet.");
    }
}
