using System.Collections.Generic;
using System.Linq;
using System;

namespace circuit_sim
{
    public class ComplexMatrix
    {
        protected int rows;
        protected int colums;
        List<ComplexNumber> values;
        public int Columns => colums;
        public int Rows => rows;
        public virtual IEnumerable<ComplexNumber> Values => values;

        protected ComplexMatrix() { }
        public ComplexMatrix(int row) : this(row, row)
        {
        }

        public int CellOffsetInValue(int rowIndex, int columnIndex)
        {
            return rowIndex * colums + columnIndex;
        }
        public ComplexMatrix(int row, int colum)
        {
            this.rows = row;
            this.colums = colum;
            values = new List<ComplexNumber>(row * colum);
            int size = row * colum;
            for (int i = 0; i < size; i++)
            {
                values.Add(0);
            }
        }
        public ComplexMatrix(int row, int colum, IEnumerable<ComplexNumber> values) : this(row, colum)
        {
            int i = 0;
            int maxSize = row * colum;
            foreach (var v in values)
            {
                if (i < maxSize)
                {
                    this.values[i] = v;
                }
                i++;
            }
        }

        public ComplexMatrix Clone()
        {
            return new ComplexMatrix(rows, colums, values);
        }

        public ComplexMatrix CroComplexMatrix(int offsetR, int offsetC)
        {
            var tmp = new ComplexMatrix(this.rows - offsetR, this.colums - offsetC);
            foreach (var (r, c) in tmp.Indexes())
            {
                tmp[r, c] = this[r + offsetR, c + offsetC];
            }
            return tmp;
        }

        public ComplexMatrix SubMatrix(int Rs, int Cs, int offsetR = 0, int offsetC = 0)
        {
            var tmp = new ComplexMatrix(Rs, Cs);
            foreach (var (r, c) in tmp.Indexes())
            {
                tmp[r, c] = this[offsetR + r, c + offsetC];
            }
            return tmp;
        }

        public override string ToString()
        {
            return $"(R:{rows}, C:{colums}, V:[{string.Join(",", Enumerable.Range(0, Math.Min(20, values.Count)).Select(m => values[m]))}])";
        }


        public string AsMatrixFormString(int round = -1)
        {
            var N = values.Select(k => k.ToString().Split('.').First().Length).Max();
            var D = round != -1 ? round : values.Select(k => k.ToString().Split('.').Last().Length).Max();
            var format = GetStringFormat(N, D);
            return string.Join(Environment.NewLine, Enumerable.Range(0, Rows).Select(r => string.Join("  ", Enumerable.Range(0, Columns).Select(c => round > 0 ? this[r, c].Round(round) : this[r, c]).Select(k => k.ToString(format)))));
        }

        public string MatrixString5 => AsMatrixFormString(5);

        public string MatrixFullString => string.Join(",", values);
        public string GetStringFormat(int N, int D)
        {
            return "{0," + (N + D) + ":" + new String('0', N) + "." + new String('0', D) + "}";
        }
        public ComplexVector Diagonal()
        {
            var result = ComplexVector.Zero(Columns);
            for (int i = 0; i < Columns; i++)
            {
                result[i] = this[i, i];
            }
            return result;
        }

        public ComplexNumber this[int index]
        {
            get
            {
                return values[index];
            }
            set
            {
                values[index] = value;
            }
        }

        public ComplexMatrix ResetValues(double value = 0)
        {
            for (int i = 0; i < values.Count; i++)
            {
                values[i] = value;
            }
            return this;
        }
        public ComplexMatrix ResetLowerTriangular(double value = 0)
        {
            if (IsSquare)
            {
                foreach (var (r, c) in Indexes())
                {
                    if (r > c)
                    {
                        this[r, c] = value;
                    }
                }
            }
            else
            {
                throw new Exception("Is Not Square Matrix");

            }

            return this;
        }

        public IEnumerable<int> UpperTriangleIndex()
        {
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < Columns; c++)
                {
                    if (r <= c)
                    {
                        yield return CellOffsetInValue(r, c);
                    }
                }
            }
        }

        public IEnumerable<int> TransposeIndexes()
        {
            foreach (var (r, c) in Indexes())
            {
                yield return CellOffsetInValue(c, r);
            }

        }
        public static int UpTriangleSize(int rows, int cols)
        {
            int i = 0;
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < cols; c++)
                {
                    if (r <= c)
                    {
                        i++;
                    }
                }

            }
            return i;
        }

        public virtual ComplexNumber this[int rowIndex, int columnIndex]
        {
            get
            {
                return values[CellOffsetInValue(rowIndex, columnIndex)];
            }
            set
            {
                values[CellOffsetInValue(rowIndex, columnIndex)] = value;
            }
        }
        public bool IsSquare => rows == colums;
        public static ComplexMatrix IdentityMatrix(int iRows, int iCols)   // Function generates the identity matrix
        {
            ComplexMatrix matrix = new ComplexMatrix(iRows, iCols);
            for (int i = 0; i < Math.Min(iRows, iCols); i++)
            {
                matrix[i, i] = 1;
            }
            return matrix;
        }
        public bool IsIdentityMatrix
        {
            get
            {
                foreach (var (r, c) in Indexes())
                {
                    if (r != c)
                    {
                        if (this[r, c].Absolute > 1e-8)
                        {
                            return false;
                        }
                    }
                    else
                    {
                        if (Math.Abs(this[r, c].Absolute - 1) > 1e-8)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }
        }
        public static ComplexMatrix IdentityMatrix(int iRows)   // Function generates the identity matrix
        {
            return IdentityMatrix(iRows, iRows);
        }

        public static ComplexMatrix ColumnwiseMatrix(ComplexMatrix main, params ComplexMatrix[] matrices)
        {
            var columns = matrices.Select(mat => mat.Columns).Sum() + main.Columns;
            var final = new ComplexMatrix(main.Rows, columns);
            final.CopyFrom(main);
            int offsetC = main.colums;
            foreach (var a in matrices)
            {
                final.CopyFrom(a, 0, offsetC);
            }
            return final;
        }
        public static ComplexMatrix RowwiseMatrix(ComplexMatrix main, params ComplexMatrix[] matrices)
        {
            var rows = matrices.Select(mat => mat.rows).Sum() + main.rows;
            var final = new ComplexMatrix(rows, main.colums);
            final.CopyFrom(main);
            int offsetR = main.rows;
            foreach (var a in matrices)
            {
                final.CopyFrom(a, offsetR, 0);
            }
            return final;
        }
        public ComplexMatrix Duplicate()                   // Function returns the copy of this matrix
        {
            ComplexMatrix matrix = new ComplexMatrix(rows, colums);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < colums; j++)
                    matrix[i, j] = this[i, j];
            return matrix;
        }

        public static ComplexMatrix FromRowiseVector(List<ComplexVector> values)
        {
            var final = new ComplexMatrix(values.Count, values[0].Size, Flatten(values.Select(k => k.Select(u => (ComplexNumber)u))));
            return final;
        }

        private static IEnumerable<T> Flatten<T>(IEnumerable<IEnumerable<T>> values)
        {
            foreach (var a in values)
            {
                foreach (var b in a)
                {
                    yield return b;
                }
            }
        }

        class SignularException : Exception
        {

        }

        public const double RelativeTolerance = 1e-6;
        public const double AbsoluteTolerance = 1e-7;

        /// <summary>
        /// Find d, such that L * d = b
        /// </summary>
        /// <param name="lowTriangular"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static ComplexVector SolveLower(ComplexMatrix L, ComplexVector b)
        {
            var d = ComplexVector.Constant(b.Size, 1);
            for (int i = 0; i < b.Size; i++)
            {
                ComplexNumber sum = 0.0;
                for (int k = 0; k < i; k++)
                {
                    sum += L[i, k] * d[k];
                }
                d[i] = b[i] - sum;
            }
            return d;
        }

        /// <summary>
        /// Find x, such that U * x = d
        /// </summary>
        /// <param name="U"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        public static ComplexVector SolveUpper(ComplexMatrix U, ComplexVector d)
        {
            var x = ComplexVector.Constant(d.Size, 1);
            for (int i = d.Size - 1; i >= 0; i--)
            {
                ComplexNumber sum = 0.0;
                for (int j = d.Size - 1; j > i; j--)
                {
                    sum += U[i, j] * x[j];
                }
                x[i] = (d[i] - sum) / U[i, i];
            }
            return x;
        }

        private bool IsDoubleEqual(double x, double y, double error = 1e-4)
        {
            return Math.Abs(x - y) <= error * Math.Sqrt(Math.Abs(x * y));
        }


        private static Exception IsNotSquareMatrixException => new Exception("Is Not Square Matrix");
        public IEnumerable<(int r, int c)> Indexes()
        {
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < colums; c++)
                {
                    yield return (r, c);
                }
            }
        }
        public IEnumerable<(int r, int c, int i)> WholeIndexes()
        {
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < colums; c++)
                {
                    yield return (r, c, r * colums + c);
                }
            }
        }

        private void InterchangeRow(int RowA, int RowB)
        {
            for (int c = 0; c < this.colums; c++)
            {
                var prev = this[RowB, c];
                this[RowB, c] = this[RowA, c];
                this[RowA, c] = prev;
            }
        }
        private void AddRowBWithScalarRowA(double scalar, int RowA, int RowB)
        {
            for (int c = 0; c < this.colums; c++)
            {
                this[RowB, c] += scalar * this[RowA, c];
            }
        }


        private double RoundToZero(double value)
        {
            return value > 0 ? value : (Math.Abs(value) < 1e-9 ? 0 : value);
        }

        public static ComplexMatrix Diagonal(int row, int column, double[] s, int size = -1)
        {
            var final = new ComplexMatrix(row, column);
            size = size == -1 ? s.Length : size;
            for (int i = 0; i < size; i++)
            {
                final[i, i] = s[i];
            }
            return final;
        }

        public static ComplexMatrix Diagonal(IEnumerable<double> v)
        {
            var count = v.Count();
            var final = new ComplexMatrix(count);
            for (int i = 0; i < count; i++)
            {
                final[i, i] = v.ElementAt(i);
            }
            return final;
        }

        public ComplexMatrix Mul(ComplexMatrix B)
        {
            if (this.Columns == B.Rows)
            {
                var result = new ComplexMatrix(rows, B.colums);
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < B.colums; j++)
                    {
                        ComplexNumber sum = 0.0;
                        for (int k = 0; k < B.rows; k++)
                        {
                            sum += this[i, k] * B[k, j];
                        }
                        result[i, j] = sum;
                    }
                }
                return result;
            }
            else
            {
                throw new Exception($"Matrix({this.Rows} x {this.Columns}) * Matrix({B.Rows} x {B.Columns})");
            }
        }

        public ComplexVector Mul(ComplexVector B)
        {
            if (this.colums == B.Size)
            {
                var result = ComplexVector.Zero(Rows);
                for (int i = 0; i < rows; i++)
                {
                    ComplexNumber sum = 0.0;
                    for (int j = 0; j < colums; j++)
                    {
                        sum += this[i, j] * B[j];
                    }
                    result[i] = sum;
                }
                return result;
            }
            else
            {
                throw new Exception($"Matrix({this.rows} x {this.colums}) * Vec({B.Size}) Is Not Allowed");
            }
        }

        public bool Equals(ComplexMatrix other)
        {
            if (rows == other.rows && colums == other.colums)
            {
                return values.SequenceEqual(other.values);
            }
            else
            {
                return false;
            }
        }

        public void AssertApproximateEquals(ComplexMatrix other)
        {
            if (!ApproximateEquals(other))
            {
                throw new Exception($"{this.MatrixString5}\n Not Approximate Equal \n{other.MatrixString5}");
            }
        }

        public bool ApproximateEquals(ComplexMatrix other, double error = AbsoluteTolerance)
        {
            if (rows == other.rows && colums == other.colums)
            {
                foreach (var (r, c) in this.Indexes())
                {
                    var a = this[r, c];
                    var b = other[r, c];
                    if (!pUtils.DoubleApproximateEqual((a - b).Absolute, RelativeTolerance, error))
                    {
                        return false;
                    }
                }
                return true;
            }
            else
            {
                return false;
            }
        }


        public ComplexMatrix MatrixExcludeRC(int rEx, int cEx)
        {
            var p = new ComplexMatrix(this.rows - 1, this.colums - 1);
            int index = 0;
            for (int r = 0; r < this.rows; r++)
            {
                if (r != rEx)
                {
                    for (int c = 0; c < this.colums; c++)
                    {
                        if (c != cEx)
                        {
                            p.values[index++] = this[r, c];
                        }
                    }
                }
            }
            return p;
        }


        public ComplexMatrix Transpose()
        {
            var final = new ComplexMatrix(colums, rows);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < colums; j++)
                {
                    final[j, i] = this[i, j];
                }
            }
            return final;
        }

        public ComplexMatrix TransposeWithInPlace(bool inPlace)
        {
            if (inPlace)
            {
                if (rows != colums)
                {
                    throw new ArgumentException("Only square matrices can be transposed in place.", "matrix");
                }
                for (int i = 0; i < rows; i++)
                {
                    for (int j = i; j < colums; j++)
                    {
                        var element = this[j, i];
                        this[j, i] = this[i, j];
                        this[i, j] = element;
                    }
                }
                return this;
            }
            else
            {
                return Transpose();
            }
        }

        public void CopyFrom(ComplexMatrix matrix, int offsetR = 0, int offsetC = 0)
        {
            foreach (var (r, c) in matrix.Indexes())
            {
                this[r + offsetR, c + offsetC] = matrix[r, c];
            }
        }
        public void AddFrom(ComplexMatrix matrix, int offsetR = 0, int offsetC = 0)
        {
            foreach (var (r, c) in matrix.Indexes())
            {
                this[r + offsetR, c + offsetC] += matrix[r, c];
            }
        }

        public enum CopyType { Row, Column }
        public void CopyFrom(ComplexVector vector, int index = 0, CopyType copyType = CopyType.Column)
        {
            for (int i = 0; i < vector.Size; i++)
            {
                switch (copyType)
                {
                    case CopyType.Row:
                        this[index, i] = vector[i];
                        break;
                    case CopyType.Column:
                        this[i, index] = vector[i];
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }
        public void AddFrom(ComplexVector vector, int index = 0, CopyType copyType = CopyType.Column)
        {
            for (int i = 0; i < vector.Size; i++)
            {
                switch (copyType)
                {
                    case CopyType.Row:
                        this[index, i] += vector[i];
                        break;
                    case CopyType.Column:
                        this[i, index] += vector[i];
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }
        public ComplexMatrix KroneckerProduct(ComplexVector vector)
        {
            return KroneckerProduct(this, ComplexMatrix.NewMatrixFromColumns(vector));
        }
        public static ComplexMatrix NewMatrixFromColumns(List<ComplexVector> columns)
        {
            var m = new ComplexMatrix(columns[0].Size, columns.Count);
            foreach (var (r, c) in m.Indexes())
            {
                m[r, c] = columns[c][r];
            }
            return m;

        }
        public static ComplexMatrix NewMatrixFromColumns(ComplexVector a, params ComplexVector[] rest)
        {
            return NewMatrixFromColumns(new[] { a }.Concat(rest).ToList());
        }
        public ComplexMatrix KroneckerProduct(ComplexMatrix matrix)
        {
            return KroneckerProduct(this, matrix);
        }
        public static ComplexMatrix KroneckerProduct(ComplexMatrix A, ComplexMatrix B)
        {
            //A : m x n,   B : p x q
            var final = new ComplexMatrix(A.rows * B.rows, A.colums * B.colums);
            foreach (var (m, n) in A.Indexes())
            {
                var a_mn = A[m, n];
                foreach (var (p, q) in B.Indexes())
                {
                    final[m * B.rows + p, n * B.colums + q] = a_mn * B[p, q];
                }
            }
            return final;
        }

        #region opertor  
        public static ComplexMatrix operator *(ComplexMatrix m, double scalar)
        {
            var final = new ComplexMatrix(m.Rows, m.Columns, m.values.Select(k => k * scalar).ToArray());
            return final;
        }
        public static ComplexMatrix operator /(ComplexMatrix m, double scalar)
        {
            var final = new ComplexMatrix(m.Rows, m.Columns, m.values.Select(k => k / scalar).ToArray());
            return final;
        }
        public static ComplexMatrix operator *(ComplexVector vector, ComplexMatrix m)
        {
            if (m.Rows == 1)
            {
                var final = new ComplexMatrix(vector.Size, m.Columns);
                foreach (var item in final.Indexes())
                {
                    final[item.r, item.c] = vector[item.r] * m[0, item.c];
                }
                return final;
            }
            else
            {
                throw new Exception($"Not Able Perfom ({vector.Size} x 1) Mul ({m.rows} x {m.colums})");
            }
        }

        public static ComplexMatrix operator *(double scalar, ComplexMatrix m)
        {
            return m * scalar;
        }
        public static ComplexMatrix operator +(ComplexMatrix m, double scalar)
        {
            var final = new ComplexMatrix(m.Rows, m.Columns, m.values.Select(k => k + scalar));
            return final;
        }
        public static ComplexMatrix operator -(ComplexMatrix m, double scalar)
        {
            var final = new ComplexMatrix(m.Rows, m.Columns, m.values.Select(k => k - scalar));
            return final;
        }
        public static ComplexMatrix operator +(ComplexMatrix a, ComplexMatrix b)
        {
            if (a.Rows == b.Rows && b.Columns == a.Columns)
            {
                var final = new ComplexMatrix(a.Rows, a.Columns, a.values.Zip(b.values, (d, t) => d + t));
                return final;
            }
            else
            {
                throw new Exception("Shape Not Equal");
            }
        }

        /// <summary>
        ///  A.Mul(B)
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static ComplexMatrix operator *(ComplexMatrix A, ComplexMatrix B)
        {
            return A.Mul(B);
        }

        /// <summary>
        ///  A.Mul(B)
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static ComplexVector operator *(ComplexMatrix A, ComplexVector B)
        {
            return A.Mul(B);
        }
        public static ComplexMatrix operator -(ComplexMatrix a)
        {
            var final = new ComplexMatrix(a.Rows, a.Columns, a.values.Select(k => -k));
            return final;
        }
        public static ComplexMatrix operator -(ComplexMatrix a, ComplexMatrix b)
        {
            if (a.Rows == b.Rows && b.Columns == a.Columns)
            {
                var final = new ComplexMatrix(a.Rows, a.Columns, a.values.Zip(b.values, (d, t) => d - t).ToArray());
                return final;
            }
            else
            {
                throw new Exception("Shape Not Equal");
            }
        }

        #endregion

    }
}
