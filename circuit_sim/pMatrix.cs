using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;

namespace circuit_sim
{
    public class pMatrix : IEquatable<pMatrix>
    {
        protected int rows;
        protected int colums;
        List<double> values;
        public int Columns => colums;
        public int Rows => rows;
        public virtual IEnumerable<double> Values => values;

        protected pMatrix() { }
        public pMatrix(int row) : this(row, row)
        {
        }

        public int CellOffsetInValue(int rowIndex, int columnIndex)
        {
            return rowIndex * colums + columnIndex;
        }
        public pMatrix(int row, int colum)
        {
            this.rows = row;
            this.colums = colum;
            values = new List<double>(row * colum);
            int size = row * colum;
            for (int i = 0; i < size; i++)
            {
                values.Add(0);
            }
        }
        public pMatrix(int row, int colum, IEnumerable<double> values) : this(row, colum)
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

        public pMatrix Clone()
        {
            return new pMatrix(rows, colums, values);
        }

        public pMatrix Sqrt(double error = 1e-9)
        {
            return new pMatrix(rows, colums, values.Select(k => Math.Abs(k) < error ? 0 : Math.Sqrt(k)));
        }

        public double Sum()
        {
            return values.Sum();
        }

        public double Rank()
        {
            if (this.rows >= this.colums)
            {
                var (Q, R) = this.HouseholderQRDecomposition();
                return R.Diagonal().Where(k => Math.Abs(k) > 1e-8).Count();
            }
            else
            {
                return this.Transpose().Rank();
            }
        }

        public pMatrix CropMatrix(int offsetR, int offsetC)
        {
            var tmp = new pMatrix(this.rows - offsetR, this.colums - offsetC);
            foreach (var (r, c) in tmp.Indexes())
            {
                tmp[r, c] = this[r + offsetR, c + offsetC];
            }
            return tmp;
        }

        public pMatrix SubMatrix(int Rs, int Cs, int offsetR = 0, int offsetC = 0)
        {
            var tmp = new pMatrix(Rs, Cs);
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
            return string.Join(Environment.NewLine, Enumerable.Range(0, Rows).Select(r => string.Join("  ", Enumerable.Range(0, Columns).Select(c => round > 0 ? Math.Round(this[r, c], round) : this[r, c]).Select(k => GetString(format, k)))));
        }

        public string MatrixString5 => AsMatrixFormString(5);

        public string MatrixFullString => string.Join(",", values);

        private string GetString(string format, double value)
        {
            if (value < 0)
            {
                return string.Format(format, value);
            }
            else
            {
                return " " + string.Format(format, value);
            }
        }
        public string GetStringFormat(int N, int D)
        {
            return "{0," + (N + D) + ":" + new String('0', N) + "." + new String('0', D) + "}";
        }
        public pVector Diagonal()
        {
            var result = pVector.Zero(Columns);
            for (int i = 0; i < Columns; i++)
            {
                result[i] = this[i, i];
            }
            return result;
        }

        public double this[int index]
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

        public pMatrix ResetValues(double value = 0)
        {
            for (int i = 0; i < values.Count; i++)
            {
                values[i] = value;
            }
            return this;
        }
        public pMatrix ResetLowerTriangular(double value = 0)
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
        public pMatrix Inverse()
        {
            if (!IsSquare)
            {
                throw new Exception("Is Not Square Matrix");
            }
            else
            {
                if (this.colums == 1)
                {
                    if (this[0] != 0)
                    {
                        return new pMatrix(1, 1, this.values.Select(k => 1.0 / k));
                    }
                    else
                    {
                        throw new Exception("Singular Matrix");
                    }
                }
                else if (this.colums == 2)
                {
                    var determinate = this.values[0] * this.values[3] - this.values[1] * this.values[2];
                    var detRecipical = 1.0 / determinate;
                    if (determinate != 0)
                    {
                        var tmp = new pMatrix(2, 2);
                        tmp.values[0] = detRecipical * this.values[3];
                        tmp.values[1] = -detRecipical * this.values[1];
                        tmp.values[2] = -detRecipical * this.values[2];
                        tmp.values[3] = detRecipical * this.values[0];
                        return tmp;
                    }
                    else
                    {
                        throw new Exception("Singular Matrix");
                    }
                }
                else
                {
                    try
                    {
                        var lup = LUPDecomposition();
                        var vectors = Enumerable.Range(0, rows).Select(m => lup.Solve(pVector.OneHot(m, rows))).ToArray();
                        var final = new pMatrix(rows, colums);
                        for (int i = 0; i < rows; i++)
                        {
                            for (int j = 0; j < colums; j++)
                            {
                                final[i, j] = vectors[j][i];
                            }
                        }
                        return final;
                    }
                    catch (SignularException)
                    {
                        throw new Exception("Singular Matrix");
                    }
                }
            }
        }
        public static pMatrix ExtractFrom(pVector vector, int rows, int cols, int offset = 0)
        {
            return new pMatrix(rows, cols, vector.Skip(offset).Take(rows * cols));
        }
        public static pMatrix UpTriangleExtractFrom(pVector vector, int rows, int cols, int offset = 0)
        {
            var final = new pMatrix(rows, cols);
            int i = 0;
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < cols; c++)
                {
                    if (r <= c)
                    {
                        final[r, c] = vector[i + offset];
                        i++;
                    }
                }

            }
            return final;
        }

        public pVector ValueAtIndexesAsVector(int[] indexes)
        {
            var a = pVector.NewWithLen(indexes.Length);
            for (int i = 0; i < indexes.Length; i++)
            {
                a[i] = this.values[indexes[i]];
            }
            return a;
        }

        public pMatrix SelectedIndexesAsMatrix(int[] rowIndexes, int[] columnIndexes)
        {
            var final = new pMatrix(rowIndexes.Length, columnIndexes.Length);
            foreach (var (r, c) in final.Indexes())
            {
                final[r, c] = this[rowIndexes[r], columnIndexes[c]];
            }
            return final;
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
        public pVector AsVector()
        {
            return new pVector(this.values);
        }

        public virtual double this[int rowIndex, int columnIndex]
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
        public static pMatrix IdentityMatrix(int iRows, int iCols)   // Function generates the identity matrix
        {
            pMatrix matrix = new pMatrix(iRows, iCols);
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
                        if (Math.Abs(this[r, c]) > 1e-8)
                        {
                            return false;
                        }
                    }
                    else
                    {
                        if (Math.Abs(this[r, c] - 1) > 1e-8)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }
        }
        public static pMatrix IdentityMatrix(int iRows)   // Function generates the identity matrix
        {
            return IdentityMatrix(iRows, iRows);
        }

        public static pMatrix ColumnwiseMatrix(pMatrix main, params pMatrix[] matrices)
        {
            var columns = matrices.Select(mat => mat.Columns).Sum() + main.Columns;
            var final = new pMatrix(main.Rows, columns);
            final.CopyFrom(main);
            int offsetC = main.colums;
            foreach (var a in matrices)
            {
                final.CopyFrom(a, 0, offsetC);
            }
            return final;
        }
        public static pMatrix RowwiseMatrix(pMatrix main, params pMatrix[] matrices)
        {
            var rows = matrices.Select(mat => mat.rows).Sum() + main.rows;
            var final = new pMatrix(rows, main.colums);
            final.CopyFrom(main);
            int offsetR = main.rows;
            foreach (var a in matrices)
            {
                final.CopyFrom(a, offsetR, 0);
            }
            return final;
        }
        public pMatrix Duplicate()                   // Function returns the copy of this matrix
        {
            pMatrix matrix = new pMatrix(rows, colums);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < colums; j++)
                    matrix[i, j] = this[i, j];
            return matrix;
        }
        public class LUPMatrix
        {
            public pMatrix P;
            public pMatrix L;
            public pMatrix U;
            public int Exchanges = 0;
            /// <summary>
            ///  L * U * V = P*B, Solve V
            /// </summary>
            /// <param name="values"></param>
            /// <returns></returns>
            public pVector Solve(pVector B)
            {
                // let Y = U* V
                // L * Y = P*B
                // 
                return SolveUpper(SolveLower(P.Mul(B)));
            }

            /// <summary>
            /// L * V = B
            /// </summary>
            /// <param name="Bs"></param>
            /// <returns></returns>
            private pVector SolveLower(pVector B)
            {
                var final = pVector.Constant(B.Size, 1);
                for (int i = 0; i < B.Size; i++)
                {
                    var sum = 0.0;
                    for (int j = 0; j < i; j++)
                    {
                        sum += this.L[i, j] * final[j];
                    }
                    final[i] = B[i] - sum;
                }
                var Bs = L.Mul(final);
                return final;
            }

            /// <summary>
            /// U * V = B
            /// </summary>
            /// <param name="B"></param>
            /// <returns></returns>
            private pVector SolveUpper(pVector B)
            {
                var final = pVector.Constant(B.Size, 1);
                for (int i = B.Size - 1; i >= 0; i--)
                {
                    var sum = 0.0;
                    for (int j = B.Size - 1; j > i; j--)
                    {
                        sum += this.U[i, j] * final[j];
                    }
                    final[i] = (B[i] - sum) / U[i, i];
                }
                return final;
            }
            public double Determinant()
            {
                var D = 1.0;
                for (int i = 0; i < U.rows; i++)
                {//product the diagonal elements of the LU matrix
                    D = D * U[i, i];
                }
                return D * Math.Pow(-1, Exchanges);
            }
        }

        public static pMatrix FromRowiseVector(List<pVector> values)
        {
            var final = new pMatrix(values.Count, values[0].Size, Flatten(values.Select(k => k.Select(u => u))));
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
        public LUPMatrix LUPDecomposition()
        {
            int n = rows;
            var A = Duplicate();
            if (!IsSquare)
            {
                throw new Exception("Is Not Square Matrix");
            }
            var Pi = Enumerable.Range(0, n).Select(m => m).ToArray();
            int exchanges = 0;
            for (int k = 0; k < n - 1; k++)
            {
                double p = 0;
                int kquote = -1;
                for (int i = k; i < n; i++)
                {
                    if (Math.Abs(A[i, k]) > p)
                    {
                        p = Math.Abs(A[i, k]);
                        kquote = i;
                    }
                }
                if (p == 0)
                {
                    throw new SignularException();
                }

                if (kquote != k)
                {//swap rows
                    exchanges++;
                }
                //exchange pi[k] with pi[k']
                int p0 = Pi[k];
                Pi[k] = Pi[kquote];
                Pi[kquote] = p0;

                //exechange aki with ak'i
                for (int i = 0; i < n; i++)
                {
                    double p1 = A[k, i];
                    A[k, i] = A[kquote, i];
                    A[kquote, i] = p1;
                }

                for (int i = k + 1; i < n; i++)
                {
                    A[i, k] = A[i, k] / A[k, k];
                    for (int j = k + 1; j < n; j++)
                    {
                        A[i, j] = A[i, j] - A[i, k] * A[k, j];
                    }
                }
            }
            var L = new pMatrix(rows, colums);
            var U = new pMatrix(rows, colums);
            var P = new pMatrix(rows, colums);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < colums; j++)
                {
                    if (i > j)
                    {
                        L[i, j] = A[i, j];
                        U[i, j] = 0;
                    }
                    else
                    {
                        L[i, j] = i == j ? 1 : 0;
                        U[i, j] = A[i, j];
                    }
                    P[i, j] = Pi[i] == j ? 1 : 0;
                }
            }
            return new LUPMatrix { L = L, P = P, U = U, Exchanges = exchanges };
        }

        public const double RelativeTolerance = 1e-6;
        public const double AbsoluteTolerance = 1e-7;
        public bool IsSymmetric
        {
            get
            {
                if (IsSquare)
                {
                    foreach (var item in Indexes())
                    {
                        var a = this[item.r, item.c];
                        var b = this[item.c, item.r];
                        if (Math.Abs(a - b) > Math.Max(RelativeTolerance * Math.Max(Math.Abs(a), Math.Abs(b)), AbsoluteTolerance))
                        {
                            return false;
                        }
                    }
                    return true;
                }
                return false;
            }
        }


        /// <summary>
        /// Find d, such that L * d = b
        /// </summary>
        /// <param name="lowTriangular"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static pVector SolveLower(pMatrix L, pVector b)
        {
            var d = pVector.Constant(b.Size, 1);
            for (int i = 0; i < b.Size; i++)
            {
                var sum = 0.0;
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
        public static pVector SolveUpper(pMatrix U, pVector d)
        {
            var x = pVector.Constant(d.Size, 1);
            for (int i = d.Size - 1; i >= 0; i--)
            {
                var sum = 0.0;
                for (int j = d.Size - 1; j > i; j--)
                {
                    sum += U[i, j] * x[j];
                }
                x[i] = (d[i] - sum) / U[i, i];
            }
            return x;
        }

        /// <summary>
        /// Find V such that, Hessian * V = Gradient
        /// 
        /// </summary>
        /// <param name="Hessian"></param>
        /// <param name="Gradient"></param>
        /// <returns></returns>
        public static pVector SolveViaLUDecomposition(pMatrix Hessian, pVector Gradient)
        {
            var (L, U) = Hessian.LUDecomposition();
            var d = pMatrix.SolveLower(L, Gradient);
            var x = pMatrix.SolveUpper(U, d);
            return x;
        }

        private bool IsDoubleEqual(double x, double y, double error = 1e-4)
        {
            return Math.Abs(x - y) <= error * Math.Sqrt(Math.Abs(x * y));
        }

        public (pMatrix Q, pMatrix R) QRDecomposition()
        {
            if (this.IsSquare)
            {
                //Using the Gram–Schmidt process
                List<pVector> eks = new List<pVector>();
                var R = new pMatrix(this.Columns);
                for (int c = 0; c < this.colums; c++)
                {
                    var ak = ColumnVectorAt(c);
                    var uk = ak;
                    for (int j = 0; j < c; j++)
                    {
                        var ak_dot_ej = ak.Dot(eks[j]);
                        R[j, c] = ak_dot_ej;
                        uk -= ak_dot_ej * eks[j];
                    }
                    var ek = uk.Normalized();
                    eks.Add(ek);
                    R[c, c] = uk.TwoNorm(); //looks like more stable 
                    //R[c, c] = ak.Dot(eks[c]);
                }
                var Q = NewMatrixFromColumns(eks);
                return (Q, R);
            }
            else
            {
                throw IsNotSquareMatrixException;
            }
        }
        public (pMatrix Q, pMatrix R) HouseholderQRDecomposition()
        {
            if (this.Rows >= this.Columns)
            {
                var Q = pMatrix.IdentityMatrix(this.rows);
                var R = this.Clone();

                for (int j = 0; j < this.colums; j++)
                {
                    var Rj = R.ColumnVectorAt(j, j);
                    var normlX = Rj.TwoNorm();
                    if (normlX != 0)
                    {
                        var s = -(R[j, j] >= 0 ? 1 : -1);
                        var u1 = R[j, j] - s * normlX;
                        if (u1 != 0)
                        {
                            var w = Rj / u1;
                            w[0] = 1;
                            var tau = -s * u1 / normlX;
                            var RDiff = (tau * w) * w.Transpose().Mul(R.CropMatrix(j, 0));
                            var QDiff = Q.CropMatrix(0, j).Mul(w) * (tau * w).Transpose();
                            R.AddFrom(-RDiff, j);
                            Q.AddFrom(-QDiff, 0, j);
                        }
                        else
                        {
                            throw new Exception("HouseholderQRDecomposition U is Zero");
                        }
                    }
                }
                return (Q, R);
            }
            else
            {
                throw new Exception("Condition M >= n not satisfied");
            }
        }

        public pVector ColumnVectorAt(int c, int offset = 0)
        {
            return new pVector(Enumerable.Range(offset, this.rows - offset).Select(r => this[r, c]));
        }

        public pVector RowVectorAt(int r, int offset = 0)
        {
            return new pVector(Enumerable.Range(offset, this.colums - offset).Select(c => this[r, c]));
        }

        public double ValueAtPosition(int positionR, int positionC)
        {
            return this[positionR - 1, positionC - 1];
        }
        public void SetValueAtPosition(int positionR, int positionC, double value)
        {
            this[positionR - 1, positionC - 1] = value;
        }
        public double ValueAtPositionAbs(int positionR, int positionC)
        {
            return Math.Abs(this[positionR - 1, positionC - 1]);
        }

        private void AssertTrue(bool condition)
        {
            if (!condition)
            {
                throw new Exception("Condition Not True");
            }
        }

        public (pMatrix EigenVector, pVector EigenValues) EigenDecompositionIII(double error = 1e-5)
        {
            bool converge = false;

            var Ak = this;

            var Identity = pMatrix.IdentityMatrix(this.rows);
            int t = 0;

            var Qk = IdentityMatrix(this.rows);
            while (!converge)
            {
                var sk = Ak.values.Last() * Identity;
                var (Q, R) = (Ak - sk).QRDecomposition();
                Ak = R.Mul(Q) + sk;

                Qk = Qk * Q;
                //var (Q, R) = (Ak).QRDecomposition();
                //Ak = R.Mul(Q);
                if (Ak.IsUpperTriangular(error))
                {
                    converge = true;
                    break;
                }
                t++;
            }

            return (Qk, Ak.Diagonal());
        }




        private (pMatrix U, pMatrix D) EigenVectorOf(pMatrix eigenValues, pMatrix U)
        {
            List<pVector> eigenVectorList = new List<pVector>();
            var Inn = pMatrix.IdentityMatrix(eigenValues.rows);
            for (int i = 0; i < eigenValues.rows; i++)
            {
                var lamda = eigenValues[i, i];
                var eigenv = SolveHomogenousSystem(this - lamda * Inn).Normalized();
                //var eigenv2 = (this * eigenv) / lamda;
                //eigenv2.AssertApproximateEqual(eigenv);
                eigenVectorList.Add(eigenv);
            }


            var eigenVectors = NewMatrixFromColumns(eigenVectorList);

            (eigenVectors * eigenValues * eigenVectors.Inverse()).AssertApproximateEquals(this);

            return (eigenVectors, eigenValues);
        }

        public (pMatrix EigenVector, pVector EigenValues) EigenDecompositionII(double error = 1e-5)
        {
            bool converge = false;
            var Ak = this;

            var Identity = pMatrix.IdentityMatrix(this.rows);
            int t = 0;
            var eigenValues = pVector.NewWithLen(this.rows);
            int n = Ak.rows;
            while (!converge)
            {
                var sk = Ak[n - 1, n - 1] * Identity;
                var (Q, R) = (Ak - sk).HouseholderQRDecomposition();
                Ak = R.Mul(Q) + sk;
                t++;
                if (Enumerable.Range(0, n - 1).All(k => Math.Abs(Ak[n - 1, k]) < error))
                {
                    eigenValues[n - 1] = Ak[n - 1, n - 1];
                    if (n > 2)
                    {
                        n -= 1;
                        Ak = Ak.SubMatrix(n, n);
                        Identity = Identity.SubMatrix(n, n);
                    }
                    else
                    {
                        eigenValues[0] = Ak[0, 0];
                        break;
                    }
                }
            }
            eigenValues = new pVector(eigenValues.Size, eigenValues.OrderByDescending(k => k));
            List<pVector> eignVectors = new List<pVector>();
            var Inn = pMatrix.IdentityMatrix(eigenValues.Size);
            for (int i = 0; i < eigenValues.Size; i++)
            {
                var lamda = eigenValues[i];
                var eigenv = SolveHomogenousSystem(this - lamda * Inn);
                eignVectors.Add(eigenv.Normalized());
            }
            return (pMatrix.NewMatrixFromColumns(eignVectors), eigenValues);

        }

        public static pVector SolveHomogenousSystem(pMatrix a)
        {
            if (a.IsSquare)
            {
                var (Q, R) = a.HouseholderQRDecomposition();
                var final = pVector.Constant(a.rows, 1);
                for (int r = final.Size - 2; r >= 0; r--)
                {
                    if (R[r, r] != 0)
                    {
                        var sum = 0.0;
                        for (int j = r + 1; j < final.Size; j++)
                        {
                            sum += R[r, j] * final[j];
                        }
                        final[r] = -sum / R[r, r];
                    }
                }
                return final;
            }
            else
            {
                throw IsNotSquareMatrixException;
            }
        }
        public bool IsDiagonalMatrix(double error)
        {
            if (IsSquare)
            {
                foreach (var (r, c) in Indexes())
                {
                    if (r != c && Math.Abs(this[r, c]) > error)
                    {
                        return false;
                    }
                }
                return true;
            }
            else
            {
                throw IsNotSquareMatrixException;
            }
        }
        public bool IsUpperTriangular(double error)
        {
            foreach (var (r, c) in Indexes())
            {
                if (r > c && Math.Abs(this[r, c]) > error)
                {
                    return false;
                }
            }
            return true;
        }

        public static pMatrix NewMatrixFromColumns(List<pVector> columns)
        {
            var m = new pMatrix(columns[0].Size, columns.Count);
            foreach (var (r, c) in m.Indexes())
            {
                m[r, c] = columns[c][r];
            }
            return m;

        }

        public static pMatrix NewMatrixFromColumns(pVector a, params pVector[] rest)
        {
            return NewMatrixFromColumns(new[] { a }.Concat(rest).ToList());
        }

        public static pMatrix New(int rows, int colums = -1)
        {
            return new pMatrix(rows, colums == -1 ? rows : colums);
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

        public bool IsPositiveDefinite
        {
            get
            {
                if (IsSymmetric)
                {
                    var pivote = Pivots();
                    return pivote.All(k => k > 0 || Math.Abs(k) < 1e-11);
                }
                return false;
            }
        }

        public bool IsNegativeDefinite
        {
            get
            {
                if (IsSymmetric)
                {
                    var pivote = Pivots();
                    return pivote.All(k => k < 0 || Math.Abs(k) < 1e-11);
                }
                return false;
            }
        }
        public bool IsPositiveSemiDefinite
        {
            get
            {
                if (IsSquare)
                {
                    var pivote = Pivots();
                    return pivote.All(k => k >= 0);
                }
                return false;
            }
        }

        public pMatrix GaussianElimation()
        {
            var D = this.Clone();
            for (int c = 0; c < this.colums - 1; c++)
            {
                //try make entry[c, c] not zero by interchange row
                if (D[c, c] == 0)
                {
                    var noZeroPr = pUtils.LoopFromAToB(c, this.rows - 1).Where(k => D[k, c] != 0).Select(k => (int?)k).FirstOrDefault();
                    if (noZeroPr.HasValue)
                    {
                        D.InterchangeRow(c, noZeroPr.Value);
                    }
                }
                for (int r = c + 1; r < this.rows; r++)
                {
                    if (D[r, c] != 0)
                    {
                        var scalar = -D[r, c] / D[c, c];
                        D.AddRowBWithScalarRowA(scalar, c, r);
                    }
                }

            }
            return D;
        }

        public double[] Pivots()
        {
            if (IsSquare)
            {
                var D = GaussianElimation();
                return Enumerable.Range(0, this.rows).Select(k => D[k, k]).ToArray();
            }
            else
            {
                throw new Exception("Pivots Only Support Square matrices");
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

        public pMatrix CholeskyDecomposition()
        {
            //LLDecomposition previous name
            //Decomposition A matrix to L.Mul(L.Transpose()) , return L
            if (IsPositiveDefinite)
            {
                var L = new pMatrix(rows, colums);
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j <= i; j++)
                    {
                        double sum = 0;
                        for (int k = 0; k < j; k++)
                        {
                            sum += L[i, k] * L[j, k];
                        }

                        L[i, j] = i == j ? Math.Sqrt(RoundToZero(this[i, i] - sum)) : (1.0 / L[j, j] * (this[i, j] - sum));
                    }
                }
                return L;
            }
            else
            {
                throw new Exception("The Cholesky factorization is only for symmetric positive definite matrices");
            }
        }

        public (pMatrix L, pMatrix U) LUDecomposition()
        {
            if (!IsSquare)
            {
                throw new Exception("Is Not Square Matrix");
            }
            else
            {
                var L = pMatrix.IdentityMatrix(this.rows);
                var U = pMatrix.IdentityMatrix(this.rows);
                int n = rows;
                for (int i = 0; i < n; i++)
                {
                    // Upper Triangular
                    for (int k = i; k < n; k++)
                    {
                        // Summation of L(i, j) * U(j, k)
                        double sum = 0.0;
                        for (int j = 0; j < i; j++)
                        {
                            sum += (L[i, j] * U[j, k]);
                        }
                        // Evaluating U(i, k)
                        U[i, k] = this[i, k] - sum;
                    }

                    // Lower Triangular
                    for (int k = i; k < n; k++)
                    {
                        if (i == k)
                        {
                            L[i, i] = 1;
                        }
                        else
                        {
                            // Summation of L(k, j) * U(j, i)
                            double sum = 0.0;
                            for (int j = 0; j < i; j++)
                                sum += (L[k, j] * U[j, i]);

                            // Evaluating L(k, i)
                            L[k, i] = (this[k, i] - sum) / U[i, i];
                        }
                    }
                }
                return (L, U);
            }
        }

        private double RoundToZero(double value)
        {
            return value > 0 ? value : (Math.Abs(value) < 1e-9 ? 0 : value);
        }

        public static pMatrix Diagonal(int row, int column, double[] s, int size = -1)
        {
            var final = new pMatrix(row, column);
            size = size == -1 ? s.Length : size;
            for (int i = 0; i < size; i++)
            {
                final[i, i] = s[i];
            }
            return final;
        }

        public static pMatrix Diagonal(IEnumerable<double> v)
        {
            var count = v.Count();
            var final = new pMatrix(count);
            for (int i = 0; i < count; i++)
            {
                final[i, i] = v.ElementAt(i);
            }
            return final;
        }
        public static pMatrix Diagonal(pMatrix A, params pMatrix[] rests)
        {
            var all = new[] { A }.Concat(rests);
            var final = new pMatrix(all.Select(k => k.rows).Sum(), all.Select(k => k.colums).Sum());
            int offsetc = 0;
            int offsetr = 0;
            foreach (var m in all)
            {
                final.CopyFrom(m, offsetr, offsetc);
                offsetc += m.colums;
                offsetr += m.rows;
            }
            return final;
        }
        public static pMatrix Diagonal(IEnumerable<pMatrix> all)
        {
            var final = new pMatrix(all.Select(k => k.rows).Sum(), all.Select(k => k.colums).Sum());
            int offsetc = 0;
            int offsetr = 0;
            foreach (var m in all)
            {
                final.CopyFrom(m, offsetr, offsetc);
                offsetc += m.colums;
                offsetr += m.rows;
            }
            return final;
        }

        public pMatrix Mul(pMatrix B)
        {
            if (this.Columns == B.Rows)
            {
                var result = new pMatrix(rows, B.colums);
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < B.colums; j++)
                    {
                        var sum = 0.0;
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

        public pVector Mul(pVector B)
        {
            if (this.colums == B.Size)
            {
                var result = pVector.Zero(Rows);
                for (int i = 0; i < rows; i++)
                {
                    var sum = 0.0;
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

        public bool Equals(pMatrix other)
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

        public void AssertApproximateEquals(pMatrix other)
        {
            if (!ApproximateEquals(other))
            {
                throw new Exception($"{this.MatrixString5}\n Not Approximate Equal \n{other.MatrixString5}");
            }
        }

        public bool ApproximateEquals(pMatrix other, double error = AbsoluteTolerance)
        {
            if (rows == other.rows && colums == other.colums)
            {
                foreach (var (r, c) in this.Indexes())
                {
                    var a = this[r, c];
                    var b = other[r, c];
                    if (!pUtils.DoubleApproximateEqual(a, b, RelativeTolerance, error))
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

        public double Determinate()
        {
            if (this.rows == 2)
            {
                return this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0];
            }
            else
            {
                //double total = 0;
                //for (int c = 0; c < this.colums; c++)
                //{
                //    var sign = c % 2 == 0 ? 1 : -1;
                //    var subMatrix = MatrixExclude(0, c);
                //    total += this[0, c] * subMatrix.Determinate() * sign;
                //}
                var (Q, R) = QRDecomposition();
                return R.Diagonal().Aggregate((a, b) => a * b);
            }
        }

        public pMatrix MatrixExcludeRC(int rEx, int cEx)
        {
            var p = new pMatrix(this.rows - 1, this.colums - 1);
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


        public double Trace()
        {
            //tr
            if (IsSquare)
            {
                return Enumerable.Range(0, this.colums).Sum(k => this[k, k]);
            }
            else
            {
                throw new Exception("Only Support Square Matrix");
            }
        }
        public pMatrix Transpose()
        {
            var final = new pMatrix(colums, rows);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < colums; j++)
                {
                    final[j, i] = this[i, j];
                }
            }
            return final;
        }

        public pMatrix TransposeWithInPlace(bool inPlace)
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

        public void CopyFrom(pMatrix matrix, int offsetR = 0, int offsetC = 0)
        {
            foreach (var (r, c) in matrix.Indexes())
            {
                this[r + offsetR, c + offsetC] = matrix[r, c];
            }
        }
        public void AddFrom(pMatrix matrix, int offsetR = 0, int offsetC = 0)
        {
            foreach (var (r, c) in matrix.Indexes())
            {
                this[r + offsetR, c + offsetC] += matrix[r, c];
            }
        }

        public enum CopyType { Row, Column }
        public void CopyFrom(pVector vector, int index = 0, CopyType copyType = CopyType.Column)
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
        public void AddFrom(pVector vector, int index = 0, CopyType copyType = CopyType.Column)
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
        public pMatrix KroneckerProduct(pVector vector)
        {
            return KroneckerProduct(this, pMatrix.NewMatrixFromColumns(vector));
        }

        public pMatrix KroneckerProduct(pMatrix matrix)
        {
            return KroneckerProduct(this, matrix);
        }
        public static pMatrix KroneckerProduct(pMatrix A, pMatrix B)
        {
            //A : m x n,   B : p x q
            var final = new pMatrix(A.rows * B.rows, A.colums * B.colums);
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
        public static pMatrix operator *(pMatrix m, double scalar)
        {
            var final = new pMatrix(m.Rows, m.Columns, m.values.Select(k => k * scalar).ToArray());
            return final;
        }
        public static pMatrix operator /(pMatrix m, double scalar)
        {
            var final = new pMatrix(m.Rows, m.Columns, m.values.Select(k => k / scalar).ToArray());
            return final;
        }
        public static pMatrix operator *(pVector vector, pMatrix m)
        {
            if (m.Rows == 1)
            {
                var final = new pMatrix(vector.Size, m.Columns);
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

        public static pMatrix operator *(double scalar, pMatrix m)
        {
            return m * scalar;
        }
        public static pMatrix operator +(pMatrix m, double scalar)
        {
            var final = new pMatrix(m.Rows, m.Columns, m.values.Select(k => k + scalar));
            return final;
        }
        public static pMatrix operator -(pMatrix m, double scalar)
        {
            var final = new pMatrix(m.Rows, m.Columns, m.values.Select(k => k - scalar));
            return final;
        }
        public static pMatrix operator +(pMatrix a, pMatrix b)
        {
            if (a.Rows == b.Rows && b.Columns == a.Columns)
            {
                var final = new pMatrix(a.Rows, a.Columns, a.values.Zip(b.values, (d, t) => d + t));
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
        public static pMatrix operator *(pMatrix A, pMatrix B)
        {
            return A.Mul(B);
        }

        /// <summary>
        ///  A.Mul(B)
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        public static pVector operator *(pMatrix A, pVector B)
        {
            return A.Mul(B);
        }
        public static pMatrix operator -(pMatrix a)
        {
            var final = new pMatrix(a.Rows, a.Columns, a.values.Select(k => -k));
            return final;
        }
        public static pMatrix operator -(pMatrix a, pMatrix b)
        {
            if (a.Rows == b.Rows && b.Columns == a.Columns)
            {
                var final = new pMatrix(a.Rows, a.Columns, a.values.Zip(b.values, (d, t) => d - t).ToArray());
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
