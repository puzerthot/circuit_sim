using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Linq;

namespace circuit_sim
{
    public class pVector : IEnumerable<double>
    {
        List<double> values;
        public int Size => values.Count;
        public IEnumerable<double> Values => values;
        public pVector(IEnumerable<double> values)
        {
            this.values = values.ToList();
        }

        public pVector Expand(params double[] new_values)
        {
            var tmp = pVector.NewWithLen(this.Size + new_values.Length);
            for (int i = 0; i < Size; i++)
            {
                tmp[i] = this.values[i];
            }

            for (int i = 0; i < new_values.Length; i++)
            {
                tmp[i + Size] = new_values[i];
            }
            return tmp;
        }

        public pVector ResetValues(double value = 0)
        {
            for (int i = 0; i < Size; i++)
            {
                values[i] = value;
            }
            return this;
        }

        private pVector() { values = new List<double>(); }
        public static pVector Combine(pVector A, params pVector[] Rests)
        {
            var tmp = new pVector();
            tmp.values.AddRange(A);
            foreach (var v in Rests)
            {
                tmp.values.AddRange(v);
            }
            return tmp;
        }
        public void ExpandWithValue(double value)
        {
            this.values.Add(value);
        }
        public void InitFrom(pVector b, int offset = 0)
        {
            for (int i = 0; i < Math.Min(this.Size, b.Size); i++)
            {
                this.values[i] = b[i + offset];
            }
        }
        public void CopyFrom(pVector b, int offset = 0)
        {
            for (int i = 0; i < b.Size; i++)
            {
                this.values[i + offset] = b[i];
            }
        }
        public void AddFrom(pVector b, int offset = 0)
        {
            for (int i = 0; i < b.Size; i++)
            {
                this.values[i + offset] += b[i];
            }
        }
        public void CopyFrom(pMatrix b, int offset = 0)
        {
            int i = offset;
            foreach (var v in b.Values)
            {
                this.values[offset] = v;
                offset++;
            }
        }
        public void UpdateAdd(pVector b)
        {
            for (int i = 0; i < Size; i++)
            {
                values[i] += b[i];
            }
        }

        public (pVector, pVector) Split(int firstSize)
        {
            if (firstSize > 0 && Size > firstSize)
            {
                var a = pVector.NewWithLen(firstSize);
                a.InitFrom(this);
                var b = pVector.NewWithLen(Size - firstSize);
                b.InitFrom(this, offset: firstSize);
                return (a, b);
            }
            else
            {
                throw new Exception("Split Position Out of Range");
            }
        }
        public pVector(int length, IEnumerable<double> initValues)
        {
            this.values = new List<double>(length);
            if (initValues != null)
            {
                int index = 0;
                foreach (var item in initValues)
                {
                    if (index < length)
                    {
                        values.Add(item);
                    }
                    else
                    {
                        break;
                    }
                    index += 1;
                }
            }
        }
        public override string ToString()
        {
            return values.Count == 1 ? values[0] + "" : $"[{string.Join(",", values)}]";
        }

        public string VectorString5 => AsVectorString(5);

        public string AsVectorString(int decimalNum = -1)
        {
            return $"[{string.Join(",", values.Select(k => decimalNum > 0 ? Math.Round(k, decimalNum) : k))}]";
        }

        public double Dot(pVector B)
        {
            return Dot(values, B.values);
        }

        public double DotWithItself()
        {
            return Dot(this);
        }

        public pMatrix Transpose()
        {
            return new pMatrix(1, Size, values);
        }

        //public pMatrix AsMatrix()
        //{
        //    return new pMatrix(Length, 1, values);
        //}
        public pMatrix ToUpperTriangleMatrix(int size)
        {
            var sizeNeeded = size * (size + 1) / 2;
            if (sizeNeeded == Size)
            {
                return pMatrix.UpTriangleExtractFrom(this, size, size);
            }
            else
            {
                throw new Exception($"Not Able Convert To Upper Triangular Matrix({size} x {size}) from Vector({Size})");
            }
        }
        public pMatrix ToMatrix(int rows, int columns)
        {
            if (rows * columns == Size)
            {
                return new pMatrix(rows, columns, this);
            }
            else
            {
                throw new Exception($"Not Able Convert To Matrix({rows} x {columns}) from Vector({Size})");
            }
        }

        public static double Dot(List<double> A, List<double> B)
        {
            if (A.Count == B.Count)
            {
                return A.Zip(B, (a, b) => a * b).Sum();
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
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

        public pVector ExcludeAtIndex(int excludeIndex)
        {
            var value = pVector.NewWithLen(Size - 1);
            int tj = 0;
            for (int i = 0; i < Size; i++)
            {
                if (i != excludeIndex)
                {
                    value[tj++] = this[i];
                }
            }
            return value;
        }

        public static double Dot(double[] A, double[] B)
        {
            if (A.Length == B.Length)
            {
                return A.Zip(B, (a, b) => a * b).Sum();
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
        }
        public double TwoNorm()
        {
            return Math.Sqrt(SquaredSum());
        }
        public double SquaredSum()
        {
            return values.Select(k => k * k).Sum();
        }

        public pVector Normalized()
        {
            var twoNorm = TwoNorm();
            if (Math.Abs(twoNorm) < 1e-10)
            {
                return pVector.Constant(this.Size, 1) / Math.Sqrt(this.Size);
            }
            else
            {
                return new pVector(values.Select(k => k / twoNorm).ToArray());
            }
        }

        public void NormalizeItSelf()
        {
            var twoNorm = TwoNorm();
            if (Math.Abs(twoNorm) < 1e-10)
            {
                double constant = 1 / Math.Sqrt(this.Size);
                for (int i = 0; i < this.Size; i++)
                {
                    this[i] = constant;
                }
            }
            else
            {
                for (int i = 0; i < this.Size; i++)
                {
                    this[i] /= twoNorm;
                }
            }
        }

        public pMatrix OuterProduct(pVector b)
        {
            var final = new pMatrix(this.Size, b.Size);
            for (int m = 0; m < this.Size; m++)
            {
                for (int n = 0; n < b.Size; n++)
                {
                    final[m, n] = this[m] * b[n];
                }
            }
            return final;
        }

        public static pVector Zero(int length)
        {
            return Constant(length, 0);
        }

        public static pVector Random(int length, Random random)
        {
            return new pVector(Enumerable.Range(0, length).Select(k => random.NextDouble()).ToArray());
        }

        public static pVector Constant(int length, double value)
        {
            return new pVector(Enumerable.Range(0, length).Select(m => value).ToArray());
        }

        public static pVector OneHot(int index, int length, double value = 1)
        {
            return new pVector(Enumerable.Range(0, length).Select(m => m == index ? value : 0.0).ToArray());
        }

        public static pVector ExtractFrom(pVector pVector, int size, int offset = 0)
        {
            return new pVector(pVector.Skip(offset).Take(size));
        }
        public static pVector From(params double[] values)
        {
            return new pVector(values);
        }

        public IEnumerator<double> GetEnumerator()
        {
            foreach (var item in values)
            {
                yield return item;
            }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public pVector Sqrt()
        {
            return pVector.From(this.Select(m => Math.Sqrt(m)).ToArray());
        }
        public static pVector operator |(pVector A, pVector B)
        {
            return pVector.Combine(A, B);
        }

        public static pVector operator +(double value, pVector vector)
        {
            return new pVector(vector.values.Select(m => m + value).ToArray());
        }
        public static pVector operator *(pVector vector, double value)
        {
            return new pVector(vector.values.Select(m => m * value).ToArray());
        }
        public static pVector operator *(double value, pVector vector)
        {
            return vector * value;
        }
        public static implicit operator pVector(double[] values)
        {
            return new pVector(values);
        }
        public static implicit operator double[](pVector vector)
        {
            return vector.values.ToArray();
        }
        public static pVector operator /(pVector A, pVector B)
        {
            if (A.Size == B.Size)
            {
                return new pVector(Enumerable.Range(0, A.Size).Select(m => A[m] / B[m]).ToArray());
            }
            else
            {
                throw new Exception("Vector Length Not Matched");
            }
        }
        public static pVector operator *(pVector A, pVector B)
        {
            if (A.Size == B.Size)
            {
                return new pVector(Enumerable.Range(0, A.Size).Select(m => A[m] * B[m]).ToArray());
            }
            else
            {
                throw new Exception("Vector Length Not Matched");
            }
        }
        public static pVector operator +(pVector vector, double value)
        {
            return new pVector(vector.values.Select(m => m + value).ToArray());
        }
        public static pVector operator -(pVector vector)
        {
            return new pVector(vector.values.Select(m => -m).ToArray());
        }
        public static pVector operator -(pVector vector, double value)
        {
            return new pVector(vector.values.Select(m => m - value).ToArray());
        }
        public static pVector operator /(pVector vector, double value)
        {
            return new pVector(vector.values.Select(m => m / value).ToArray());
        }
        public static pVector operator +(pVector A, pVector B)
        {
            if (A.values.Count == B.values.Count)
            {
                return new pVector(A.values.Zip(B.values, (a, b) => a + b).ToArray());
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
        }

        public static pVector operator -(pVector A, pVector B)
        {
            if (A.values.Count == B.values.Count)
            {
                return new pVector(A.values.Zip(B.values, (a, b) => a - b).ToArray());
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
        }

        public static pVector operator /(double A, pVector B)
        {
            return new pVector(B.values.Select(k => A / k).ToArray());
        }

        public int MaxIndex()
        {
            int mI = 0;
            double max = double.MinValue;
            for (int i = 0; i < values.Count; i++)
            {
                if (values[i] > max)
                {
                    mI = i;
                    max = values[i];
                }
            }
            return mI;
        }

        public static pVector NewWithLen(int size)
        {
            return new pVector(size, Enumerable.Range(0, size).Select(m => 0.0));
        }

        public pVector pVectorZeroMark(params bool[] marskArray)
        {
            return this.values.Select((v, i) => marskArray.Length > i ? (marskArray[i] ? 0.0 : v) : v).ToArray();
        }
        public void AssertApproximateEqual(pVector other, double error = 1e-8)
        {
            if (!ApproximateEquals(other, error))
            {
                throw new Exception($"Vector({this.VectorString5}) != Vector({other.VectorString5})");
            }
        }

        public bool ApproximateEquals(pVector other, double error = 1e-8)
        {
            if (Size == other.Size)
            {
                return values.Zip(other.values, (a, b) => (a, b)).All(k => Math.Abs(k.a - k.b) < error);
            }
            else
            {
                return false;
            }
        }

        public pVector Clone()
        {
            return new pVector(this.values);
        }

        public pMatrix KroneckerProduct(pMatrix pMatrix)
        {
            return pMatrix.KroneckerProduct(pMatrix.NewMatrixFromColumns(this), pMatrix);
        }

        public pVector CrossProduct(pVector B)
        {
            if (this.Size == B.Size && B.Size == 3)
            {
                //| i   j   k  |
                //| ax  ay  az |
                //| bx  by  bz |
                //i * (ay* bz - az * by) - j*(ax * bz - az* bx) + k*(ax* by - ay*bx)
                var tmp = pVector.NewWithLen(3);
                tmp[0] = this[1] * B[2] - this[2] * B[1];
                tmp[1] = this[2] * B[0] - this[0] * B[2];
                tmp[2] = this[0] * B[1] - this[1] * B[0];
                return tmp;

            }
            else
            {
                throw new Exception($"Cross Product Only Define for Vector Size 3 X 3, But {this.Size} x {B.Size} is Provided");
            }
        }

        public pVector SubVector(int startIndex, int length)
        {
            var aa = pVector.NewWithLen(length);
            for (int i = 0; i < length; i++)
            {
                aa[i] = this[i + startIndex];
            }
            return aa;
        }
        public static pVector Sum(IEnumerable<pVector> vs)
        {
            if (vs.Any())
            {
                return vs.Aggregate((a, b) => a + b);
            }
            else
            {
                return null;
            }
        }
    }
}
