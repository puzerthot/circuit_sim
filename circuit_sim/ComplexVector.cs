using System.Collections;
using System;
using System.Collections.Generic;
using System.Linq;

namespace circuit_sim
{
    public class ComplexVector : IEnumerable<ComplexNumber>
    {
        List<ComplexNumber> values;
        public int Size => values.Count;
        public IEnumerable<ComplexNumber> Values => values;
        public ComplexVector(IEnumerable<ComplexNumber> values)
        {
            this.values = values.ToList();
        }
        public ComplexVector(IEnumerable<double> values)
        {
            this.values = values.Select(k => (ComplexNumber)k).ToList();
        }

        public ComplexVector Expand(params ComplexNumber[] new_values)
        {
            var tmp = ComplexVector.NewWithLen(this.Size + new_values.Length);
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

        public ComplexVector ResetValues(double value = 0)
        {
            for (int i = 0; i < Size; i++)
            {
                values[i] = value;
            }
            return this;
        }

        private ComplexVector() { values = new List<ComplexNumber>(); }
        public static ComplexVector Combine(ComplexVector A, params ComplexVector[] Rests)
        {
            var tmp = new ComplexVector();
            tmp.values.AddRange(A);
            foreach (var v in Rests)
            {
                tmp.values.AddRange(v);
            }
            return tmp;
        }
        public void ExpandWithValue(ComplexNumber value)
        {
            this.values.Add(value);
        }
        public void InitFrom(ComplexVector b, int offset = 0)
        {
            for (int i = 0; i < Math.Min(this.Size, b.Size); i++)
            {
                this.values[i] = b[i + offset];
            }
        }
        public void CopyFrom(ComplexVector b, int offset = 0)
        {
            for (int i = 0; i < b.Size; i++)
            {
                this.values[i + offset] = b[i];
            }
        }
        public void AddFrom(ComplexVector b, int offset = 0)
        {
            for (int i = 0; i < b.Size; i++)
            {
                this.values[i + offset] += b[i];
            }
        }
        public void CopyFrom(ComplexMatrix b, int offset = 0)
        {
            int i = offset;
            foreach (var v in b.Values)
            {
                this.values[offset] = v;
                offset++;
            }
        }
        public void UpdateAdd(ComplexVector b)
        {
            for (int i = 0; i < Size; i++)
            {
                values[i] += b[i];
            }
        }

        public (ComplexVector, ComplexVector) Split(int firstSize)
        {
            if (firstSize > 0 && Size > firstSize)
            {
                var a = ComplexVector.NewWithLen(firstSize);
                a.InitFrom(this);
                var b = ComplexVector.NewWithLen(Size - firstSize);
                b.InitFrom(this, offset: firstSize);
                return (a, b);
            }
            else
            {
                throw new Exception("Split Position Out of Range");
            }
        }
        public ComplexVector(int length, IEnumerable<double> initValues) : this(length, initValues.Select(k => (ComplexNumber)k))
        {
        }
        public ComplexVector(int length, IEnumerable<ComplexNumber> initValues)
        {
            this.values = new List<ComplexNumber>(length);
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
            return $"[{string.Join(",", values.Select(k => decimalNum > 0 ? k.Round(decimalNum) : k))}]";
        }

        public ComplexNumber Dot(ComplexVector B)
        {
            return Dot(values, B.values);
        }

        public ComplexNumber DotWithItself()
        {
            return Dot(this);
        }

        public ComplexMatrix Transpose()
        {
            return new ComplexMatrix(1, Size, values);
        }

        //public ComplexMatrix AsMatrix()
        //{
        //    return new ComplexMatrix(Length, 1, values);
        //}

        public ComplexMatrix ToMatrix(int rows, int columns)
        {
            if (rows * columns == Size)
            {
                return new ComplexMatrix(rows, columns, this);
            }
            else
            {
                throw new Exception($"Not Able Convert To Matrix({rows} x {columns}) from Vector({Size})");
            }
        }

        public static ComplexNumber Dot(List<ComplexNumber> A, List<ComplexNumber> B)
        {
            if (A.Count == B.Count)
            {
                return A.Zip(B, (a, b) => a * b).Aggregate((a, b) => a + b);
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
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

        public ComplexVector ExcludeAtIndex(int excludeIndex)
        {
            var value = ComplexVector.NewWithLen(Size - 1);
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

        public static ComplexNumber Dot(ComplexNumber[] A, ComplexNumber[] B)
        {
            if (A.Length == B.Length)
            {
                return A.Zip(B, (a, b) => a * b).Aggregate((a, b) => a + b);
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
        }
        public ComplexNumber SquaredSum()
        {
            return values.Select(k => k * k).Aggregate((a, b) => a + b);
        }

        public ComplexMatrix OuterProduct(ComplexVector b)
        {
            var final = new ComplexMatrix(this.Size, b.Size);
            for (int m = 0; m < this.Size; m++)
            {
                for (int n = 0; n < b.Size; n++)
                {
                    final[m, n] = this[m] * b[n];
                }
            }
            return final;
        }

        public static ComplexVector Zero(int length)
        {
            return Constant(length, 0);
        }

        public static ComplexVector Random(int length, Random random)
        {
            return new ComplexVector(Enumerable.Range(0, length).Select(k => random.NextDouble()).ToArray());
        }

        public static ComplexVector Constant(int length, double value)
        {
            return new ComplexVector(Enumerable.Range(0, length).Select(m => value).ToArray());
        }

        public static ComplexVector OneHot(int index, int length, double value = 1)
        {
            return new ComplexVector(Enumerable.Range(0, length).Select(m => m == index ? value : 0.0).ToArray());
        }

        public static ComplexVector ExtractFrom(ComplexVector ComplexVector, int size, int offset = 0)
        {
            return new ComplexVector(ComplexVector.Skip(offset).Take(size));
        }
        public static ComplexVector From(params ComplexNumber[] values)
        {
            return new ComplexVector(values);
        }

        public IEnumerator<ComplexNumber> GetEnumerator()
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

        public static ComplexVector operator |(ComplexVector A, ComplexVector B)
        {
            return ComplexVector.Combine(A, B);
        }

        public static ComplexVector operator +(double value, ComplexVector vector)
        {
            return new ComplexVector(vector.values.Select(m => m + value).ToArray());
        }
        public static ComplexVector operator *(ComplexVector vector, double value)
        {
            return new ComplexVector(vector.values.Select(m => m * value).ToArray());
        }
        public static ComplexVector operator *(double value, ComplexVector vector)
        {
            return vector * value;
        }
        public static implicit operator ComplexVector(double[] values)
        {
            return new ComplexVector(values.Select(k => (ComplexNumber)k));
        }
        public static implicit operator ComplexVector(ComplexNumber[] values)
        {
            return new ComplexVector(values);
        }
        public static implicit operator ComplexNumber[](ComplexVector vector)
        {
            return vector.values.ToArray();
        }
        public static ComplexVector operator *(ComplexVector A, ComplexVector B)
        {
            if (A.Size == B.Size)
            {
                return new ComplexVector(Enumerable.Range(0, A.Size).Select(m => A[m] * B[m]).ToArray());
            }
            else
            {
                throw new Exception("Vector Length Not Matched");
            }
        }
        public static ComplexVector operator +(ComplexVector vector, double value)
        {
            return new ComplexVector(vector.values.Select(m => m + value).ToArray());
        }
        public static ComplexVector operator -(ComplexVector vector)
        {
            return new ComplexVector(vector.values.Select(m => -m).ToArray());
        }
        public static ComplexVector operator -(ComplexVector vector, double value)
        {
            return new ComplexVector(vector.values.Select(m => m - value).ToArray());
        }
        public static ComplexVector operator /(ComplexVector vector, double value)
        {
            return new ComplexVector(vector.values.Select(m => m / value).ToArray());
        }
        public static ComplexVector operator +(ComplexVector A, ComplexVector B)
        {
            if (A.values.Count == B.values.Count)
            {
                return new ComplexVector(A.values.Zip(B.values, (a, b) => a + b).ToArray());
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
        }

        public static ComplexVector operator -(ComplexVector A, ComplexVector B)
        {
            if (A.values.Count == B.values.Count)
            {
                return new ComplexVector(A.values.Zip(B.values, (a, b) => a - b).ToArray());
            }
            else
            {
                throw new Exception("Size Not Equal");
            }
        }

        public static ComplexVector operator /(double A, ComplexVector B)
        {
            return new ComplexVector(B.values.Select(k => A / k).ToArray());
        }


        public static ComplexVector NewWithLen(int size)
        {
            return new ComplexVector(size, Enumerable.Range(0, size).Select(m => 0.0));
        }

        public ComplexVector ComplexVectorZeroMark(params bool[] marskArray)
        {
            return this.values.Select((v, i) => marskArray.Length > i ? (marskArray[i] ? 0.0 : v) : v).ToArray();
        }
        public void AssertApproximateEqual(ComplexVector other, double error = 1e-8)
        {
            if (!ApproximateEquals(other, error))
            {
                throw new Exception($"Vector({this.VectorString5}) != Vector({other.VectorString5})");
            }
        }

        public bool ApproximateEquals(ComplexVector other, double error = 1e-8)
        {
            if (Size == other.Size)
            {
                return values.Zip(other.values, (a, b) => (a, b)).All(k => (k.a - k.b).Absolute < error);
            }
            else
            {
                return false;
            }
        }

        public ComplexVector Clone()
        {
            return new ComplexVector(this.values);
        }

        public ComplexMatrix KroneckerProduct(ComplexMatrix ComplexMatrix)
        {
            return ComplexMatrix.KroneckerProduct(ComplexMatrix.NewMatrixFromColumns(this), ComplexMatrix);
        }

        public ComplexVector CrossProduct(ComplexVector B)
        {
            if (this.Size == B.Size && B.Size == 3)
            {
                //| i   j   k  |
                //| ax  ay  az |
                //| bx  by  bz |
                //i * (ay* bz - az * by) - j*(ax * bz - az* bx) + k*(ax* by - ay*bx)
                var tmp = ComplexVector.NewWithLen(3);
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

        public ComplexVector SubVector(int startIndex, int length)
        {
            var aa = ComplexVector.NewWithLen(length);
            for (int i = 0; i < length; i++)
            {
                aa[i] = this[i + startIndex];
            }
            return aa;
        }
        public static ComplexVector Sum(IEnumerable<ComplexVector> vs)
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
