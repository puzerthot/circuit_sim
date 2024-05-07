using System;

namespace circuit_sim
{
    public struct ComplexNumber
    {
        public double Real;
        public double Imag;
        public ComplexNumber(double real, double img)
        {
            Real = real;
            Imag = img;
        }
        public double Absolute => Math.Sqrt(Real * Real + Imag * Imag);

        public static implicit operator ComplexNumber(double value)
        {
            return new ComplexNumber(value, 0);
        }

        public ComplexNumber Round(int round)
        {
            return new ComplexNumber(Math.Round(Real, round), Math.Round(Imag, round));
        }
        public static ComplexNumber operator +(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.Real + b.Real, a.Imag + b.Imag);
        }
        public static ComplexNumber operator -(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.Real - b.Real, a.Imag - b.Imag);
        }
        public static ComplexNumber operator *(ComplexNumber a, ComplexNumber b)
        {
            return new ComplexNumber(a.Real * b.Real - a.Imag * b.Imag, a.Real * b.Imag + b.Real * a.Imag);
        }
        public static ComplexNumber operator /(ComplexNumber a, ComplexNumber b)
        {
            var scalar = 1 / (b.Real * b.Real + b.Imag * b.Imag);
            return new ComplexNumber(scalar * (a.Real * b.Real + a.Imag * b.Imag), scalar * (a.Imag * b.Real - a.Real * b.Imag));
        }
        public static ComplexNumber operator *(ComplexNumber a, double scalar)
        {
            return new ComplexNumber(a.Real * scalar, a.Real * scalar);
        }
        public static ComplexNumber operator -(ComplexNumber a)
        {
            return new ComplexNumber(-a.Real, -a.Real);
        }
        public static ComplexNumber operator /(ComplexNumber a, double scalar)
        {
            return new ComplexNumber(a.Real / scalar, a.Real / scalar);
        }
        public static ComplexNumber operator /(double scalar, ComplexNumber a)
        {
            return new ComplexNumber(scalar / a.Real, scalar / a.Real);
        }
        public override string ToString()
        {
            return $"{Real} + j*{Imag}";
        }

        public string ToString(string format)
        {
            return $"{GetString(format, Real)} + j*{GetString(format, Imag)}";
        }
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
    }
}
