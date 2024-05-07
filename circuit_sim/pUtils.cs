using System;
using System.Collections.Generic;

namespace circuit_sim
{
    public class pUtils
    {
        public static bool DoubleApproximateEqual(double a, double b, double relativeError = 1e-6, double absoluteError = 1e-8)
        {
            return !(double.IsNaN(a) || double.IsNaN(b) || (Math.Abs(a - b) > Math.Max(relativeError * Math.Max(Math.Abs(a), Math.Abs(b)), absoluteError)));
        }
        public static IEnumerable<int> LoopFromAToB(int start, int end)
        {
            for (int i = start; i <= end; i++)
            {
                yield return i;
            }
        }
    }
}
