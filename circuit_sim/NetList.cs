using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;
using System.Windows.Forms;


namespace circuit_sim
{
    class NetList
    {
        public class DiodeBranch : Branch
        {
            public override bool IsNonlinear => true;
            public const double q = -1.60217663e-19;// charge on an electron: 1.60217663 × 10-19 coulombs
            public const double k = 1.380649e-23; //Boltzman's constant: 1.380649 × 10-23 m2 kg s-2 K-1
            public double T; //absolute temperature in degree Kelvin
            public double Eta; //non-ideality factor
            public double Isat; //reverse saturation current;

            public double finalC;
            public bool finalCLock = false;
            public DiodeBranch(string id) : base(id)
            {
            }
            public double CurrentGivenVd(double vd)
            {
                // ref: circuit and system simulation methods.pdf, page 24
                var current = Isat * (Math.Exp(finalC * vd) - 1);
                return current;
            }
            public double CurrentDerivativeGivenVd(double vd)
            {
                return Isat * Math.Exp(finalC * vd) * finalC;
            }
            public override void PrepareFor(SimulationConfig simulationConfig)
            {
                T = simulationConfig.TemperatureInDegreeKelvin;
                if (!finalCLock)
                {
                    finalC = q / (Eta * k * T);
                }
            }

            public override (double Ieq, double Geq) LinearizedModel(double vd)
            {
                var id = CurrentGivenVd(vd);
                var Geq = CurrentDerivativeGivenVd(vd);
                var Ieq = id - vd * Geq;
                return (Ieq, Geq);
            }
            public override void SetValue(string valueString)
            {
                var values = valueString.Split(';');
                foreach (var value in values)
                {
                    if (!string.IsNullOrEmpty(value.Trim()))
                    {
                        var items = value.Split(':');
                        var key = items[0];
                        var rkey = key.ToLower().Trim();
                        switch (rkey)
                        {
                            case "isat":
                                Isat = double.Parse(items[1]);
                                break;
                            case "eta":
                                Eta = double.Parse(items[1]);
                                break;
                            case "finalc":
                                finalC = double.Parse(items[1]);
                                finalCLock = true;
                                break;
                            default:
                                throw new InvalidDataException($"Diode Dont have a value for key: {key}");
                        }
                    }
                }
            }
        }
        public class Branch
        {
            public virtual bool IsNonlinear => false;
            public Branch(string id)
            {
                ID = id;
            }
            public virtual void PrepareFor(SimulationConfig simulationConfig) { }
            public enum BranchTypeEnum { Current, Resistor, Voltage, Diode }
            public BranchTypeEnum BranchType;
            public string ID { get; private set; }
            public int FromNode;
            public int ToNode;
            public double Value;
            public virtual void SetValue(string valueString)
            {
                Value = double.Parse(valueString);
            }
            public virtual (double Ieq, double Geq) LinearizedModel(double Vd)
            {
                throw new NotSupportedException($"Branch({BranchType}) Not Support Linearization");
            }
            public static Branch NewBranchOf(string nodeId)
            {
                switch (nodeId[0].ToString().ToUpper())
                {
                    case "I": return new Branch(nodeId) { BranchType = BranchTypeEnum.Current };
                    case "R": return new Branch(nodeId) { BranchType = BranchTypeEnum.Resistor };
                    case "V": return new Branch(nodeId) { BranchType = BranchTypeEnum.Voltage };
                    case "D": return new DiodeBranch(nodeId) { BranchType = BranchTypeEnum.Diode };
                    default:
                        throw new NotSupportedException("Not Support Node:" + nodeId);

                }
            }
            public static Branch FromLine(string[] headers, string line)
            {
                var values = line.Split(new string[] { "," }, StringSplitOptions.None);
                var branch = NewBranchOf(values[0]);
                branch.FromNode = int.Parse(values[1]);
                branch.ToNode = int.Parse(values[2]);
                branch.SetValue(values[3]);
                return branch;
            }

            public override string ToString()
            {
                return $"{ID}: {Value}, {FromNode} -> {ToNode}";
            }
        }
        List<Branch> Branches;
        Dictionary<int, NodeInfo> NodeIDDict = new Dictionary<int, NodeInfo>();
        class NodeInfo
        {
            public int Index;
        }
        public NetList(string filePath)
        {
            string[] headers = null;
            List<Branch> branches = new List<Branch>();
            foreach (var line in File.ReadAllLines(filePath))
            {
                if (headers == null)
                {
                    headers = line.Split(new string[] { "," }, StringSplitOptions.None);
                }
                else
                {
                    branches.Add(Branch.FromLine(headers, line));
                }
            }
            Branches = branches;

            for (int i = 0; i < branches.Count; i++)
            {
                var branch = branches[i];
                BindNodeAndBranch(branch.FromNode);
                BindNodeAndBranch(branch.ToNode);
            }
        }
        public class SimulationConfig
        {
            public double TemperatureInDegreeKelvin;
            public void SetTemperatureInDegree(double degree)
            {
                TemperatureInDegreeKelvin = 273.15 + degree;
            }
        }

        public (pMatrix A, pVector J) Full(pVector Vi)
        {
            var Y = new pMatrix(NodeIDDict.Count);
            var J = pVector.NewWithLen(NodeIDDict.Count);
            foreach (var b in Branches)
            {
                //ref: circuit and system simulation methods.pdf,  page 19, 
                var i = NodeIDDict[b.FromNode].Index;
                var j = NodeIDDict[b.ToNode].Index;
                if (b.BranchType == Branch.BranchTypeEnum.Resistor)
                {
                    var conductance = 1 / b.Value;
                    Y[i, i] += conductance;
                    Y[i, j] -= conductance;
                    Y[j, i] -= conductance;
                    Y[j, j] += conductance;
                }
                else if (b.BranchType == Branch.BranchTypeEnum.Current)
                {
                    J[i] -= b.Value;
                    J[j] += b.Value;
                }
                else if (b.IsNonlinear)
                {
                    var (Ieq, Geq) = b.LinearizedModel(Vi[i] - Vi[j]);
                    if (Ieq != 0.0)
                    {
                        J[i] -= Ieq;
                        J[j] += Ieq;
                    }

                    if (Geq != 0.0)
                    {
                        Y[i, i] += Geq;
                        Y[i, j] -= Geq;
                        Y[j, i] -= Geq;
                        Y[j, j] += Geq;
                    }
                }
            }
            return (Y, J);
        }

        public (pMatrix A, pVector J) Partial(pVector Vi, int excludedIndex)
        {
            //to definiate admitance, ref: circuit and system simulation methods.pdf,  page 21
            var Y = new pMatrix(NodeIDDict.Count);
            var J = pVector.NewWithLen(NodeIDDict.Count);
            foreach (var b in Branches)
            {
                //ref: circuit and system simulation methods.pdf,  page 19, 
                var i = NodeIDDict[b.FromNode].Index;
                var j = NodeIDDict[b.ToNode].Index;
                if (b.BranchType == Branch.BranchTypeEnum.Resistor)
                {
                    var conductance = 1 / b.Value;
                    Y[i, i] += conductance;
                    Y[i, j] -= conductance;
                    Y[j, i] -= conductance;
                    Y[j, j] += conductance;
                }
                else if (b.BranchType == Branch.BranchTypeEnum.Current)
                {
                    J[i] -= b.Value;
                    J[j] += b.Value;
                }
                else if (b.IsNonlinear)
                {
                    var (Ieq, Geq) = b.LinearizedModel(Vi[i] - Vi[j]);
                    if (Ieq != 0.0)
                    {
                        J[i] -= Ieq;
                        J[j] += Ieq;
                    }

                    if (Geq != 0.0)
                    {
                        Y[i, i] += Geq;
                        Y[i, j] -= Geq;
                        Y[j, i] -= Geq;
                        Y[j, j] += Geq;
                    }
                }
            }
            var rY = Y.MatrixExcludeRC(excludedIndex, excludedIndex);
            var rJ = J.ExcludeAtIndex(excludedIndex);

            return (rY, rJ);
        }
        public Dictionary<int, double> Simulation(SimulationConfig config)
        {
            foreach (var b in Branches)
            {
                b.PrepareFor(config);
            }
            int selectGroudNodeIndex = 0;
            bool isAnyNonlinear = Branches.Any(k => k.IsNonlinear);
            bool isConverge = false;
            var Vi = pVector.NewWithLen(NodeIDDict.Count);
            var random = new Random(0);
            do
            {
                var (rY, rJ) = Partial(Vi, selectGroudNodeIndex);
                var rV = pMatrix.SolveViaLUDecomposition(rY, rJ);
                double totalDiff = 0.0;
                for (int i = 0; i < Vi.Size; i++)
                {
                    if (i == selectGroudNodeIndex)
                    {
                        Vi[i] = 0;
                    }
                    else
                    {
                        var tj = i > selectGroudNodeIndex ? i - 1 : i;
                        totalDiff += Math.Pow(Vi[i] - rV[tj], 2);
                        Vi[i] = rV[tj];

                    }
                }

                isConverge = Math.Sqrt(totalDiff) < 1e-8;

            } while (isAnyNonlinear && !isConverge);

            var tmp = new Dictionary<int, double>();
            for (int i = 0; i < Vi.Size; i++)
            {
                var value = Vi[i];
                var nodeId = NodeIDDict.Where(k => k.Value.Index == i).Single().Value.Index;
                tmp[nodeId] = value;
            }
            return tmp;
        }

        private void BindNodeAndBranch(int nodeId)
        {
            if (!NodeIDDict.ContainsKey(nodeId))
            {
                NodeIDDict[nodeId] = new NodeInfo() { Index = NodeIDDict.Count() };
            }
        }
    }
}
