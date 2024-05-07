using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;


namespace circuit_sim
{
    class NetList
    {
        class DiodeBranch : Branch
        {
            public const double q = 1.60217663e-19;// charge on an electron: 1.60217663 × 10-19 coulombs
            public const double k = 1.380649e-23; //Boltzman's constant: 1.380649 × 10-23 m2 kg s-2 K-1
            public double T; //absolute temperature in degree Kelvin
            public double Eta; //non-ideality factor
            public double Isat; //reverse saturation current;
            public DiodeBranch(string id) : base(id)
            {
            }
            public double CurrentGivenVd(double vd)
            {
                // ref: circuit and system simulation methods.pdf, page 24
                var current = Isat * (Math.Exp(q * vd / (Eta * k * T)) - 1);
                return current;
            }
            public override void SetValue(string valueString)
            {
                base.SetValue(valueString);
            }
        }
        class Branch
        {
            public Branch(string id)
            {
                ID = id;
            }
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
        public Dictionary<int, double> Simulation(SimulationConfig config)
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
            }

            //to definiate admitance, ref: circuit and system simulation methods.pdf,  page 21
            int selectGroudNodeIndex = 0;
            var rY = new pMatrix(Y.Rows - 1, Y.Columns - 1);
            foreach (var (r, c) in Y.Indexes())
            {
                if (r != selectGroudNodeIndex && c != selectGroudNodeIndex)
                {
                    var tr = r > selectGroudNodeIndex ? r - 1 : r;
                    var tc = c > selectGroudNodeIndex ? c - 1 : c;
                    rY[tr, tc] = Y[r, c];
                }
            }

            var rJ = pVector.NewWithLen(J.Size - 1);
            for (int i = 0; i < J.Size; i++)
            {
                var tj = i > selectGroudNodeIndex ? i - 1 : i;
                if (i != selectGroudNodeIndex)
                {
                    rJ[tj] = J[i];
                }
            }

            var rV = pMatrix.SolveViaLUDecomposition(rY, rJ);

            var test = rY * rV;

            var v = pVector.NewWithLen(J.Size);
            for (int i = 0; i < J.Size; i++)
            {
                if (i == selectGroudNodeIndex)
                {
                    v[i] = 0;
                }
                else
                {
                    var tj = i > selectGroudNodeIndex ? i - 1 : i;
                    v[i] = rV[tj];
                }
            }

            var tmp = new Dictionary<int, double>();
            for (int i = 0; i < v.Size; i++)
            {
                var value = v[i];
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
