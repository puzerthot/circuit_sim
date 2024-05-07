using pLib;
using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;

namespace circuit_sim
{
    class NetList
    {
        class Branch
        {
            public enum BranchTypeEnum { Current, Resistor, Voltage }
            public BranchTypeEnum BranchType;
            public string ID;
            public int FromNode;
            public int ToNode;
            public double Value;
            public static Branch FromLine(string[] headers, string line)
            {
                var values = line.Split(new string[] { "," }, StringSplitOptions.None);
                var branch = new Branch();
                branch.ID = values[0];
                branch.FromNode = int.Parse(values[1]);
                branch.ToNode = int.Parse(values[2]);
                branch.Value = double.Parse(values[3]);
                if (branch.ID.StartsWith("I", StringComparison.OrdinalIgnoreCase))
                {
                    branch.BranchType = BranchTypeEnum.Current;
                }
                else if (branch.ID.StartsWith("R", StringComparison.OrdinalIgnoreCase))
                {
                    branch.BranchType = BranchTypeEnum.Resistor;
                }
                else if (branch.ID.StartsWith("V", StringComparison.OrdinalIgnoreCase))
                {
                    branch.BranchType = BranchTypeEnum.Resistor;
                }
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

        public Dictionary<int, double> CaculateNodeVolatage()
        {
            var Y = new pLib.pMatrix(NodeIDDict.Count);
            var J = pLib.pVector.NewWithLen(NodeIDDict.Count);
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
