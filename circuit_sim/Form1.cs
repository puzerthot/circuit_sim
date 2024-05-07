using System;
using System.ComponentModel;
using System.Drawing;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Markup;

namespace circuit_sim
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();

            Work();
        }

        private void Work()
        {
            //var tmp = new NetList("C:\\Users\\paulcheuk\\Desktop\\Playground\\circuit_sim\\circuit_sim\\bin\\Debug\\net1.txt");
            //var config = new NetList.SimulationConfig();
            //var c = tmp.Simulation(config);

            var tmp = new NetList("C:\\Users\\paulcheuk\\Desktop\\Playground\\circuit_sim\\circuit_sim\\bin\\Debug\\net2.txt");
            var config = new NetList.SimulationConfig();
            config.SetTemperatureInDegree(27);
            var c = tmp.Simulation(config);


            //var dio = new NetList.DiodeBranch("D3");
            //var config = new NetList.SimulationConfig();
            //config.SetTemperatureInDegree(27);
            //dio.PrepareFor(config);
            //dio.SetValue("Eta:0.9999;Isat:1E-14;");
            //dio.finalC = 12.025;


            //double r = 100;
            //var g = 1 / r;
            //var vs = 55;

            //var errorFunc = new Func<double, double>((vn) => vs * g - vn * g - dio.CurrentGivenVd(vn));
            //var dErrorFunc = new Func<double, double>((vn) => -g - dio.CurrentDerivativeGivenVd(vn));

            //var finalValue = NewtonMethod(0, errorFunc, dErrorFunc);

            //Console.WriteLine($"final:{finalValue}, error:{errorFunc(finalValue)}, current:{dio.CurrentGivenVd(finalValue)}");

            foreach(var item in c)
            {
                Console.WriteLine($"{item.Key}:{item.Value}");
            }
        }

        public double NewtonMethod(double init, Func<double, double> error, Func<double, double> derivative, double eplison = 1e-9)
        {
            double ak = init;
            for (int i = 0; i < 10000; i++)
            {
                var errorValue = error(ak);
                var derivativeValue = derivative(ak);
                ak -= errorValue / derivativeValue;

                if (Math.Abs(errorValue) < eplison)
                {
                    break;
                }
            }
            return ak;
        }
    }
}
