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
            var tmp = new NetList("C:\\Users\\paulcheuk\\Desktop\\Playground\\circuit_sim\\circuit_sim\\bin\\Debug\\net1.txt");
            var config = new NetList.SimulationConfig();
            var c = tmp.Simulation(config);

            //var tmp = new NetList("C:\\Users\\paulcheuk\\Desktop\\Playground\\circuit_sim\\circuit_sim\\bin\\Debug\\net2.txt");
            //var config = new NetList.SimulationConfig();
            //config.SetTemperatureInDegree(27);
            //var c = tmp.Simulation(config);
        }

    }
}
