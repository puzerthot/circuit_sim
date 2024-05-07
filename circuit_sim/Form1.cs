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
            var tmp = new NetList("C:\\Users\\paulcheuk\\Desktop\\Playground\\circuit_sim\\circuit_sim\\bin\\Debug\\text.txt");

            var c = tmp.CaculateNodeVolatage();
        }

    }
}
