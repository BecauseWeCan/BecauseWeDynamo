using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Topology;
using Text;

namespace Panelization
{
    public class EdgeConnector
    {
        //FIELDS
        public Edge Edge { get; private set; }

        internal EdgeConnector(Edge Edge, double MinInset, int ConnectorCount)
        {
            //initialize
            this.Edge = Edge;

        }

    }
}
