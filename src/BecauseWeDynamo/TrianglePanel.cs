using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using TriangleRigging;

namespace Fabrication
{
    class TrianglePanel: List<Point>
    {
        public double[] Angle { get; set; }
        public double[] Radius { get; set; } 

        internal void TrianglePanel(IEnumerable<Point> Points, double EdgeOffset, double CornerOffset, double Thickness, int direction=0)
        {
            Triangle t = Triangle.BySurfaceAndPoints()
            this.AddRange(Points);

        }

    }



}
