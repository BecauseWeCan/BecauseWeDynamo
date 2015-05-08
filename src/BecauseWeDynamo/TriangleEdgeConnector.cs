using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Fabrication
{
    public class TriangleEdgeConnector
    {
        double radius;
        double spacing;
        public int[] numHoles { get; set; }
        public TriangleEdge edge { get; set; }

        internal TriangleEdgeConnector(TriangleEdge tEdge, double radiusHole, double spacingHole)
        {
            spacing = spacingHole;
            radius = radiusHole;
            edge = tEdge;
            numHoles = new int[3] { 0, 0, 0 };
            for (int i = 0; i < edge.triangles.Count; i++)
            {
                numHoles[i] = (int)(Math.Ceiling((edge.triangles[i].surface.ClosestPointTo(edge.midpoint).DistanceTo(edge.midpoint) + 3 * radius) / spacing)) + edge.triangles[i].numHoles - 1;
            }
            numHoles[2] = numHoles[0] + numHoles[1] + 1;
        }

        public static TriangleEdgeConnector ByEdge(TriangleEdge edge, double radiusHole, double spacingHole) { return new TriangleEdgeConnector(edge, radiusHole, spacingHole); }

        public List<Circle> PlaceHoles(double DisplayFactor = 1)
        {
            if (numHoles[2] < 2 || edge.isOuterEdge) return null;
            List<Circle> result = new List<Circle>(numHoles[2]);

            for (int i = 0; i < edge.triangles.Count; i++)
            {
                for (int k = 0; k < numHoles[i]; k++)
                {
                    result.Add(Circle.ByCenterPointRadiusNormal(
                        (Point)edge.midpoint.Translate(Vector.ByTwoPoints(edge.midpoint, edge.triangles[i].surface.ClosestPointTo(edge.midpoint)), spacing * (k + 1)),
                        DisplayFactor * radius,
                        edge.triangles[i].Normal
                        ));
                }
            }

            result.Add(Circle.ByCenterPointRadiusNormal(
                            (Point)edge.midpoint,
                            DisplayFactor * radius,
                            edge.Normal
                            ));
            return result;
        }
        public List<Circle> PlaceHolesByHoleCount(int numHoles, double DisplayFactor = 1)
        {
            if (this.numHoles[2] == numHoles) return PlaceHoles(DisplayFactor);
            return null;
        }
    }
}
