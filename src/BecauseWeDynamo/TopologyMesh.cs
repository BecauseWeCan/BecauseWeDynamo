using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;

namespace Topology
{
    public class TriangleMesh : IDisposable
    {
        //**FIELDS**
        private bool disposed = false;
        internal Dictionary<Point, Vertex> V;
        internal Dictionary<string, Edge> E;
        internal List<Edge> E1;
        internal List<Edge> E2;
        internal List<Edge> E3;
        internal List<Edge> E0;

        //**PROPERTIES**QUERY
        public List<Triangle> Faces { get; set; }
        public List<Edge> Edges { get; set; }
        public List<Edge> EdgesOuter { get { return E1; } }
        public List<Vertex> Vertices { get { return V.Values.ToList(); } }
        public List<Spline> Splines { get; set; }
        public Point[] Points { get; set; }

        internal TriangleMesh(Surface[] Surfaces, Point[] Points)
        {
            // initialize
            Faces = new List<Triangle>(Surfaces.Length);
            V = new Dictionary<Point, Vertex>(Points.Length);
            E = new Dictionary<string, Edge>();
            Edges = new List<Edge>();
            E1 = new List<Edge>();
            E2 = new List<Edge>();
            E3 = new List<Edge>();
            E0 = new List<Edge>();
            Splines = new List<Spline>();

            // store input points
            this.Points = Points;
            // create vertex lookup table from points
            for (int i = 0; i < Points.Length; i++) V.Add(Points[i], new Vertex(Points[i]));
            // create faces from surfaces
            for (int i = 0; i < Surfaces.Length; i++)
            {
                // find face vertices in lookup table
                Autodesk.DesignScript.Geometry.Vertex[] vtx = Surfaces[i].Vertices;
                List<Vertex> v = new List<Vertex>(vtx.Length);
                for (int j = 0; j < vtx.Length; j++)
                {
                    Point pt = vtx[j].PointGeometry;
                    for (int k = 0; k < Points.Length; k++)
                        if (pt.IsAlmostEqualTo(Points[k]))
                        { v.Add(V[Points[k]]); break; }
                    pt.Dispose();
                }
                vtx.ForEach(x => x.Dispose());
                // create face based on vertices
                Triangle t = new Triangle(v);
                Faces.Add(t);
                // create or find edges
                for (int j = 0; j < t.E.Count; j++)
                {
                    bool edgeFound = false;
                    if (Edges.Count > 0)
                    {
                        for (int k = 0; k < Edges.Count; k++)
                        {
                            if (Edges[k].Vertices.Contains(t.E[j].V[0]) && Edges[k].Vertices.Contains(t.E[j].V[1]))
                            {
                                Edges[k].E.Add(t.E[j]);
                                t.E[j].AddEdge(Edges[k]);
                                edgeFound = true;
                                break;
                            }
                        }
                    }
                    if (!edgeFound)
                    {
                        Edge e = new Edge();
                        e.E.Add(t.E[j]);
                        t.E[j].AddEdge(e);
                        Edges.Add(e);
                    }
                }
            }
            for (int i = 0; i < Edges.Count; i++)
            {
                if (Edges[i].E.Count == 1) E1.Add(Edges[i]);
                else if (Edges[i].E.Count == 2)
                {
                    double[] AngleNormal = Edges[i].GetAngleNormal(Edges[i].E.ElementAt(0), Edges[i].E.ElementAt(1));
                    Edges[i].Angle = AngleNormal[0];
                    Edges[i].N = new double[] { AngleNormal[1], AngleNormal[2], AngleNormal[3] };
                    E2.Add(Edges[i]);
                }
                else if (Edges[i].E.Count == 3)
                {
                    Vertex A = Edges[i].Vertices[0];
                    Vertex B = Edges[i].Vertices[1];
                    List<HalfEdge> eA = new List<HalfEdge>(2);
                    List<HalfEdge> eB = new List<HalfEdge>(2);
                    for (int j = 0; j < 3; j++)
                    {
                        if (Edges[i].E.ElementAt(j).V[0].Equals(A)) eA.Add(Edges[i].E.ElementAt(j));
                        else if (Edges[i].E.ElementAt(j).V[0].Equals(B)) eB.Add(Edges[i].E.ElementAt(j));
                    }
                    if (eA.Count + eB.Count != 3) continue;

                    HalfEdge e0, e1, e2;
                    List<HalfEdge> e;
                    if (eA.Count == 1) { e0 = eA[0]; e = eB; } else { e0 = eB[0]; e = eA; }
                    double[] n0 = Edges[i].GetAngleNormal(e0, e[0]);
                    double[] n1 = Edges[i].GetAngleNormal(e0, e[1]);
                    double[] a1, a2;
                    if (n0[0] < n1[0]) { e1 = e[0]; e2 = e[1]; a1 = n0; a2 = n1; }
                    else { e1 = e[1]; e2 = e[0]; a1 = n1; a2 = n0; }
                    Edges[i].E = new HashSet<HalfEdge> { e0, e1, e2 };
                    Edges[i].Angle = new double[] { a1[0], a2[0] };
                    Edges[i].N = new double[] { a1[1], a1[2], a1[3], a2[1], a2[2], a2[3] };
                    E3.Add(Edges[i]);
                }
                else E0.Add(Edges[i]);
            }
        }

        public static TriangleMesh BySurfacesPoints(Surface[] Surfaces, Point[] Points) { return new TriangleMesh(Surfaces, Points); }

        public TriangleMesh AddEdgeNames(Point Point, int SplineCount = 2, int EdgeCount = 3)
        {
            // find start point
            Vertex v = Vertices.Find(vtx => vtx.X == Point.X && vtx.Y == Point.Y && vtx.Z == Point.Z);
            if (v.Equals(null)) return null;
            List<Edge> E0 = new List<Edge>(Edges);
            Spline S0 = new Spline(); S0.Vertices.Add(v);
            HashSet<Edge> vE = new HashSet<Edge>(v.Edges); vE.IntersectWith(E1);
            HalfEdge e = null;
            for (int i = 0; i < vE.Count; i++) if (vE.ElementAt(i).E.ElementAt(0).V[0].Equals(v)) e = vE.ElementAt(i).E.ElementAt(0);
            if (e.Equals(null)) return null;
            int s = 0;
            int sN = 1;
            while (vE.Count > 0)
            {
                if (sN > 1) e = vE.ElementAt(0).E.ElementAt(0);
                e.Face.ReOrderVertices(S0.Vertices[sN - 1]);
                e.Edge.Name = "s" + s.ToString("D" + SplineCount) + "-" + sN.ToString("D" + EdgeCount);
                E.Add(e.Edge.Name, e.Edge);
                E0.Remove(e.Edge);
                S0.Edges.Add(e.Edge);
                Vertex vS = e.Edge.GetOtherVertex(S0.Vertices[sN - 1]);
                if (!S0.Vertices.Contains(vS)) S0.Vertices.Add(vS);
                vE = new HashSet<Edge>(S0.Vertices[sN].Edges);
                vE.IntersectWith(E1);
                vE.IntersectWith(E0);
                sN++;
            }
            Splines.Add(S0);
            s++;
            // splines
            while (E0.Count > 0 && this.E.Count < Edges.Count)
            {
                for (int i = 0; i < Splines[s].Edges.Count; i++)
                {
                    Spline S = new Spline();
                    sN = 1;
                    int eN = 1;
                    int tN = 1;
                    List<Face> T0 = Splines[s].Edges[i].Faces.FindAll(f => f.Name.Equals(""));
                    if (T0.Count > 1) { }
                    S.Vertices.Add(T0[0].Vertices[2]);
                    for (int j = T0[0].E.Count - 1; j > 0; j--) if (T0[0].E[i].Edge.Name.Equals(""))
                    {
                        T0[0].E[i].Edge.Name = "e" + s.ToString("D" + SplineCount) + (s + 1).ToString("D" + SplineCount) + "-" + eN.ToString("D" + EdgeCount);
                        E.Add(T0[0].E[i].Edge.Name, T0[0].E[i].Edge);
                        E0.Remove(T0[0].E[i].Edge);
                        eN++;
                    }
                    T0[0].Name = "t" + s.ToString("D" + SplineCount) + (s + 1).ToString("D" + SplineCount) + "-" + tN.ToString("D" + EdgeCount);
                    tN++;
                    List<Face> T1 = T0[0].E[1].Edge.Faces.FindAll(f => f.Name.Equals(""));
                    if (T0.Count > 1) { }
                    Vertex v0 = (T1[0] as Triangle).GetOtherVertex(T0[0].E[1].Edge);
                    T1[0].ReOrderVertices(v0);
                    T1[0].E[0].Edge.Name = "s" + s.ToString("D" + SplineCount) + "-" + sN.ToString("D" + EdgeCount);
                    E.Add(T1[0].E[0].Edge.Name, T1[0].E[0].Edge);
                    E0.Remove(T1[0].E[0].Edge);
                    S.Edges.Add(T1[0].E[0].Edge);
                    S.Vertices.Add(e.Edge.GetOtherVertex(S.Vertices[sN - 1]));
                    sN++;
                }

            }
            return this;
        }

        public void Dispose() { Dispose(true); GC.SuppressFinalize(this); }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                Faces.ForEach(f => f.Dispose());
            }
            disposed = true;
        }
    }
}
