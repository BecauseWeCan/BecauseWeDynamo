using System;
using System.Collections.Generic;
using System.Linq;
using Autodesk.DesignScript.Geometry;
using Topology;

namespace Topology
{
    /// <summary>
    /// Triangle Mesh Object
    /// </summary>
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
        /// <summary>
        /// Face List
        /// </summary>
        public List<Triangle> Faces { get; set; }
        /// <summary>
        /// Edge List
        /// </summary>
        public List<Edge> Edges { get; set; }
        /// <summary>
        /// Outer Edge List
        /// </summary>
        public List<Edge> EdgesOuter { get { return E1; } }
        /// <summary>
        /// Vertex List
        /// </summary>
        public List<Vertex> Vertices { get { return V.Values.ToList(); } }
        /// <summary>
        /// Spline List
        /// </summary>
        public List<Spline> Splines { get; set; }
        /// <summary>
        /// Point List (Vertex Geometry List)
        /// </summary>
        public Point[] Points { get; set; }

        //**CONSTRUCTOR
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

            //**CONSTRUCT MESH
            // create vertex lookup table from points
            for (int i = 0; i < Points.Length; i++)
            {
                if (V.ContainsKey(Points[i])) continue;
                V.Add(Points[i], new Vertex(Points[i]));
            }
            this.Points = V.Keys.ToArray();
            // create faces from surfaces
            int eCount = 1;
            for (int i = 0; i < Surfaces.Length; i++)
            {
                // 

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
                Triangle t = new Triangle(v, Surfaces[i].NormalAtParameter(0.5, 0.5));
                t.Name = "t" + (i + 1).ToString("D" + 4);
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
                        e.Name = eCount.ToString("D" + 3);
                        e.E.Add(t.E[j]);
                        t.E[j].AddEdge(e);
                        Edges.Add(e);
                        eCount++;
                    }
                }
            }
            //**ITERATE THROUGH EDGES FOR ANGLE CALCS
            CalculateAngles();
        }
        internal TriangleMesh(Autodesk.DesignScript.Geometry.Mesh Mesh)
        {
            // initialize
            Faces = new List<Triangle>(Mesh.FaceIndices.Length);
            V = new Dictionary<Point, Vertex>(Mesh.VertexPositions.Length);
            E = new Dictionary<string, Edge>();
            Edges = new List<Edge>();
            E1 = new List<Edge>();
            E2 = new List<Edge>();
            E3 = new List<Edge>();
            E0 = new List<Edge>();
            Splines = new List<Spline>();

            //**CONSTRUCT MESH
            Points = Mesh.VertexPositions;
            // create vertex lookup table from points
            for (int i = 0; i < Points.Length; i++)
            {
                if (V.ContainsKey(Mesh.VertexPositions[i])) continue;
                V.Add(Points[i], new Vertex(Points[i]));
            }
            Points = V.Keys.ToArray();
            // create faces from surfaces
            for (int i = 0; i < Mesh.FaceIndices.Length; i++)
            {
                // find face vertices in lookup table
                IndexGroup iF = Mesh.FaceIndices[i];
                if (iF.Count == 3)
                {
                    List<Vertex> v = new List<Vertex>(3);
                    v.Add(Vertices[(int)iF.A]);
                    v.Add(Vertices[(int)iF.B]);
                    v.Add(Vertices[(int)iF.C]);
                    // create face based on vertices
                    Triangle t = Triangle.ByVertices(v);
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
                else
                {
                    bool flipDiagonal = false;
                    if (Vertices[(int)iF.A].DistanceTo(Vertices[(int)iF.C]) > Vertices[(int)iF.B].DistanceTo(Vertices[(int)iF.D]))
                        flipDiagonal = true;
                    List<Vertex> v1 = new List<Vertex>(3);
                    List<Vertex> v2 = new List<Vertex>(3);
                    if (!flipDiagonal)
                    {
                        v1.Add(Vertices[(int)iF.A]);
                        v1.Add(Vertices[(int)iF.B]);
                        v1.Add(Vertices[(int)iF.C]);
                        v2.Add(Vertices[(int)iF.C]);
                        v2.Add(Vertices[(int)iF.D]);
                        v2.Add(Vertices[(int)iF.A]);
                    }
                    else
                    {
                        v1.Add(Vertices[(int)iF.B]);
                        v1.Add(Vertices[(int)iF.C]);
                        v1.Add(Vertices[(int)iF.D]);
                        v2.Add(Vertices[(int)iF.D]);
                        v2.Add(Vertices[(int)iF.A]);
                        v2.Add(Vertices[(int)iF.B]);
                    }
                    // create face based on vertices
                    Triangle t1 = Triangle.ByVertices(v1);
                    Triangle t2 = Triangle.ByVertices(v2);
                    Faces.Add(t1);
                    Faces.Add(t2);
                    for (int j = 0; j < t1.E.Count; j++)
                    {
                        bool edgeFound = false;
                        if (Edges.Count > 0)
                        {
                            for (int k = 0; k < Edges.Count; k++)
                            {
                                if (Edges[k].Vertices.Contains(t1.E[j].V[0]) && Edges[k].Vertices.Contains(t1.E[j].V[1]))
                                {
                                    Edges[k].E.Add(t1.E[j]);
                                    t1.E[j].AddEdge(Edges[k]);
                                    edgeFound = true;
                                    break;
                                }
                            }
                        }
                        if (!edgeFound)
                        {
                            Edge e = new Edge();
                            e.E.Add(t1.E[j]);
                            t1.E[j].AddEdge(e);
                            Edges.Add(e);
                        }
                    }
                    for (int j = 0; j < t2.E.Count; j++)
                    {
                        bool edgeFound = false;
                        if (Edges.Count > 0)
                        {
                            for (int k = 0; k < Edges.Count; k++)
                            {
                                if (Edges[k].Vertices.Contains(t2.E[j].V[0]) && Edges[k].Vertices.Contains(t2.E[j].V[1]))
                                {
                                    Edges[k].E.Add(t2.E[j]);
                                    t2.E[j].AddEdge(Edges[k]);
                                    edgeFound = true;
                                    break;
                                }
                            }
                        }
                        if (!edgeFound)
                        {
                            Edge e = new Edge();
                            e.E.Add(t2.E[j]);
                            t2.E[j].AddEdge(e);
                            Edges.Add(e);
                        }
                    }
                }
            }
            //**ITERATE THROUGH EDGES FOR ANGLE CALCS
            CalculateAngles();
        }

        //**METHODS**CREATE
        /// <summary>
        /// creates a triangle mesh object based on given surfaces and points
        /// </summary>
        /// <param name="Surfaces">Surface Array</param>
        /// <param name="Points">Point Array</param>
        /// <returns>TriangleMesh Object</returns>
        public static TriangleMesh BySurfacesPoints(Surface[] Surfaces, Point[] Points) { return new TriangleMesh(Surfaces, Points); }
        /// <summary>
        /// creates a triangle mesh object based on given native mesh object
        /// </summary>
        /// <param name="Mesh">Mesh</param>
        /// <returns>TriangleMesh Object</returns>
        public static TriangleMesh ByMesh(Autodesk.DesignScript.Geometry.Mesh Mesh) { return new TriangleMesh(Mesh); }

        //**METHODS**ACTIONS
        /// <summary>
        /// finds vertex indexed by given point
        /// </summary>
        /// <param name="Point">Point Geometry</param>
        /// <returns>Vertex</returns>
        public Vertex GetVertexAtPoint(Point Point)
        {
            Vertex Result = null;
            for (int k = 0; k < V.Keys.Count; k++) if (Point.IsAlmostEqualTo(Points[k])) { Result = V[Points[k]]; break; }
            return Result;
        }
        /// <summary>
        /// finds edge indexed by given line
        /// </summary>
        /// <param name="Line">Line Geometry</param>
        /// <returns>Edge</returns>
        public Edge GetEdgeAtLine(Curve Line)
        {
            Edge Result = null;
            for (int k = 0; k < Edges.Count; k++) if (Edges[k].IsAtCurve(Line)) { Result = Edges[k]; break; }
            return Result;
        }
        /// <summary>
        /// analyzes mesh to find angles
        /// </summary>
        public void CalculateAngles()
        {
            //**ITERATE THROUGH EDGES FOR ANGLE CALCS
            for (int i = 0; i < Edges.Count; i++)
            {
                if (Edges[i].E.Count == 1) E1.Add(Edges[i]);
                else if (Edges[i].E.Count == 2)
                {
                    double[] AngleNormal = Edges[i].GetAngleNormal(Edges[i].E.ElementAt(0), Edges[i].E.ElementAt(1));
                    Edges[i].Angle = new double[] { AngleNormal[0] };
                    Edges[i].N = new double[] { AngleNormal[1], AngleNormal[2], AngleNormal[3] };
                    Edges[i].E.ElementAt(0).Angle = AngleNormal[0];
                    Edges[i].E.ElementAt(1).Angle = AngleNormal[0];
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
                    Vector e1N = e1.Face.Normal, e2N = e2.Face.Normal;
                    double[] a3 = { e1N.X - e2N.X, e1N.Y - e2N.Y, e1N.Z - e2N.Z };
                    e1N.Dispose(); e2N.Dispose();
                    Edges[i].E = new List<HalfEdge> { e0, e1, e2 };
                    Edges[i].Angle = new double[] { a1[0], a2[0] - a1[0], 360 - a2[0] };
                    Edges[i].N = new double[] { -a1[1], -a1[2], -a1[3], -a3[0], -a3[1], -a3[2], a2[1], a2[2], a2[3] };
                    Edges[i].E.ElementAt(0).Angle = a1[0];
                    Edges[i].E.ElementAt(1).Angle = Math.Min(a1[0], a2[0] - a1[0]);
                    Edges[i].E.ElementAt(2).Angle = a2[0] - a1[0];
                    E3.Add(Edges[i]);
                }
                else E0.Add(Edges[i]);
            }
        }

        //**METHODS**IN PROGRESS
        internal TriangleMesh AddEdgeNames(Point Point, int SplineCount = 2, int EdgeCount = 3)
        {
            // find start point
            Vertex v = GetVertexAtPoint(Point);
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
        internal TriangleMesh GetSplineData(Curve[] Start, Curve[] End, Curve[] Lines)
        {
            List<Edge> E = new List<Edge>(Start.Length);
            List<Curve> Search = new List<Curve>(Lines);
            Search.RemoveAll(c => Start.Contains(c));
            Search.RemoveAll(c => End.Contains(c));
            for (int i = 0; i < Search.Count; i++)
            {
                Edge e = GetEdgeAtLine(Search[i]);
                HalfEdge h;
                if (e.Equals(null))
                {
                    Vertex v1 = GetVertexAtPoint(Lines[i].StartPoint);
                    Vertex v2 = GetVertexAtPoint(Lines[i].EndPoint);
                    e = new Edge(new HashSet<HalfEdge> { new HalfEdge(v1, v2), new HalfEdge(v2, v1) });
                    e.E[0].AddEdge(e);
                    e.E[1].AddEdge(e);
                    v1.Edges.Add(e);
                    v2.Edges.Add(e);
                }
                E.Add(e);

            }
            // case zero
            for (int i = 0; i < Start.Length; i++)
            {
                List<Vertex> V = new List<Vertex>();
                Edge e0 = GetEdgeAtLine(Start[i]);
                V.Add(e0.E[0].V[0]);
                V.Add(e0.E[0].V[1]);
                HashSet<Edge> H = V[V.Count - 1].Edges;
                H.IntersectWith(E);
                V.Add(H.ElementAt(0).GetOtherVertex(V[V.Count - 1]));
                E.Remove(H.ElementAt(0));



            }
            return this;
        }

        //**METHODS**DISPOSE
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

    /// <summary>
    /// Mesh Object with Topology and Geometry
    /// </summary>
    public class Mesh : IDisposable
    {
        //**FIELDS**
        private bool disposed = false;
        internal Dictionary<Point, Vertex> V;
        internal Dictionary<string, Edge> E;
        internal List<Edge> E1;
        internal List<Edge> E2;
        internal List<Edge> E3;
        internal List<Edge> E0;
        internal int Df = 0;
        internal int De = 0;
        internal double MinFaceAngle = 180;

        //**PROPERTIES**QUERY
        /// <summary>
        /// Face List
        /// </summary>
        public List<Face> Faces { get; set; }
        /// <summary>
        /// Edge List
        /// </summary>
        public List<Edge> Edges { get; set; }
        /// <summary>
        /// Naked Edge List
        /// </summary>
        public List<Edge> EdgesOuter { get { return E1; } }
        /// <summary>
        /// Vertex List
        /// </summary>
        public List<Vertex> Vertices { get { return V.Values.ToList(); } }
        /// <summary>
        /// Spline List
        /// </summary>
        public List<Spline> Splines { get; set; }
        /// <summary>
        /// Point List (Vertex Geometry)
        /// </summary>
        public Point[] Points { get; set; }

        //**CONSTRUCTOR**
        internal Mesh()
        {
            V = new Dictionary<Point, Vertex>();
            E = new Dictionary<string, Edge>();
            Edges = new List<Edge>();
            E1 = new List<Edge>();
            E2 = new List<Edge>();
            E3 = new List<Edge>();
            E0 = new List<Edge>();
            Splines = new List<Spline>();
        }
        internal Mesh(Point[] Points)
            : this()
        {
            De = Points.Length;
            V = new Dictionary<Point, Vertex>(Points.Length);
            // create vertex lookup table from points
            for (int i = 0; i < Points.Length; i++)
            {
                if (V.ContainsKey(Points[i])) continue;
                V.Add(Points[i], new Vertex(Points[i]));
            }
            this.Points = V.Keys.ToArray();
        }
        internal Mesh(Point[] Points, Surface[] Surfaces)
            : this(Points)
        {
            Df = Surfaces.Length.ToString().Length;
            De = (Surfaces.Length + Points.Length).ToString().Length;
            // initialize
            Faces = new List<Face>(Surfaces.Length);
            //**CONSTRUCT MESH
            // create faces from surfaces
            int eCount = 1;
            for (int i = 0; i < Surfaces.Length; i++)
            {
                // find face vertices in lookup table
                // create face based on vertices
                if (Surfaces[i].Vertices.Length == 3)
                {
                    Triangle t = new Triangle(FindFaceVertices(Surfaces[i]), Surfaces[i].NormalAtParameter(0.5, 0.5));
                    t.Name = "t" + (i + 1).ToString("D" + Df);
                    Faces.Add(t);
                    eCount = FindFaceEdges(t, eCount);
                }
                else
                {
                    Face f = new Face(FindFaceVertices(Surfaces[i]), Surfaces[i].NormalAtParameter(0.5,0.5));
                    f.Name = "f" + (i + 1).ToString("D" + Df);
                    Faces.Add(f);
                    eCount = FindFaceEdges(f, eCount);
                }
            }
            //**ITERATE THROUGH EDGES FOR ANGLE CALCS
            CalculateAngles();
        }

        //**METHODS**CREATE
        /// <summary>
        /// creates mesh object with given points and surfaces
        /// surfaces must be planar but can be concave as well as polygonal
        /// </summary>
        /// <param name="Points">Vertex Point Array</param>
        /// <param name="Surfaces">Surface Array</param>
        /// <returns>Mesh Object</returns>
        public static Mesh BySurfacesPoints(Point[] Points, Surface[] Surfaces) { return new Mesh(Points, Surfaces); }
        /// <summary>
        /// creates mesh object with given surfaces and derives vertices
        /// surfaces must be planar but can be concave as well as polygonal
        /// </summary>
        /// <param name="Surfaces">Surface Array</param>
        /// <returns>Mesh Object</returns>
        public static Mesh BySurfaces(Surface[] Surfaces) 
        {
            List<Point> Points = new List<Point>();
            for (int i = 0; i < Surfaces.Length; i++ )
            {
                for (int j=0; j < Surfaces[i].Vertices.Length; j++)
                {
                    if (Points.Contains(Surfaces[i].Vertices[j].PointGeometry)) continue;
                    Points.Add(Surfaces[i].Vertices[j].PointGeometry);
                }
            }
            return new Mesh(Points.ToArray(), Surfaces); 
        }

        //**METHODS**ACTIONS
        /// <summary>
        /// returns topological entity (vertex) at point
        /// </summary>
        /// <param name="Point">Vertex Geometry</param>
        /// <returns>Vertex</returns>
        public Vertex GetVertexAtPoint(Point Point)
        {
            Vertex Result = null;
            for (int k = 0; k < V.Keys.Count; k++) if (Point.IsAlmostEqualTo(Points[k])) { Result = V[Points[k]]; break; }
            return Result;
        }
        /// <summary>
        /// returns topological entity (edge) at line
        /// </summary>
        /// <param name="Line">Edge Geometry</param>
        /// <returns>Edge</returns>
        public Edge GetEdgeAtLine(Curve Line)
        {
            Edge Result = null;
            for (int k = 0; k < Edges.Count; k++) if (Edges[k].IsAtCurve(Line)) { Result = Edges[k]; break; }
            return Result;
        }
        /// <summary>
        /// calculates edge angles based on input geometry
        /// </summary>
        public void CalculateAngles()
        {
            //**ITERATE THROUGH EDGES FOR ANGLE CALCS
            for (int i = 0; i < Edges.Count; i++)
            {
                // case 0: is OuterEdge
                if (Edges[i].E.Count == 1) E1.Add(Edges[i]);
                // case 1: is toplogically consistent edge
                else if (Edges[i].E.Count == 2)
                {
                    double[] AngleNormal = Edges[i].GetAngleNormal(Edges[i].E.ElementAt(0), Edges[i].E.ElementAt(1));
                    Edges[i].Angle = new double[] { AngleNormal[0] };
                    Edges[i].N = new double[] { AngleNormal[1], AngleNormal[2], AngleNormal[3] };
                    Edges[i].E.ElementAt(0).Angle = AngleNormal[0];
                    Edges[i].E.ElementAt(1).Angle = AngleNormal[0];
                    E2.Add(Edges[i]);
                }
                // case 2: topologically inconsistent
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
                    Vector e1N = e1.Face.Normal, e2N = e2.Face.Normal;
                    double[] a3 = { e1N.X - e2N.X, e1N.Y - e2N.Y, e1N.Z - e2N.Z };
                    e1N.Dispose(); e2N.Dispose();
                    Edges[i].E = new List<HalfEdge> { e0, e1, e2 };
                    Edges[i].Angle = new double[] { a1[0], a2[0] - a1[0], 360 - a2[0] };
                    Edges[i].N = new double[] { -a1[1], -a1[2], -a1[3], -a3[0], -a3[1], -a3[2], a2[1], a2[2], a2[3] };
                    Edges[i].E.ElementAt(0).Angle = a1[0];
                    Edges[i].E.ElementAt(1).Angle = Math.Min(a1[0], a2[0] - a1[0]);
                    Edges[i].E.ElementAt(2).Angle = a2[0] - a1[0];
                    E3.Add(Edges[i]);
                }
                else E0.Add(Edges[i]);
            }
        }
        /// <summary>
        /// looks for Face Edge in Mesh Edge List,
        /// adds topological information if found
        /// otherwise creates Edge and adds to Mesh Edge List
        /// </summary>
        /// <param name="Face">Face</param>
        /// <param name="EdgeCount">ExistingEdgeCount</param>
        /// <returns>NewEdgeCount</returns>
        public int FindFaceEdges(Face Face, int EdgeCount)
        {
            int eCount = EdgeCount;
            for (int j = 0; j < Face.E.Count; j++)
            {
                bool edgeFound = false;
                if (Edges.Count > 0)
                {
                    for (int k = 0; k < Edges.Count; k++)
                    {
                        if (Edges[k].Vertices.Contains(Face.E[j].V[0]) && Edges[k].Vertices.Contains(Face.E[j].V[1]))
                        {
                            Edges[k].E.Add(Face.E[j]);
                            Face.E[j].AddEdge(Edges[k]);
                            edgeFound = true;
                            break;
                        }
                    }
                }
                if (!edgeFound)
                {
                    Edge e = new Edge();
                    e.Name = eCount.ToString("D" + De);
                    e.E.Add(Face.E[j]);
                    Face.E[j].AddEdge(e);
                    Edges.Add(e);
                    eCount++;
                }
            }
            return eCount;
        }
        /// <summary>
        /// looks for Face Vertex in Mesh Vertex List,
        /// adds topological information if found
        /// otherwise creates Vertex and adds to Mesh Vertex List
        /// </summary>
        /// <param name="Surface"></param>
        /// <returns></returns>
        public List<Vertex> FindFaceVertices(Surface Surface)
        {
            Autodesk.DesignScript.Geometry.Vertex[] vtx = Surface.Vertices;
            List<Vertex> v = new List<Vertex>(vtx.Length);
            for (int j = 0; j < vtx.Length; j++)
            {
                Point pt = vtx[j].PointGeometry;
                Vertex search = GetVertexAtPoint(pt);
                if (!search.Equals(null))
                {
                    v.Add(search);
                    pt.Dispose();
                }
                else
                {
                    V.Add(pt, new Vertex(pt));
                    v.Add(GetVertexAtPoint(pt));
                }
            }
            vtx.ForEach(x => x.Dispose());
            return v;
        }
        /// <summary>
        /// flips face normals
        /// </summary>
        /// <returns>Mesh</returns>
        public Mesh FlipFaceNormals()
        {
            Faces.ForEach(f => f.E.ForEach(e => e.FlipDirection()));
            return this;
        }

        //**METHODS**IN PROGRESS
        internal Mesh AddEdgeNames(Point Point, int SplineCount = 2, int EdgeCount = 3)
        {
            // find start point
            Vertex v = GetVertexAtPoint(Point);
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
        internal Mesh GetSplineData(Curve[] Start, Curve[] End, Curve[] Lines)
        {
            List<Edge> E = new List<Edge>(Start.Length);
            List<Curve> Search = new List<Curve>(Lines);
            Search.RemoveAll(c => Start.Contains(c));
            Search.RemoveAll(c => End.Contains(c));
            for (int i = 0; i < Search.Count; i++)
            {
                Edge e = GetEdgeAtLine(Search[i]);
                HalfEdge h;
                if (e.Equals(null))
                {
                    Vertex v1 = GetVertexAtPoint(Lines[i].StartPoint);
                    Vertex v2 = GetVertexAtPoint(Lines[i].EndPoint);
                    e = new Edge(new HashSet<HalfEdge> { new HalfEdge(v1, v2), new HalfEdge(v2, v1) });
                    e.E[0].AddEdge(e);
                    e.E[1].AddEdge(e);
                    v1.Edges.Add(e);
                    v2.Edges.Add(e);
                }
                E.Add(e);

            }
            // case zero
            for (int i = 0; i < Start.Length; i++)
            {
                List<Vertex> V = new List<Vertex>();
                Edge e0 = GetEdgeAtLine(Start[i]);
                V.Add(e0.E[0].V[0]);
                V.Add(e0.E[0].V[1]);
                HashSet<Edge> H = V[V.Count - 1].Edges;
                H.IntersectWith(E);
                V.Add(H.ElementAt(0).GetOtherVertex(V[V.Count - 1]));
                E.Remove(H.ElementAt(0));



            }
            return this;
        }

        //**METHODS**DISPOSE
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
