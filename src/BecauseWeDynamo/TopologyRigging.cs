using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;
using Text;

namespace Topology
{
    public class TriangleVertex: Vertex
    {
        //**QUERY**PROPERTIES
        /// <summary>
        /// Point of vertex
        /// </summary>
        public List<TrianglePanel> Triangles { get; private set; }
        public new List<TriangleEdge> Edges { get; private set; }
        public int SplineId { get; set; }

        //*CONSTRUCTOR
        internal TriangleVertex(Point pt): base(pt)
        {
            Triangles = new List<TrianglePanel>();
            Edges = new List<TriangleEdge>();
            SplineId = -1;
        }

        //**CREATE
        /// <summary>
        /// make vertex from Point
        /// </summary>
        /// <param name="pt">Point</param>
        /// <returns></returns>
        public static TriangleVertex ByPoint(Point pt) { return new TriangleVertex(pt); }

        //**ACTIONS**METHODS
        /// <summary>
        /// adds Triangle information to vertex
        /// </summary>
        /// <param name="triangle">Triangle object to be added</param>
        public void AddTriangle(TrianglePanel triangle)
        {
            if (!Triangles.Contains(triangle)) Triangles.Add(triangle);
        }
        public void AddEdge(TriangleEdge edge)
        {
            if (!Edges.Contains(edge)) Edges.Add(edge);
        }
    }

    public class TriangleEdge: IEquatable<TriangleEdge>
    {
        //**QUERY**PROPERTIES
        public List<TriangleVertex> Vertices { get; private set; }
        public TriangleVertex A { get { return Vertices[0]; } }
        public TriangleVertex B { get { return Vertices[1]; } }
        public List<TrianglePanel> Triangles { get; private set; }
        public Point MidPoint { get { return Point.ByCoordinates((A.X + B.X) / 2, (A.Y + B.Y) / 2, (A.Z + B.Z) / 2); } }
        public string Name { get; set; }
        public Boolean IsOuterEdge { get { return (Triangles.Count == 1); } }

        //*CONSTRUCTOR
        internal TriangleEdge(TriangleVertex pta, TriangleVertex ptb, string name, List<TrianglePanel> Triangles)
        {
            Vertices = new List<TriangleVertex> { pta, ptb };
            Name = name;
            this.Triangles = Triangles;
        }

        //**CREATE
        public static TriangleEdge ByPoints(TriangleVertex pta, TriangleVertex ptb, string name, List<TrianglePanel> Triangles) { return new TriangleEdge(pta, ptb, name, Triangles); }

        //**ACTIONS**METHODS
        public void AddTriangle(TrianglePanel t) { if (!Triangles.Contains(t)) Triangles.Add(t); }   
        public Boolean HasVertices(TriangleVertex vtx1, TriangleVertex vtx2) { return ((vtx1.Equals(A) && vtx2.Equals(B)) || (vtx1.Equals(B) && vtx2.Equals(A))); }
        public Line GetLine()  
        {
            Point a = A.Point;
            Point b = B.Point;
            Line l = Line.ByStartPointEndPoint(a, b);
            a.Dispose(); b.Dispose();
            return l;
        }
        public Vector GetNormal()
        {
            Vector N = Vector.ByCoordinates(0,0,0);
            for (int i = 0; i < Triangles.Count; i++)
            {
                Vector v = Triangles[i].Normal;
                N = N.Add(v);
                v.Dispose();
            }
            N = N.Normalized();
            return N;
        }
        public List<PolyCurve> GetEdgeLabels(double factor)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            int a, b;
            for (int i = 0; i < Triangles.Count; i++)
            {
                Point c = Triangles[i].Center;
                for (int j = 0; j < Triangles[i].Edges.Count; j++)
                {
                    if (Triangles[i].Edges[j].Name.Equals(this.Name))
                    {
                        Point m = Triangles[i].Edges[j].MidPoint;
                        Point ta = Triangles[i].Edges[j].A.Point;
                        Point tb = Triangles[i].Edges[j].B.Point;
                        a = Triangles[i].Vertices.IndexOf(Triangles[i].Edges[j].A);
                        b = Triangles[i].Vertices.IndexOf(Triangles[i].Edges[j].B);
                        Vector A = Vector.ByTwoPoints(m, ta);
                        Vector B = Vector.ByTwoPoints(m, tb);
                        Vector M = Vector.ByTwoPoints(c, m);
                        Point p = Triangles[i].Geometry.ClosestPointTo(m);
                        Word w;
                        if (((a - b) % 3 + 3) % 3 == 1)
                            w = Word.ByStringOriginVectors(Triangles[i].Edges[j].Name, p, B, M);
                        else if (((a - b) % 3 + 3) % 3 == 2)
                            w = Word.ByStringOriginVectors(Triangles[i].Edges[j].Name, p, A, M);
                        else
                            w = null;
                        labels.AddRange(w.display(factor));

                        m.Dispose(); ta.Dispose(); tb.Dispose(); A.Dispose(); B.Dispose(); M.Dispose(); p.Dispose(); w.Dispose();
                    }
                }
                c.Dispose();
            }
            return labels;
        }
        public TriangleVertex GetOtherVertex(TriangleVertex vtx)
        {
            if (vtx.Equals(A)) return B;
            else if (vtx.Equals(B)) return A;
            return null;
        }
        public override bool Equals(Object Object) { return this.Equals(Object as TriangleEdge); }
        public bool Equals(TriangleEdge Edge)
        {
            if (Object.ReferenceEquals(Edge, null)) return false;
            if (Object.ReferenceEquals(this, Edge)) return true;
            if (this.GetType() != Edge.GetType()) return false;
            return ((A.Equals(Edge.A) && B.Equals(Edge.B)) || (A.Equals(Edge.B) && B.Equals(Edge.A)));
        }
        public int GetHashCode(TriangleEdge e)
        {
            return e.GetHashCode();
        }
    }

    public class TriangleEdgeConnector
    {
        private double radius;
        private double spacing;
        private bool disposed = false;
        public int[] numHoles { get; set; }
        public TriangleEdge Edge { get; set; }

        internal TriangleEdgeConnector(TriangleEdge tEdge, double radiusHole, double spacingHole)
        {
            spacing = spacingHole;
            radius = radiusHole;
            Edge = tEdge;
            numHoles = new int[3] { 0, 0, 0 };
            for (int i = 0; i < Edge.Triangles.Count; i++)
            {
                numHoles[i] = (int)(Math.Ceiling((Edge.Triangles[i].Geometry.ClosestPointTo(Edge.MidPoint).DistanceTo(Edge.MidPoint) + 3 * radius) / spacing)) + Edge.Triangles[i].numHoles - 1;
            }
            numHoles[2] = numHoles[0] + numHoles[1] + 1;
        }

        public static TriangleEdgeConnector ByEdge(TriangleEdge edge, double radiusHole, double spacingHole) { return new TriangleEdgeConnector(edge, radiusHole, spacingHole); }

        public List<Circle> PlaceHoles(double DisplayFactor = 1)
        {
            if (numHoles[2] < 2 || Edge.IsOuterEdge) return null;
            List<Circle> result = new List<Circle>(numHoles[2]);
            Point m = Point.ByCoordinates(0, 0, 0);
            Point c = Point.ByCoordinates(0, 0, 0);
            Vector T = Vector.ByCoordinates(0, 0, 0);
            for (int i = 0; i < Edge.Triangles.Count; i++)
            {
                m = Edge.Triangles[i].Geometry.ClosestPointTo(Edge.MidPoint);
                for (int k = 0; k < numHoles[i]; k++)
                {
                    T = Vector.ByTwoPoints(Edge.MidPoint, m);
                    c = (Point)Edge.MidPoint.Translate(T, spacing * (k + 1));
                    result.Add(Circle.ByCenterPointRadiusNormal(c, DisplayFactor * radius, Edge.Triangles[i].Normal));
                }
            }
            result.Add(Circle.ByCenterPointRadiusNormal(Edge.MidPoint, DisplayFactor * radius, Edge.GetNormal()));
            m.Dispose();
            c.Dispose();
            T.Dispose();
            return result;
        }
        public List<Circle> PlaceHolesByHoleCount(int numHoles, double DisplayFactor = 1)
        {
            if (this.numHoles[2] == numHoles) return PlaceHoles(DisplayFactor);
            return null;
        }
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                this.Edge = null;
            }
            disposed = true;
        }

    }



    public class TrianglePanel: IDisposable
    {
        //*PRIVATE**PROPERTIES
        public Dictionary<string, Object> Parameters;
        public Boolean HasName;
        private double Height;

        bool disposed = false;

        //**PROPERTIES - **QUERY**
        public CoordinateSystem CS { get; private set; }
        public List<TriangleEdge> Edges { get; private set; }
        public List<Point> Points { get; private set; }
        public List<TriangleVertex> Vertices { get; private set; }
        public Point Circumcenter { get; private set; }
        public Point Center { get; private set; }
        public Vector Normal { get; set; }
        public Geometry Geometry { get; set; }
        public string Name { get; set; }
        public int splineId { get; set; }
        public int numHoles { get; set; }
        public Curve[] PerimeterCurves { get {
            if (!(Geometry is Solid) && !(Geometry is Surface)) return null;
            if (Geometry is Surface) return ((Surface)Geometry).PerimeterCurves();
            if (Geometry is Solid)
            {
                Geometry[] g = Geometry.Explode();
                Curve[] result = ((Surface) g[0]).PerimeterCurves();
                g.ForEach(i => i.Dispose());
                return result;
            }
            return null;
        }}
        public Boolean HasEdges { get {
                if (Edges.Count == 3) return true;
                else return false;
        }}

        //**CONSTRUCTOR
        internal TrianglePanel(Geometry inputGeometry, List<Point> Points)
        {
            Geometry = inputGeometry;
            Points = new List<Point>(3);
            Edges = new List<TriangleEdge>(3);
            Vertices = new List<TriangleVertex>(3);
            Name = "";
            HasName = false;
            Height = 0;
            numHoles = 1;

            Points.AddRange(Points);
            Normal = Plane.ByBestFitThroughPoints(Points).Normal;
            Circle c = Circle.ByThreePoints(Points[0], Points[1], Points[2]);
            Circumcenter = c.CenterPoint;
            c.Dispose();
            Center = Point.ByCoordinates((Points[0].X + Points[1].X + Points[2].X) / 3.0, (Points[0].Y + Points[1].Y + Points[2].Y) / 3.0, (Points[0].Z + Points[1].Z + Points[2].Z) / 3.0);
            CS = CoordinateSystem.ByOriginVectors(Center, Vector.ByTwoPoints(Points[0], Points[1]), Vector.ByTwoPoints(Points[0], Points[2]));
        }

        //**CREATE**
        public static TrianglePanel BySurfaceAndPoints(Surface srf, List<Point> Points) { return new TrianglePanel(srf, Points); }


        public void BuildGeometricProperties(TriangleEdge edge)
        {
            int a = Vertices.IndexOf(edge.A);
            int b = Vertices.IndexOf(edge.B);
            int c = 3 - (a + b);
            Vector X;
            Vector Y;
            if (((a - b) % 3 + 3) % 3 == 1)
            {
                X = Vector.ByTwoPoints(Points[b], Points[a]);
                Y = Vector.ByTwoPoints(Points[b], Points[c]);
            }
            else if (((a - b) % 3 + 3) % 3 == 2)
            {
                X = Vector.ByTwoPoints(Points[a], Points[b]);
                Y = Vector.ByTwoPoints(Points[a], Points[c]);
            }
            else
            {
                X = Vector.XAxis();
                Y = Vector.YAxis();
            }
            CS = CoordinateSystem.ByOriginVectors(Center, X, Y);
            X.Dispose(); Y.Dispose();
            double w = Points[a].DistanceTo(Points[b]);
            double sa = Points[c].DistanceTo(Points[a]);
            double sb = Points[c].DistanceTo(Points[b]);
            double s = (w + sa + sb) / 2;
            Height = 2 * Math.Sqrt(s * (s - w) * (s - sa) * (s - sb)) / w;
        }
        /// <summary>
        /// adjusts CS
        /// </summary>
        public void BuildCS()
        {
            for (int i = 0; i < Edges.Count; i++)
            {
                if (Edges[i].Name.Substring(1, 2).Equals(Name.Substring(3, 2)))
                {
                    CS = CS.Rotate(CS.Origin, CS.ZAxis, 180);
                }
            }
        }
        /// <summary>
        /// add edge reference to triangle
        /// </summary>
        /// <param name="edge">triangle edge</param>
        public void AddEdge(TriangleEdge edge) { if (!Edges.Contains(edge)) Edges.Add(edge); }
        /// <summary>
        /// add vertex reference to triangle
        /// </summary>
        /// <param name="vtx">triangle vertex</param>
        public void AddVertex(TriangleVertex vtx) { if (!Vertices.Contains(vtx)) Vertices.Add(vtx); }
        /// <summary>
        /// adds circles to Edges
        /// </summary>
        /// <param name="factor">text scale</param>
        /// <returns></returns>
        public List<PolyCurve> GetEdgeLabels(double factor)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < Edges.Count; j++)
            {
                if (Edges[j].IsOuterEdge) continue;

                int a = Vertices.IndexOf(Edges[j].A);
                int b = Vertices.IndexOf(Edges[j].B);

                Point p = Geometry.ClosestPointTo(Edges[j].MidPoint);
                Point m = Edges[j].MidPoint;
                Point A = Edges[j].A.Point;
                Point B = Edges[j].B.Point;
                Vector mA = Vector.ByTwoPoints(m, A);
                Vector mB = Vector.ByTwoPoints(m, B);
                Vector cm = Vector.ByTwoPoints(Center,m);
                Word w;
                if (((a - b) % 3 + 3) % 3 == 1)  w = Word.ByStringOriginVectors(Edges[j].Name, p, mB, cm);
                else if (((a - b) % 3 + 3) % 3 == 2) w = Word.ByStringOriginVectors(Edges[j].Name, p, mA, cm);
                else w=null;
                labels.AddRange(w.display(factor));
                p.Dispose(); m.Dispose(); A.Dispose(); B.Dispose(); mA.Dispose(); mB.Dispose(); cm.Dispose(); w.Dispose();
            } // end Edges
            return labels;
        }// end action
        /// <summary>
        /// get profiles and letters
        /// </summary>
        /// <param name="factor">letter size</param>
        /// <returns>curves</returns>
        public List<Curve> GetTriangleCurves(double factor)
        {
            List<Curve> curves = new List<Curve>();
            curves.AddRange(GetEdgeLabels(factor));
            curves.AddRange(PerimeterCurves);
            return curves;
        }
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                Circumcenter.Dispose();
                Center.Dispose();
                CS.Dispose();
                Normal.Dispose();
                Geometry.Dispose();
            }
            disposed = true;
        }
    }

    public class TriangleRigging: IDisposable
    {
        //**GLOBAL**VARIABLES
        private Dictionary<Point, TriangleVertex> V;
        bool disposed = false;
        //**PROPERTIES - **QUERY**
        public List<Point> VertexPoints { get; set; }
        public List<TrianglePanel> Triangles { get; set; }
        public List<List<TriangleVertex>> Splines { get; set; }
        public List<TriangleEdge> Edges { get; set; }
        public List<TriangleVertex> Vertices { get { return V.Values.ToList(); } }

        //**CONSTRUCTORS
        internal TriangleRigging(Solid[][] Solids, Point[][] SolidPoints, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            V = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<TrianglePanel>(Solids.Length);
            Edges = new List<TriangleEdge>();
            Splines = new List<List<TriangleVertex>>(OrderedSplines.Length);

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < Solids.Length; i++)
            {
                Triangles.Add(new TrianglePanel(Solids[i][0], SolidPoints[i].ToList()));
                for (int j = 0; j < SolidPoints[i].Length; j++)
                {
                    if (!V.ContainsKey(SolidPoints[i][j]))
                    {
                        V.Add(SolidPoints[i][j], new TriangleVertex(SolidPoints[i][j]));
                    }
                    V[SolidPoints[i][j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(V[SolidPoints[i][j]]);
                }
            }

            VertexPoints = V.Keys.ToList();

            BuildEdges(OrderedSplines);
        }
        internal TriangleRigging(Surface[] TriangleSurfaces, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            V = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<TrianglePanel>(TriangleSurfaces.Length);
            Edges = new List<TriangleEdge>();
            Splines = new List<List<TriangleVertex>>(OrderedSplines.Length);

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < TriangleSurfaces.Length; i++)
            {
                List<Point> pts = new List<Point>(3);
                TriangleSurfaces[i].Vertices.ForEach(vtx => pts.Add(vtx.PointGeometry));
                Triangles.Add(new TrianglePanel(TriangleSurfaces[i], pts));
                for (int j = 0; j < pts.Count; j++)
                {
                    if (!V.ContainsKey(pts[j]))
                    {
                        V.Add(pts[j], new TriangleVertex(pts[j]));
                    }
                    V[pts[j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(V[pts[j]]);
                }
            }

            VertexPoints = V.Keys.ToList();

            BuildEdges(OrderedSplines);
        }
        

        //**PRIVATE**METHODS
        internal void BuildEdges(Point[][] OrderedSplines)
        {
            //build splines
            int numD = OrderedSplines.Length.ToString().Length;
            for (int i = 0; i < OrderedSplines.Length; i++)
            {
                List<TriangleVertex> Spline = new List<TriangleVertex>();
                int z = 0;
                int edgeCount = 1;
                int numDs = OrderedSplines[i].Length.ToString().Length;
                // iterate through spline Vertices
                for (int j = 0; j < OrderedSplines[i].Length; j++)
                {
                    // search vertex Points for match
                    for (int k = 0; k < VertexPoints.Count; k++)
                    {
                        if (VertexPoints[k].IsAlmostEqualTo(OrderedSplines[i][j]))
                        {
                            Spline.Add(V[VertexPoints[k]]);
                            V[VertexPoints[k]].SplineId = i;

                            if (j > 0)
                            {
                                // check Triangles
                                // find Triangles on edge
                                List<TrianglePanel> listT = new List<TrianglePanel>();
                                for (int t = 0; t < V[VertexPoints[k]].Triangles.Count; t++)
                                {
                                    if (V[VertexPoints[z]].Triangles.Contains(V[VertexPoints[k]].Triangles[t]))
                                        listT.Add(V[VertexPoints[k]].Triangles[t]);
                                }

                                if (listT.Count > 0)
                                {
                                    string edgeId = "s" + i.ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDs);
                                    TriangleEdge e = new TriangleEdge(V[VertexPoints[z]], V[VertexPoints[k]], edgeId, listT);
                                    Edges.Add(e);
                                    // add edge key to Vertices
                                    for (int t = 0; t < listT.Count; t++)
                                    {
                                        listT[t].AddEdge(e); 
                                        listT[t].BuildGeometricProperties(e);
                                    }

                                    V[VertexPoints[z]].AddEdge(e);
                                    V[VertexPoints[k]].AddEdge(e);
                                    edgeCount++;
                                }
                            }
                            z = k;
                            break;
                        }
                    } // end search
                } // end spline Point iteration
                // check Triangles
                // find Triangles on edge
                if (OrderedSplines[i].Length > 2)
                {
                    List<TrianglePanel> listTlast = new List<TrianglePanel>();
                    for (int t = 0; t < Spline[0].Triangles.Count; t++)
                    {
                        if (V[VertexPoints[z]].Triangles.Contains(Spline[0].Triangles[t]))
                            listTlast.Add(Spline[0].Triangles[t]);
                    }

                    if (listTlast.Count > 0)
                    {
                        string edgeId = "s" + i.ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDs);
                        TriangleEdge e = new TriangleEdge(V[VertexPoints[z]], Spline[0], edgeId, listTlast);
                        Edges.Add(e);
                        // add edge key to Vertices
                        for (int t = 0; t < listTlast.Count; t++)
                        {
                            listTlast[t].AddEdge(e);
                        }
                        V[VertexPoints[z]].AddEdge(e);
                        Spline[0].AddEdge(e);
                        edgeCount++;
                    }
                }
                Splines.Add(Spline);
            } // end build spline


            //build Edges
            for (int s = 0; s < Splines.Count - 1; s++)
            {
                // for each spline
                int edgeCount = 0;
                int triangleCount = 1;
                int numDe = Splines[s].Count.ToString().Length;
                for (int v = 0; v < Splines[s].Count; v++)
                {
                    List<TriangleVertex> Search = new List<TriangleVertex>();
                    for (int t = 0; t < Splines[s][v].Triangles.Count; t++)
                    {
                        if (!Splines[s][v].Triangles[t].HasEdges)
                        {
                            for (int vt = 0; vt < Splines[s][v].Triangles[t].Vertices.Count; vt++)
                            {
                                if (!Search.Contains(Splines[s][v].Triangles[t].Vertices[vt]) && Splines[s][v].Triangles[t].Vertices[vt].SplineId == s + 1)
                                    Search.Add(Splines[s][v].Triangles[t].Vertices[vt]);
                            }
                        }
                    }
                    Search.Remove(Splines[s][v]);

                    for (int sv = 0; sv < Search.Count; sv++)
                    {
                        // check Triangles
                        // find Triangles on edge
                        List<TrianglePanel> listT = new List<TrianglePanel>();
                        for (int t = 0; t < Search[sv].Triangles.Count; t++)
                        {
                            if (Splines[s][v].Triangles.Contains(Search[sv].Triangles[t]))
                                listT.Add(Search[sv].Triangles[t]);
                        }

                        if (listT.Count > 0)
                        {
                            string edgeId = "e" + s.ToString("D" + numD) + (s + 1).ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDe);
                            TriangleEdge e = new TriangleEdge(Splines[s][v], Search[sv], edgeId, listT);
                            Edges.Add(e);
                            // add edge key to Vertices
                            for (int t = 0; t < listT.Count; t++)
                            {
                                listT[t].AddEdge(e);
                                if (!listT[t].HasName)
                                {
                                    listT[t].Name = "t" + s.ToString("D" + numD) + (s + 1).ToString("D" + numD) + "-" + triangleCount.ToString("D" + numDe);
                                    listT[t].HasName = true;
                                    listT[t].splineId = s;
                                    listT[t].BuildCS();
                                    triangleCount++;
                                }
                            }
                            Splines[s][v].AddEdge(e);
                            Search[sv].AddEdge(e);
                            edgeCount++;
                        }
                    }


                } // end spline

            } // end build Edges 
        }

        //**STATIC**METHODS - **CREATE**
        public static TriangleRigging BySolidsPointsSplines(Solid[][] Solids, Point[][] SolidsPoints, Point[][] Splines) { return new TriangleRigging(Solids, SolidsPoints, Splines); }

        //**METHODS - **ACTION**
        public List<List<PolyCurve>> GetEdgeLabels(double factor)
        {
            List<List<PolyCurve>> labels = new List<List<PolyCurve>>(Triangles.Count);
            for (int i = 0; i < Triangles.Count; i++) { labels.Add(Triangles[i].GetEdgeLabels(factor)); }
            return labels;
        }
        /// <summary>
        /// gets index set of Triangles sorted by angles
        /// </summary>
        /// <returns></returns>
        public List<int> GetSortedTriangleIndexByName()
        {
            List<int> index = new List<int>(Triangles.Count);
            Dictionary<int, string> sorted = new Dictionary<int, string>();
            for (int i = 0; i < Triangles.Count; i++) { sorted.Add(i, Triangles[i].Name); }
            foreach (var item in sorted.OrderBy(i => i.Value)) { index.Add(item.Key); }

            return index;
        }

        public void SetHolesBySpline(int[] holes)
        {
            for (int t = 0; t < Triangles.Count; t++) Triangles[t].numHoles = holes[Triangles[t].splineId];
        }
        public List<TriangleEdgeConnector> PlaceEdgeConnectorHolesBySpline(int[] holes, double radius, double spacing)
        {
            SetHolesBySpline(holes);
            List<TriangleEdgeConnector> result = new List<TriangleEdgeConnector>();
            for (int i = 0; i < Edges.Count; i++)
            {
                if (Edges[i].IsOuterEdge) continue;
                result.Add(new TriangleEdgeConnector(Edges[i], radius, spacing));
            }
            return result;
        }
        public List<TriangleEdgeConnector> GetEdgeConnectorByHoleCount(int HoleCount, List<TriangleEdgeConnector> input)
        {
            List<TriangleEdgeConnector> result = new List<TriangleEdgeConnector>();
            for (int i = 0; i < input.Count; i++) { if (input[i].numHoles[2] == HoleCount) result.Add(input[i]); }
            return result;
        }
        public void AddParamterToTriangles(string name, Object[] data)
        {
            if (data.Length == Triangles.Count) for (int i = 0; i < Triangles.Count; i++) { Triangles[i].Parameters.Add(name, data[i]); }
        }
        public List<int> GetTriangleIndexByParameterValue(string parameter, Object value)
        {
            if (!Triangles[0].Parameters.Keys.Contains(parameter)) return null;
            List<int> result = new List<int>();
            for (int i = 0; i < Triangles.Count; i++)
            {
                if (Triangles[i].Parameters[parameter].Equals(value)) result.Add(i);
            }
            return result;
        }

        public List<List<int>> GetTriangleRowIndexBySpline()
        {
            int numRow = Splines.Count - 1;
            int numD = Splines.Count.ToString().Length;
            List<int> Index = GetSortedTriangleIndexByName();
            List<int> Row = new List<int>(Index.Count);
            List<List<int>> result = new List<List<int>>(numRow);

            for (int i = 0; i < Index.Count; i++) Row.Add(Convert.ToInt32(Triangles[Index[i]].Name.Substring(1, numD)));
            for (int i = 0; i < numRow; i++)
            {
                List<int> list = new List<int>();
                for (int j = Row.IndexOf(i); j < Row.Count; j++)
                {
                    if (Row[j] > i) break;
                    list.Add(j);
                }
                result.Add(list);
            }
            return result;
        }

        [MultiReturn(new[] { "VertexPoints", "Index" })]
        public static Dictionary<string, object> ConvertToMeshComponent(TriangleRigging triangleRigging)
        {
            List<int> indexGroup = new List<int>(triangleRigging.Triangles.Count);
            for (int t = 0; t < triangleRigging.Triangles.Count; t++) for (int p = 0; p < 3; p++) indexGroup.Add(triangleRigging.VertexPoints.IndexOf(triangleRigging.Triangles[t].Points[p]));
            
            return new Dictionary<string, object>
            {
                { "VertexPoints",  triangleRigging.VertexPoints},
                { "Index", indexGroup }
            };
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                Triangles.ForEach(t => t.Dispose());
                VertexPoints.ForEach(v => v.Dispose());
            }
            disposed = true;
        }



    } // end class


}