using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;

/*
namespace TriangleRiggingTest
{
    public class TriangleVertex
    {
        private List<Triangle> Triangles;
        private List<TriangleEdge> Edges;
        private Point id;
        private int Spline;
        private TriangleRigging Rig;

        public Point point { get { return id; } }
        public List<Triangle> triangles { get { return Triangles; } }
        public List<TriangleEdge> edges { get { return Edges; } }
        public int splineId
        {
            get { return Spline; }
            set { Spline = value; }
        }

        internal TriangleVertex(Point pt, TriangleRigging rig)
        {
            id = pt;
            Triangles = new List<Triangle>();
            Edges = new List<TriangleEdge>();
            Spline = -1;
            Rig = rig;
        }

        public static TriangleVertex ByPoint(Point pt, TriangleRigging rig)
        {
            return new TriangleVertex(pt, rig);
        }

        public void AddTriangle(Triangle triangle)
        {
            if (!Triangles.Contains(triangle)) Triangles.Add(triangle);
        }

        public void AddEdge(TriangleEdge edge)
        {
            if (!Edges.Contains(edge)) Edges.Add(edge);
        }

        public Boolean AlmostEquals(TriangleVertex vtx)
        {
            return (vtx.point.IsAlmostEqualTo(id));
        }
    }

    public class TriangleEdge
    {
        private TriangleVertex a;
        private TriangleVertex b;
        private List<Triangle> Triangles;
        private string id;
        private Boolean IsOuterEdge;

        public TriangleVertex A { get { return a; } }
        public TriangleVertex B { get { return b; } }
        public string name
        {
            get { return id; }
            set { id = value; }
        }
        public Boolean isOuterEdge
        {
            get { return IsOuterEdge; }
            set { IsOuterEdge = value; }
        }
        public List<Triangle> triangles { get { return Triangles; } }

        internal TriangleEdge(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles, Boolean isOuterEdge)
        {
            a = pta;
            b = ptb;
            id = name;
            Triangles = triangles;
            IsOuterEdge = isOuterEdge;
        }

        internal TriangleEdge(TriangleVertex pta, TriangleVertex ptb)
        {
            a = pta;
            b = ptb;
            id = "";
            Triangles = new List<Triangle>();
            IsOuterEdge = false;
        }

        public static TriangleEdge ByPoints(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles)
        {
            return new TriangleEdge(pta, ptb, name, triangles, false);
        }

        public static TriangleEdge ByPoints(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles, Boolean boo)
        {
            return new TriangleEdge(pta, ptb, name, triangles, boo);
        }

        public void AddTriangle(Triangle t)
        {
            if (!Triangles.Contains(t)) Triangles.Add(t);
        }

        public Boolean AlmostEquals(TriangleEdge e)
        {
            return (e.A.AlmostEquals(a) && e.B.AlmostEquals(b));
        }

        public Boolean HasVertices(TriangleVertex vtx1, TriangleVertex vtx2)
        {
            return ((vtx1.AlmostEquals(a) && vtx2.AlmostEquals(b)) || (vtx1.AlmostEquals(b) && vtx2.AlmostEquals(a)));
        }

    }

    public class Triangle
    {
        private List<Point> Points;
        private List<TriangleVertex> Vertices;
        private List<TriangleEdge> Edges;
        private Surface shape;
        public Boolean HasName;
        private string Name;

        public List<Point> points { get { return Points; } }
        public List<TriangleVertex> vertices { get { return Vertices; } }
        public Surface surface
        {
            get { return shape; }
            set { shape = value; }
        }
        public string name
        {
            get { return Name; }
            set { Name = value; }
        }

        internal Triangle(Surface srf, List<Point> points)
        {
            shape = srf;
            Points = new List<Point>(3);
            Edges = new List<TriangleEdge>(3);
            Vertices = new List<TriangleVertex>(3);
            Name = "";
            HasName = false;
            Points.AddRange(points);
        }

        public static Triangle BySurfaceAndPoints(Surface srf, List<Point> points)
        {
            return new Triangle(srf, points);
        }

        public void AddEdges(TriangleEdge edge)
        {
            if (!Edges.Contains(edge)) Edges.Add(edge);
        }

        public void AddVertex(TriangleVertex vtx)
        {
            if (!Vertices.Contains(vtx)) Vertices.Add(vtx);
        }

        public Boolean HasEdges()
        {
            if (Edges.Count == 3) return true;
            else return false;
        }
    }

    public class TriangleRigging
    {


        //**GLOBAL VARIABLES
        private Dictionary<Point, TriangleVertex> Vertices;
        private List<Point> VertexPoints;
        private List<Triangle> Triangles;
        private List<List<TriangleVertex>> Splines;
        private List<TriangleEdge> Edges;

        //**QUERY
        public List<Triangle> triangles { get { return Triangles; } }
        public List<TriangleEdge> edges { get { return Edges; } }
        public List<TriangleVertex> vertices { get { return Vertices.Values.ToList(); } }
        public List<List<TriangleVertex>> splines { get { return Splines; } }


        //**CONSTRUCTORS
        internal TriangleRigging(Solid[][] Solids, Point[][] SolidsPoint, Point[][] OrderedSplines)
        {
            // initialize Vertices and Triangles
            Vertices = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<Triangle>(Solids.Length);
            Edges = new List<TriangleEdge>();
            Splines = new List<List<TriangleVertex>>(OrderedSplines.Length);

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            for (int i = 0; i < Solids.Length; i++)
            {
                Triangles.Add(new Triangle((Surface)Solids[i][0].Explode()[0], SolidsPoint[i].ToList()));
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    if (!Vertices.ContainsKey(SolidsPoint[i][j]))
                    {
                        Vertices.Add(SolidsPoint[i][j], new TriangleVertex(SolidsPoint[i][j], this));
                    }
                    Vertices[SolidsPoint[i][j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(Vertices[SolidsPoint[i][j]]);
                }
            }

            VertexPoints = Vertices.Keys.ToList();

            //build splines
            for (int i = 0; i < OrderedSplines.Length; i++)
            {
                List<TriangleVertex> Spline = new List<TriangleVertex>();
                int z = 0;
                int edgeCount = 1;
                // iterate through spline vertices
                for (int j = 0; j < OrderedSplines[i].Length; j++)
                {
                    // search vertex points for match
                    for (int k = 0; k < VertexPoints.Count; k++)
                    {
                        if (VertexPoints[k].IsAlmostEqualTo(OrderedSplines[i][j]))
                        {
                            Spline.Add(Vertices[VertexPoints[k]]);
                            Vertices[VertexPoints[k]].splineId = i;

                            if (j > 0)
                            {
                                // check triangles
                                // find triangles on edge
                                List<Triangle> listT = new List<Triangle>();
                                for (int t = 0; t < Vertices[VertexPoints[k]].triangles.Count; t++)
                                {
                                    if (Vertices[VertexPoints[z]].triangles.Contains(Vertices[VertexPoints[k]].triangles[t]))
                                        listT.Add(Vertices[VertexPoints[k]].triangles[t]);
                                }

                                if (listT.Count > 0)
                                {
                                    string edgeId = "s" + i + "-" + edgeCount;
                                    TriangleEdge e = new TriangleEdge(Vertices[VertexPoints[z]], Vertices[VertexPoints[k]], edgeId, listT, false);
                                    if (listT.Count == 1) e.isOuterEdge = true;
                                    Edges.Add(e);
                                    // add edge key to vertices
                                    for (int t = 0; t < listT.Count; t++) { listT[t].AddEdges(e); }
                                    Vertices[VertexPoints[z]].AddEdge(e);
                                    Vertices[VertexPoints[k]].AddEdge(e);
                                    edgeCount++;
                                }
                            }
                            z = k;
                            break;
                        }
                    } // end search
                } // end spline point iteration
                // check triangles
                // find triangles on edge
                List<Triangle> listTlast = new List<Triangle>();
                for (int t = 0; t < Spline[0].triangles.Count; t++)
                {
                    if (Vertices[VertexPoints[z]].triangles.Contains(Spline[0].triangles[t]))
                        listTlast.Add(Spline[0].triangles[t]);
                }

                if (listTlast.Count > 0)
                {
                    string edgeId = "s" + i + "-" + edgeCount;
                    TriangleEdge e = new TriangleEdge(Vertices[VertexPoints[z]], Spline[0], edgeId, listTlast, false);
                    if (listTlast.Count == 1) e.isOuterEdge = true;
                    Edges.Add(e);
                    // add edge key to vertices
                    for (int t = 0; t < listTlast.Count; t++) { listTlast[t].AddEdges(e); }
                    Vertices[VertexPoints[z]].AddEdge(e);
                    Spline[0].AddEdge(e);
                    edgeCount++;
                }
                Splines.Add(Spline);
            } // end build spline


            //build edges


            for (int s = 0; s < Splines.Count - 1; s++)
            {
                // for each spline
                int edgeCount = 0;
                for (int v = 0; v < Splines[s].Count; v++)
                {
                    List<TriangleVertex> Search = new List<TriangleVertex>();
                    for (int t = 0; t < Splines[s][v].triangles.Count; t++)
                    {
                        if (!Splines[s][v].triangles[t].HasEdges())
                        {
                            for (int vt = 0; vt < Splines[s][v].triangles[t].vertices.Count; vt++)
                            {
                                if (!Search.Contains(Splines[s][v].triangles[t].vertices[vt]) && Splines[s][v].triangles[t].vertices[vt].splineId == s + 1)
                                    Search.Add(Splines[s][v].triangles[t].vertices[vt]);
                            }
                        }
                    }
                    Search.Remove(Splines[s][v]);

                    for (int sv = 0; sv < Search.Count; sv++)
                    {
                        // check triangles
                        // find triangles on edge
                        List<Triangle> listT = new List<Triangle>();
                        for (int t = 0; t < Search[sv].triangles.Count; t++)
                        {
                            if (Splines[s][v].triangles.Contains(Search[sv].triangles[t]))
                                listT.Add(Search[sv].triangles[t]);
                        }

                        if (listT.Count > 0)
                        {
                            string edgeId = "e" + s.ToString() + (s + 1).ToString() + "-" + edgeCount;
                            TriangleEdge e = new TriangleEdge(Splines[s][v], Search[sv], edgeId, listT, false);
                            if (listT.Count == 1) e.isOuterEdge = true;
                            Edges.Add(e);
                            // add edge key to vertices
                            for (int t = 0; t < listT.Count; t++) { listT[t].AddEdges(e); }
                            Splines[s][v].AddEdge(e);
                            Search[sv].AddEdge(e);
                            edgeCount++;
                        }
                    }


                } // end spline

            } // end build edges 



        }

        //**CREATE
        public static TriangleRigging BySolidsPointsPoints(Solid[][] Solids, Point[][] SolidsPoints, Point[][] Splines)
        {
            return new TriangleRigging(Solids, SolidsPoints, Splines);
        }

        //**ACTION

    }
}
*/