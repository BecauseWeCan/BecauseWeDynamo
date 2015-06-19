﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;

namespace Fabrication
{
    public class TriangleEdgeConnector
    {
        double radius;
        double spacing;
        public int[] numHoles { get; set;}
        public TriangleEdge edge { get; set;}

        internal TriangleEdgeConnector(TriangleEdge tEdge, double radiusHole, double spacingHole)
        {
            spacing = spacingHole;
            radius = radiusHole;
            edge = tEdge;
            numHoles = new int[3]{0,0,0};
            for (int i = 0; i < edge.triangles.Count; i++)
            {
                numHoles[i] = (int) (Math.Ceiling((edge.triangles[i].surface.ClosestPointTo(edge.midpoint).DistanceTo(edge.midpoint) + 3 * radius) / spacing)) + edge.triangles[i].numHoles - 1;
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

    public class TriangleSpline
    {
        //**PROPERTIES
        List<TriangleEdge> edges;
        Dictionary<Point, TriangleVertex> vtxs;
        List<Point> pts;
        Point[] stpts;

        //**QUERY**PROPERTIES
        public List<TriangleEdge> Edges { get { return edges; } }
        public List<Point> Points { get { return pts; } }

        //**CONSTRUCTOR
        internal TriangleSpline(Curve[] Curves, Point[] StartPoints)
        {
            edges = new List<TriangleEdge>(Curves.Length);
            vtxs = new Dictionary<Point, TriangleVertex>();
            stpts = StartPoints;

            for (int i =0; i<Curves.Length; i++)
            {
                Point stpt = Curves[i].StartPoint;
                Point endpt = Curves[i].EndPoint;
                bool hasStart = false;
                bool hasEnd = false;
                for (int j = 0; j < vtxs.Keys.ToList().Count; j++)
                {
                    if (vtxs.Keys.ToList()[j].IsAlmostEqualTo(Curves[i].StartPoint))
                    {
                        stpt = vtxs.Keys.ToList()[j];
                        hasStart = true;
                    }
                    if (vtxs.Keys.ToList()[j].IsAlmostEqualTo(Curves[i].EndPoint))
                    {
                        endpt = vtxs.Keys.ToList()[j];
                        hasEnd = true;
                    }
                    if (hasEnd && hasStart) break;
                }
                if (!hasStart) vtxs.Add(stpt, new TriangleVertex(stpt));
                if (!hasEnd) vtxs.Add(endpt, new TriangleVertex(endpt));
                TriangleEdge e = new TriangleEdge(vtxs[stpt], vtxs[endpt], "e" + i.ToString("D3"), new List<Triangle>(), true);
                vtxs[stpt].AddEdge(e);
                vtxs[endpt].AddEdge(e);
                edges.Add(e);
            }
            for (int i = 0; i < StartPoints.Length; i++)
            {
                bool inList = false;
                for (int j = 0; j < vtxs.Keys.ToList().Count; j++)
                {
                    if (vtxs.Keys.ToList()[j].IsAlmostEqualTo(StartPoints[i]))
                    {
                        inList = true;
                        break;
                    }
                }
                if (!inList) vtxs.Add(StartPoints[i], new TriangleVertex(StartPoints[i]));

            }
            pts = vtxs.Keys.ToList();
        }

        //**CREATE
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Circles"></param>
        /// <param name="StartPoints"></param>
        /// <returns></returns>
        public static TriangleSpline ByCurvesAndPoints(Curve[] Curves, Point[] StartPoints)
        {
            return new TriangleSpline(Curves, StartPoints);
        }

        //**ACTIONS*METHODS
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public List<List<Point>> GetSplines()
        {
            List<List<Point>> result = new List<List<Point>>(stpts.Length);
            for (int i = 0; i < stpts.Length; i++)
            {
                List<Point> spline = new List<Point>();
                for (int n = 0; n < pts.Count; n++)
                {
                    if (pts[n].IsAlmostEqualTo(stpts[i]))
                    {
                        spline.Add(pts[n]);
                        break;
                    }
                }
                if (vtxs[spline[0]].edges.Count == 0) 
                {
                    result.Add(spline);
                    continue;
                }

                bool added = true;
                while (added)
                {
                    if (vtxs[spline[spline.Count - 1]].edges.Count < 1) break;
                    added = false;
                    for (int k = 0; k < vtxs[spline[spline.Count - 1]].edges.Count; k++)
                    {
                        if (!spline.Contains(vtxs[spline[spline.Count - 1]].edges[k].A.point))
                        {
                            spline.Add(vtxs[spline[spline.Count - 1]].edges[k].A.point);
                            added = true;
                            break;
                        }
                        if (!spline.Contains(vtxs[spline[spline.Count - 1]].edges[k].B.point))
                        {
                            spline.Add(vtxs[spline[spline.Count - 1]].edges[k].B.point);
                            added = true;
                            break;
                        }
                    }
                }
                result.Add(spline);
            }
            return result;
        }

    }

    public class TriangleVertex
    {
        private List<Triangle> Triangles;
        private List<TriangleEdge> Edges;
        private Point id;
        private int Spline;

        //**QUERY**PROPERTIES
        /// <summary>
        /// point of vertex
        /// </summary>
        public Point point {get { return id; } }
        /// <summary>
        /// triangle with this vertex
        /// </summary>
        public List<Triangle> triangles { get { return Triangles; } }
        /// <summary>
        /// Edges with this vertex
        /// </summary>
        public List<TriangleEdge> edges { get { return Edges; } }
        /// <summary>
        /// spline location of vertex
        /// </summary>
        public int splineId { 
            get { return Spline; }
            set { Spline = value; }
        }

        //*CONSTRUCTOR
        internal TriangleVertex(Point pt)
        {
            id = pt;
            Triangles = new List<Triangle>();
            Edges = new List<TriangleEdge>();
            Spline = -1;
        }

        //**CREATE
        /// <summary>
        /// make vertex from point
        /// </summary>
        /// <param name="pt">point</param>
        /// <returns></returns>
        public static TriangleVertex ByPoint(Point pt) { return new TriangleVertex(pt); }

        //**ACTIONS**METHODS
        /// <summary>
        /// adds Triangle information to vertex
        /// </summary>
        /// <param name="triangle">Triangle object to be added</param>
        public void AddTriangle(Triangle triangle)
        {
            if (!Triangles.Contains(triangle)) Triangles.Add(triangle);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="edge"></param>
        public void AddEdge(TriangleEdge edge)
        {
            if (!Edges.Contains(edge)) Edges.Add(edge);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="vtx"></param>
        /// <returns></returns>
        public Boolean AlmostEquals(TriangleVertex vtx)
        {
            return (vtx.point.IsAlmostEqualTo(id));
        }
    }

    public class TriangleEdge
    {
        private TriangleVertex a;
        private TriangleVertex b;
        private Point Midpoint;
        private List<Triangle> Triangles;
        private string id;
        private Boolean IsOuterEdge;
        private List<TriangleVertex> Vertices;

        //**QUERY**PROPERTIES
        public Vector Normal { get; set; }
        public TriangleVertex A { get { return a; } }
        public TriangleVertex B { get { return b; } }
        public List<Triangle> triangles { get { return Triangles; } }
        public Point midpoint { get { return Midpoint;  } }
        public string name { 
            get { return id; }
            set { id = value; }
        }
        public Boolean isOuterEdge
        {
            get { return IsOuterEdge; }
            set { IsOuterEdge = value; }
        }

        //*CONSTRUCTOR
        internal TriangleEdge(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles, Boolean isOuterEdge)
        {
            Vertices = new List<TriangleVertex>(2);
            a = pta;
            b = ptb;
            Midpoint = Point.ByCoordinates((pta.point.X + ptb.point.X) / 2, (pta.point.Y + ptb.point.Y) / 2, (pta.point.Z + ptb.point.Z) / 2);
            id = name;
            Triangles = triangles;
            IsOuterEdge = isOuterEdge;
            Vertices.Add(pta);
            Vertices.Add(ptb);
            Normal = Vector.ByCoordinates(0,0,0);
        }

        //**CREATE
        /// <summary>
        /// 
        /// </summary>
        /// <param name="pta"></param>
        /// <param name="ptb"></param>
        /// <param name="name"></param>
        /// <param name="Triangles"></param>
        /// <returns></returns>
        public static TriangleEdge ByPoints(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles) { return new TriangleEdge(pta, ptb, name, triangles, false); }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="pta"></param>
        /// <param name="ptb"></param>
        /// <param name="name"></param>
        /// <param name="Triangles"></param>
        /// <param name="boo"></param>
        /// <returns></returns>
        public static TriangleEdge ByPoints(TriangleVertex pta, TriangleVertex ptb, string name, List<Triangle> triangles, Boolean boo) { return new TriangleEdge(pta, ptb, name, triangles, boo); }

        //**ACTIONS**METHODS
        /// <summary>
        /// adds triangle reference data to internal list
        /// </summary>
        /// <param name="t"></param>
        public void AddTriangle(Triangle t)
        {
            if (!Triangles.Contains(t)) Triangles.Add(t);
        }
        /// <summary>
        /// compare this edge to another
        /// </summary>
        /// <param name="e">edge</param>
        /// <returns></returns>
        public Boolean AlmostEquals(TriangleEdge e) { return (e.A.AlmostEquals(a) && e.B.AlmostEquals(b)); }
        /// <summary>
        /// sees if input vertices are edge endpoints
        /// </summary>
        /// <param name="vtx1">first vertex</param>
        /// <param name="vtx2">second vertex</param>
        /// <returns></returns>
        public Boolean HasVertices(TriangleVertex vtx1, TriangleVertex vtx2) { return ((vtx1.AlmostEquals(a) && vtx2.AlmostEquals(b)) || (vtx1.AlmostEquals(b) && vtx2.AlmostEquals(a))); }
        public Line GetEdgeGeometry()
        {
            return Line.ByStartPointEndPoint(A.point, B.point);
        }
    }

    public class Triangle
    {
        //*PRIVATE PROPERTIES
        public Dictionary<string, Object> Parameters;
        private int holes;
        private double Height;
        private int SplineId;
        private List<Point> Points;
        private List<TriangleVertex> Vertices;
        private List<TriangleEdge> Edges;
        private Surface shape;
        public Boolean HasName;
        private string Name;
        private Point Circumcenter;
        private Point Center;
        private CoordinateSystem CS; 

        //**QUERY**PROPERTIES
        public CoordinateSystem contextCoordinateSystem { get { return CS; } }
        public List<TriangleEdge> edges { get { return Edges; } }
        public List<Point> points { get { return Points; } }
        public List<TriangleVertex> vertices { get { return Vertices; } }
        public Point circumcenter { get { return Circumcenter; } }
        public Point center { get { return Center; } }
        public Curve[] PerimeterCurves{get {return shape.PerimeterCurves(); } }
        public Vector Normal { get; set; }
        public Surface surface { 
            get { return shape; }
            set { shape = value; }
        }
        public string name
        {
            get { return Name; }
            set { Name = value; }
        }
        public int splineId { 
            get { return SplineId; } 
            set { SplineId = value; }
        }
        public int numHoles
        {
            get { return holes; }
            set { holes = value; }
        }

        //**CONSTRUCTOR
        internal Triangle(Surface srf, List<Point> points)
        {
            shape = srf;
            Points = new List<Point>(3);
            Edges = new List<TriangleEdge>(3);
            Vertices = new List<TriangleVertex>(3);
            Name = "";
            HasName = false;
            Height = 0;
            holes = 1;

            Points.AddRange(points);
            Normal = srf.NormalAtParameter();
            Circumcenter = Circle.ByThreePoints(points[0], points[1], points[2]).CenterPoint;
            Center = Point.ByCoordinates( (points[0].X + points[1].X + points[2].X)/3.0, (points[0].Y + points[1].Y + points[2].Y)/3.0, (points[0].Z + points[1].Z + points[2].Z)/3.0);
            CS = CoordinateSystem.ByOriginVectors(Center, Vector.ByTwoPoints(points[0], points[1]), Vector.ByTwoPoints(points[0], points[2]));
        }

        //**CREATE**
        /// <summary>
        /// 
        /// </summary>
        /// <param name="srf"></param>
        /// <param name="points"></param>
        /// <returns></returns>
        public static Triangle BySurfaceAndPoints(Surface srf, List<Point> points) { return new Triangle(srf,points); }

        //**ACTIONS**METHODS
        /// <summary>
        /// build internal geometry
        /// </summary>
        /// <param name="edge">base edge</param>
        public void BuildGeometricProperties(TriangleEdge edge)
        {
            int a = vertices.IndexOf(edge.A);
            int b = vertices.IndexOf(edge.B);
            int c = 3 - (a + b);
            if (((a - b) % 3 + 3) % 3 == 1)  CS = CoordinateSystem.ByOriginVectors(center, Vector.ByTwoPoints(Points[b], Points[a]), Vector.ByTwoPoints(Points[b], Points[c]));
            else if (((a - b) % 3 + 3) % 3 == 2) CS = CoordinateSystem.ByOriginVectors(center, Vector.ByTwoPoints(Points[a], Points[b]), Vector.ByTwoPoints(Points[a], Points[c]));
            double w = points[a].DistanceTo(points[b]);
            double sa = points[c].DistanceTo(points[a]);
            double sb = points[c].DistanceTo(points[b]);
            double s = (w + sa + sb) / 2;
            Height = 2 * Math.Sqrt(s * (s - w) * (s - sa) * (s - sb)) / w;
        }
        /// <summary>
        /// adjusts CS
        /// </summary>
        public void BuildCS()
        {
            for (int i = 0; i < edges.Count; i++)
            {
                if (edges[i].name.Substring(1, 2).Equals(Name.Substring(3, 2)))
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
        /// gets edge
        /// </summary>
        /// <returns>boolean</returns>
        public Boolean HasEdges()
        {
            if (Edges.Count == 3) return true;
            else return false;
        }
        /// <summary>
        /// adds circles to Edges
        /// </summary>
        /// <param name="factor">text scale</param>
        /// <returns></returns>
        public List<PolyCurve> GetEdgeLabels(double factor)
        {
            List<PolyCurve> labels = new List<PolyCurve>();
            for (int j = 0; j < edges.Count; j++)
                {
                    if (edges[j].isOuterEdge) continue;

                    int a = vertices.IndexOf(edges[j].A);
                    int b = vertices.IndexOf(edges[j].B);

                    if (((a - b) % 3 + 3) % 3 == 1)
                    {
                        labels.AddRange(Word.ByString(
                                    edges[j].name,
                                    surface.ClosestPointTo(edges[j].midpoint),
                                    Vector.ByTwoPoints(edges[j].midpoint, edges[j].B.point),
                                    Vector.ByTwoPoints(center, edges[j].midpoint)
                                    ).display(factor));
                    }
                    else if (((a - b) % 3 + 3) % 3 == 2)
                    {
                        labels.AddRange(Word.ByString(
                                    edges[j].name, //string
                                    surface.ClosestPointTo(edges[j].midpoint), //cs point
                                    Vector.ByTwoPoints(edges[j].midpoint, edges[j].A.point), // X-axis
                                    Vector.ByTwoPoints(center, edges[j].midpoint) // Y-axis
                                    ).display(factor));
                    } // end switch
                } // end Edges

            return labels;
        }// end action
        /// <summary>
        /// adds circles to Edges
        /// </summary>
        /// <param name="factor">text scale</param>
        /// <returns></returns>
        public List<Circle> PlaceHoles(int holes, double radius, double spacing)
        {
            List<Circle> circles = new List<Circle>();
            for (int j = 0; j < edges.Count; j++)
            {
                if (edges[j].isOuterEdge) continue;

                int a = vertices.IndexOf(edges[j].A);
                int b = vertices.IndexOf(edges[j].B);

                if (((a - b) % 3 + 3) % 3 == 1)
                {
                    double offset = spacing * Math.Ceiling( (surface.ClosestPointTo(edges[j].midpoint).DistanceTo(edges[j].midpoint) + 3*radius)/spacing );
                    for (int k = 0; k < holes; k++)
                    {
                        circles.Add( Circle.ByCenterPointRadiusNormal(
                            (Point)edges[j].midpoint.Translate(Vector.ByTwoPoints(edges[j].midpoint, surface.ClosestPointTo(edges[j].midpoint)), offset), 
                            radius,
                            Normal
                            ) );
                        offset += spacing;
                    }
                }
                else if (((a - b) % 3 + 3) % 3 == 2)
                {
                    double offset = spacing * Math.Ceiling((surface.ClosestPointTo(edges[j].midpoint).DistanceTo(edges[j].midpoint) + 3 * radius) / spacing);
                    for (int k = 0; k < holes; k++)
                    {
                        circles.Add( Circle.ByCenterPointRadiusNormal(
                            (Point)edges[j].midpoint.Translate(Vector.ByTwoPoints(edges[j].midpoint, surface.ClosestPointTo(edges[j].midpoint)), offset), 
                            radius,
                            Normal
                            ) );
                        offset += spacing;
                    }
                } // end switch
            } // end Edges
            return circles;
        }// end action

        /// <summary>
        /// get profiles and letters
        /// </summary>
        /// <param name="factor">letter size</param>
        /// <returns>curves</returns>
        public List<Curve> GetTriangleCurves(double factor)
        {
            List<Curve> curves = new List<Curve>();
            List<PolyCurve> letters =  GetEdgeLabels(factor);
            curves.AddRange(letters);
            curves.AddRange(surface.PerimeterCurves());
            return curves;
        }


    }

    public class TriangleRigging
    {
        //**CLASS VARIABLES
        private Dictionary<Point, TriangleVertex> Vertices;
        public List<Point> VertexPoints { get; set; }
        public List<Triangle> Triangles { get; set; }
        public List<List<TriangleVertex>> Splines { get; set; }
        public List<TriangleEdge> Edges {get; set;}

        public List<TriangleVertex> vertices{ get { return Vertices.Values.ToList(); } }



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
                Triangles.Add(new Triangle((Surface) Solids[i][0].Explode()[0], SolidsPoint[i].ToList()));
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    if (!Vertices.ContainsKey(SolidsPoint[i][j]))
                    {
                        Vertices.Add(SolidsPoint[i][j], new TriangleVertex(SolidsPoint[i][j]));
                    }
                    Vertices[SolidsPoint[i][j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(Vertices[SolidsPoint[i][j]]);
                }
            }

            VertexPoints = Vertices.Keys.ToList();

            //build splines
            int numD = OrderedSplines.Length.ToString().Length;
            for (int i = 0; i < OrderedSplines.Length; i++)
            {
                List<TriangleVertex> Spline = new List<TriangleVertex>();
                int z = 0;
                int edgeCount = 1;
                int numDs = OrderedSplines[i].Length.ToString().Length;
                // iterate through spline vertices
                for (int j = 0; j < OrderedSplines[i].Length; j++ )
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
                                // check Triangles
                                // find Triangles on edge
                                List<Triangle> listT = new List<Triangle>();
                                for (int t = 0; t < Vertices[VertexPoints[k]].triangles.Count; t++)
                                {
                                    if (Vertices[VertexPoints[z]].triangles.Contains(Vertices[VertexPoints[k]].triangles[t]))
                                        listT.Add(Vertices[VertexPoints[k]].triangles[t]);
                                }

                                if (listT.Count > 0)
                                {
                                    string edgeId = "s" + i.ToString("D"+ numD) + "-" + edgeCount.ToString("D"+numDs);
                                    TriangleEdge e = new TriangleEdge(Vertices[VertexPoints[z]], Vertices[VertexPoints[k]], edgeId, listT, false);
                                    if (listT.Count == 1) e.isOuterEdge = true;
                                    Edges.Add(e);
                                    // add edge key to vertices
                                    for (int t = 0; t < listT.Count; t++) { 
                                        listT[t].AddEdge(e); listT[t].BuildGeometricProperties(e);
                                        e.Normal = e.Normal.Add(listT[t].Normal);
                                    }
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
                // check Triangles
                // find Triangles on edge
                if (OrderedSplines[i].Length > 2)
                {
                    List<Triangle> listTlast = new List<Triangle>();
                    for (int t = 0; t < Spline[0].triangles.Count; t++)
                    {
                        if (Vertices[VertexPoints[z]].triangles.Contains(Spline[0].triangles[t]))
                            listTlast.Add(Spline[0].triangles[t]);
                    }

                    if (listTlast.Count > 0)
                    {
                        string edgeId = "s" + i.ToString("D" + numD) + "-" + edgeCount.ToString("D" + numDs);
                        TriangleEdge e = new TriangleEdge(Vertices[VertexPoints[z]], Spline[0], edgeId, listTlast, false);
                        if (listTlast.Count == 1) e.isOuterEdge = true;
                        Edges.Add(e);
                        // add edge key to vertices
                        for (int t = 0; t < listTlast.Count; t++) { 
                            listTlast[t].AddEdge(e);
                            e.Normal = e.Normal.Add(listTlast[t].Normal);
                        }
                        Vertices[VertexPoints[z]].AddEdge(e);
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
                for (int v = 0; v < Splines[s].Count; v++ )
                {
                    List<TriangleVertex> Search = new List<TriangleVertex>();
                    for (int t = 0; t < Splines[s][v].triangles.Count; t++ )
                    {
                        if (!Splines[s][v].triangles[t].HasEdges())
                        {
                            for (int vt = 0; vt < Splines[s][v].triangles[t].vertices.Count; vt++ )
                            {
                                if ( !Search.Contains(Splines[s][v].triangles[t].vertices[vt]) && Splines[s][v].triangles[t].vertices[vt].splineId == s+1 ) 
                                    Search.Add(Splines[s][v].triangles[t].vertices[vt]);
                            }
                        }
                    }
                    Search.Remove(Splines[s][v]);

                    for (int sv = 0; sv < Search.Count; sv++ )
                    {
                        // check Triangles
                        // find Triangles on edge
                        List<Triangle> listT = new List<Triangle>();
                        for (int t = 0; t < Search[sv].triangles.Count; t++)
                        {
                            if (Splines[s][v].triangles.Contains(Search[sv].triangles[t]))
                                listT.Add(Search[sv].triangles[t]);
                        }

                        if (listT.Count > 0)
                        {
                            string edgeId = "e" + s.ToString("D" + numD) + (s + 1).ToString("D" + numD) + "-" + edgeCount.ToString("D"+numDe);
                            TriangleEdge e = new TriangleEdge(Splines[s][v], Search[sv], edgeId, listT, false);
                            if (listT.Count == 1) e.isOuterEdge = true;
                            Edges.Add(e);
                            // add edge key to vertices
                            for (int t = 0; t < listT.Count; t++) { 
                                listT[t].AddEdge(e);
                                e.Normal = e.Normal.Add(listT[t].Normal);
                                if (!listT[t].HasName)
                                {
                                    listT[t].name = "t" + s.ToString("D" + numD) + (s + 1).ToString("D" + numD) + "-" + triangleCount.ToString("D" + numDe);
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
        internal TriangleRigging(Solid[][] Solids, Point[][] SolidsPoint, Point start)
        {
             // initialize Vertices and Triangles
            Vertices = new Dictionary<Point, TriangleVertex>();
            Triangles = new List<Triangle>(Solids.Length);
            Edges = new List<TriangleEdge>();

            // iterate over Solids and assign topology to Vertices
            // and create Triangles based on Solids
            int length = Solids.Length.ToString().Length;
            int eCount = 1;
            List<TriangleEdge> HalfEdges = new List<TriangleEdge>();
            for (int i = 0; i < Solids.Length; i++)
            {
                Triangles.Add(new Triangle((Surface) Solids[i][0].Explode()[0], SolidsPoint[i].ToList()));
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    if (!Vertices.ContainsKey(SolidsPoint[i][j]))
                    {
                        Vertices.Add(SolidsPoint[i][j], new TriangleVertex(SolidsPoint[i][j]));
                    }
                    Vertices[SolidsPoint[i][j]].AddTriangle(Triangles[i]);
                    Triangles[i].AddVertex(Vertices[SolidsPoint[i][j]]);
                }
                for (int j = 0; j < SolidsPoint[i].Length; j++)
                {
                    bool addEdge = true;
                    int k = (j+1)%(SolidsPoint[i].Length);
                    TriangleEdge e = new TriangleEdge(Vertices[SolidsPoint[i][j]], Vertices[SolidsPoint[i][k]], "e" + eCount.ToString("D" + length), new List<Triangle> { Triangles[i] }, true);
                    int remove = -1;
                    for (int n = 0; n < HalfEdges.Count; n ++)
                    {
                        if (HalfEdges[n].HasVertices(e.A, e.B))
                        {
                            Triangles[i].AddEdge(HalfEdges[n]);
                            HalfEdges[n].AddTriangle(Triangles[i]);
                            HalfEdges[n].isOuterEdge = false;
                            addEdge = false;
                            remove = n;
                            break;
                        }
                    }
                    if (addEdge)
                    {
                        Edges.Add(e);
                        HalfEdges.Add(e);
                        Triangles[i].AddEdge(e);
                        Vertices[SolidsPoint[i][j]].AddEdge(e);
                        Vertices[SolidsPoint[i][k]].AddEdge(e);
                        eCount++;
                    }
                    else
                    {
                        HalfEdges.RemoveAt(remove);
                    }

                }
            }
        }
       
        //**CREATE
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Solids"></param>
        /// <param name="SolidsPoints"></param>
        /// <param name="Splines"></param>
        /// <returns></returns>
        public static TriangleRigging BySolidsPointsPoints(Solid[][] Solids, Point[][] SolidsPoints, Point[][] Splines)
        {
            return new TriangleRigging(Solids, SolidsPoints, Splines);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Solids"></param>
        /// <param name="SolidsPoints"></param>
        /// <param name="start"></param>
        /// <returns></returns>
        public static TriangleRigging BySolidsPointsPoints(Solid[][] Solids, Point[][] SolidsPoints, Point start)
        {
            return new TriangleRigging(Solids, SolidsPoints, start);
        }

        //**ACTION
        /// <summary>
        /// circles Edges
        /// </summary>
        /// <param name="factor">text scale</param>
        /// <returns>list of circles ordered by triangle then edge then letter</returns>
        public List<List<PolyCurve>> GetEdgeLabels(double factor)
        {
            List<List<PolyCurve>> labels = new List<List<PolyCurve>>(Triangles.Count);
            for (int i = 0; i < Triangles.Count; i++ ) { labels.Add(Triangles[i].GetEdgeLabels(factor));  }
            return labels;
        }
        /// <summary>
        /// gets index set of Triangles sorted by angles
        /// </summary>
        /// <returns></returns>
        public List<int> GetSortedTriangleIndicesByName()
        {
            List<int> index = new List<int>(Triangles.Count);
            Dictionary<int, string> sorted = new Dictionary<int, string>();
            for (int i = 0; i < Triangles.Count; i++) { sorted.Add(i, Triangles[i].name); }
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
                if (Edges[i].isOuterEdge) continue;
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
        public void AddParamterToTriangles(string name, Object[] data )
        {
            if (data.Length == Triangles.Count) for (int i = 0; i < Triangles.Count; i++) { Triangles[i].Parameters.Add(name, data[i]); }
        }
        public List<int> GetTriangleIndexByParameterValue( string parameter, Object value )
        {
            if (!Triangles[0].Parameters.Keys.Contains(parameter)) return null;
            List<int> result = new List<int>();
            for (int i = 0; i < Triangles.Count; i++ )
            {
                if ( Triangles[i].Parameters[parameter].Equals(value) )  ;
            }
                return result;
        }
    } // end class


}
