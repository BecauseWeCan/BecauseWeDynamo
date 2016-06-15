using System;
using System.IO;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;
using Geometry;

namespace Fabrication
{
    /// <summary>
    /// Exports vectors to dxf;
    /// </summary>
    public class DXF
    {
        //**PROPERTIES
        internal Fabrication.DXFLibrary.Document dxf;
        internal Fabrication.DXFLibrary.Tables tables;
        internal Fabrication.DXFLibrary.Table layers;
        internal List<string> layerNames;

        //**CONSTRUCTOR
        internal DXF()
        {
            // dxflibrary document implementation
            dxf = new Fabrication.DXFLibrary.Document();
            tables = new Fabrication.DXFLibrary.Tables();
            dxf.SetTables(tables);
            layers = new Fabrication.DXFLibrary.Table("LAYER");
            tables.addTable(layers);
            layerNames = new List<string>();
            Fabrication.DXFLibrary.Layer layer = new Fabrication.DXFLibrary.Layer("0", 7, "CONTINUOUS");
            layers.AddTableEntry(layer);
            layerNames.Add("0");
        }

        //**CREATE
        /// <summary>
        /// Creats dxf document setup with no entities.
        /// </summary>
        /// <returns>default dxf document with no entities</returns>
        public static DXF ByDefault() { return new DXF(); }

        //**METHODS - **ACTIONS
        /// <summary>
        /// adds planar curve on XY plane to dxf data as arcs and lines
        /// </summary>
        /// <param name="Curve">Curve</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCurveAsArcLineXY(Curve Curve, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            Curve[] Curves = Curve.ApproximateWithArcAndLineSegments();
            AddArcLineXY(Curves, LayerName);
            Curves.ForEach(c => c.Dispose());
        }
        /// <summary>
        /// adds planar curves on XY plane to dxf data as arcs and lines
        /// </summary>
        /// <param name="Curves">Curve Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCurvesAsArcLineXY(Curve[] Curves, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            List<Curve> ArcLineCurves = new List<Curve>();
            Curves.ForEach(c => ArcLineCurves.AddRange(c.ApproximateWithArcAndLineSegments()));
            AddArcLineXY(ArcLineCurves.ToArray(), LayerName);
            ArcLineCurves.ForEach(c => c.Dispose());
        }
        /// <summary>
        /// adds planar curves on XY plane to dxf data as arcs and lines
        /// </summary>
        /// <param name="Curves">Curve Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCurvesAsArcLineXY(Curve[][] Curves, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            List<Curve> ArcLineCurves = new List<Curve>();
            Curves.ForEach(c => ArcLineCurves.AddRange(c.ApproximateWithArcAndLineSegments()));
            AddArcLineXY(ArcLineCurves.ToArray(), LayerName);
            ArcLineCurves.ForEach(c => c.Dispose());
        }
        /// <summary>
        /// adds curves to dxf data as arcs and lines
        /// </summary>
        /// <param name="Curves">Curve Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCurvesAsArcLine(Curve[][] Curves, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            List<Curve> ArcLineCurves = new List<Curve>();
            Curves.ForEach(c => ArcLineCurves.AddRange(c.ApproximateWithArcAndLineSegments()));
            AddArcLineCurves(ArcLineCurves.ToArray(), LayerName);
            ArcLineCurves.ForEach(c => c.Dispose());
        }
        /// <summary>
        /// adds curves to dxf data as arcs and lines
        /// </summary>
        /// <param name="Curves">Curve Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCurvesAsArcLine(Curve[] Curves, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            List<Curve> ArcLineCurves = new List<Curve>();
            Curves.ForEach(c => ArcLineCurves.AddRange(c.ApproximateWithArcAndLineSegments()));
            AddArcLineCurves(ArcLineCurves.ToArray(), LayerName);
            ArcLineCurves.ForEach(c => c.Dispose());
        }
        /// <summary>
        /// adds curve to dxf data as arcs and lines
        /// </summary>
        /// <param name="Curve">Curve</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCurveAsArcLine(Curve Curve, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            Curve[] Curves = Curve.ApproximateWithArcAndLineSegments();
            AddArcLineCurves(Curves, LayerName);
            Curves.ForEach(c => c.Dispose());
        }
        /// <summary>
        /// adds arc to dxf data
        /// </summary>
        /// <param name="Arc">Arc</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddArc(Arc Arc, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            AddArc(Arc, LayerName);
        }
        /// <summary>
        /// adds circles to dxf data
        /// </summary>
        /// <param name="Circles">Circle Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCircles(Circle[] Circles, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            Circles.ForEach(c => AddCircle(c, LayerName));
        }
        /// <summary>
        /// adds circles to dxf data
        /// </summary>
        /// <param name="Circles">Circle Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCircles(Circle[][] Circles, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            Circles.ForEach(c => AddCircle(c, LayerName));
        }
        /// <summary>
        /// adds circle to dxf data
        /// </summary>
        /// <param name="Circle">Circle</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddCircle(Circle Circle, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            AddCircle(Circle, LayerName);
        }
        /// <summary>
        /// adds lines to dxf data
        /// </summary>
        /// <param name="Lines">Line Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddLines(Curve[][] Lines, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            Lines.ForEach(ln => AddLine(ln));
        }
        /// <summary>
        /// adds lines to dxf data
        /// </summary>
        /// <param name="Lines">Line Array</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddLines(Curve[] Lines, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            Lines.ForEach(ln => AddLine(ln));
        }
        /// <summary>
        /// adds lines to dxf data
        /// </summary>
        /// <param name="Line">Line</param>
        /// <param name="LayerName">Layer Name (defaults to 0)</param>
        /// <param name="ACADcolor">Layer Color (defaults to white)</param>
        public void AddLine(Curve Line, string LayerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(LayerName, ACADcolor);
            AddLine(Line, LayerName);
        }
        /// <summary>
        /// writes dxf data to specified file path
        /// will overwrite existing file
        /// </summary>
        /// <param name="FilePath">FilePath (w/o .dxf extension)</param>
        public void WriteFile(string FilePath)
        {
            FileStream f1 = new FileStream(FilePath + ".dxf", FileMode.Create);
            Fabrication.DXFLibrary.Writer.Write(dxf, f1);
            f1.Close();
        }

        //**INTERNAL**METHODS
        internal void CheckLayer(string layerName, short ACADcolor)
        {
            if (!layerNames.Contains(layerName))
            {
                Fabrication.DXFLibrary.Layer layer = new Fabrication.DXFLibrary.Layer(layerName, ACADcolor, "CONTINUOUS");
                layers.AddTableEntry(layer);
                layerNames.Add(layerName);
            }
        }
        internal void AddLine(Curve Line, string LayerName)
        {
            Fabrication.DXFLibrary.Line line = new Fabrication.DXFLibrary.Line(LayerName, Line.StartPoint.X, Line.StartPoint.Y, Line.StartPoint.Z, Line.EndPoint.X, Line.EndPoint.Y, Line.EndPoint.Z);
            dxf.add(line);
        }
        internal void AddCircle(Circle Circle, string LayerName)
        {
            Fabrication.DXFLibrary.Circle circle = new Fabrication.DXFLibrary.Circle(Circle.CenterPoint.X, Circle.CenterPoint.Y, Circle.CenterPoint.Z, Circle.Radius, Circle.Normal.X, Circle.Normal.Y, Circle.Normal.Z, LayerName);
            dxf.add(circle);
        }
        internal void AddArcXY(Arc Arc, string LayerName)
        {
            point c = point.ByPoint(Arc.CenterPoint);
            point s = point.ByPoint(Arc.StartPoint);
            vector A = vector.ByTwoPoints(c, s);
            double startAngle = math.Mod(360 + Math.Sign(A.Y) * math.toDegrees(vector.Xaxis().AngleBetween(A)), 360);
            Fabrication.DXFLibrary.Arc arc = new Fabrication.DXFLibrary.Arc(Arc.CenterPoint.X, Arc.CenterPoint.Y, Arc.Radius, startAngle, startAngle + Arc.SweepAngle, LayerName);
            dxf.add(arc);
        }
        internal void AddArcLineXY(Curve[] Curves, string LayerName)
        {
            for (int i = 0; i < Curves.Length; i++)
            {
                if (Curves[i].Equals(null)) continue;
                Vector stTan = Curves[i].TangentAtParameter(0).Normalized();
                Vector endTan = Curves[i].TangentAtParameter(Curves[i].EndParameter()).Normalized();
                if (stTan.IsParallel(endTan) && stTan.IsAlmostEqualTo(endTan)) AddLine(Curves[i], LayerName);
                else
                {
                    Arc Arc = Arc.ByThreePoints(Curves[i].StartPoint, Curves[i].PointAtParameter(0.5), Curves[i].EndPoint);
                    AddArcXY(Arc, LayerName);
                    Arc.Dispose();
                }
                stTan.Dispose();
                endTan.Dispose();
            }
        }
        internal void AddArc(Arc Arc, string LayerName)
        {
            vector n = vector.ByCoordinates(Arc.Normal.X, Arc.Normal.Y, Arc.Normal.Z);
            point c = point.ByPoint(Arc.CenterPoint);
            point s = point.ByPoint(Arc.StartPoint);
            vector A = vector.ByTwoPoints(c,s);
            vector[] XYZ = vector.NormalizedBasisDXF(n);
            double[] a = A.ChangeBasis(XYZ);
            double[] center = c.ChangeBasis(XYZ);
            double startAngle = math.Mod(360 + Math.Sign(a[1]) * math.toDegrees(XYZ[0].AngleBetween(A)), 360);
            if (A.IsParallel(vector.Xaxis()) && A.X < 0) startAngle += 180;
            double endAngle = math.Mod(startAngle + Arc.SweepAngle, 360);
            Fabrication.DXFLibrary.Arc arc = new Fabrication.DXFLibrary.Arc(center[0], center[1], center[2], Arc.Radius, startAngle, endAngle, Arc.Normal.X, Arc.Normal.Y, Arc.Normal.Z, LayerName);
            dxf.add(arc);
        }
        internal void AddArcLineCurves(Curve[] Curves, string LayerName)
        {
            for (int i = 0; i < Curves.Length; i++)
            {
                if (Curves[i].Equals(null)) continue;
                Vector stTan = Curves[i].TangentAtParameter(0).Normalized();
                Vector endTan = Curves[i].TangentAtParameter(Curves[i].EndParameter()).Normalized();
                if (stTan.IsParallel(endTan) && stTan.IsAlmostEqualTo(endTan)) AddLine(Curves[i], LayerName);
                else
                {
                    Arc Arc = Arc.ByThreePoints(Curves[i].StartPoint, Curves[i].PointAtParameter(0.5), Curves[i].EndPoint);
                    AddArc(Arc, LayerName);
                    Arc.Dispose();
                }
                stTan.Dispose();
                endTan.Dispose();
            }
        }
    }

}