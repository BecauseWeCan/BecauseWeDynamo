using System;
using System.IO;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;

namespace Utilities
{
    /// <summary>
    /// Exports vectors to dxf;
    /// XY: all vectorwork.
    /// XYZ: vectorwork composed of lines only.
    /// </summary>
    public class ExportToDXF
    {
        //**PROPERTIES
        internal DXFLibrary.Document dxf;
        internal DXFLibrary.Tables tables;
        internal DXFLibrary.Table layers;
        internal List<string> layerNames;

        //**CONSTRUCTOR
        internal ExportToDXF()
        {
            // dxflibrary document implementation
            dxf = new DXFLibrary.Document();
            tables = new DXFLibrary.Tables();
            dxf.SetTables(tables);
            layers = new DXFLibrary.Table("LAYER");
            tables.addTable(layers);
            layerNames = new List<string>();
            DXFLibrary.Layer layer = new DXFLibrary.Layer("0", 7, "CONTINUOUS");
            layers.AddTableEntry(layer);
            layerNames.Add("0");
        }

        //**INTERNAL**METHODS
        internal void CheckLayer(string layerName, short ACADcolor)
        {
            if (!layerNames.Contains(layerName))
            {
                DXFLibrary.Layer layer = new DXFLibrary.Layer(layerName, ACADcolor, "CONTINUOUS");
                layers.AddTableEntry(layer);
                layerNames.Add(layerName);
            }
        }

        //**CREATE
        public static ExportToDXF ByDefault() { return new ExportToDXF(); }

        //**METHODS - **ACTIONS
        public ExportToDXF AddPolyCurvesAsArcsLines(PolyCurve[][] polycurves, string layerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(layerName, ACADcolor);
            List<Curve> Curves = new List<Curve>();
            for (int i = 0; i < polycurves.Length; i++) for (int j = 0; j < polycurves[i].Length; j++) for (int k = 0; k < polycurves[i][j].Curves().Length; k++)
                        Curves.AddRange(polycurves[i][j].Curves()[k].ApproximateWithArcAndLineSegments());
            Vector stTan = null;
            Vector endTan = null;
            Point stpt = null;
            Point endpt = null;
            for (int i = 0; i < Curves.Count; i++)
            {
                stTan = Curves[i].TangentAtParameter(0);
                endTan = Curves[i].TangentAtParameter(Curves[i].EndParameter());
                stpt = Curves[i].StartPoint;
                endpt = Curves[i].EndPoint;

                if (stTan.IsParallel(endTan))
                {
                    DXFLibrary.Line line = new DXFLibrary.Line(layerName, stpt.X, stpt.Y, endpt.X, endpt.Y);
                    dxf.add(line);
                }
                else
                {
                    Arc arc = Arc.ByStartPointEndPointStartTangent(stpt, endpt, stTan);
                    arc = Arc.ByCenterPointStartPointEndPoint(arc.CenterPoint, endpt, stpt);
                    CoordinateSystem cs = arc.CoordinateSystemAtParameter(0);
                    double rotate = cs.XAxis.AngleBetween(Vector.XAxis());
                    double signY = cs.XAxis.Y;
                    double signX = cs.XAxis.X;
                    double startAngle = rotate;
                    if (signY < 0) startAngle = -rotate;
                    if (signX > 0) startAngle = 180 - startAngle;
                    DXFLibrary.Arc arcs = new DXFLibrary.Arc(arc.CenterPoint.X, arc.CenterPoint.Y, arc.Radius, -startAngle - arc.SweepAngle, -startAngle, layerName);
                    dxf.add(arcs);
                    arc.Dispose();
                    cs.Dispose();
                }
            }
            if (stTan != null) stTan.Dispose();
            if (endTan != null) endTan.Dispose();
            if (stpt != null) stpt.Dispose();
            if (endpt != null) endpt.Dispose();
            Curves.ForEach(c => c.Dispose());
            return this;
        }
        public ExportToDXF AddCircles(Circle[][] circles, string layerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(layerName, ACADcolor);
            for (int i = 0; i < circles.Length; i++) for (int j = 0; j < circles[i].Length; j++)
            {
                DXFLibrary.Circle circle = new DXFLibrary.Circle(circles[i][j].CenterPoint.X, circles[i][j].CenterPoint.Y, circles[i][j].Radius, layerName);
                dxf.add(circle);
            }
            return this;
        }
        public ExportToDXF AddPolyCurvesAsLines(PolyCurve[][] polycurves, string layerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(layerName, ACADcolor);
            List<Curve> Curves = new List<Curve>();
            for (int i = 0; i < polycurves.Length; i++) for (int j = 0; j < polycurves[i].Length; j++) for (int k = 0; k < polycurves[i][j].Curves().Length; k++)
                        Curves.AddRange(polycurves[i][j].Curves()[k].ApproximateWithArcAndLineSegments());
            for (int i = 0; i < Curves.Count; i++)
            {
                DXFLibrary.Line line = new DXFLibrary.Line("PROFILES", Curves[i].StartPoint.X, Curves[i].StartPoint.Y, Curves[i].StartPoint.Z, Curves[i].EndPoint.X, Curves[i].EndPoint.Y, Curves[i].EndPoint.Z);
                dxf.add(line);
            }
            Curves.ForEach(c => c.Dispose());
            return this;
        }
        public ExportToDXF AddLines(Line[] lines, string layerName = "0", short ACADcolor = (short) 7)
        {
            CheckLayer(layerName, ACADcolor);
            for (int i = 0; i < lines.Length; i++)
            {
                DXFLibrary.Line line = new DXFLibrary.Line("PROFILES", lines[i].StartPoint.X, lines[i].StartPoint.Y, lines[i].StartPoint.Z, lines[i].EndPoint.X, lines[i].EndPoint.Y, lines[i].EndPoint.Z);
                dxf.add(line);
            }
            return this;
        }
        public void WriteFile(string filename)
        {
            using (FileStream f1 = new FileStream(filename + ".dxf", FileMode.Create))
            {
                DXFLibrary.Writer.Write(dxf, f1);
                f1.Close();
            }
        }
    }

    /// <summary>
    /// exports all curves as lines
    /// </summary>
    public class ExportLinesToDXF
    {
        //**CONSTRUCTOR
        internal ExportLinesToDXF(List<PolyCurve[]> polycurves, string filename)
        {
            DXFLibrary.Document dxf = new DXFLibrary.Document();
            DXFLibrary.Tables tables = new DXFLibrary.Tables();
            dxf.SetTables(tables);
            DXFLibrary.Table layers = new DXFLibrary.Table("LAYER");
            tables.addTable(layers);
            DXFLibrary.Layer layer = new DXFLibrary.Layer("PROFILES", 30, "CONTINUOUS");
            layers.AddTableEntry(layer);

            List<Curve> Curves = new List<Curve>();
            for (int i = 0; i < polycurves.Count; i++) for (int j = 0; j < polycurves[i].Length; j++) for (int k = 0; k < polycurves[i][j].Curves().Length; k++) 
                Curves.AddRange(polycurves[i][j].Curves()[k].ApproximateWithArcAndLineSegments());
            for (int i = 0; i < Curves.Count; i++)
            {
                DXFLibrary.Line line = new DXFLibrary.Line(
                    "PROFILES", 
                    Curves[i].StartPoint.X, 
                    Curves[i].StartPoint.Y, 
                    Curves[i].StartPoint.Z, 
                    Curves[i].EndPoint.X, 
                    Curves[i].EndPoint.Y, 
                    Curves[i].EndPoint.Z);
                dxf.add(line);
            }
            Curves.ForEach(c => c.Dispose());
            using (FileStream f1 = new FileStream(filename, FileMode.Create))
            {
                DXFLibrary.Writer.Write(dxf, f1);
                f1.Close();
            }
        }

        //**CREATE
        public static ExportLinesToDXF ByPolyCurves(List<PolyCurve[]> polycurves, string filename) { return new ExportLinesToDXF(polycurves, filename); }
    }

 }
