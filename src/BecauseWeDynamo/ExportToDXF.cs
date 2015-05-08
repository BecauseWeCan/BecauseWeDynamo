using System;
using System.IO;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;

namespace Utilities
{
    /// <summary>
    /// exports circles as lines and circles
    /// using: dxflibrary by Cosmin Radoi from GitHub
    /// </summary>
    public class ExportCurvesAsArcsToDXF
    {
        //**Properties
        private DXFLibrary.Document dxf;
        private DXFLibrary.Tables tables;
        private DXFLibrary.Table layers;
        private List<string> layerNames;

        //**CONSTRUCTOR
        internal ExportCurvesAsArcsToDXF()
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
        
        //**CREATE
        /// <summary>
        /// initializes dxfEporter with default layer "0" with color white
        /// </summary>
        /// <returns>initialized dxfEporter</returns>
        public static ExportCurvesAsArcsToDXF ByDefault() { return new ExportCurvesAsArcsToDXF(); }     
        /// <summary>
        /// initializes dxfEporter with polycurves added
        /// </summary>
        /// <param name="polycurves">polycurves to add</param>
        /// <param name="layerName">layer name (default is "0")</param>
        /// <param name="ACADcolor">autocad color as short (defaualt is whtie)</param>
        /// <returns></returns>
        public static ExportCurvesAsArcsToDXF ByPolyCurves(List<PolyCurve[]> polycurves, string layerName = "0", short ACADcolor = (short) 7) 
        {
            ExportCurvesAsArcsToDXF dxfWriter = new ExportCurvesAsArcsToDXF();
            dxfWriter.AddPolyCurves(polycurves.ToArray(), layerName, ACADcolor);
            return dxfWriter; 
        }

        //**ACTIONS
        /// <summary>
        /// writes the dxf file
        /// </summary>
        /// <param name="filename">string for file location including filename, extensions is added</param>
        public void WriteFile(string filename)
        {
            FileStream f1 = new FileStream(filename + ".dxf", FileMode.Create);
            DXFLibrary.Writer.Write(dxf, f1);
            f1.Close();
        }
        /// <summary>
        /// writes polycurve data to dxfExporter
        /// </summary>
        /// <param name="circles">array of circles</param>
        /// <param name="layerName">layer name (default is "0")</param>
        /// <param name="ACADcolor">autocad color as short (default is white)</param>
        /// <returns>dxfExporter</returns>
        public ExportCurvesAsArcsToDXF AddCircles(Circle[][] circles, string layerName = "0", short ACADcolor = (short) 7)
        {
            if (!layerNames.Contains(layerName))
            {
                DXFLibrary.Layer layer = new DXFLibrary.Layer(layerName, ACADcolor, "CONTINUOUS");
                layers.AddTableEntry(layer);
                layerNames.Add(layerName);
            }
            for (int i = 0; i < circles.Length; i++)
            {
                for (int j = 0; j < circles[i].Length; j++)
                {
                    DXFLibrary.Circle circle = new DXFLibrary.Circle(circles[i][j].CenterPoint.X, circles[i][j].CenterPoint.Y, circles[i][j].Radius, layerName);
                    dxf.add(circle);
                }
            }
            return this;
        }
        /// <summary>
        /// writes polycurve data to dxfExporter
        /// </summary>
        /// <param name="polycurves">array of polycurves</param>
        /// <param name="layerName">layer name (default is "0")</param>
        /// <param name="ACADcolor">autocad color as short (default = white)</param>
        /// <returns>dxfExporter</returns>
        public ExportCurvesAsArcsToDXF AddPolyCurves(PolyCurve[][] polycurves, string layerName = "0", short ACADcolor = (short) 7)
        {
            if (!layerNames.Contains(layerName))
            {
                DXFLibrary.Layer layer = new DXFLibrary.Layer(layerName, ACADcolor, "CONTINUOUS");
                layers.AddTableEntry(layer);
                layerNames.Add(layerName);
            }

            List<Curve> Curves = new List<Curve>();
            for (int i = 0; i < polycurves.Length; i++)
            {
                for (int j = 0; j < polycurves[i].Length; j++)
                {
                    for (int k = 0; k < polycurves[i][j].Curves().Length; k++)
                    {
                        Curves.AddRange(polycurves[i][j].Curves()[k].ApproximateWithArcAndLineSegments());
                    }
                }
            }

            for (int i = 0; i < Curves.Count; i++)
            {
                    Vector stTan = Curves[i].TangentAtParameter(0);
                    Vector endTan = Curves[i].TangentAtParameter(Curves[i].EndParameter());
                    Point stpt = Curves[i].StartPoint;
                    Point endpt = Curves[i].EndPoint;

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
                    }
                }
            return this;
            }
    } // end Class


    /// <summary>
    /// exports circles to lines
    /// using: dxflibrary by Cosmin Radoi from GitHub
    /// </summary>
    public class ExportCurvesAsLinesToDXF
    {
        //**CONSTRUCTOR
        /// <summary>
        /// exports 2D circles to a dxf file.
        /// </summary>
        /// <param name="circles">input circles to be exported</param>
        /// <param name="filename">target file location and filename including extension</param>
        internal ExportCurvesAsLinesToDXF(List<PolyCurve[]> polycurves, string filename)
        {
            DXFLibrary.Document dxf = new DXFLibrary.Document();
            DXFLibrary.Tables tables = new DXFLibrary.Tables();
            dxf.SetTables(tables);
            DXFLibrary.Table layers = new DXFLibrary.Table("LAYER");
            tables.addTable(layers);
            DXFLibrary.Layer layer = new DXFLibrary.Layer("PROFILES", 30, "CONTINUOUS");
            layers.AddTableEntry(layer);

            List<Curve> Curves = new List<Curve>();
            for (int i = 0; i < polycurves.Count; i++)
            {
                for (int j = 0; j < polycurves[i].Length; j++)
                {
                    for (int k = 0; k < polycurves[i][j].Curves().Length; k++)
                    {
                        Curves.AddRange(polycurves[i][j].Curves()[k].ApproximateWithArcAndLineSegments());
                    }

                }
            }

            for (int i = 0; i < Curves.Count; i++)
            {
                Point stpt = Curves[i].StartPoint;
                Point endpt = Curves[i].EndPoint;
                DXFLibrary.Line line = new DXFLibrary.Line("PROFILES", stpt.X, stpt.Y, stpt.Z, endpt.X, endpt.Y, endpt.Z);
                dxf.add(line);
            }

            FileStream f1 = new FileStream(filename, FileMode.Create);
            DXFLibrary.Writer.Write(dxf, f1);
            f1.Close();

            
        }

        //**CREATE
        public static ExportCurvesAsLinesToDXF ByPolyCurves(List<PolyCurve[]> polycurves, string filename)
        {
            return new ExportCurvesAsLinesToDXF(polycurves, filename);
        }

    }

 }
