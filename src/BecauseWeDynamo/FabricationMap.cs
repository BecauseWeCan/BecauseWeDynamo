using System;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;

namespace Fabrication
{
    /// <summary>
    /// map geometry to XY plane
    /// </summary>
    public static class MaptoXY
    {
        /// <summary>
        /// map polycurves in index to XY plane
        /// with given spacing and at given y-coordinate
        /// </summary>
        public static List<List<PolyCurve>> MapPolyCurvesToRow(PolyCurve[][] Curves, CoordinateSystem[] CS, List<int> Index, double Xspacing, double Ycoordinate)
        {
            List<List<PolyCurve>> result = new List<List<PolyCurve>>(Index.Count);
            for (int i = 0; i < Index.Count; i++)
            {
                List<PolyCurve> temp = new List<PolyCurve>();
                for (int j = 0; j < Curves[Index[i]].Length; j++)
                {
                    CoordinateSystem TargetCS = CoordinateSystem.ByOrigin(Xspacing * i, Ycoordinate);
                    temp.Add((PolyCurve)Curves[Index[i]][j].Transform(CS[Index[i]], TargetCS));
                    TargetCS.Dispose();
                }
                result.Add(temp);
            }
            return result;
        }
        /// <summary>
        /// map circles in index to XY plane
        /// with given spacing and at given y-coordinate
        /// </summary>
        public static List<List<Circle>> MapCirclesToRow(Circle[][] Circles, CoordinateSystem[] CS, List<int> Index, double Xspacing, double Ycoordinate)
        {
            List<List<Circle>> result = new List<List<Circle>>(Index.Count);
            for (int i = 0; i < Index.Count; i++)
            {
                List<Circle> temp = new List<Circle>();
                for (int j = 0; j < Circles[Index[i]].Length; j++)
                {
                    CoordinateSystem TargetCS = CoordinateSystem.ByOrigin(Xspacing * i, Ycoordinate);
                    temp.Add((Circle)Circles[Index[i]][j].Transform(CS[Index[i]], TargetCS));
                    TargetCS.Dispose();
                }
                result.Add(temp);
            }
            return result;
        }
        /// <summary>
        /// creates a grid of polycurves on XY-plane 
        /// from an array of curves with a designated coordinate system
        /// </summary>
        public static List<PolyCurve> MapPolyCurve(PolyCurve[] PolyCurves, CoordinateSystem[] CS, int[] Index, double Xspacing = 0, int XmaxCount = 0, double Yspacing = 0, int YmaxCount = 0)
        {
            int Length = Index.Length;
            if (XmaxCount * YmaxCount > 0 && XmaxCount * YmaxCount < Index.Length) Length = XmaxCount * YmaxCount;
            double[,] position = new double[Length, 2];
            for (int i = 0; i < Length; i++)
            {
                if (XmaxCount + YmaxCount == 0)
                {
                    position[i, 0] = i;
                    position[i, 1] = 0;
                }
                else if (XmaxCount == 0)
                {
                    position[i, 1] = i % YmaxCount;
                    position[i, 0] = Math.Floor((double)i / YmaxCount);
                }
                else
                {
                    position[i, 0] = i % XmaxCount;
                    position[i, 1] = Math.Floor((double)i / XmaxCount);
                }
            }
            List<PolyCurve> result = new List<PolyCurve>();
            for (int i = 0; i < Length; i++)
            {
                CoordinateSystem targetCS = CoordinateSystem.ByOrigin(Xspacing * position[i, 0], Yspacing * position[i, 1]);
                result.Add((PolyCurve)PolyCurves[Index[i]].Transform(CS[i], targetCS));
                targetCS.Dispose();
            }
            return result;
        }
        /// <summary>
        /// creates a grid of curves on XY-plane 
        /// from an array of curves with a designated coordinate system
        /// </summary>
        public static List<Curve> MapCurve(Curve[] Curves, CoordinateSystem[] CS, int[] Index, double Xspacing = 0, int XmaxCount = 0, double Yspacing = 0, int YmaxCount = 0)
        {
            int Length = Index.Length;
            if (XmaxCount * YmaxCount > 0 && XmaxCount * YmaxCount < Index.Length) Length = XmaxCount * YmaxCount;
            double[,] position = new double[Length, 2];
            for (int i = 0; i < Length; i++)
            {
                if (XmaxCount + YmaxCount == 0)
                {
                    position[i, 0] = i;
                    position[i, 1] = 0;
                }
                else if (XmaxCount == 0)
                {
                    position[i, 1] = i % YmaxCount;
                    position[i, 0] = Math.Floor((double)i / YmaxCount);
                }
                else
                {
                    position[i, 0] = i % XmaxCount;
                    position[i, 1] = Math.Floor((double)i / XmaxCount);
                }
            }
            List<Curve> result = new List<Curve>();
            for (int i = 0; i < Length; i++)
            {
                CoordinateSystem targetCS = CoordinateSystem.ByOrigin(Xspacing * position[i, 0], Yspacing * position[i, 1]);
                result.Add((Curve)Curves[Index[i]].Transform(CS[i], targetCS));
                targetCS.Dispose();
            }
            return result;
        }
        /// <summary>
        /// creates a grid of geometry on XY-plane 
        /// from an array of curves with a designated coordinate system
        /// </summary>
        public static List<Autodesk.DesignScript.Geometry.Geometry> MapGeometry(Autodesk.DesignScript.Geometry.Geometry[] Geometry, CoordinateSystem[] CS, int[] Index, double Xspacing = 0, int XmaxCount = 0, double Yspacing = 0, int YmaxCount = 0)
        {
            int Length = Index.Length;
            if (XmaxCount * YmaxCount > 0 && XmaxCount * YmaxCount < Index.Length) Length = XmaxCount * YmaxCount;
            double[,] position = new double[Length, 2];
            for (int i = 0; i < Length; i++)
            {
                if (XmaxCount + YmaxCount == 0)
                {
                    position[i, 0] = i;
                    position[i, 1] = 0;
                }
                else if (XmaxCount == 0)
                {
                    position[i, 1] = i % YmaxCount;
                    position[i, 0] = Math.Floor((double)i / YmaxCount);
                }
                else
                {
                    position[i, 0] = i % XmaxCount;
                    position[i, 1] = Math.Floor((double)i / XmaxCount);
                }
            }
            List<Autodesk.DesignScript.Geometry.Geometry> result = new List<Autodesk.DesignScript.Geometry.Geometry>();
            for (int i = 0; i < Length; i++)
            {
                CoordinateSystem targetCS = CoordinateSystem.ByOrigin(Xspacing * position[i, 0], Yspacing * position[i, 1]);
                result.Add(Geometry[Index[i]].Transform(CS[i], targetCS));
                targetCS.Dispose();
            }
            return result;
        }
        /// <summary>
        /// creates a grid of polycurves on XY-plane 
        /// from an array of curves with a designated coordinate system
        /// </summary>
        public static List<List<PolyCurve>> MapPolyCurves(PolyCurve[][] PolyCurves, CoordinateSystem[] CS, int[] Index, double Xspacing = 0, int XmaxCount = 0, double Yspacing = 0, int YmaxCount = 0)
        {
            int Length = Index.Length;
            if (XmaxCount * YmaxCount > 0 && XmaxCount * YmaxCount < Index.Length) Length = XmaxCount * YmaxCount;
            double[,] position = new double[Length, 2];
            for (int i = 0; i < Length; i++)
            {
                if (XmaxCount + YmaxCount == 0)
                {
                    position[i, 0] = i;
                    position[i, 1] = 0;
                }
                else if (XmaxCount == 0)
                {
                    position[i, 1] = i % YmaxCount;
                    position[i, 0] = Math.Floor((double)i / YmaxCount);
                }
                else
                {
                    position[i, 0] = i % XmaxCount;
                    position[i, 1] = Math.Floor((double)i / XmaxCount);
                }
            }
            List<List<PolyCurve>> result = new List<List<PolyCurve>>();
            for (int i = 0; i < Length; i++)
            {
                CoordinateSystem targetCS = CoordinateSystem.ByOrigin(Xspacing * position[i, 0], Yspacing * position[i, 1]);
                List<PolyCurve> temp = new List<PolyCurve>(PolyCurves[Index[i]].Length);
                for (int j = 0; j < PolyCurves[Index[i]].Length; j++)
                {
                    temp.Add((PolyCurve)PolyCurves[Index[i]][j].Transform(CS[Index[i]], targetCS));
                }
                targetCS.Dispose();
                result.Add(temp);
            }
            return result;
        }
        /// <summary>
        /// creates a grid of curves on XY-plane 
        /// from an array of curves with a designated coordinate system
        /// </summary>
        public static List<List<Curve>> MapCurves(Curve[][] Curves, CoordinateSystem[] CS, int[] Index, double Xspacing = 0, int XmaxCount = 0, double Yspacing = 0, int YmaxCount = 0)
        {
            int Length = Index.Length;
            if (XmaxCount * YmaxCount > 0 && XmaxCount * YmaxCount < Index.Length) Length = XmaxCount * YmaxCount;
            double[,] position = new double[Length, 2];
            for (int i = 0; i < Length; i++)
            {
                if (XmaxCount + YmaxCount == 0)
                {
                    position[i, 0] = i;
                    position[i, 1] = 0;
                }
                else if (XmaxCount == 0)
                {
                    position[i, 1] = i % YmaxCount;
                    position[i, 0] = Math.Floor((double)i / YmaxCount);
                }
                else
                {
                    position[i, 0] = i % XmaxCount;
                    position[i, 1] = Math.Floor((double)i / XmaxCount);
                }
            }
            List<List<Curve>> result = new List<List<Curve>>();
            for (int i = 0; i < Length; i++)
            {
                CoordinateSystem targetCS = CoordinateSystem.ByOrigin(Xspacing * position[i, 0], Yspacing * position[i, 1]);
                List<Curve> temp = new List<Curve>(Curves[Index[i]].Length);
                for (int j = 0; j < Curves[Index[i]].Length; j++)
                {
                    temp.Add((Curve)Curves[Index[i]][j].Transform(CS[Index[i]], targetCS));
                }
                targetCS.Dispose();
                result.Add(temp);
            }
            return result;
        }
        /// <summary>
        /// creates a grid of circles on XY-plane 
        /// from an array of curves with a designated coordinate system
        /// </summary>
        public static List<List<Circle>> MapCircles(Circle[][] Circles, CoordinateSystem[] CS, int[] Index, double Xspacing = 0, int XmaxCount = 0, double Yspacing = 0, int YmaxCount = 0)
        {
            int Length = Index.Length;
            if (XmaxCount * YmaxCount > 0 && XmaxCount * YmaxCount < Index.Length) Length = XmaxCount * YmaxCount;
            double[,] position = new double[Length, 2];
            for (int i = 0; i < Length; i++)
            {
                if (XmaxCount + YmaxCount == 0)
                {
                    position[i, 0] = i;
                    position[i, 1] = 0;
                }
                else if (XmaxCount == 0)
                {
                    position[i, 1] = i % YmaxCount;
                    position[i, 0] = Math.Floor((double)i / YmaxCount);
                }
                else
                {
                    position[i, 0] = i % XmaxCount;
                    position[i, 1] = Math.Floor((double)i / XmaxCount);
                }
            }
            List<List<Circle>> result = new List<List<Circle>>();
            for (int i = 0; i < Length; i++)
            {
                CoordinateSystem targetCS = CoordinateSystem.ByOrigin(Xspacing * position[i, 0], Yspacing * position[i, 1]);
                List<Circle> temp = new List<Circle>(Circles[Index[i]].Length);
                for (int j = 0; j < Circles[Index[i]].Length; j++)
                {
                    temp.Add((Circle)Circles[Index[i]][j].Transform(CS[Index[i]], targetCS));
                }
                targetCS.Dispose();
                result.Add(temp);
            }
            return result;
        }
        /// <summary>
        /// creates a grid of geometry on XY-plane 
        /// from an array of curves with a designated coordinate system
        /// </summary>
        public static List<List<Autodesk.DesignScript.Geometry.Geometry>> MapGeometries(Autodesk.DesignScript.Geometry.Geometry[][] Geometry, CoordinateSystem[] CS, int[] Index, double Xspacing = 0, int XmaxCount = 0, double Yspacing = 0, int YmaxCount = 0)
        {
            int Length = Index.Length;
            if (XmaxCount * YmaxCount > 0 && XmaxCount * YmaxCount < Index.Length) Length = XmaxCount * YmaxCount;
            double[,] position = new double[Length, 2];
            for (int i = 0; i < Length; i++)
            {
                if (XmaxCount + YmaxCount == 0)
                {
                    position[i, 0] = i;
                    position[i, 1] = 0;
                }
                else if (XmaxCount == 0)
                {
                    position[i, 1] = i % YmaxCount;
                    position[i, 0] = Math.Floor((double)i / YmaxCount);
                }
                else
                {
                    position[i, 0] = i % XmaxCount;
                    position[i, 1] = Math.Floor((double)i / XmaxCount);
                }
            }
            List<List<Autodesk.DesignScript.Geometry.Geometry>> result = new List<List<Autodesk.DesignScript.Geometry.Geometry>>();
            for (int i = 0; i < Length; i++)
            {
                CoordinateSystem targetCS = CoordinateSystem.ByOrigin(Xspacing * position[i, 0], Yspacing * position[i, 1]);
                List<Autodesk.DesignScript.Geometry.Geometry> temp = new List<Autodesk.DesignScript.Geometry.Geometry>(Geometry[Index[i]].Length);
                for (int j = 0; j < Geometry[Index[i]].Length; j++)
                {
                    temp.Add(Geometry[Index[i]][j].Transform(CS[Index[i]], targetCS));
                }
                targetCS.Dispose();
                result.Add(temp);
            }
            return result;
        }

    }
}
