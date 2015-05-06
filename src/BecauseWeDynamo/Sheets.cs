using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Autodesk.DesignScript.Geometry;
using Autodesk.DesignScript.Interfaces;
using Autodesk.DesignScript.Runtime;

namespace Fabrication
{
    public class Sheets
    {
        List<List<PolyCurve>> curves;
        List<List<Circle>> circles;
        List<CoordinateSystem> cs;


        internal Sheets(List<List<PolyCurve>> Curves, List<CoordinateSystem> CS)
        {
            cs = CS;
            curves = Curves;
        }
        internal Sheets(List<List<Circle>> Circles, List<CoordinateSystem> CS)
        {
            cs = CS;
            circles = Circles;
        }
        
        //**CREATE
        /// <summary>
        /// creates a sheet layout of an array of Polycurves with designated coordinate system
        /// </summary>
        /// <param name="Curves">list of polycurves</param>
        /// <param name="CS">list of polycurves coordinate system</param>
        /// <returns>sheet object</returns>
        public static Sheets ByPolyCurvesAndCS(List<List<PolyCurve>> Curves, List<CoordinateSystem> CS)
        {
            return new Sheets(Curves, CS);
        }
        /// <summary>
        /// creates a sheet layout of an array of Circles with designated coordinate system
        /// </summary>
        /// <param name="Circles">list of circles</param>
        /// <param name="CS">list of circles coordinate system</param>
        /// <returns>sheet object</returns>
        public static Sheets ByCirclesAndCS(List<List<Circle>> Circles, List<CoordinateSystem> CS)
        {
            return new Sheets(Circles, CS);
        }
        /// <summary>
        /// creates a sheet layout of an array of Curves with designated coordinate system
        /// </summary>
        /// <param name="Curves">list of curves</param>
        /// <param name="CS">list of curves coordinate system</param>
        /// <returns>sheet object</returns>
        public static Sheets ByCurvesAndCS(List<List<Curve>> Curves, List<CoordinateSystem> CS)
        {
            List<List<PolyCurve>> result = new List<List<PolyCurve>>(Curves.Count);
            for (int i = 0; i < Curves.Count; i++)
            {
                List<PolyCurve> shape = new List<PolyCurve>(Curves[i].Count);
                for (int j = 0; j < Curves[i].Count; j++)
                {
                    shape.Add( PolyCurve.ByJoinedCurves( new List<Curve>{Curves[i][j]}) );
                }
                result.Add(shape);
            }
            return new Sheets(result, CS);
        }

        //**ACTIONS
        /// <summary>
        /// returns a row containing given index list
        /// </summary>
        /// <param name="index">index list</param>
        /// <param name="X">spacing in X direction</param>
        /// <param name="Y">Y-coordinate</param>
        /// <returns>row of polycurves as list of polycurves</returns>
        public List<List<PolyCurve>> GetPolyCurveIndexAsRow(List<int> index, double X, double Y)
        {
            List<List<PolyCurve>> result = new List<List<PolyCurve>>(index.Count);
            for (int i = 0; i < index.Count; i++ )
            {
                List<PolyCurve> shape = new List<PolyCurve>(curves[index[i]].Count);
                for (int j = 0; j < curves[index[i]].Count; j++ )
                {
                    shape.Add( (PolyCurve) curves[index[i]][j].Transform( cs[index[i]], CoordinateSystem.ByOrigin(X * i, Y)));
                }
                result.Add(shape);
            }
            return result;
        }
        /// <summary>
        /// returns a row containing given index list
        /// </summary>
        /// <param name="index">index list</param>
        /// <param name="X">spacing in X direction</param>
        /// <param name="Y">Y-coordinate</param>
        /// <returns>row of circles as list of circles</returns>
        public List<List<Circle>> GetCircleIndexAsRow(List<int> index, double X, double Y)
        {
            List<List<Circle>> result = new List<List<Circle>>(index.Count);
            for (int i = 0; i < index.Count; i++)
            {
                List<Circle> shape = new List<Circle>(circles[index[i]].Count);
                for (int j = 0; j < circles[index[i]].Count; j++)
                {
                    shape.Add((Circle)circles[index[i]][j].Transform(cs[index[i]], CoordinateSystem.ByOrigin(X * i, Y)));
                }
                result.Add(shape);
            }
            return result;
        }
    }
}
