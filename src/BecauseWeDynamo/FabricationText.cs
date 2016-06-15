using System;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;

namespace Fabrication
{

    internal class letter : List<PolyCurve>, IDisposable
    {
        private char ltr;
        bool disposed = false;

        internal letter(char c)
        {
            Point[,] points = new Point[3, 5];
            for (int i = 0; i < 3; i++) for (int j = 0; j < 5; j++) points[i, j] = Point.ByCoordinates(i - 1, -j - 1);
            ltr = c;
            switch (c)
            {
                case 'a':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[0, 2], points[1, 0], points[2, 2], points[2, 4] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 2], points[2, 2] }));
                    break;
                case 'b':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[1, 2], points[2, 3], points[1, 4], points[0, 4] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 2], points[1, 2] }));
                    break;
                case 'c':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 1], points[1, 0], points[0, 1], points[0, 3], points[1, 4], points[2, 3] }));
                    break;
                case 'd':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 4] }));
                    break;
                case 'e':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 0], points[1, 0], points[0, 1], points[0, 3], points[1, 4], points[2, 4] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 2], points[2, 2] }));
                    break;
                case 'f':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 0], points[1, 0], points[0, 1], points[0, 4], points[0, 2], points[2, 2] }));
                    break;
                case 'g':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 1], points[1, 0], points[0, 1], points[0, 3], points[1, 4], points[2, 3], points[2, 2], points[1, 2] }));
                    break;
                case 'h':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[0, 4], points[0, 2], points[2, 2], points[2, 0], points[2, 4] }));
                    break;
                case 'i':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[2, 0], points[1, 0], points[1, 4], points[0, 4], points[2, 4] }));
                    break;
                case 'j':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[2, 0], points[1, 0], points[1, 3], points[0, 4] }));
                    break;
                case 'k':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[0, 4] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 4], points[0, 2], points[2, 0] }));
                    break;
                case 'l':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[0, 4], points[2, 4] }));
                    break;
                case 'm':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[0, 0], points[1, 2], points[2, 0], points[2, 4] }));
                    break;
                case 'n':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[0, 0], points[2, 4], points[2, 0] }));
                    break;
                case 'o':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 3], points[0, 1], points[1, 0] }));
                    break;
                case 'p':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[1, 2], points[0, 2] }));
                    break;
                case 'q':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 3], points[0, 1], points[1, 0] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[1, 3], points[2, 4] }));
                    break;
                case 'r':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[1, 2], points[0, 2], points[2, 4] }));
                    break;
                case 's':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 1], points[1, 0], points[0, 1], points[2, 3], points[1, 4], points[0, 3] }));
                    break;
                case 't':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[2, 0], points[1, 0], points[1, 4] }));
                    break;
                case 'u':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[0, 4], points[2, 4], points[2, 0] }));
                    break;
                case 'v':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[0, 2], points[1, 4], points[2, 2], points[2, 0] }));
                    break;
                case 'w':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[0, 4], points[1, 2], points[2, 4], points[2, 0] }));
                    break;
                case 'x':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[2, 4] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[2, 0] }));
                    break;
                case 'y':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[1, 2] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[1, 4], points[1, 2], points[2, 0] }));
                    break;
                case 'z':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[2, 0], points[0, 4], points[2, 4] }));
                    break;
                case '-':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 2], points[2, 2] }));
                    break;
                case '0':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 3], points[0, 1], points[1, 0] }));
                    this.Add(PolyCurve.ByPoints(new Point[] { points[1, 3], points[2, 1] }));
                    break;
                case '1':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 1], points[1, 0], points[1, 4], points[0, 4], points[2, 4] }));
                    break;
                case '2':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[2, 0], points[2, 1], points[0, 4], points[2, 4] }));
                    break;
                case '3':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 0], points[2, 0], points[1, 1], points[2, 1], points[2, 3], points[1, 4], points[0, 4] }));
                    break;
                case '4':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 4], points[2, 0], points[0, 3], points[2, 3] }));
                    break;
                case '5':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 0], points[0, 0], points[0, 1], points[2, 1], points[2, 3], points[1, 4], points[0, 4] }));
                    break;
                case '6':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 0], points[1, 0], points[0, 1], points[0, 4], points[1, 4], points[2, 3], points[2, 1], points[1, 1], points[0, 2] }));
                    break;
                case '7':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 1], points[0, 0], points[2, 0], points[2, 1], points[1, 4] }));
                    break;
                case '8':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[2, 0], points[2, 1], points[0, 2], points[0, 4], points[2, 4], points[2, 2], points[0, 1], points[0, 0], points[2, 0] }));
                    break;
                case '9':
                    this.Add(PolyCurve.ByPoints(new Point[] { points[0, 4], points[1, 4], points[2, 3], points[2, 1], points[1, 0], points[0, 1], points[0, 2], points[1, 3], points[2, 2] }));
                    break;
                default:
                    break;
            } // end switch

            for (int i = 0; i < 3; i++) for (int j = 0; j < 5; j++) points[i, j].Dispose();
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing) { this.ForEach(p => p.Dispose()); this.Clear(); }
            disposed = true;
        }
    }// end letter

    /// <summary>
    /// Word object composed on polycurve letters
    /// </summary>
    public class word : IDisposable
    {
        internal string Word;
        internal CoordinateSystem cs;
        bool disposed = false;

        /// <summary>
        /// letter array in word
        /// </summary>
        internal letter[] Letters { get; private set; }

        internal word(string String, Point Origin, Vector X, Vector Y)
        {
            cs = CoordinateSystem.ByOriginVectors(Origin, X, Y);
            Word = String;
            char[] chars = String.ToLower().ToCharArray();
            Letters = new letter[chars.Length];
            double dblLtr = (chars.Length - 1) * (-1.5);
            for (int i = 0; i < chars.Length; i++)
            {
                Letters[i] = new letter(chars[i]);
                for (int j = 0; j < Letters[i].Count; j++) Letters[i][j] = (PolyCurve)Letters[i][j].Transform(cs).Translate(cs.XAxis, dblLtr);
                dblLtr = dblLtr + 3.0;
            }
        }

        /// <summary>
        /// creates word object by string, origin point, orthonormal vectors X, Y
        /// </summary>
        public static word ByStringOriginVectors(String Word, Point Origin, Vector X, Vector Y) { return new word(Word, Origin, X, Y); }

        /// <summary>
        /// creates word object by string and coordinate system
        /// </summary>
        public static word ByStringCS(String Word, CoordinateSystem CS) { return new word(Word, CS.Origin, CS.XAxis, CS.YAxis); }

        /// <summary>
        /// transform word object to new coordinate system
        /// </summary>
        /// <param name="newCS">New Coordinate System</param>
        /// <returns>Transformed Word</returns>
        public word Transform(CoordinateSystem newCS)
        {
            for (int i = 0; i < Letters.Length; i++) for (int j = 0; j < Letters[i].Count; j++)
                { Letters[i][j] = Letters[i][j].Transform(newCS) as PolyCurve; }
            return this;
        }

        /// <summary>
        /// translates word object to new coordinates
        /// based on vector direction and scalar distance
        /// </summary>
        /// <param name="Direction">Direction Vector</param>
        /// <param name="Distance">Distance</param>
        /// <returns>Translated Object</returns>
        public word Translate(Vector Direction, double Distance)
        {
            for (int i = 0; i < Letters.Length; i++) for (int j = 0; j < Letters[i].Count; j++)
                { Letters[i][j] = Letters[i][j].Translate(Direction, Distance) as PolyCurve; }
            return this;
        }

        /// <summary>
        /// scales word about object origin by given scale factor
        /// </summary>
        /// <param name="Factor">Scale Factor</param>
        /// <returns>Scaled Object</returns>
        public word Scale(double Factor=1)
        {
            for (int i = 0; i < Letters.Length; i++) for (int j = 0; j < Letters[i].Count; j++)
                { Letters[i][j] = Letters[i][j].Scale(cs.Origin, cs.Origin.Translate(cs.XAxis, 1) as Point, cs.Origin.Translate(cs.XAxis, Factor) as Point) as PolyCurve; }
            return this;
        }

        /// <summary>
        /// creates polycurves for text
        /// </summary>
        /// <param name="Scale"></param>
        /// <returns></returns>
        public List<PolyCurve> Display(double Scale=1)
        {
            List<PolyCurve> lines = new List<PolyCurve>();
            for (int i = 0; i < Letters.Length; i++) for (int j = 0; j < Letters[i].Count; j++)
                { lines.Add(Letters[i][j].Scale(cs.Origin, cs.Origin.Translate(cs.XAxis, 1) as Point, cs.Origin.Translate(cs.XAxis, Scale) as Point) as PolyCurve); }
            return lines;
        }
        /// <summary>
        /// dispose function
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }
        /// <summary>
        /// protected dispose function
        /// </summary>
        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing)
            {
                Letters.ForEach(l => l.Dispose());
                cs.Dispose();
            }
            disposed = true;
        }

    }

}
