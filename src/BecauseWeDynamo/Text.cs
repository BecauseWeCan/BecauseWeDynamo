using System;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;

namespace Text
{

    public class Letter : IDisposable
    {
        private PolyCurve[] pl;
        private char ltr;
        bool disposed = false;

        public PolyCurve[] drawn { get { return pl; } }

        internal Letter(char c)
        {
            Point[,] points = new Point[3, 5];
            for (int i = 0; i < 3; i++) for (int j = 0; j < 5; j++) points[i, j] = Point.ByCoordinates(i - 1, -j - 1);
            ltr = c;
            switch (c)
            {
                case 'a':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[0, 2], points[1, 0], points[2, 2], points[2, 4] }),
                        PolyCurve.ByPoints( new Point[] { points[0, 2], points[2, 2] })};
                    break;
                case 'b':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[1, 3], points[2, 4], points[1, 5], points[0, 5] }),
                        PolyCurve.ByPoints( new Point[] { points[0, 3], points[1, 3] })};
                    break;
                case 'c':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[2, 1], points[1, 0], points[0, 1], points[0, 3], points[1, 4], points[2, 3] })};
                    break;
                case 'd':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 4] })};
                    break;
                case 'e':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[2, 0], points[1, 0], points[0, 1], points[0, 3], points[1, 4], points[2, 4] }),
                        PolyCurve.ByPoints( new Point[] { points[0, 2], points[2, 2] })};
                    break;
                case 'f':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[2, 0], points[1, 0], points[0, 1], points[0, 4], points[0, 2], points[2, 2] })};
                    break;
                case 'g':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[2, 1], points[1, 0], points[0, 1], points[0, 3], points[1, 4], points[2, 3], points[2, 2], points[1,2] })};
                    break;
                case 'h':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[0, 4], points[0, 2], points[2, 2], points[2, 0], points[2, 4] })};
                    break;
                case 'i':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[2, 0], points[1, 0], points[1, 4], points[0, 4], points[2, 4] })};
                    break;
                case 'j':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[2, 0], points[1, 0], points[1, 3], points[0, 4] })};
                    break;
                case 'k':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[0, 4] }),
                        PolyCurve.ByPoints( new Point[] { points[2, 4], points[0, 2], points[2, 0] })};
                    break;
                case 'l':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[0, 4], points[2, 4] })};
                    break;
                case 'm':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[0, 0], points[1, 2], points[2, 0], points[2, 4] })};
                    break;
                case 'n':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[0, 0], points[2, 4], points[2, 0] })};
                    break;
                case 'o':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 3], points[0, 1], points[1, 0] })};
                    break;
                case 'p':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[1, 2], points[0, 2] })};
                    break;
                case 'q':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 3], points[0, 1], points[1, 0] }),
                        PolyCurve.ByPoints( new Point[] { points[1, 3], points[1, 4] })};
                    break;
                case 'r':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[0, 0], points[1, 0], points[2, 1], points[1, 2], points[0, 2], points[2,4] })};
                    break;
                case 's':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[2, 1], points[1, 0], points[0, 1], points[2, 3], points[1, 4], points[0, 3] })};
                    break;
                case 't':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[2, 0], points[1, 0], points[1, 4] })};
                    break;
                case 'u':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[0, 4], points[2, 4], points[2, 0] })};
                    break;
                case 'v':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[0, 2], points[1, 4], points[2, 2], points[2, 0] })};
                    break;
                case 'w':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[0, 4], points[1, 2], points[2, 4], points[2, 0] })};
                    break;
                case 'x':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[2, 4] }),
                        PolyCurve.ByPoints( new Point[] { points[0, 4], points[2, 0] })};
                    break;
                case 'y':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[1, 2] }),
                        PolyCurve.ByPoints( new Point[] { points[1, 4], points[1, 2], points[2, 0] })};
                    break;
                case 'z':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 0], points[2, 0], points[0, 4], points[2, 4] })};
                    break;

                case '-':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints( new Point[] { points[0, 2], points[2, 2] })};
                    break;

                case '0':
                    pl = new PolyCurve[2] {
                        PolyCurve.ByPoints( new Point[] { points[1, 0], points[2, 1], points[2, 3], points[1, 4], points[0, 3], points[0, 1], points[1, 0] }),
                        PolyCurve.ByPoints( new Point[] { points[1, 3], points[2, 1] })};
                    break;
                case '1':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[0, 1], points[1, 0], points[1, 4], points[0, 4], points[2, 4] })};
                    break;
                case '2':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[0, 0], points[2, 0], points[2, 1], points[0, 4], points[2, 4] })};
                    break;
                case '3':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[0, 0], points[2, 0], points[1, 1], points[2, 1], points[2, 3], points[1,4], points[0, 4] })};
                    break;
                case '4':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[2, 4], points[2, 0], points[0, 3], points[2, 3] })};
                    break;
                case '5':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[2, 0], points[0, 0], points[0, 1], points[2, 1], points[2, 3], points[1, 4], points[0, 4] })};
                    break;
                case '6':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[2, 0], points[1, 0], points[0, 1], points[0, 4], points[1, 4], points[2, 3], points[2, 1], points[1, 1], points[0, 2] })};
                    break;
                case '7':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[0, 1], points[0, 0], points[2, 0], points[2, 1], points[1, 4] })};
                    break;
                case '8':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[2, 0], points[2, 1], points[0, 2], points[0, 4], points[2, 4], points[2, 2], points[0, 1], points[0, 0], points[2, 0] })};
                    break;
                case '9':
                    pl = new PolyCurve[1] {
                        PolyCurve.ByPoints(
                        new Point[] { points[0, 4], points[1, 4], points[2, 3], points[2, 1], points[1, 0], points[0, 1], points[0, 2], points[1, 3], points[2, 2] })};
                    break;

                default:
                    pl = null;
                    break;
            } // end switch

            for (int i = 0; i < 3; i++) for (int j = 0; j < 5; j++) points[i, j].Dispose();
        }

        public Letter Transform(CoordinateSystem newCS)
        {
            pl.ForEach(p => p = p.Transform(newCS) as PolyCurve);
            return this;
        }

        public Letter Translate(Vector vec, double distance)
        {
            pl.ForEach(p => p = p.Translate(vec, distance) as PolyCurve);
            return this;
        }

        public Letter Scale(CoordinateSystem cs, double factor)
        {
            Point from = cs.Origin.Translate(cs.XAxis, 1) as Point;
            Point to = cs.Origin.Translate(cs.XAxis, factor) as Point;
            pl.ForEach(p => p = p.Scale(cs.Origin, from, to) as PolyCurve);
            from.Dispose();
            to.Dispose();
            return this;
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing) pl.ForEach(p => p.Dispose());
            disposed = true;
        }
    }// end letter

    public class Word : IDisposable
    {
        private string word;
        private Letter[] Letters;
        private CoordinateSystem cs;
        bool disposed = false;

        public Letter[] letters { get { return Letters; } }

        internal Word(string str, Point origin, Vector X, Vector Y)
        {
            cs = CoordinateSystem.ByOriginVectors(origin, X, Y);
            word = str;
            char[] ltrs = str.ToCharArray();
            Letters = new Letter[ltrs.Length];
            double dblLtr = (ltrs.Length - 1) * (-1.5);
            for (int i = 0; i < ltrs.Length; i++)
            {
                Letters[i] = new Letter(ltrs[i]).Transform(cs).Translate(cs.XAxis, dblLtr);
                dblLtr = dblLtr + 3.0;
            }
        }

        public static Word ByString(String word, Point origin, Vector X, Vector Y)
        {
            return new Word(word, origin, X, Y);
        }


        public Word Transform(CoordinateSystem newCS)
        {
            for (int i = 0; i < Letters.Length; i++)
            {
                Letters[i] = Letters[i].Transform(newCS);
            }
            return this;
        }

        public Word Translate(Vector vec, double distance)
        {
            for (int i = 0; i < Letters.Length; i++)
            {
                Letters[i] = Letters[i].Translate(vec, distance);
            }
            return this;
        }

        public Word Scale(double factor)
        {
            for (int i = 0; i < Letters.Length; i++)
            {
                Letters[i] = Letters[i].Scale(cs, factor);
            }
            return this;
        }

        public List<PolyCurve> display(double factor)
        {
            List<PolyCurve> lines = new List<PolyCurve>();
            for (int i = 0; i < Letters.Length; i++) for (int j = 0; j < Letters[i].drawn.Length; j++)
            {
                lines.Add( Letters[i].drawn[j].Scale(cs.Origin, cs.Origin.Translate(cs.XAxis, 1) as Point, cs.Origin.Translate(cs.XAxis, factor) as Point) as PolyCurve );
            }
            return lines;
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposed) return;
            if (disposing) for (int i = 0; i < Letters.Length; i++) Letters[i].Dispose();
            disposed = true;
        }

    }

}
