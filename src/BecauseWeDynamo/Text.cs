using System;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;

namespace Text
{

    public class Letter : List<PolyCurve>, IDisposable
    {
        private char ltr;
        bool disposed = false;

        internal Letter(char c)
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

    public class Word : IDisposable
    {
        private string word;
        private Letter[] letters;
        private CoordinateSystem cs;
        bool disposed = false;

        public Letter[] Letters { get { return letters; } }

        internal Word(string str, Point origin, Vector X, Vector Y)
        {
            cs = CoordinateSystem.ByOriginVectors(origin, X, Y);
            word = str;
            char[] ltrs = str.ToCharArray();
            letters = new Letter[ltrs.Length];
            double dblLtr = (ltrs.Length - 1) * (-1.5);
            for (int i = 0; i < ltrs.Length; i++)
            {
                letters[i] = new Letter(ltrs[i]);
                for (int j = 0; j < letters[i].Count; j++) letters[i][j] = (PolyCurve) letters[i][j].Transform(cs).Translate(cs.XAxis, dblLtr);
                dblLtr = dblLtr + 3.0;
            }
        }

        public static Word ByStringOriginVectors(String word, Point origin, Vector X, Vector Y) { return new Word(word, origin, X, Y); }
        public static Word ByStringCS(String word, CoordinateSystem CS) { return new Word(word, CS.Origin, CS.XAxis, CS.YAxis); }


        public Word Transform(CoordinateSystem newCS)
        {
            for (int i = 0; i < letters.Length; i++) for (int j = 0; j < letters[i].Count; j++)
                { letters[i][j] = letters[i][j].Transform(newCS) as PolyCurve; }
            return this;
        }

        public Word Translate(Vector vec, double distance)
        {
            for (int i = 0; i < letters.Length; i++) for (int j = 0; j < letters[i].Count; j++)
                { letters[i][j] = letters[i][j].Translate(vec, distance) as PolyCurve; }
            return this;
        }

        public Word Scale(double factor)
        {
            for (int i = 0; i < letters.Length; i++) for (int j = 0; j < letters[i].Count; j++)
                { letters[i][j] = letters[i][j].Scale(cs.Origin, cs.Origin.Translate(cs.XAxis, 1) as Point, cs.Origin.Translate(cs.XAxis, factor) as Point) as PolyCurve; }
            return this;
        }

        public List<PolyCurve> display(double factor)
        {
            List<PolyCurve> lines = new List<PolyCurve>();
            for (int i = 0; i < letters.Length; i++) for (int j = 0; j < letters[i].Count; j++)
                { lines.Add(letters[i][j].Scale(cs.Origin, cs.Origin.Translate(cs.XAxis, 1) as Point, cs.Origin.Translate(cs.XAxis, factor) as Point) as PolyCurve); }
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
            if (disposing) for (int i = 0; i < letters.Length; i++) letters[i].Dispose();
            disposed = true;
        }

    }

}
