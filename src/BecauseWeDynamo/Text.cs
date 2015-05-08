using System;
using System.Collections.Generic;
using Autodesk.DesignScript.Geometry;

namespace Fabrication
{
    public class Letter
    {
        private PolyCurve pl;
        private char ltr;

        public PolyCurve drawn { get { return pl; } }

        internal Letter(char c)
        {
            ltr = c;
            switch (c)
            {
                case 't':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1,  1, -1), 
                        Point.ByCoordinates(1, -1), 
                        Point.ByCoordinates(0, -1), 
                        Point.ByCoordinates(0, -5)
                    });
                    break;
                case 's':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(1, -2), 
                        Point.ByCoordinates(0, -1), 
                        Point.ByCoordinates(-1, -2), 
                        Point.ByCoordinates(1, -4),
                        Point.ByCoordinates(0, -5), 
                        Point.ByCoordinates(-1, -4)
                    });
                    break;
                case 'e':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -3), 
                        Point.ByCoordinates(1, -3), 
                        Point.ByCoordinates(1, -2), 
                        Point.ByCoordinates(0, -1),
                        Point.ByCoordinates(-1, -2), 
                        Point.ByCoordinates(-1, -4),
                        Point.ByCoordinates(0, -5),
                        Point.ByCoordinates(1, -5)
                    });
                    break;
                case '-':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -3), 
                        Point.ByCoordinates(1, -3)
                    });
                    break;
                case '0':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(0, -1), 
                        Point.ByCoordinates(1, -2), 
                        Point.ByCoordinates(1, -4), 
                        Point.ByCoordinates(0, -5),
                        Point.ByCoordinates(-1, -4), 
                        Point.ByCoordinates(-1, -2),
                        Point.ByCoordinates(0, -1)
                    });
                    break;
                case '1':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -2), 
                        Point.ByCoordinates(0, -1), 
                        Point.ByCoordinates(0, -5), 
                        Point.ByCoordinates(1, -5),
                        Point.ByCoordinates(-1, -5)
                    });
                    break;
                case '2':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -1), 
                        Point.ByCoordinates(1, -1), 
                        Point.ByCoordinates(1, -2), 
                        Point.ByCoordinates(-1, -5),
                        Point.ByCoordinates(1, -5)
                    });
                    break;
                case '3':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -1), 
                        Point.ByCoordinates(1, -1), 
                        Point.ByCoordinates(0, -2),
                        Point.ByCoordinates(1, -2),
                        Point.ByCoordinates(1, -4),
                        Point.ByCoordinates(0, -5),
                        Point.ByCoordinates(-1, -5)
                    });
                    break;
                case '4':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(1, -4), 
                        Point.ByCoordinates(-1, -4), 
                        Point.ByCoordinates(1, -1),
                        Point.ByCoordinates(1, -5)
                    });
                    break;
                case '5':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(1, -1), 
                        Point.ByCoordinates(-1, -1), 
                        Point.ByCoordinates(-1, -2),
                        Point.ByCoordinates(1, -2),
                        Point.ByCoordinates(1, -4),
                        Point.ByCoordinates(0, -5),
                        Point.ByCoordinates(-1, -5)
                    });
                    break;
                case '6':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -3), 
                        Point.ByCoordinates(0, -2),
                        Point.ByCoordinates(1, -2),
                        Point.ByCoordinates(1, -4),
                        Point.ByCoordinates(0, -5),
                        Point.ByCoordinates(-1, -5),
                        Point.ByCoordinates(-1, -2),
                        Point.ByCoordinates(0, -1),
                        Point.ByCoordinates(1, -1)
                    });
                    break;
                case '7':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -2), 
                        Point.ByCoordinates(-1, -1),
                        Point.ByCoordinates(1, -1),
                        Point.ByCoordinates(1, -2),
                        Point.ByCoordinates(0, -5)
                    });
                    break;
                case '8':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(-1, -1),
                        Point.ByCoordinates(1, -1),
                        Point.ByCoordinates(1, -2),
                        Point.ByCoordinates(-1, -3),
                        Point.ByCoordinates(-1, -5),
                        Point.ByCoordinates(1, -5),
                        Point.ByCoordinates(1, -3),
                        Point.ByCoordinates(-1, -2),
                        Point.ByCoordinates(-1, -1)
                    });
                    break;
                case '9':
                    pl = PolyCurve.ByPoints(new List<Point>
                    {
                        Point.ByCoordinates(1, -3), 
                        Point.ByCoordinates(0, -4),
                        Point.ByCoordinates(-1, -3),
                        Point.ByCoordinates(-1, -2),
                        Point.ByCoordinates(0, -1),
                        Point.ByCoordinates(1, -2),
                        Point.ByCoordinates(1, -4),
                        Point.ByCoordinates(0, -5),
                        Point.ByCoordinates(-1, -5)
                    });
                    break;

                default:
                    pl = null;
                    break;
            } // end switch
        }

        public Letter Transform(CoordinateSystem newCS)
        {
            pl = pl.Transform(newCS) as PolyCurve;
            return this;
        }

        public Letter Translate(Vector vec, double distance)
        {
            pl = pl.Translate(vec, distance) as PolyCurve;
            return this;
        }

        public Letter Scale(CoordinateSystem cs, double factor)
        {
            Point from = cs.Origin.Translate(cs.XAxis, 1) as Point;
            Point to = cs.Origin.Translate(cs.XAxis, factor) as Point;
            pl = pl.Scale(cs.Origin, from, to) as PolyCurve;
            return this;
        }


        

    }// end letter

    public class Word
    {
        private string word;
        private Letter[] Letters;
        private CoordinateSystem cs;

        public Letter[] letters { get { return Letters; } }


        internal Word(string str, Point origin, Vector X, Vector Y)
        {
            cs = CoordinateSystem.ByOriginVectors(origin, X, Y);
            word = str;
            char[] ltrs = str.ToCharArray();
            Letters = new Letter[ltrs.Length];
            double dblLtr = (ltrs.Length-1)*(-1.5);
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

        public PolyCurve[] display(double factor)
        {
            PolyCurve[] lines = new PolyCurve[Letters.Length];
            for (int i = 0; i < Letters.Length; i++)
            {
                lines[i] = Letters[i].drawn.Scale(cs.Origin, cs.Origin.Translate(cs.XAxis, 1) as Point, cs.Origin.Translate(cs.XAxis, factor) as Point) as PolyCurve;
            }

            return lines;
        }

    }

}
