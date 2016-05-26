using System;

namespace Fabrication.DXFLibrary 
{
    class Arc : Entity
    {
		public Arc(double x, double y, double radius, double startAngle, double endAngle, string layer):base("ARC",layer)
		{
			this.dataAcceptanceList.AddRange(new int[] { 39, 10, 20, 30, 40, 100, 50, 51, 210, 220, 230 });
			this.AddReplace(10, x);
			this.AddReplace(20, y);
			this.AddReplace(40, radius);
            this.AddReplace(50, startAngle);
            this.AddReplace(51, endAngle);
		}

        /// <summary>
        /// right hand rule for angle direction
        /// for normal vector (nx,ny,nz)
        /// zero vector = (-sgn(nz)*ny, sgn(nz)*nx, 0)
        /// </summary>
        /// <param name="x">Center X-coordinate</param>
        /// <param name="y">Center Y-coordinate</param>
        /// <param name="z">Center Z-coordinate</param>
        /// <param name="radius">Radius</param>
        /// <param name="startAngle">start angle</param>
        /// <param name="endAngle">end angle</param>
        /// <param name="nx">Normal Vector X-coordinate</param>
        /// <param name="ny">Normal Vector Y-coordinate</param>
        /// <param name="nz">Normal Vector Z-coordinate</param>
        /// <param name="layer">layer</param>
        public Arc(double x, double y, double z, double radius, double startAngle, double endAngle, double nx, double ny, double nz, string layer):base("ARC", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 39, 10, 20, 30, 40, 100, 50, 51, 210, 220, 230 });
            this.AddReplace(10, x);
            this.AddReplace(20, y);
            this.AddReplace(30, z);
            this.AddReplace(40, radius);
            this.AddReplace(50, startAngle);
            this.AddReplace(51, endAngle);
            this.AddReplace(210, nx);
            this.AddReplace(220, ny);
            this.AddReplace(230, nz);
        }
    }

    class Circle : Entity
    {
        public Circle(double x, double y, double radius, string layer)
            : base("CIRCLE", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 39, 10, 20, 30, 40, 210, 220, 230 });
            this.AddReplace(10, x);
            this.AddReplace(20, y);
            this.AddReplace(40, radius);
        }
        public Circle(double x, double y, double z, double radius, double nx, double ny, double nz, string layer)
            : base("CIRCLE", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 39, 10, 20, 30, 40, 210, 220, 230 });
            this.AddReplace(10, x);
            this.AddReplace(20, y);
            this.AddReplace(30, z);
            this.AddReplace(40, radius);
            this.AddReplace(210, nx);
            this.AddReplace(220, ny);
            this.AddReplace(230, nz);
        }
    }

    class Line : Entity
    {
        public Line(string layer, double xi, double yi, double zi, double xf, double yf, double zf)
            : base("LINE", layer)
        {
            this.dataAcceptanceList.AddRange(new int[10] { 39, 10, 20, 30, 11, 21, 31, 210, 220, 230 });
            this.AddReplace(10, xi);
            this.AddReplace(20, yi);
            this.AddReplace(30, zi);
            this.AddReplace(11, xf);
            this.AddReplace(21, yf);
            this.AddReplace(31, zf);
        }

        public Line(string layer, double xi, double yi, double xf, double yf)
            : base("LINE", layer)
        {
            this.dataAcceptanceList.AddRange(new int[10] { 39, 10, 20, 30, 11, 21, 31, 210, 220, 230 });
            this.AddReplace(10, xi);
            this.AddReplace(20, yi);
            this.AddReplace(11, xf);
            this.AddReplace(21, yf);
        }

        public void setInitialPoint(double x, double y, double z =0)
        {
            this.AddReplace(10, x);
            this.AddReplace(20, y);
            if (z > 0) this.AddReplace(30, z);
        }
        public void setFinalPoint(double x, double y, double z = 0)
        {
            this.AddReplace(11, x);
            this.AddReplace(21, y);
            if (z > 0) this.AddReplace(31, z);
        }
    }
    class PolyLine : Entity
    {
        public PolyLine(string layer)
            : base("POLYLINE", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 66, 10, 20, 30, 39, 70, 40, 41, 71, 72, 73, 74, 75, 210, 220, 230, 66 });
            this.AddElement(new SeqEnd(layer));
            this.AddReplace(66, (short)1);
        }
        public void AddVertex(Vertex v) { this.InsertElement(this.ElementCount() - 1, v); }
        public void AddVertex(double x, double y) { this.AddVertex(new Vertex(x, y, this.layer)); }
        public void AddVertex(double x, double y, double z) { this.AddVertex(new Vertex(x, y, z, this.layer)); }
    }

    class Text : Entity
    {
        public Text(string text, double x, double y, double height, string layer)
            : base("TEXT", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 39, 10, 20, 30, 40, 1, 50, 41, 51, 7, 71, 72, 11, 21, 31, 210, 220, 230, 73 });
            this.AddReplace(10, x);
            this.AddReplace(20, y);
            this.AddReplace(1, text);
            this.AddReplace(40, height);
        }
        public Text(string text, double x, double y, double z, double height, string layer)
            : base("TEXT", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 39, 10, 20, 30, 40, 1, 50, 41, 51, 7, 71, 72, 11, 21, 31, 210, 220, 230, 73 });
            this.AddReplace(10, x);
            this.AddReplace(20, y);
            this.AddReplace(30, z);
            this.AddReplace(1, text);
            this.AddReplace(40, height);
        }
        /// <summary>
        /// Horizontal text justification type (optional, default = 0) integer codes (not bit-coded)
        /// 0 = Left; 1= Center; 2 = Right
        /// 3 = Aligned (if vertical alignment = 0)
        /// 4 = Middle (if vertical alignment = 0)
        /// 5 = Fit (if vertical alignment = 0)
        /// </summary>
        public void SetHorizontalJustification(short flag)
        {
            this.AddReplace(72, flag);
        }
        ///<summary>
        /// Vertical text justification type (optional, default = 0): integer codes (not bit-coded):
        /// 0 = Baseline; 1 = Bottom; 2 = Middle; 3 = Top
        /// </summary>
        public void SetVerticalJustification(short flag) { this.AddReplace(73, flag); }
        /// <summary>
        /// The secound alignament coords.
        /// </summary>
        /// <param name="x">x coord (horizontal)</param>
        /// <param name="y">y coord (vertical)</param>
        public void SetSecoundAlignament(double x, double y, double z = 0)
        {
            this.AddReplace(11, x);
            this.AddReplace(21, y);
            this.AddReplace(31, z);
        }

        public void Rotate(double angle) { this.AddReplace(50, angle); }
    }

    class Vertex : Entity
    {
        public Vertex(double x, double y, string layer)
            : base("VERTEX", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 10, 20, 30, 70, 40, 41, 42, 50, 71, 72, 73, 74, });
            this.AddReplace(10, x);
            this.AddReplace(20, y);
        }
        public Vertex(double x, double y, double z, string layer)
            : base("VERTEX", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 10, 20, 30, 70, 40, 41, 42, 50, 71, 72, 73, 74, });
            this.AddReplace(10, x);
            this.AddReplace(20, y);
            this.AddReplace(30, z);
        }
    }

    class Insert : Entity
    {
        public Insert(string block, double x, double y, string layer)
            : base("INSERT", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 66, 2, 10, 20, 30, 41, 42, 43, 50, 70, 71, 44, 45, 210, 220, 230 });
            this.AddData(new Data(2, block));
            this.AddData(new Data(10, x));
            this.AddData(new Data(20, y));
        }
        public Insert(string block, double x, double y, double z, string layer)
            : base("INSERT", layer)
        {
            this.dataAcceptanceList.AddRange(new int[] { 66, 2, 10, 20, 30, 41, 42, 43, 50, 70, 71, 44, 45, 210, 220, 230 });
            this.AddData(new Data(2, block));
            this.AddData(new Data(10, x));
            this.AddData(new Data(20, y));
            this.AddData(new Data(30, z));
        }
    }
    class SeqEnd : Entity
    {
        public SeqEnd(string layer) : base("SEQEND", layer) { }
    }
}