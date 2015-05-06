using System;

namespace DXFLibrary
{
	/// <summary>
	/// The representation of LINE entity.
	/// </summary>
	class Line : DXFLibrary.Entity
	{
		public Line(string layer,double xi, double yi,double zi, double xf, double yf, double zf):base("LINE",layer)
		{
			this.dataAcceptanceList.AddRange(new int[10] { 39, 10, 20, 30, 11, 21, 31, 210, 220, 230 });
			this.AddReplace(10,xi);
			this.AddReplace(20,yi);
            this.AddReplace(30, zi);
			this.AddReplace(11,xf);
			this.AddReplace(21,yf);
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

		public void setInitialPoint(double x, double y)
		{
			this.AddReplace(10,x);
			this.AddReplace(20,y);
		}
		public void setFinalPoint(double x, double y)
		{
			this.AddReplace(11,x);
			this.AddReplace(21,x);
            this.AddReplace(31,x);
		}
	}
}
