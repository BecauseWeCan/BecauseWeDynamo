using System;

namespace DXFLibrary
{
    abstract class Entity : Element
    {
        //**CONSTRUCTOR
        public Entity(string name, string layer)
            : base()
        {
            this.dataAcceptanceList.AddRange(new int[] { 0, 5, 102, 330, 360, 100, 67, 410, 8, 6, 62, 370, 48, 60, 92, 310 });
            this.startTag = new Data(0, name);
            this.AddData(new Data(8, layer));
        }

        //PROPERTY
        public string layer
        {
            get { return (string)(this.GetDataFor(8)).data; }
            set { this.AddReplace(8, s); }
        }
    }

	class Block : Element
	{
		public Block():base()
		{
			this.dataAcceptanceList.AddRange(new int[14]{ 0, 5, 102, 330, 100, 8, 70,10 ,20 ,30 ,3, 1, 4, 2});
			this.startTag = new Data(0,"BLOCK");
		}
		public Block(string s):base(s){}
		
        public void SetEndBlk(EndBlk eb) 
		{
			if(this.elements.Count > 0 && ((Element)this.elements[this.elements.Count-1]).Name == "ENDBLK") this.elements.RemoveAt(this.elements.Count-1);
			this.elements.Add(eb);
		}
		public void AddEntity(Entity e)
		{
			if(this.data.Count==0) this.elements.Add(e);
			else this.elements.Insert(this.elements.Count-1,e);
		}
		public void SetLayer(string l)
		{
			int ind = this.GetIndexFor(8);
			if(ind>-1)
			{
				this.data.RemoveAt(ind);
				this.data.Insert(ind,new Data(8,l));
			}
			else this.data.Add(new Data(8,l));
		}
		public void SetPosition(double x, double y, double z)
		{
			Data dx,dy,dz;
			bool swx = false,swy = false,swz = false;
			foreach(Data d in this.data)
			{
				if(d.code==10) 
				{
					dx = d;
					swx = true;
				}
				if(d.code==20) 
				{
					dy = d;
					swy=true;
				}
				if(d.code==30) 
				{
					dz = d;
					swz = true;
				}
			}
			if(swx) dx.data = x;
			else
			{
				dx.code = 10;
				dx.data = x;
				this.data.Add(dx);
			}
			if(swy) dy.data = y;
			else
			{
				dy.code = 20;
				dy.data = y;
				this.data.Add(dy);
			}
			if(swz) dz.data = z;
			else
			{
				dy.code = 30;
				dy.data = y;
				this.data.Add(dy);
			}
		}
		public void SetName(string name) { this.AddReplace(2,name); }
		public void SetHandle(string handle) { this.AddReplace(5,handle); }
		public void SetFlag(short flag) { this.AddReplace(70,flag); }
	}

    class EndBlk : Element
    {
        public EndBlk()
        {
            this.dataAcceptanceList.AddRange(new int[6] { 0, 5, 8, 102, 330, 100 });
            this.startTag = new Data(0, "ENDBLK");
        }
        public EndBlk(Block b)
            : this()
        {
            if (b.GetIndexFor(5) != -1) this.AddData(b.GetDataFor(5));
            if (b.GetIndexFor(8) != -1) this.AddData(b.GetDataFor(8));
        }
    }

    class Table : Element
    {
        public Table(string type)
            : base()
        {
            this.dataAcceptanceList.AddRange(new int[] { 0, 2, 5, 102, 360, 102, 330, 100, 70 });
            this.AddReplace(2, type);
            this.startTag = new Data(0, "TABLE");
            this.endTag = new Data(0, "ENDTAB");
        }
        public void AddTableEntry(TableEntry te)
        {
            this.AddElement(te);
        }
    }

    class TableEntry : Element
    {
        internal TableEntry(string type)
            : base()
        {
            this.dataAcceptanceList.AddRange(new int[] { 0, 5, 105, 102, 330, 360, 100 });
            this.startTag = new Data(0, type);
        }
        internal TableEntry(string dxfcode, bool withDxfCode) : base(dxfcode) { }
    }

    class Variable : Element
    {
        public Variable(string nume, int dataType, object val)
        {
            startTag = new Data(0, 0);
            endTag = new Data(0, 0);
            this.data.Add(new Data(9, nume));
            this.data.Add(new Data(dataType, val));
        }
        public string VarName { get { return (string)((Data)this.data[0]).data; } }
        public object Value
        { get { return ((Data)this.data[1]).data; } }
    }
}
