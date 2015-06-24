using System;

namespace DXFLibrary
{
	class Document
	{
		internal Header header;
		internal Entities entities;
		internal Blocks blocks;
		internal Tables tables;
		public Document()
		{
			this.entities = new Entities();
		}

		public void SetHeader(Header h)
		{
			this.header = h;
		}
		public void SetEntities(Entities e)
		{
			this.entities = e;
		}
		public void SetBlocks(Blocks b)
		{
			this.blocks = b;
		}
		public void SetTables(Tables t)
		{
			this.tables = t;
		}
		public void add(Entity e)
		{
			this.entities.AddEntity(e);
		}
		public void add(Block b)
		{
			this.blocks.addBlock(b);
		}
	}

    class UnexpectedElement : Exception
    {
        internal UnexpectedElement() : base("Unrecognized DXF element") { }
    }
}
