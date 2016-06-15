using System;

namespace Fabrication.DXFLibrary
{

    /// <summary>
    /// The representation of a dxf document.
    /// </summary>
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
        internal void SetHeader(Header h)
        {
            this.header = h;
        }
        internal void SetEntities(Entities e)
        {
            this.entities = e;
        }
        internal void SetBlocks(Blocks b)
        {
            this.blocks = b;
        }
        internal void SetTables(Tables t)
        {
            this.tables = t;
        }
        internal void add(Entity e)
        {
            this.entities.AddEntity(e);
        }
        internal void add(Block b)
        {
            this.blocks.addBlock(b);
        }
    }

    class UnexpectedElement : Exception
    {
        internal UnexpectedElement() : base("Unrecognized DXF element") { }
    }
}
