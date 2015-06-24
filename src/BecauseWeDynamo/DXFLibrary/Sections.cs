using System;

namespace DXFLibrary
{
    class Header : Section
    {
        //**CONSTRUCTOR
        internal Header() : base("HEADER") { }

        //**METHODS
        public int VariableCount() { return this.elements.Count; }
        public Variable getVariable(int index) { return (Variable)this.elements[index]; }
        public object valueOf(string varName)
        {
            foreach (Variable v in this.elements) if (v.VarName == varName) return v.Value;
            return null;
        }
        public void addVariable(Variable v) { this.elements.Add(v); }
    }
    class Blocks : Section
	{
		public Blocks():base("BLOCKS") {}
		public Blocks(string s):base(s,true) {}
		public void addBlock(Block b) { this.elements.Add(b); }
	}
    class Entities : Section
    {
        public Entities() : base("ENTITIES") { }
        public void AddEntity(Entity e) { this.AddElement(e); }
    }
    class Tables : Section
    {
        public Tables() : base("TABLES") { }
        public void addTable(Table t) { this.AddElement(t); }
    }
}

