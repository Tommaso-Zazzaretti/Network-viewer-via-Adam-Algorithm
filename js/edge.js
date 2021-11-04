class Edge {
	constructor(node1,node2,graph) {
		this.s = node1;      //Source node obj reference
		this.t = node2;      //Target node obj reference
		this.id = "edge"+node1.id+"-"+node2.id;
		this.PATH = null;    //d3 path reference
		this.s.edgeList.push(this);     //Push this edge in s edge list (to drag)
		this.t.edgeList.push(this);		//Push this edge in t edge list (to drag)
		this.isDrawn == false;
		this.G = graph;
	}

	draw(edgeColor,edgeWidth,edgeOpacity) {
		let SVG = this.G.SVG_BOX;
		let WIDTH = this.G.W;
		let HEIGHT = this.G.H;
		if(!this.isDrawn){
			this.isDrawn = true;
			let path = d3.path();
			let PIXEL_UNIT = GRAPHIC_PROPERTIES.PIXEL_UNIT;
			path.moveTo(this.s.x*PIXEL_UNIT+WIDTH/2, this.s.y*PIXEL_UNIT+HEIGHT/2);
			path.lineTo(this.t.x*PIXEL_UNIT+WIDTH/2, this.t.y*PIXEL_UNIT+HEIGHT/2);
			this.PATH = SVG.append("path")
			.attr("d", path.toString())
			.attr("id",this.id)
			.attr("stroke",edgeColor)
			.attr("fill","none")
			.attr("stroke-width",edgeWidth)
			.attr("opacity",edgeOpacity); 
		}
	}

	delete() {
		if(this.PATH!=null) { 
			this.PATH.remove(); 
			this.PATH = null;
		}
	}

	//Update drawing with a new d3 path passed as argument. If is NULL, compute default path between S and T
	updateDrawing(newPath) {
		if(this.PATH != null){  //Update only if exist d3 element 
			if(newPath == null){
				let width = this.G.W;
				let height = this.G.H;
				newPath = d3.path();
				let PIXEL_UNIT = GRAPHIC_PROPERTIES.PIXEL_UNIT;
				newPath.moveTo(this.s.x*PIXEL_UNIT+width/2, this.s.y*PIXEL_UNIT+height/2);
				newPath.lineTo(this.t.x*PIXEL_UNIT+width/2, this.t.y*PIXEL_UNIT+height/2);
			}
			this.PATH.attr("d",newPath.toString());
		}
	}
}