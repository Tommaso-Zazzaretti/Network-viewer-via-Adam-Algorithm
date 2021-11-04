class Node {
	constructor(id,posx,posy,graph) {
		this.id = id;            //Node id in G(V,E)
		this.x = posx;           //Current x position
		this.y = posy;           //Current y position
		this.CIRCLE = null;      //d3 circle reference
		this.TEXT   = null;      //d3 text reference
		this.isDragged = false;  //Drag/drop boolean status variable
		this.edgeList = new Array();
		this.isDrawn == false;
		this.G = graph;
	}

	draw(nodeColor,textColor,text) {
		if(!this.isDrawn){
			let SVG = this.G.SVG_BOX;
			let WIDTH = this.G.W;
			let HEIGHT = this.G.H;
			this.isDrawn=true;
			let PIXEL_UNIT = GRAPHIC_PROPERTIES.PIXEL_UNIT;
			//Circle elem
			this.CIRCLE = SVG.append("circle")
			.attr("id","circle"+this.id)
			.attr("cx",this.x*PIXEL_UNIT+WIDTH/2)
			.attr("cy",this.y*PIXEL_UNIT+HEIGHT/2)
			.attr("r",GRAPHIC_PROPERTIES.NODES_RADIUS)
			.attr("fill",nodeColor)
			.attr("stroke",GRAPHIC_PROPERTIES.NODES_BORDER_COLOR)
			.attr("stroke-width",2);
			//Text elem
			this.TEXT = SVG.append("text")
			.attr("id","text"+this.id)
			.attr("x",this.x*PIXEL_UNIT+WIDTH/2)
			.attr("y", this.y*PIXEL_UNIT+HEIGHT/2)
			.attr("fill",textColor)
			.text(text.toString())
			.style("font-size", GRAPHIC_PROPERTIES.NODES_RADIUS*1.15)
			.attr("text-anchor", "middle")
			.attr("font-weigth", "bold")
			.attr("dominant-baseline", "central") 
			.style("font-family","Arial");

			//Mouse over event setting:
			this.CIRCLE.on("mouseover", function(d) { d3.select(this).style("cursor", "pointer"); });
			this.TEXT.on("mouseover",   function(d) { d3.select(this).style("cursor", "pointer"); });

			//Mouse out event setting:
			this.CIRCLE.on("mouseout", function(d) { d3.select(this).style("cursor", "default"); });
			this.TEXT.on("mouseout",   function(d) { d3.select(this).style("cursor", "default"); });

			let dragHandler = d3.drag()
				.on('start',function() {  	this.isDragged=true;  }.bind(this)  )
				.on('end'  ,function() {  	this.isDragged=false; }.bind(this)  )
				.on('drag', function() {    let WIDTH = this.G.W;
											let HEIGHT = this.G.H;
											//Node x,y position update
											this.x = (d3.event.x-WIDTH/2)/PIXEL_UNIT; 
											this.y = (d3.event.y-HEIGHT/2)/PIXEL_UNIT;
											let X = this.x*PIXEL_UNIT+WIDTH/2; 
											let Y = this.y*PIXEL_UNIT+HEIGHT/2;
											//Circle x,y update
											this.CIRCLE.attr("cx",X).attr("cy",Y);
											//Text x,y update
											this.TEXT.attr("x",X).attr("y",Y); 
											//Update edge with dragged node
											for(let i=0;i<this.edgeList.length;i++){
												this.edgeList[i].updateDrawing(null);
											}
				  	   				   }.bind(this));
			//Drag event setting:
			this.CIRCLE.call(dragHandler);
			this.TEXT.call(dragHandler);
		}
	}

	removeEvents(){
		this.CIRCLE.on("mouseover",null);
		this.CIRCLE.on("mouseout",null);
		this.CIRCLE.on('mousedown.drag', null);
		
		this.TEXT.on("mouseover",null);
		this.TEXT.on("mouseout",null);
		this.TEXT.on('mousedown.drag', null);
	}

	delete() {
		if(this.CIRCLE!=null) { 
			this.CIRCLE.remove(); 
			this.CIRCLE = null;
		}
		if(this.TEXT!=null) { 
			this.TEXT.remove(); 
			this.TEXT = null;
		}
		if(this.isDragged) {
			this.isDragged = false;
		}
		this.isDrawn = false;
	}

	updateDrawing() {
		if((!this.isDragged) && (this.CIRCLE!=null) && (this.TEXT!=null)){
			let width = this.G.W;
			let height = this.G.H;
			let PIXEL_UNIT = GRAPHIC_PROPERTIES.PIXEL_UNIT;
			this.CIRCLE.attr("cx",this.x*PIXEL_UNIT+width/2).attr("cy",this.y*PIXEL_UNIT+height/2);
			this.TEXT.attr("x",this.x*PIXEL_UNIT+width/2).attr("y",this.y*PIXEL_UNIT+height/2)	
		}
	}
}