class Graph {
	constructor(V,E,name,svgBox,w,h) {
		this.V   = [...V]
		this.E   = [...E]
		this.ADJ = new Array(V.length);
		this.vIndex = this.V.length;
		this.eIndex = this.E.length;
		this.indexMapV = {};  //Dictionary key: NodeId -> NodeIndex in this.V | ids: [3,67,4,82] => [0,1,2,3];
		this.indexMapE = {};  //Dictionary key: EdgeId -> EdgeIndex in this.E | ids: [[2,23],[1,43]] => [0,1];
		this.edgeList = new Array(E.length);
		this.isDrawn = false;
		this.graphName = name;
		
		//Nodes objects init
		for(let i=0;i<V.length;i++){
			let nodeId = V[i];
			let x = Math.random();	let y = Math.random();
			this.V[i] = new Node(nodeId,x,y,this);
			this.indexMapV[nodeId] = i;
			this.ADJ[i] = new Array();
		}
		//Edges objects and adj lists init
		for(let i=0;i<E.length;i++){
			let source = E[i][0];
			let target = E[i][1];
			if(source>target){ let tmp=source; source=target; target=tmp; }
			let sourceIndex = this.indexMapV[source];
			let targetIndex = this.indexMapV[target];
			this.E[i] = new Edge(this.V[sourceIndex],this.V[targetIndex],this);
			this.edgeList[i] = [sourceIndex,targetIndex];
			this.indexMapE[source+"-"+target] = i;
			this.ADJ[sourceIndex].push(targetIndex);
			this.ADJ[targetIndex].push(sourceIndex);
		}
		//console.log(this.indexMapE)
		this.SVG_BOX = svgBox;
		this.W = w;
		this.H = h;
	}
	

	draw(){
		if(!this.isDrawn){
			this.isDrawn = true;
			for(let i=0;i<this.E.length;i++){
				this.E[i].draw(GRAPHIC_PROPERTIES.EDGES_COLOR, GRAPHIC_PROPERTIES.EDGES_WIDTH);
			}
			this.center = new Node("graphCenter",0,0,this);
			for(let i=0;i<this.V.length;i++){
				this.V[i].draw(GRAPHIC_PROPERTIES.GRAPH_NODES_COLOR, GRAPHIC_PROPERTIES.NODES_TEXT_COLOR, this.V[i].id);
			}
		}
	}

	updateDrawing(){
		this.updateDrawingE();
		this.updateDrawingV();
	}

	updateDrawingE(){
		if(this.isDrawn){
			for(let i=0;i<this.E.length;i++){
				this.E[i].updateDrawing(null);
			}
		}
	}

	updateDrawingV(){
		if(this.isDrawn){
			for(let i=0;i<this.V.length;i++){
				this.V[i].updateDrawing();
			}
		}
	}

	toString(color){
		let first = "    "+this.graphName+" INFO:";
		// V
		let Vline = new Array();
		let str="    V = [";  
		for(let i=0;i<this.vIndex-1;i++){
			if(i%20==0 && i!=0){Vline.push(str); str="         ";}
			str+=this.V[i].id.toString().split("-")[0]+",";
		}
		str+=this.V[this.vIndex-1].id.toString().split("-")[0]+"]";
		Vline.push(str);
		// E
		let Eline = new Array();
		str="    E = ["; 
		for(let i=0;i<this.eIndex-1;i++){
			if(i%8==0 && i!=0){Eline.push(str); str="         ";}
			str+="["+this.E[i].s.id.toString().split("-")[0]+","+this.E[i].t.id.toString().split("-")[0]+"],";
		}
		str+="["+this.E[this.eIndex-1].s.id.toString().split("-")[0]+","+this.E[this.eIndex-1].t.id.toString().split("-")[0]+"]]";
		Eline.push(str);
		let maxLen = Math.max(Math.max(...Vline.map(e=>e.length)),Math.max(...Eline.map(e=>e.length)),first.length)+4;
		let padding = " ".repeat(maxLen);
		for(let i=0;i<Vline.length;i++){Vline[i]+=" ".repeat(maxLen-Vline[i].length); }
		for(let i=0;i<Eline.length;i++){Eline[i]+=" ".repeat(maxLen-Eline[i].length); }
		first+=" ".repeat(maxLen-first.length);
		let outputString=""+padding+"\n"+first+"\n"+padding+"\n";
		for(let i=0;i<Vline.length;i++){
			outputString+=Vline[i]+"\n";
		}
		outputString+=""+padding+"\n";
		for(let i=0;i<Eline.length;i++){
			outputString+=Eline[i]+"\n";
		}
		outputString+=""+padding+"\n";
		if(color!=undefined){console.log("%c"+outputString,'background: #000; color:'+color)};
		return outputString;
	}
}