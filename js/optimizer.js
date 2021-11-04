class Optimizer {
	constructor(graph,updateStartIdx,updateEndIdx,visualCriteria,hyperParameters,max_iter_num,FRAME_MSEC,f){
		this.nodes = graph.V;
		this.edges = graph.edgeList;
		this.graph = graph;
		this.adj   = graph.ADJ;
		this.startIndex = updateStartIdx;
		this.endIndex   = updateEndIdx;
		this.visualCriteria  = [];
		this.hyperParameters  = [];
		this.lossFunctionsValues = {};
		//Init data structures according to the required visual criteria
		this.initDataStructures(visualCriteria,hyperParameters);
	}

	initDataStructures(visualCriteria,hyperParameters){
		for(let i=0;i<visualCriteria.length;i++){
			let name = visualCriteria[i].name;
			if(name=="STRESS"){
				this.addStress(visualCriteria[i].learningRate)
			}
			if(name=="IDEAL_EDGE_LENGTH"){
				this.addEdgeLenght(visualCriteria[i].learningRate,hyperParameters.EDGES_LENGTH);
			}
			if(name=="ANGLE_RESOLUTION"){
				this.addAngleResolution(visualCriteria[i].learningRate,hyperParameters.ANR_SENSITIVITY)
			}
			if(name=="NEIGHBORS_PRESERVATION"){
				this.addNeighborPreservation(visualCriteria[i].learningRate,hyperParameters.NEIGH_HALF_MARGIN)
			}
			if(name=="VERTEX_RESOLUTION"){
				this.addVertexResolution(visualCriteria[i].learningRate,(this.nodes.length*10/50),hyperParameters.VERTEX_SCALING)
			}
		}
		//Momentum struct init
		if(hyperParameters.MOMENTUM.enable){
			this.MOMENTUM_ENABLE = true;
			this.MOMENTUM = hyperParameters.MOMENTUM.value;
			this.hyperParameters.push({name:"MOMENTUM", value: hyperParameters.MOMENTUM.value})
			this.initExpWeightedAveragesVelocity();
			if(this.MOMENTUM<0) { this.MOMENTUM = 0; }
			if(this.MOMENTUM>=1){ this.MOMENTUM = 0.999 }
		} else {
			this.MOMENTUM_ENABLE = false;
		}
		//ADAM second moment struct init
		if(hyperParameters.MOMENTUM2.enable){
			this.ADAM_ENABLE = true; 
			this.MOMENTUM2 = hyperParameters.MOMENTUM2.value;
			this.hyperParameters.push({name:"2TH_MOMENTUM", value: hyperParameters.MOMENTUM2.value})
			this.initExpWeightedAveragesVelocity2();
			if(this.MOMENTUM2<0) { this.MOMENTUM2 = 0; }
			if(this.MOMENTUM2>=1){ this.MOMENTUM2 = 0.999 }
		} else {
			this.MOMENTUM_ENABLE = false;
		}
		//Start learning rate init
		this.INIT_LR =hyperParameters.START_LEARNING_RATE;
		this.LEARNING_RATE = this.INIT_LR;
		this.hyperParameters.push({name:"LEARNING_RATE", value: this.LEARNING_RATE})
		//Exponential decay init
		if(hyperParameters.EXP_DECAY_LEARING_RATE.enable){
			this.EXP_DECAY_ENABLE = true;
			this.DECAY = hyperParameters.EXP_DECAY_LEARING_RATE.value;
		}
		//Gradients array init
		this.NODES_GRADIENTS = new Array(this.nodes.length);
		//Gradient module tresh init (stop condition)
		this.GRADIENT_MIN_MODULE = hyperParameters.GRADIENT_MIN_MODULE;
	}

	//--------------------------------------------STRESS---------------------------------------------------------
	addStress(learningRate){
		this.visualCriteria.push({name:"STRESS", learningRate: learningRate})
		this.initTheoreticDistMatrix()
		this.LEARNING_RATE_STRESS = learningRate;
		this.lossFunctionsValues.stress = 0.0;
	}
	removeStress(){
		this.visualCriteria = this.visualCriteria.filter(e=>e.name != "STRESS");
		this.theoreticDistMatrix = null;
		this.LEARNING_RATE_STRESS = null
		this.lossFunctionsValues.stress = null;
	}
	initTheoreticDistMatrix(){
		this.theoreticDistMatrix = new Array(this.nodes.length);
		for(let sourceNode=0;sourceNode<this.nodes.length;sourceNode++){
			this.theoreticDistMatrix[sourceNode] = this.dijsktra(this.adj,sourceNode);
		}
	}

	//--------------------------------------------IDEAL_EDGE_LENGTH-----------------------------------------------
	addEdgeLenght(learningRate,EDGES_LENGTH){
		this.visualCriteria.push({name:"IDEAL_EDGE_LENGTH", learningRate: learningRate})
		this.hyperParameters.push({name:"EDGES_LEN", value: EDGES_LENGTH})
		this.LEARNING_RATE_EDGELEN = learningRate;
		this.edgesLength = new Array(this.edges.lenght);

		for(let i=0;i<this.edges.length;i++){
			this.edgesLength[i] = EDGES_LENGTH;
		}
		this.lossFunctionsValues.edgeLen = 0.0;
	}
	removeEdgeLenght(){
		this.visualCriteria = this.visualCriteria.filter(e=>e.name != "IDEAL_EDGE_LENGTH");
		this.hyperParameters = this.hyperParameters.filter(e=>e.name != "EDGES_LEN");
		this.edgesLength = null;
		this.LEARNING_RATE_EDGELEN = null
		this.lossFunctionsValues.edgeLen = null;
	}
	setEdgesLengthArray(lenArray){
		this.edgesLength = lenArray;
	}
	scaleEdgesLengthArray(alpha){
		for(let i=0;i<this.edgesLength.length;i++){
			this.edgesLength[i]*=alpha;
		}
	}

	//--------------------------------------------ANGLE RESOLUTION-----------------------------------------------
	addAngleResolution(learningRate,ANR_SENSITIVITY){
		this.visualCriteria.push({name:"ANGLE_RESOLUTION", learningRate: learningRate})
		this.hyperParameters.push({name:"ANGLE_SENSITIVITY", value: ANR_SENSITIVITY})
		this.LEARNING_RATE_ANGLE = learningRate;
		this.ANR_SENSITIVITY = ANR_SENSITIVITY;
		this.lossFunctionsValues.angleRes = 0.0;
	}
	removeAngleResolution(learningRate,ANR_SENSITIVITY){
		this.visualCriteria = this.visualCriteria.filter(e=>e.name != "ANGLE_RESOLUTION");
		this.hyperParameters = this.hyperParameters.filter(e=>e.name != "ANGLE_SENSITIVITY");
		this.LEARNING_RATE_ANGLE = null;
		this.ANR_SENSITIVITY = null;
		this.lossFunctionsValues.angleRes = 0.0;
	}

	//--------------------------------------------NEIGHBORS PRESERVATION---------------------------------------
	addNeighborPreservation(learningRate,NEIGH_HALF_MARGIN){
		this.visualCriteria.push({name:"NEIGHBORS_PRESERVATION", learningRate: learningRate})
		this.hyperParameters.push({name:"NEIGHBORS_MARGIN", value: NEIGH_HALF_MARGIN})
		this.initAdjMatrix();
		this.initNeighborsNumber();
		this.initCumulativeJaccardIndex();
		this.LEARNING_RATE_NEIGH = learningRate;
		this.NEIGH_HALF_MARGIN   = NEIGH_HALF_MARGIN;
		this.lossFunctionsValues.neighPres = 0.0;
	}
	removeNeighborPreservation(){
		this.visualCriteria = this.visualCriteria.filter(e=>e.name != "NEIGHBORS_PRESERVATION");
		this.hyperParameters = this.hyperParameters.filter(e=>e.name != "NEIGHBORS_MARGIN");
		this.realAdjMatrix = null;
		this.nodesNeighborsNumber = null;
		this.cumulativeJaccardIndex = null;
		this.LEARNING_RATE_NEIGH = null;
		this.NEIGH_HALF_MARGIN = null;
		this.lossFunctionsValues.neighPres = 0.0;
	}
	initAdjMatrix(){ //O(V^2)
		this.realAdjMatrix = new Array(this.nodes.length);
		for(let i=0;i<this.nodes.length;i++){
			this.realAdjMatrix[i] = new Array(this.nodes.length);
			for(let j=0;j<this.nodes.length;j++){
				let k = this.adj[i].indexOf(j);
				if(k!=-1){ this.realAdjMatrix[i][j]=1; } 
				else {     this.realAdjMatrix[i][j]=0; }
			}
		}
	}
	initNeighborsNumber(){ //O(V)
		this.nodesNeighborsNumber = new Array(this.nodes.length);
		for(let i=0;i<this.nodes.length;i++){
			this.nodesNeighborsNumber[i] = this.adj[i].length;
		}
	}
	initCumulativeJaccardIndex(){ //O(V^2)
		// [B] Calcolo valori del Jaccard index in base a quanti elementi sono sbagliati
		this.cumulativeJaccardIndex = new Array();
		let positivesCount = [...this.realAdjMatrix].flat().reduce((a,b)=>a+b,0);
		let len = this.nodes.length;
		for(let i=0;i<(len*len);i++){
			let loss = i/(i+positivesCount);
			this.cumulativeJaccardIndex.push(loss);
		}
	}
	//--------------------------------------------VERTEX_RESOLUTION---------------------------------------
	addVertexResolution(learningRate,VERTEX_MAX_DISTANCE,VERTEX_SCALING){
		this.visualCriteria.push({name:"VERTEX_RESOLUTION", learningRate: learningRate})
		this.hyperParameters.push({name:"VERTEX_RESOLUTION_DIST", value: VERTEX_MAX_DISTANCE})
		this.LEARNING_RATE_VERTEX = learningRate;
		this.VERTEX_MAX_DISTANCE = VERTEX_MAX_DISTANCE;
		if(VERTEX_SCALING != null ){
			this.VERTEX_SCALING = VERTEX_SCALING;
		} else {
			this.VERTEX_SCALING = 1 / Math.sqrt(this.nodes.length);
		}
		this.hyperParameters.push({name:"VERTEX_RESOLUTION_SCALE", value: this.VERTEX_SCALING})

		this.lossFunctionsValues.vertexRes = 0.0;
	}
	removeVertexResolution(){
		this.visualCriteria = this.visualCriteria.filter(e=>e.name != "VERTEX_RESOLUTION");
		this.hyperParameters = this.hyperParameters.filter(e=>e.name != "VERTEX_RESOLUTION_DIST");
		this.hyperParameters = this.hyperParameters.filter(e=>e.name != "VERTEX_RESOLUTION_SCALE");
		this.LEARNING_RATE_VERTEX = null;
		this.VERTEX_MAX_DISTANCE = null;
		this.VERTEX_SCALING = null;
		this.lossFunctionsValues.vertexRes = 0.0;
	}
	//------------------------------------------MOMENTUM,GRADIENT-------------------------------------------
	initExpWeightedAveragesVelocity(){
		this.velocity = new Array(this.nodes.length);
		for(let i=0;i<this.nodes.length;i++){
			this.velocity[i]  = {dx: 0, dy: 0};
		}
	}
	initExpWeightedAveragesVelocity2(){
		this.velocity2 = new Array(this.nodes.length);
		for(let i=0;i<this.nodes.length;i++){
			this.velocity2[i]  = {dx: 0, dy: 0};
		}
	}
	initGradientToZero(){
		for(let i=0;i<this.nodes.length;i++){
			this.NODES_GRADIENTS[i] = {dx: 0, dy: 0};
		}
	}

	//--------------------------------------UTILITY FUNCTIONS-----------------------------------------------------------------
	euclideanDistance(p1,p2){
		let distx=p1.x-p2.x;
		let disty=p1.y-p2.y;
		return Math.sqrt(distx*distx + disty*disty);
	}
	scaleTheoreticDistMatrix(alpha){
		for(let i=0;i<this.nodes.length;i++){
			for(let j=0;j<this.nodes.length;j++){
				this.theoreticDistMatrix[i][j]*=alpha;
			}
		}
	}
	dijsktra(adj,source){
		let pred = new Array(adj.length);
		let dist = new Array(adj.length)
		let flag = new Array(adj.length)
		for(let i=0; i<adj.length; i++){ //Inizializzo distanze, nodiEsplorabili e predecessori
			dist[i] = 10000000;	flag[i] = true; pred[i] = undefined;
		}
		dist[source] = 0; 
		let explored = 0;
		while(explored<adj.length){ //Finchè non ho finito di fare tutti i nodi..
			let distMin = 10000000;	let nodeU = -1;
			for(let i=0;i<adj.length;i++){ //Seleziono nodo con distanza minore da esplorare
				if(dist[i]<=distMin && flag[i]==true){ nodeU=i; distMin=dist[i]; }
			}
			flag[nodeU]=false; 
			explored=explored+1;
			if(dist[nodeU]==10000000){ break; } //Se il nodo con dist minore ha dist infinita significa che non posso raggiungere niente
			for(let i=0;i<adj[nodeU].length;i++){ //Aggiorno distanze di tutti i successori del nodo esplorato
				let pathLen = dist[nodeU] + 1;
				if(pathLen < dist[adj[nodeU][i]]){
					dist[adj[nodeU][i]] = pathLen; pred[adj[nodeU][i]] = nodeU
				}
			} 
		}
		return dist;
	}

	stressIteration(){
		let loss=0.0
		//O((V^2)/2) Per ogni coppia (u,v) di nodi calcolo le derivate parziali xu yu xv yv e le sommo rispettivamente 
		for(let u=0;u<this.nodes.length-1;u++){ 
			for(let v=u+1;v<this.nodes.length;v++){
				let theoreticD = this.theoreticDistMatrix[u][v];
				let distx = this.nodes[u].x - this.nodes[v].x;
				let disty = this.nodes[u].y - this.nodes[v].y;
				let distUV = Math.sqrt((distx*distx)+(disty*disty));
				let alpha = ((2*(distUV-theoreticD))/(theoreticD*theoreticD)); //Per non calcolare la parte uguale del gradiente 4 volte
				if(u>=this.startIndex && u<this.endIndex){
					this.NODES_GRADIENTS[u].dx -= this.LEARNING_RATE_STRESS*alpha*(distx/distUV);
					this.NODES_GRADIENTS[u].dy -= this.LEARNING_RATE_STRESS*alpha*(disty/distUV);
				}
				if(v>=this.startIndex && v<this.endIndex){
					this.NODES_GRADIENTS[v].dx -= this.LEARNING_RATE_STRESS*alpha*(distx/distUV)*(-1);
					this.NODES_GRADIENTS[v].dy -= this.LEARNING_RATE_STRESS*alpha*(disty/distUV)*(-1);
				}
				loss += ((distUV-theoreticD)*(distUV-theoreticD))/(theoreticD*theoreticD);
			}
			
		}
		this.lossFunctionsValues.stress = loss;
	}

	idealEdgeLengthIteration(){
		let sommaQij=0.0
		let nodesGradient = new Array();
		for(let i=0;i<this.nodes.length;i++){
			nodesGradient.push({x:0.0, y:0.0});
		}
		// O(E) iteration
		for(let i=0;i<this.edges.length;i++){ 
			let u = this.edges[i][0];
			let v = this.edges[i][1];
			let l = this.edgesLength[i];
			//let l = this.EDGES_LENGTH;
			let distx = this.nodes[u].x - this.nodes[v].x;
			let disty = this.nodes[u].y - this.nodes[v].y;
			let distUV = Math.sqrt((distx*distx)+(disty*disty));
			let alpha = 2*((distUV-l)/l); 
			nodesGradient[u].x += alpha*(distx/distUV);
			nodesGradient[u].y += alpha*(disty/distUV);
			nodesGradient[v].x -= alpha*(distx/distUV);
			nodesGradient[v].y -= alpha*(disty/distUV);
			sommaQij += ((distUV-l)/l)*((distUV-l)/l);
		}
		let loss=Math.sqrt((1/this.edges.length)*sommaQij);
		if(!Number.isNaN(loss) || loss>10e-5){ 
			for(let i=this.startIndex;i<this.endIndex;i++){ 
				this.NODES_GRADIENTS[i].dx -= this.LEARNING_RATE_EDGELEN*(1/(2*this.edges.length*loss))*nodesGradient[i].x; 
				this.NODES_GRADIENTS[i].dy -= this.LEARNING_RATE_EDGELEN*(1/(2*this.edges.length*loss))*nodesGradient[i].y;
			}
		}
		this.lossFunctionsValues.edgeLen = loss;
	}

	angleResolutionIteration(){
		let loss=0;
		//O(V) Per ogni nodo..
		for(let n=0;n<this.nodes.length;n++){
			let gradoNodo = this.adj[n].length;
			if(gradoNodo>1){
				//O(D) Calcolo angoli nodi attaccati al nodo n
				let angleArray = new Array();
				for(let i=0;i<this.adj[n].length;i++){ 
					let nodeIndex = this.adj[n][i];
					let vectNI = {x: this.nodes[nodeIndex].x - this.nodes[n].x, y: this.nodes[nodeIndex].y - this.nodes[n].y}
					let normNI = Math.sqrt((vectNI.x*vectNI.x) + (vectNI.y*vectNI.y));
					let angleNI = 2*Math.atan(vectNI.y / (normNI+vectNI.x)); 
					if(Number.isNaN(angleNI)) { angleNI = Math.PI;}
					angleArray.push({angle: angleNI, node: nodeIndex, vect: vectNI, norm: normNI});
				}
				//O(d*log(d)) Ordino i nodi attaccati ad n sulla base degli angoli
				angleArray.sort((a,b) => (a.angle > b.angle) ? 1 : ((a.angle < b.angle) ? -1 : 0)); //console.log(angleArray)
				//O(D) Per ogni coppia di nodi consecutiva i, i+1, calcolo gradiente rispetto alla differenza angle[i-1]-angle[i]
				for(let k=0;k<angleArray.length;k++){
					let h=0;
					(k==angleArray.length-1) ? h=0 : h=k+1;
					let angleI = angleArray[k].angle; let vectSI = angleArray[k].vect; 
					let normSI = angleArray[k].norm;  let i = angleArray[k].node;
					let angleJ = angleArray[h].angle; let vectSJ = angleArray[h].vect; 
					let normSJ = angleArray[h].norm;  let j = angleArray[h].node;
					let deltaAngle = 0;
					(k==angleArray.length-1) ? deltaAngle=(2*Math.PI)+angleJ-angleI : deltaAngle=angleJ-angleI; //BIG - LITTLE
					loss += Math.exp((-this.ANR_SENSITIVITY)*deltaAngle)
					let derivExp = Math.exp((-this.ANR_SENSITIVITY)*deltaAngle)
					let atanArgSI = (vectSI.y / (normSI+vectSI.x)); 
					let atanArgSJ = (vectSJ.y / (normSJ+vectSJ.x)); 
					let derivAtanSI = (1/(1+atanArgSI*atanArgSI));
					let derivAtanSJ = (1/(1+atanArgSJ*atanArgSJ));
					//let gradXJ = ( ((-1)*vectSJ.y*(1+((vectSJ.x)/normSJ)))/(((normSJ+vectSJ.x)*(normSJ+vectSJ.x) + (vectSJ.y*vectSJ.y) ) ));
					//let part1 = (2/(1+Math.pow((vectSJ.y/(normSJ+vectSJ.x)),2)));
					//let gradYJ = part1*((1/(normSJ+vectSJ.x)) - ((vectSJ.y*vectSJ.y)/(normSJ*(normSJ+vectSJ.x)*(normSJ+vectSJ.x)) ));
					//let gradXI = ( ((-1)*vectSI.y*(1+((vectSI.x)/normSI)))/(((normSI+vectSI.x)*(normSI+vectSI.x)) + (vectSI.y*vectSI.y)));
					//let part2 = (2/(1+Math.pow((vectSI.y/(normSI+vectSI.x)),2)));
					//let gradYI = part2*((1/(normSI+vectSI.x)) - ((vectSI.y*vectSI.y)/(normSI*(normSI+vectSI.x)*(normSI+vectSI.x)) ));
					let gradXJ = ( ((-1)*vectSJ.y*(1+((vectSJ.x)/normSJ)))/((normSJ+vectSJ.x)*(normSJ+vectSJ.x)) );
					let gradYJ = ((1/(normSJ+vectSJ.x)) - ((vectSJ.y*vectSJ.y)/(normSJ*(normSJ+vectSJ.x)*(normSJ+vectSJ.x)) ));
					let gradXI = ( ((-1)*vectSI.y*(1+((vectSI.x)/normSI)))/((normSI+vectSI.x)*(normSI+vectSI.x)) );
					let gradYI = ((1/(normSI+vectSI.x)) - ((vectSI.y*vectSI.y)/(normSI*(normSI+vectSI.x)*(normSI+vectSI.x)) ));
					if(j>=this.startIndex && j<this.endIndex){ 
						this.NODES_GRADIENTS[j].dx += this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSJ)*gradXJ;
						this.NODES_GRADIENTS[j].dy += this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSJ)*gradYJ;
					}
					if(i>=this.startIndex && i<this.endIndex){
						this.NODES_GRADIENTS[i].dx -= this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSI)*gradXI; 
						this.NODES_GRADIENTS[i].dy -= this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSI)*gradYI;
					}
					if(n>=this.startIndex && n<this.endIndex){
						this.NODES_GRADIENTS[n].dx -= this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSJ)*gradXJ;
						this.NODES_GRADIENTS[n].dy -= this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSJ)*gradYJ;
						this.NODES_GRADIENTS[n].dx += this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSI)*gradXI; 
						this.NODES_GRADIENTS[n].dy += this.LEARNING_RATE_ANGLE*(derivExp)*(derivAtanSI)*gradYI;
					}
				}
			}
		}
		this.lossFunctionsValues.angleRes = loss;
	}

	neighborsPreservationIteration(){
		let nodesNumber = this.nodes.length;
		// [1] calcolo matrice distanze
		let euclideanDistMatrix = new Array(nodesNumber);
		for(let i=0;i<nodesNumber;i++){
			euclideanDistMatrix[i] = new Array(nodesNumber);
			for(let j=0;j<nodesNumber;j++){
				if(i==j){ euclideanDistMatrix[i][j]=0; }
				else    { euclideanDistMatrix[i][j] = this.euclideanDistance(this.nodes[i],this.nodes[j]); }
			}
		}
		// [2] Per ogni nodo i computo un raggio di una circonferenza entro cui un nodo j può essere definito neighbor di i
		let nodesRadius = new Array(nodesNumber)
		for(let i=0;i<nodesNumber;i++){
		    let sortedDistRow = [...euclideanDistMatrix[i]].sort((a,b)=>a-b);
		    let k = this.nodesNeighborsNumber[i];
		    if(k==nodesNumber-1){ nodesRadius[i] = (sortedDistRow[k]  + this.NEIGH_HALF_MARGIN)/2; } 
		    else                { nodesRadius[i] = (sortedDistRow[k]+sortedDistRow[k+1])/2; }
		}
		// [3] Calcolo ogni Kij (valore continuo se j appartiene al neighborhood di i (Kij>0) oppure no (Kij<0))
		// Se il Kij è dello stesso segno della adjIJ reale (+1,-1), kij*ADJij sarà > 0. 
		// Quindi se concorde (>0), e se kij*ADJij è più grande del semi-margine rispetto al nodo i, significa che la coppia
		// <i,j> rispetta la matrice di adiacenza (OK!), e che j sta nella circonferenza del nodo i, quindi errore_IJ = 0; (ReLU)
		let error = new Array();
		let predAdjMatrix = new Array();
		let errorIndexPair = new Array();
		for(let i=0;i<nodesNumber;i++){
			error[i] = new Array(nodesNumber);
			predAdjMatrix[i] = new Array(nodesNumber);
			for(let j=0;j<nodesNumber;j++){
				if(i==j){ //UN NODO NON E' MAI VICINO DI SE STESSO!
					error[i][j] = 0; 
					errorIndexPair.push([0,i,j]); 
					predAdjMatrix[i][j] = 0; 
				} else {
					let kij = 0;  //k indica se j è interno alla circonferenza di i (>0) o viceversa
					if(nodesRadius[i]<this.NEIGH_HALF_MARGIN/5){ //DISTANZA MINIMA PER CUI CONSIDERO UN VICINO
						kij = -(euclideanDistMatrix[i][j] - this.NEIGH_HALF_MARGIN/5);
					} else {
						if(nodesRadius[i]>this.NEIGH_HALF_MARGIN){ //PONGO UN FRENO ALL'ESPANSIONE DEL GRAFO
							kij = -(euclideanDistMatrix[i][j] - this.NEIGH_HALF_MARGIN);
						} else {
							kij = -(euclideanDistMatrix[i][j] - nodesRadius[i]);
						}
					}
					//Con kij calcolo l'errore rispetto alla coppia <I,J>
					let e = (this.NEIGH_HALF_MARGIN - (kij*(2*this.realAdjMatrix[i][j]-1)));
					error[i][j] = Math.max(0.0,e);
					errorIndexPair.push([error[i][j],i,j]);
					((kij>=0) ? predAdjMatrix[i][j] = 1 : predAdjMatrix[i][j] = 0)
				}
			}
		}
		// [4] Ordino gli errori in maniera decrescente
		errorIndexPair.sort((a,b)=>b[0]-a[0]);
		// [5] Calcolo il margine del valore del JaccardIndex che si ottiene aggiungendo un nuovo elemento
		let deltaJ = new Array(nodesNumber);
		for(let i=0;i<deltaJ.length;i++){ deltaJ[i] = new Array(nodesNumber); }
		for(let z=0;z<this.cumulativeJaccardIndex.length;z++){
			let i=errorIndexPair[z][1]; let j=errorIndexPair[z][2];
			if(z>0){ deltaJ[i][j] = this.cumulativeJaccardIndex[z]-this.cumulativeJaccardIndex[z-1]; } 
			else   { deltaJ[i][j] = this.cumulativeJaccardIndex[0]; }
		}
		// [6] Calcolo gradienti della Lovasz hinge loss rispetto a ogni elemento ij della matrice di adiacenza:
		let loss = 0; //LOVASZ_HINGE_LOSS = SUM(error_vector * deltaJaccardIndexValue_vector)
		for(let i=0;i<nodesNumber;i++){
			for(let j=0;j<nodesNumber;j++){
				if(error[i][j]>0){
					//GRADIENT OF ERROR BETWEEN (I,J) (without kij gradient)
					//let gradErrorXij = deltaJ[i][j]*(1-2*this.realAdjMatrix[i][j]);
					//let gradErrorYij = deltaJ[i][j]*(1-2*this.realAdjMatrix[i][j]);
					let gradErrorXij = 2*error[i][j]*deltaJ[i][j]*(1-2*this.realAdjMatrix[i][j]);
					let gradErrorYij = 2*error[i][j]*deltaJ[i][j]*(1-2*this.realAdjMatrix[i][j]);
					//GRADIENT OF kij (EUCLIDEAN DIST(nodeI,nodeJ))
					let distxIJ = this.nodes[i].x - this.nodes[j].x;	
					let distyIJ = this.nodes[i].y - this.nodes[j].y;
					let distIJ = euclideanDistMatrix[i][j];
					// I
					if(i>=this.startIndex && i<this.endIndex){
						this.NODES_GRADIENTS[i].dx += this.LEARNING_RATE_NEIGH*gradErrorXij*(distxIJ/distIJ);
						this.NODES_GRADIENTS[i].dy += this.LEARNING_RATE_NEIGH*gradErrorYij*(distyIJ/distIJ);
					}
					// J
					if(j>=this.startIndex && j<this.endIndex){
						this.NODES_GRADIENTS[j].dx += this.LEARNING_RATE_NEIGH*gradErrorXij*(distxIJ/distIJ)*(-1);
						this.NODES_GRADIENTS[j].dy += this.LEARNING_RATE_NEIGH*gradErrorYij*(distyIJ/distIJ)*(-1);
					}
				}
				loss+=deltaJ[i][j]*error[i][j];
			}
		}

		/*if((iteration_count%200==0) && (printLoss==true)){
			//JACCARD_INDEX(A,B) = (|A && B|) / (|A || B|)
			let intersect=0; let union=0;
			for(let i=0;i<nodesNumber;i++){
				for(let j=0;j<nodesNumber;j++){
					if(this.realAdjMatrix[i][j]==1 || predAdjMatrix[i][j]==1){ union+=1;	    }
					if(this.realAdjMatrix[i][j]==1 && predAdjMatrix[i][j]==1){ intersect+=1;	}
				}
			}
			let jaccardIndexValue = undefined;
			if(union!=0){ jaccardIndexValue = (intersect/union); }
			//console.log("Lovàzs_Hinge_Loss => "+loss+"\nJaccard_Index => "+jaccardIndexValue);
		}*/
		this.lossFunctionsValues.neighPres = loss;
	}

	vertexResolutionIteration(){
		var loss=0.0
		for(let u=0;u<this.nodes.length-1;u++){ 
			for(let v=u+1;v<this.nodes.length;v++){
				let distx = this.nodes[u].x - this.nodes[v].x;
				let disty = this.nodes[u].y - this.nodes[v].y;
				let distUV = Math.sqrt((distx*distx)+(disty*disty));
				let ReLuArg = (1 - (distUV/(this.VERTEX_SCALING*this.VERTEX_MAX_DISTANCE)))
				if(ReLuArg>0){
					let alpha = (2*ReLuArg/(this.VERTEX_SCALING*this.VERTEX_MAX_DISTANCE)); 
					if(u>=this.startIndex && u<this.endIndex){
						this.NODES_GRADIENTS[u].dx -= this.LEARNING_RATE_VERTEX*alpha*(distx/distUV)*(-1);
						this.NODES_GRADIENTS[u].dy -= this.LEARNING_RATE_VERTEX*alpha*(disty/distUV)*(-1);
					}
					if(v>=this.startIndex && v<this.endIndex){
						this.NODES_GRADIENTS[v].dx -= this.LEARNING_RATE_VERTEX*alpha*(distx/distUV);
						this.NODES_GRADIENTS[v].dy -= this.LEARNING_RATE_VERTEX*alpha*(disty/distUV);
					}
					loss += (ReLuArg*ReLuArg);
				}
			}
		}
		this.lossFunctionsValues.vertexRes = loss;
	}

	updatePosition(iteration_count){
		for(let i=0;i<this.nodes.length;i++){
			if(!this.nodes[i].isDragged){  //DRAGGED NODE IS INDIPENDENT!!! event.x, event.y!!!!
				//SIMPLE GRADIENT DESCENT
				if(((!this.MOMENTUM_ENABLE) && (!this.ADAM_ENABLE)) || ((!this.MOMENTUM_ENABLE) && (this.ADAM_ENABLE))){
					this.nodes[i].x += this.LEARNING_RATE*(this.NODES_GRADIENTS[i].dx);
					this.nodes[i].y += this.LEARNING_RATE*(this.NODES_GRADIENTS[i].dy);
				} 
				//GRADIENT DESCENT WITH MOMENTUM
				if((this.MOMENTUM_ENABLE) && (!this.ADAM_ENABLE)){
					this.velocity[i].dx = (this.MOMENTUM*this.velocity[i].dx)+((1-this.MOMENTUM)*(this.NODES_GRADIENTS[i].dx));
					this.velocity[i].dy = (this.MOMENTUM*this.velocity[i].dy)+((1-this.MOMENTUM)*(this.NODES_GRADIENTS[i].dy));
					this.nodes[i].x += this.LEARNING_RATE*(this.velocity[i].dx/(1-Math.pow(this.MOMENTUM,iteration_count+1)));
					this.nodes[i].y += this.LEARNING_RATE*(this.velocity[i].dy(1-Math.pow(this.MOMENTUM,iteration_count+1)));
				}
				//GRADIENT DESCENT WITH ADAM (MOMENTUM + ADAPTIVE-LEARNING-RATES)
				if((this.MOMENTUM_ENABLE) && (this.ADAM_ENABLE)){
					this.velocity[i].dx = (this.MOMENTUM*this.velocity[i].dx)+((1-this.MOMENTUM)*(this.NODES_GRADIENTS[i].dx));
					this.velocity[i].dy = (this.MOMENTUM*this.velocity[i].dy)+((1-this.MOMENTUM)*(this.NODES_GRADIENTS[i].dy));
					this.velocity2[i].dx=(this.MOMENTUM2*this.velocity2[i].dx)+((1-this.MOMENTUM2)*(Math.pow(this.NODES_GRADIENTS[i].dx,2)));
					this.velocity2[i].dy=(this.MOMENTUM2*this.velocity2[i].dy)+((1-this.MOMENTUM2)*(Math.pow(this.NODES_GRADIENTS[i].dy,2)));
					let biasedVx  = this.velocity[i].dx/(1-Math.pow(this.MOMENTUM,iteration_count+1))
					let biasedVy  = this.velocity[i].dy/(1-Math.pow(this.MOMENTUM,iteration_count+1))
					let biasedV2x = this.velocity2[i].dx/(1-Math.pow(this.MOMENTUM2,iteration_count+1))
					let biasedV2y = this.velocity2[i].dy/(1-Math.pow(this.MOMENTUM2,iteration_count+1))
					this.nodes[i].x += this.LEARNING_RATE*(biasedVx/(Math.sqrt(biasedV2x)+10e-5));
					this.nodes[i].y += this.LEARNING_RATE*(biasedVy/(Math.sqrt(biasedV2y)+10e-5));
				}
			}
		}
		if(this.EXP_DECAY_ENABLE){
			this.LEARNING_RATE = this.INIT_LR * Math.exp(-1*this.DECAY*iteration_count);
		}
	}

	updateNodesPosition(iteration_count){
		this.initGradientToZero();
		for(let i=0;i<this.visualCriteria.length;i++){
			let name = this.visualCriteria[i].name;
			if(name == "STRESS"){
				this.stressIteration();
			}
			if(name == "IDEAL_EDGE_LENGTH"){
				this.idealEdgeLengthIteration();
			}
			if(name == "ANGLE_RESOLUTION"){
				this.angleResolutionIteration();
			}
			if(name == "NEIGHBORS_PRESERVATION"){
				this.neighborsPreservationIteration();
			}
			if(name == "CROSSING_NUMBER"){
				this.crossingNumberIteration();
			}
			if(name == "VERTEX_RESOLUTION"){
				this.vertexResolutionIteration();
			}
		}
		this.updatePosition(iteration_count);
	}

	printLosses(iteration_count){
		let name =" ".repeat(4);
		let value=" ".repeat(4);
		for (var key in this.lossFunctionsValues) {
		  	if(key=="stress")     { name+="|   STRESS   "; }
		  	if(key=="edgeLen")    { name+="|  EDGE_LEN  "; }
		  	if(key=="angleRes")   { name+="| ANGLE_RES. "; }
		  	if(key=="neighPres")  { name+="| NEIGH_PRES "; }
		  	if(key=="vertexRes")  { name+="| VERTEX_RES "; }
		  	let lossValue = this.lossFunctionsValues[key].toString().substring(0,10);
		  	if(lossValue.length<10){ lossValue+=" ".repeat(10-lossValue.length); }
		  	value+="| "+lossValue+" ";
		}
		name +="|"+" ".repeat(4);	value+="|"+" ".repeat(4);
		let line  = " ".repeat(4)+"+"; line+="-".repeat(name.length-10); line+="+"+" ".repeat(4);
		let title = " ".repeat(4)+"LOSS VALUES AT ITERATION: "+iteration_count+" ".repeat(4);;
		let lr    = " ".repeat(4)+ "Current Learning Rate : "+this.LEARNING_RATE.toString().substring(0,8)+" ".repeat(4);
		let module = " ".repeat(4)+"Current Gradient L2   : "+this.gradientModule().toString().substring(0,8)+" ".repeat(4);
		let max = 0;
		if(this.EXP_DECAY_ENABLE){
			max = Math.max(name.length,line.length,title.length,module.length,lr.length);
		} else {
			max = Math.max(name.length,line.length,title.length,module.length,);
		}
		line +=" ".repeat(max-line.length)
		name +=" ".repeat(max-name.length)
		value+=" ".repeat(max-value.length)
		title+=" ".repeat(max-title.length)
		lr+=" ".repeat(max-lr.length)
		module+=" ".repeat(max-module.length)
		let padding = " ".repeat(max);

		let outputString = "";
		if(this.EXP_DECAY_ENABLE){
			outputString = "%c"+padding+"\n"+title+"\n"+line+"\n"+name+"\n"+line+"\n"+value+"\n"+line+"\n"+lr+"\n"+module+"\n"+padding;
		} else {
			outputString = "%c"+padding+"\n"+title+"\n"+line+"\n"+name+"\n"+line+"\n"+value+"\n"+line+"\n"+module+"\n"+padding;
		}
		console.log(outputString,'background: #000; color: #badd55')
	}

	printStop(iteration_count){
		let str = "    |    STOP AT ITERATION : "+iteration_count+"    |    ";
		let line = "    +"+"-".repeat(str.length-10)+"+    ";
		console.log("%c"+line+"\n"+str+"\n"+line,'background: #000; color: orange')
	}

	gradientModule(){
		let dot=0;
		for(let i=0;i<this.NODES_GRADIENTS.length;i++){
			dot+=(this.NODES_GRADIENTS[i].dx * this.NODES_GRADIENTS[i].dx);
			dot+=(this.NODES_GRADIENTS[i].dy * this.NODES_GRADIENTS[i].dy);
		}
		return Math.sqrt(dot);
	}
}