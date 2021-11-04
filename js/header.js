//PARAMETRI ANIMAZIONE
const FRAME_MSEC = 10;            //Durata di ogni frame di aggiornamento posizione dei nodi e archi (in millisecondi)
const MAX_ITERATION_NUM = 2000;   //Massimo numero iterazioni dell'algoritmo

//PROPRIETA' GRAFICHE
let GRAPHIC_PROPERTIES = {
		PIXEL_UNIT         : 90,                    //Quanti pixel corrispondono a una distanza unitaria tra nodi (dist 1 = 90 pixel)
		GRAPH_NODES_COLOR  : "rgb(  0, 139, 139)",  //Colore interno di tutti i nodi
		NODES_BORDER_COLOR : "rgb(154, 205,  50)",  //Colore del bordo di tutti i nodi
	    NODES_TEXT_COLOR   : "rgb(154, 205,  50)",  //Colore delle etichette testuali dei nodi
	    NODES_RADIUS       : 12,                    //Raggio (in pixel) di tutti i nodi
		EDGES_COLOR        : "rgb(154, 205,  50)",  //Colore di tutti gli archi
		EDGES_WIDTH        : 2                      //Spessore archi
}

//IPERPARAMETRI DELL'OTTIMIZZATORE
let HYPER_PARAMETERS = { 
		EDGES_LENGTH           : 1.0,
	   	ANR_SENSITIVITY        : 1,
	   	NEIGH_HALF_MARGIN      : 1.0,
	   	CROSSING_NUMBER_MARGIN : 0.05,
	   	VERTEX_SCALING         : null, // Default 1/sqrt(|V|)
	   	VERTEX_MAX_DISTANCE    : 5.0,
	  	MAX_GRAVITY_RADIUS     : 30,

	   	MOMENTUM               : {enable: true, value: 0.9 },
	   	MOMENTUM2              : {enable: true, value: 0.99 },
	   	START_LEARNING_RATE    : 0.05,
	   	EXP_DECAY_LEARING_RATE : {enable: true, value: 0.0001 },
	   	GRADIENT_MIN_MODULE    : 0.001
};

//CRITERI ESTETICI DA MIGLIORARE E COEFFICIENTI DI MISCELAZIONE 
let VISUAL_CRITERIA=[];
VISUAL_CRITERIA.push({name:"STRESS"                , learningRate:0.005})
VISUAL_CRITERIA.push({name:"IDEAL_EDGE_LENGTH"     , learningRate:0.0})
VISUAL_CRITERIA.push({name:"NEIGHBORS_PRESERVATION", learningRate:0.0})
VISUAL_CRITERIA.push({name:"ANGLE_RESOLUTION"      , learningRate:0.0})
VISUAL_CRITERIA.push({name:"VERTEX_RESOLUTION"     , learningRate:0.0})