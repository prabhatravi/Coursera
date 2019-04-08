// Hex game
// 
// Compile with: g++ 4-Homework.cpp -o 4-Homework -lboost_regex -std=c++11
//
// Brief explanation:
// The board is represented by a graph.
// Each board position corresponds to one node identified by a number (from 1 to dimension^2) and marked with one color (WHITE=available, BLUE=player1 and RED=player2).
// There is an edge connecting two nodes if those nodes correspond to neighbor positions in the board. 
// There are 4 extra virtual nodes (WEST, EAST, NORTH and SOUTH) connected to nodes of the start and end posiitons in the board. WEST and EAST start BLUE. NORTH and SOUTH start RED.
// Each node is initialized in its own set.
// At each move, the selected position is marked with player's color and the node set is combined with sets of all neighbors with the same color. If start virtual node and end virtual node are in the same set, so there is a path connecting them.
//
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <ctime>
#include <cstdlib>
#include <limits>		// Max integer definition library
#include <string>		// String manipulation library
#include <fstream>		// File manipulation library
#include <boost/regex.hpp> 	// Boost regular expression library
using namespace std;

//==============================================================================
// General definitions
//==============================================================================

// INFINIT is used to represent no edge/path between two nodes 
const int INFINIT = numeric_limits<int>::max();

// Node color definition 
// WHITE means available positions
// BLACK means invalid color
// BLUE and RED represent player's colors
enum class nodeColor {BLUE, RED, WHITE, BLACK};

// Overload operator << to print a nodeColor
ostream &operator<<(ostream &output, nodeColor c) {
  switch (c) {
    case nodeColor::BLUE:	output << "X"; break;
    case nodeColor::RED:	output << "O"; break;
    case nodeColor::WHITE:	output << "."; break;
    default: 			output << "~"; break;
  }  
  return output;
}

// Game player definition 
// HUMAN means user
// COMPUTER means random choices
enum class gamePlayer {HUMAN, COMPUTER};

// Overload operator << to print gamePlayer
ostream &operator<<(ostream &output, gamePlayer p) {
  switch (p) {
    case gamePlayer::HUMAN:	output << "Human   "; break;
    case gamePlayer::COMPUTER:	output << "Computer"; break;
    default:			output << "Unknown ";
  }
  return output;
}

//==============================================================================
// Class Node
// Used to store information about nodes/edges in the adjacency list of a graph
// Adjacency lists is a list of Nodes (identified by numbers)
// Each node contains a list of neighbors and edge weights
//==============================================================================
class Node {
  public:
    Node():number(-1),weight(INFINIT),color(nodeColor::WHITE) { edges.clear(); };
    Node(int n, int w=INFINIT, nodeColor c=nodeColor::WHITE):number(n),weight(w),color(c) { edges.clear(); };
    ~Node() {};
    Node& operator=(const Node &n);
    bool operator==(const Node &n);
    void setNodeNumber(int n) { number = n; };
    int getNodeNumber() { return number; };
    void setEdgeWeight(int w=INFINIT) { weight = w; };
    int getEdgeWeight() { return weight; };
    void setNodeColor(nodeColor c) { color = c; };
    nodeColor getNodeColor() { return color; };
    void setEdge(Node N) { edges.push_back(N); };
    list<Node> getEdges() { return edges; };
    
  private:
    int number;
    int weight;
    nodeColor color;
    list<Node> edges;	
};

// Assignment operator definition for Node
Node& Node::operator=(const Node &n) {
  this->number = n.number;
  this->weight = n.weight;
  this->color = n.color;
  this->edges = n.edges;
  return *this;
}

// Return true if two nodes have the same node number and false otherwise
bool Node::operator==(const Node &n) {
  if (this->number == n.number)
    return true;
  return false;
}

// Overload operator << to print Node
ostream &operator<<(ostream &output, Node n) {
  if (n.getEdgeWeight() == INFINIT)
    output << "(" << n.getNodeNumber() << ")" ;
  else
    output << "(" << n.getNodeNumber() << ":" << n.getEdgeWeight() << ")";
  return output;
}

// Overload operator << to print list<Node>
ostream &operator<<(ostream &output, list<Node> L) {
  for(list<Node>::iterator i=L.begin(); i != L.end(); ++i)
    output << (*i) << " ";
  return output;
}

//==============================================================================
// Class Edge
// Used to store information about edges (source, destination and weight)
//==============================================================================
class Edge {
  public:
    Edge():src(-1),dst(-1),weight(INFINIT) {};
    Edge(int s, int d, int w):src(s),dst(d),weight(w) {};
    ~Edge() {};
    Edge& operator=(const Edge &e);
    bool operator==(const Edge &e);
    void set(int s, int d, int w) { src=s; dst=d; weight=w; };
    void setSrc(int s) { src = s; }; 
    void setDst(int d) { dst = d; }; 
    void setWeight(int w) { weight = w; }; 
    int getSrc() { return src; }; 
    int getDst() { return dst; }; 
    int getWeight() { return weight; }; 
    
  private:
    int src;	// Source node
    int dst;	// Destination node
    int weight;	// Weight of the edge between src and dst
};

// Assignment operator definition for Edge
Edge& Edge::operator=(const Edge &e) {
  this->src = e.src;
  this->dst = e.dst;
  this->weight = e.weight;
  return *this;
}

// Return true if two edges have the same src and dst and false otherwise
bool Edge::operator==(const Edge &e) {
  if ((this->src == e.src) && (this->dst == e.dst))
    return true;
  return false;
}

// Return true if 'e1' weight is less then 'e2' weight and false otherwise
bool operator<(Edge& e1, Edge& e2) {
  if (e1.getWeight() < e2.getWeight())
    return true;
  return false;
}

// Overload operator << to print Edge
ostream &operator<<(ostream &output, Edge e) {
  output << "(" << e.getSrc() << "," << e.getDst() << "," << e.getWeight() << ")";
  return output;
}

// Overload operator << to print list<Edge>
ostream &operator<<(ostream &output, list<Edge> L) {
  for(list<Edge>::iterator i=L.begin(); i != L.end(); ++i)
    output << (*i) << " ";
  return output;
}

//==============================================================================
// Graph Class
// Represent a Graph through an adjacency list of nodes
//==============================================================================
class Graph {
  public:
    Graph():numV(0),numE(0) { adjList.clear(); };	// Create an empty graph
    ~Graph() {};
    void init(int numVertices);
    int get_edge_value(int x, int y);
    void set_edge_value(int x, int y, int value);
    nodeColor get_node_color(int x);
    void set_node_color(int x, nodeColor c);
    bool adjacent(int x, int y);
    list<Node> neighbors(int x);
    int V();
    int E();
    list<Node> vertices();
  
  private:
    int numV;			// Number of nodes of the Graph
    int numE;			// Number of edges of the Graph
    list<Node> adjList;		// Adjacency list of nodes representing the Graph
};

// Initialize an empty graph with numVertices nodes
// Create an adjacency list with all nodes and no edges
void Graph::init(int numVertices) {
  numV = numVertices;
  numE = 0;
  adjList.clear();
  for(int i=1; i<=numVertices; ++i) {
    Node newNode(i);
    adjList.push_back(newNode);
  }
}

// Return edge weight between two nodes
// Return INFINIT if edge doesn't exist
int Graph::get_edge_value(int x, int y) {
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i)
    if ((*i).getNodeNumber() == x) {
      list<Node> iEdges = (*i).getEdges();
      for(list<Node>::iterator j=iEdges.begin(); j != iEdges.end(); ++j)
	 if ((*j).getNodeNumber() == y)
	   return (*j).getEdgeWeight(); 
    }
  return INFINIT;
}

// Set edge weight between two nodes
void Graph::set_edge_value(int x, int y, int value) {
  bool found;
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i) {
    // Add 'y' in the list of 'x' neighbors (if doesn't exist)
    // Set edge weight to value
    if ((*i).getNodeNumber() == x) {
      found = false;
      list<Node> iEdges = (*i).getEdges();
      for(list<Node>::iterator j=iEdges.begin(); j != iEdges.end(); ++j)
	 if ((*j).getNodeNumber() == y) {  
	   (*j).setEdgeWeight(value);
	   found = true;
	 }
      if (!found) {
	Node newNodeY(y,value);	
	(*i).setEdge(newNodeY);
      }
    }  
    // Add 'x' in the list of 'y' neighbors (if doesn't exist)
    // Set edge weight to value
    if ((*i).getNodeNumber() == y) {
      found = false;
      list<Node> iEdges = (*i).getEdges();
      for(list<Node>::iterator j=iEdges.begin(); j != iEdges.end(); ++j)
	 if ((*j).getNodeNumber() == x) {  
	   (*j).setEdgeWeight(value);
	   found = true;
	 }
      if (!found) {
	Node newNodeX(x,value);	
	(*i).setEdge(newNodeX);
      }
    }
  }
  if (!found)
    ++numE;	  	// Increment the number of edges in the graph
}

// Return node color
// Return BLACK if node doesn't exist
nodeColor Graph::get_node_color(int x) {
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i)
    if ((*i).getNodeNumber() == x)
      return (*i).getNodeColor(); 
  return nodeColor::BLACK;
}

// Set node color
void Graph::set_node_color(int x, nodeColor c) {
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i)
    if ((*i).getNodeNumber() == x)
      (*i).setNodeColor(c); 
}

// Return true if two nodes are neighbors and false otherwise
bool Graph::adjacent(int x, int y) {
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i)
    if ((*i).getNodeNumber() == x) {
      list<Node> iEdges = (*i).getEdges();
      for(list<Node>::iterator j=iEdges.begin(); j != iEdges.end(); ++j)
	if ((*j).getNodeNumber() == y)
	  return true;
    }
  return false;
}

// Return a list<Node> containing the list of neighbors of 'x'
list<Node> Graph::neighbors(int x) {
  list<Node> xEdges;
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i)
    if ((*i).getNodeNumber() == x)
      xEdges = (*i).getEdges();
  return xEdges;
}
   
// Return the number of nodes in the Graph
int Graph::V() {
  return numV;
}

// Return the number of edges in the Graph
int Graph::E() {
  return numE;
}

// Return a list<Node> containing all nodes in the Graph
list<Node> Graph::vertices() {
  list<Node> nodesList;
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i) {
    Node newNode((*i).getNodeNumber(),(*i).getEdgeWeight());
    nodesList.push_back(newNode);
  }
  return nodesList;
}

// Overload operator << to print a graph
ostream &operator<<(ostream &output, Graph g) {
  output << endl;
  output << "   ";
  list<Node> gNodes = g.vertices();
  for(list<Node>::iterator i=gNodes.begin(); i != gNodes.end(); ++i)
    if ((*i).getNodeNumber() < 10)
      output << "  " << (*i).getNodeNumber();
    else
      output << " " << (*i).getNodeNumber();
  output << endl;
  for(list<Node>::iterator i=gNodes.begin(); i != gNodes.end(); ++i) {
    if ((*i).getNodeNumber() < 10)
      output << "  " << (*i).getNodeNumber();
    else
      output << " " << (*i).getNodeNumber();
    int shift=0;
    list<Node> iEdges = g.neighbors((*i).getNodeNumber());
    for(list<Node>::iterator j=iEdges.begin(); j != iEdges.end(); ++j) {
      int walk=(*j).getNodeNumber()-shift;
      for(int k=0; k<walk; ++k) {
	output << "  -";
	shift++;
      }
      if ((*j).getEdgeWeight() < 10)
	output << "  " << (*j).getEdgeWeight();
      else
	output << " " << (*j).getEdgeWeight();
      shift++;
    }
    while (shift<g.V()) {
      output << "  -";
      shift++;
    }
    output << endl;
  }
  return output;
}

//==============================================================================
// Coordinate Class
// Used to represent a coordinate (i,j) 
//==============================================================================
class Coord {
  public:
    Coord():i(0),j(0) {};
    Coord(int i, int j):i(i),j(j) {};
    Coord operator=(Coord c) { i=c.i; j=c.j; return *this; };
    friend ostream &operator<<(ostream &output, Coord c);
    int geti() { return i; };
    int getj() { return j; };
    
  private:
    int i, j;
};

// Overload operator << to print a Coord element
ostream &operator<<(ostream &output, Coord c) {
  output << "(" << c.i << "," << c.j << ")";
  return output;
}

// Overload operator << to print a Coord []
ostream &operator<<(ostream &output, vector<Coord> vec) {
  for (auto i:vec) 
    output << i << " ";
  return output;
}

//==============================================================================
// Board Class
// Used to generate the board of the hex game 
//==============================================================================
class Board {
  public:
    Board() {};
    Board(int d) { init(d); };	// Constructor calls init method to initialize the board with dimension (d x d)
    ~Board() {};
    void init(int d);
    friend ostream &operator<<(ostream &output, Board b);
    int getDimension() { return dimension; };	// Return the dimension of the board
    int coordToPos(Coord c) { return (c.geti()-1)*dimension+c.getj(); };	// Transform coordinates into positions in the board
    nodeColor getPosColor(Coord c);
    void setPosColor(Coord coordinate, nodeColor color);
    bool isCoordValid(Coord c);
    bool isCoordOccupied(Coord c);
    bool hasConnectedPath(int player);
    vector<Coord> availableCoord();
    
  private:
    int dimension;	// Dimension of the board (dimension x dimension)
    Graph g;		// Graph that represents the board
    int startNode[2];	// Store start node for player1 (0) and player2 (1)
    int endNode[2];	// Store end node for player1 (0) and player2 (1)
    map<int, int> set;	// Map used to store nodes set
    void makeSet(int u);
    int findSet(int u);
    void combineSet(int u, int v);
    void combineNeighborSet(int u);
    void initSets();
    bool sameSet(int u, int v);
};

// Create a graph with d^2 nodes and 4 virtual nodes
// Each node of the graph represents one position of the board
// Each edge of the graph connects adjacent positions in the board
// Each virtual node represents a special node to represent WEST, EAST, NORTH and SOUTH 
void Board::init(int d) {
  dimension = d;
  const int edgeweight = 1;
  g.init(d*d+4);	// Initialize the graph with dimension^2 nodes + 4 virtual nodes
 
  for (int i=1; i<=dimension; ++i)
    for (int j=1; j<=dimension; ++j) {
      // neighbors is a vector containing coordinates of all position adjacent to board position (i,j)
      vector<Coord> neighbors = {{i-1,j},{i-1,j+1},{i,j-1},{i+1, j-1}, {i+1, j},{i,j+1}};
      // For each valid neighbor position, create an edge between node (i,j) and the node corresponding to this position
      for (auto k:neighbors)
	if (isCoordValid(k))
	  g.set_edge_value(coordToPos({i,j}),coordToPos(k),edgeweight);
    }
  int vnodeWest, vnodeEast, vnodeNorth, vnodeSouth;
  vnodeWest = coordToPos({dimension+1,1});	// Node (dimension+1,1) represents WEST virtual node
  g.set_node_color(vnodeWest,nodeColor::BLUE);	// Assign BLUE color to WEST virtual node
  for (int i=1; i<=dimension; ++i)		// Connect all nodes in the column 1 to WEST virtual node
    g.set_edge_value(vnodeWest,coordToPos({i,1}),edgeweight);
  vnodeEast = coordToPos({dimension+1,2});	// Node (dimension+1,2) represents EAST virtual node
  g.set_node_color(vnodeEast,nodeColor::BLUE);	// Assign BLUE color to EAST virtual node
  for (int i=1; i<=dimension; ++i)		// Connect all nodes in the column 'dimension' to EAST virtual node
    g.set_edge_value(vnodeEast,coordToPos({i,dimension}),edgeweight);
  vnodeNorth = coordToPos({dimension+1,3});	// Node (dimension+1,3) represents NORTH virtual node
  g.set_node_color(vnodeNorth,nodeColor::RED);	// Assign RED color to NORTH virtual node
  for (int j=1; j<=dimension; ++j)		// Connect all nodes in the row 1 to NORTH virtual node
    g.set_edge_value(vnodeNorth,coordToPos({1,j}),edgeweight);
  vnodeSouth = coordToPos({dimension+1,4});	// Node (dimension+1,4) represents SOUTH virtual node
  g.set_node_color(vnodeSouth,nodeColor::RED);	// Assign RED color to SOUTH virtual node
  for (int j=1; j<=dimension; ++j)		// Connect all nodes in the row 'dimension' to SOUTH virtual node
    g.set_edge_value(vnodeSouth,coordToPos({dimension,j}),edgeweight);
  startNode[0]=vnodeWest;	// Assign WEST to start node of player1
  endNode[0]=vnodeEast;		// Assign EAST to end node of player1
  startNode[1]=vnodeNorth;	// Assign NORTH to start node of player2
  endNode[1]=vnodeSouth;	// Assign SOUTH to end node of player2  
  initSets();		// Assign each node of the graph to its own set
}

// Assign node 'u' to a new set 
void Board::makeSet(int u) {
  static int setNumber = 0;
  set[u] = setNumber++;
}

// Return the set of node 'u' 
int Board::findSet(int u) {
  return set[u]; 
}

// Move all nodes in the same set of node 'v' to the set of node 'u' 
void Board::combineSet(int u, int v) {
  int uSet = set[u], vSet = set[v];
  for (map<int,int>::iterator it=set.begin(); it != set.end(); ++it)
    if (it->second == vSet)
      set[it->first] = uSet;
}

// Add each node of the graph to its own set
void Board::initSets() {
  list<Node> allNodes = g.vertices();  
  for(auto i:allNodes)
    makeSet(i.getNodeNumber());
}

// Move n and all neighbors of 'n' with n's color to the same set 
void Board::combineNeighborSet(int u) {
  for(auto i:g.neighbors(u))
    if (g.get_node_color(u)==g.get_node_color(i.getNodeNumber()))
      combineSet(u,i.getNodeNumber());
}

// Return true if two nodes are in the same set and false otherwise 
bool Board::sameSet(int u, int v) {
  if (findSet(u) == findSet(v))
    return true;
  return false;
}

// Return true if there is a path from startNode to endNode of that player 
bool Board::hasConnectedPath(int player) {
  if (sameSet(startNode[player],endNode[player]))
    return true;
  return false;
}

// Return the color of a board position 
nodeColor Board::getPosColor(Coord c) {
  return g.get_node_color(coordToPos(c));
}

// Set the color of a board position 
void Board::setPosColor(Coord coordinate, nodeColor color) {
  g.set_node_color(coordToPos(coordinate), color);
  combineNeighborSet(coordToPos(coordinate));
}

// Return true if 'c' is a valid coordinate in the board an false otherwise 
bool Board::isCoordValid(Coord c) {
  if ((c.geti()<1) || (c.geti()>getDimension()) || (c.getj()<1) || (c.getj()>getDimension()))
      return false;
  return true;
}

// Return true if 'c' is an occupied position in the board and false otherwise 
bool Board::isCoordOccupied(Coord c) {
  if ((getPosColor(c)==nodeColor::RED) || ((getPosColor(c))==nodeColor::BLUE))
      return true;
  return false;
}

// Return a vector with available (free) positions in the board 
vector<Coord> Board::availableCoord() {
  vector<Coord> available;
  for (int i=1; i<=dimension; ++i)
    for (int j=1; j<=dimension; ++j)
      if (getPosColor({i,j})==nodeColor::WHITE)
	available.push_back({i,j});
  return available;
}

// Overload operator << to print a board
ostream &operator<<(ostream &output, Board b) {
  // Show NORTH label
  for (int j=1; j<=b.dimension/2+1; ++j)
    output << "    ";
  output << "  NORTH" << endl << endl;
  // Show column names (letters)
  output << "         ";
  for (int j=1; j<=b.dimension; ++j)
    output << static_cast<char>(j+'A'-1) << "   ";
  output << endl << endl;
  // Show line numbers, board positions and line numbers again
  for (int i=1; i<=b.dimension; ++i) {
    if (i==((b.dimension/2)+1)) {
      for (int k=1; k<=(i-1)*2; ++k)
	output << " ";
      output << "WEST  ";
    }
    else
      for (int k=1; k<=(i-1)*2+6; ++k)
	output << " ";
    if (i<10)
      output << i << "   ";
    else
      output << i << "  ";      
    for (int j=1; j<=b.dimension; ++j)
      output << b.g.get_node_color(b.coordToPos({i,j})) << "   ";
    output << i ;
    if (i==((b.dimension/2)+1))
      output << "  EAST";
    output << endl << endl;
  }
  output << "     ";
  // Show column names (letters) again
  for (int j=1; j<=(b.dimension)*2-1+6; ++j)
    output << " ";
  for (int j=1; j<=b.dimension; ++j)
    output << static_cast<char>(j+'A'-1) << "   ";
  output << endl;
  // Show "SOUTH" label
  output << endl;
  for (int j=1; j<=b.dimension+1; ++j)
    output << "    ";
  output << "  SOUTH" << endl;
  return output;
}

//==============================================================================
// Game Class
// Used to control the game flow 
//==============================================================================
class Game {
  public:
    Game() { player[0]=player[1]=gamePlayer::COMPUTER; round=0; srand(time(NULL)); };
    ~Game() {};
    void init();
    
  private:
    Board b;
    gamePlayer player[2];
    int round;
    inline nodeColor currentColor() { return static_cast<nodeColor>(round%2); };	// Return the color of current player
    inline int currentMove() { return (round/2)+1; };	// Return the current move number
    inline int currentPlayer() { return (round%2); };	// Return the current player
    int playRound();
    Coord readLocation();
    Coord randomLocation();
};

// Start and control the game flow
// Init the board and select which player will start the game according to user input
void Game::init() {
  int dimension;
  char gofirst;
  cout << "Application start..." << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  do {	// Get board dimension from user and validate it
    cout << "Enter number of hexes on 1 side (min: 1 - max: 22): ";
    cin >> dimension;
    if ((dimension < 1) || (dimension > 22))
      cout << "ERROR - Invalid dimension!" << endl;
  } while ((dimension < 1) || (dimension > 22));
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Playing with " << dimension << "x" << dimension << " board." << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  do {	// Get user's choice about starting the game or not and validate it
    cout << "Human, do you want to go first (Y or N)? ";
    cin >> gofirst;
    gofirst=toupper(gofirst);
    if ((gofirst!='Y') && (gofirst!='N'))
      cout << "ERROR - Invalid answer!" << endl;
  } while ((gofirst!='Y') && (gofirst!='N'));
  if (gofirst=='Y')
    player[0] = gamePlayer::HUMAN;    
  else
    player[1] = gamePlayer::HUMAN;     
  cout << "-----------------------------------------------------------------------" << endl;
  // Show information about the game definitions
  cout << "************************************************************" << endl;
  cout << "* KEY                                                      *" << endl;
  cout << "************************************************************" << endl;
  cout << "*      ITEM      * SYMBOL  *             NOTES             *" << endl;
  cout << "************************************************************" << endl;
  cout << "* Empty Location *    .    *                               *" << endl;
  cout << "* " << player[0] << "       *    " << static_cast<nodeColor>(0) << "    * connects West-East            *" << endl;
  cout << "* " << player[1] << "       *    " << static_cast<nodeColor>(1) << "    * connects North-South          *" << endl;
  cout << "************************************************************" << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  // Create the board with users's specified dimension
  b.init(dimension);
  // Play round by round until we have a winner
  int winner =-1;
  while (winner==-1) {
    winner = playRound();
    round++;
  }
  // Show winner information
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Winner: " << player[winner] << "(" << static_cast<nodeColor>(winner) << ")" << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Board: " << endl;
  cout << b << endl;
  if (player[winner] == gamePlayer::HUMAN) {
    cout << "-----------------------------------------------------------------------" << endl;
    cout << "Argggh! You have bested me, Human!" << endl;
    cout << "-----------------------------------------------------------------------" << endl;
  }
  else {
    cout << "-----------------------------------------------------------------------" << endl;
    cout << "Loooser! Not good enough, Human!" << endl;
    cout << "-----------------------------------------------------------------------" << endl;
  }
}

// Control a round play
int Game::playRound() {
  Coord currentPos;
  cout << "Board: " << endl;
  cout << b << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Move " << currentMove() << ": " << player[currentPlayer()] << "(" << currentColor() << ")" << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  // Select if location will be entered by user input (HUMAN) or random choice (COMPUTER)
  if (player[currentPlayer()] == gamePlayer::HUMAN)  
    currentPos = readLocation();
  else
    currentPos = randomLocation();
  b.setPosColor(currentPos,currentColor());	// Set entered position with player's color
  // Check if start virtual node and end virtual node are connected (in the same set)
  if (b.hasConnectedPath(currentPlayer()))
    return currentPlayer();
  return -1;
}

// Validate user location input and return node number chosen by player
Coord Game::readLocation() {
  string userLocation;
  char col;
  int lin;
  Coord location;
  bool valid;
  boost::regex collinformat( "^([A-Za-z])([0-9]+)$", boost::regex::icase);	// regex that matches with LetterNumber format
  boost::regex lincolformat( "^([0-9]+)([A-Za-z])$", boost::regex::icase);	// regex that matches with NumberLetter format
  boost::match_results<string::const_iterator> result;
  do {
    do {
      valid = true;
      cout << "Enter location: ";
      cin >> userLocation;
      if (regex_match(userLocation, result, collinformat)) {	// Match coordinate if it is in LetterNumber format
	col = toupper(string(result[1]).at(0));		// Parse column letter
	lin = stoi(result[2]);				// Parse row number
      }
      else
	if (regex_match(userLocation, result, lincolformat)) {	// Match coordinate if it is in NumberLetter format
	  col = toupper(string(result[2]).at(0));	// Parse column letter
	  lin = stoi(result[1]);			// Parse row number
	}
	else
	  valid = false;
      if (valid) {
	location = {lin,(col-'A'+1)};
	if (!b.isCoordValid(location))	// Check if entered coordinate is valid in the board
	  valid=false;
      }
      if (!valid)	// If coordinate is not in a valid formar shows appropriate error message
	cout << "ERROR - Invalid position!" << endl;
    } while (!valid);
    if (b.isCoordOccupied(location)) {	// If coordinate is occupied in the board shows appropriate error message
      cout << "ERROR - Location already occupied!" << endl;
      valid = false;
    }
  } while (!valid);
  return location;
}

// Random generate a valid location
Coord Game::randomLocation() {
  vector<Coord> available = b.availableCoord();	// Create a vector of valid available coordinates in the board
  cout << "Enter location: ";
  fflush(stdout);
  sleep(1);
  Coord randCoord = available[rand() % available.size()];	// Select a valid available random coordinate
  cout << static_cast<char>(randCoord.getj()+'A'-1);
  fflush(stdout);
  sleep(1);
  cout << randCoord.geti() << endl;
  fflush(stdout);
  sleep(1);  
  return randCoord;
}

//==============================================================================
// Main Function
//==============================================================================
int main()
{
  Game g;
  g.init();
  return 0;  
}