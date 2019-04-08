// Compile with: g++ 3-Homework.cpp -o 3-Homework -lboost_regex -std=c++11
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <ctime>
#include <cstdlib>
#include <limits>		// Max integer definition library
#include <fstream>		// File manipulation library
#include <string>		// String manipulation library
#include <boost/regex.hpp> 	// Boost regular expression library
using namespace std;

//==============================================================================
// General definitions
//==============================================================================

// INFINIT is used to represent no edge/path between two nodes 
const int INFINIT = numeric_limits<int>::max();

//==============================================================================
// Class Node
// Used to store information about nodes/edges in the adjacency list of a graph
// Adjacency lists is a list of Nodes (identified by numbers)
// Each node contains a list of neighbors containing edge weight
//==============================================================================
class Node {
  public:
    Node():number(-1),weight(INFINIT) { edges.clear(); };
    Node(int n, int w=INFINIT):number(n),weight(w) { edges.clear(); };
    ~Node() {};
    Node& operator=(const Node &n);
    bool operator== (const Node& n) const;
    void setNodeNumber(int n) { number = n; };
    int getNodeNumber() { return number; };
    void setEdgeWeight(int w=INFINIT) { weight = w; };
    int getEdgeWeight() { return weight; };
    void setEdge(Node N) { edges.push_back(N); };
    list<Node> getEdges() { return edges; };
    
  private:
    int number;
    int weight;	
    list<Node> edges;	
};

// Assignment operator definition for Node
Node& Node::operator=(const Node &n)
{
  this->number = n.number;
  this->weight = n.weight;
  this->edges = n.edges;
  return *this;
}

// Return true if two nodes have the same node number and false otherwise
bool Node::operator== (const Node& n) const {
  if (this->number == n.number) return true;
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
   for(list<Node>::iterator i=L.begin(); i != L.end(); ++i) {
     output << (*i) << " ";
   }
   return output;
}

//==============================================================================
// Graph Class
// Represent a Graph through an adjacency list of nodes
//==============================================================================
class Graph {
  public:
    Graph():numV(0),numE(0) { adjList.clear(); };	// Create an empty graph
    Graph(int numVertices);
    Graph(string filename);
    ~Graph() {};
    int get_edge_value(int x, int y);
    void set_edge_value(int x, int y, int value);
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

// Constructor of Graph Class
// Create an adjacency list with all nodes and no edges
Graph::Graph(int numVertices) {
  numV = numVertices;
  numE = 0;
  adjList.clear();
  for(int i=0; i<numVertices; ++i) {
    Node newNode(i);
    adjList.push_back(newNode);
  }
}

// Constructor of Graph Class
// Read the graph information from file passed as argument
Graph::Graph(string filename) {
  ifstream file (filename);
  string firstline, otherlines;
  boost::regex firstlineformat( "^\\d+$", boost::regex::icase);	// regex that matches with the first line format
  boost::regex otherlinesformat( "^\\s*(\\d+)\\s+(\\d+)\\s+(\\d+)\\s*$", boost::regex::icase);	// regex that matches with edge lines format
  boost::match_results<string::const_iterator> result;
  
  // Initialize an empty graph
  numV = 0;
  numE = 0;
  adjList.clear();

  if (file.is_open()) {	// Open file 'filename'
    getline(file,firstline);	// Read the first line of file
    if (regex_match(firstline, result, firstlineformat)) {	// Check if the first line is in the apropriate format (only one number)
      // Create nodes
      numV = stoi(result[0]);	// Assign the number read from file to the number of nodes in the graph (numV)
      numE = 0;
      // Create adjacency list with all nodes and no edges
      adjList.clear();
      for(int i=0; i<numV; ++i) {
	Node newNode(i);
	adjList.push_back(newNode);
      }
    }
    
    while(getline(file,otherlines))	// Read one edge line of file
      if (regex_match(otherlines, result, otherlinesformat)) {	// Check if the line is in the apropriate format: i j cost (three numbers)
	int i = stoi(result[1]);
	int j = stoi(result[2]);
	int cost = stoi(result[3]);
	set_edge_value(i,j,cost);	// Create the desired edge in the graph
      }
    file.close();	// Close file 'filename'
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
  
  // Add 'y' in the list of 'x' neighbors (if doesn't exist)
  // Set edge weight to value
  for(list<Node>::iterator i=adjList.begin(); i != adjList.end(); ++i) {
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
// Class Edge
// Used to store information about edges (source, destination and weight)
//==============================================================================
class Edge {
  public:
    Edge():src(-1),dst(-1),weight(INFINIT) {};
    Edge(int s, int d, int w):src(s),dst(d),weight(w) {};
    ~Edge() {};
    Edge& operator=(const Edge &e);
    bool operator== (const Edge& e) const;
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
Edge& Edge::operator=(const Edge &e)
{
  this->src = e.src;
  this->dst = e.dst;
  this->weight = e.weight;
  return *this;
}

// Return true if two edges have the same src and dst and false otherwise
bool Edge::operator== (const Edge& e) const {
  if ((this->src == e.src) && (this->dst == e.dst)) return true;
  return false;
}

// Return true if 'e1' weight is less then 'e2' weight and false otherwise
bool operator< (Edge& e1, Edge& e2) {
  if (e1.getWeight() < e2.getWeight()) return true;
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
// PriorityQueue Class (Generic)
// Store ordered elements of any TYPE
//==============================================================================
template<class TYPE>
class PriorityQueue {
  public:
    PriorityQueue() { pq.clear(); };	// Create an empty list of TYPE
    ~PriorityQueue() {};
    void insert(TYPE e);
    TYPE top();
    int size();
    
  private:
    list<TYPE> pq;	// List of ordered elements from min (higher priority) to max (lower priority)
};

// Insert element into priority queue
template<class TYPE>
void PriorityQueue<TYPE>::insert(TYPE e) {
  typename list<TYPE>::iterator found = find(pq.begin(), pq.end(), e);
  if (found == pq.end()) {
    pq.push_back(e);
    pq.sort();
  }
}

// Return the element with higher priority (firts element), removing it from the queue
template<class TYPE>
TYPE PriorityQueue<TYPE>::top() {
  TYPE e;
  if (! pq.empty()) {
    typename list<TYPE>::iterator i=pq.begin();
    e = (*i);
    pq.pop_front();
  }
  return e;
}

// Return the number of elements in the priority queue
template<class TYPE>
int PriorityQueue<TYPE>::size() {
  return pq.size();
}

// Overload operator << to print a priority queue
template<class TYPE>
ostream &operator<<(ostream &output, PriorityQueue<TYPE> p) {
  while (p.size() != 0)
    output << p.top() << " ";
  return output;
}

//==============================================================================
// PrimMST Class
// Implements Prim's Algorithm to find a minimum spanning tree in a graph
//==============================================================================
class PrimMST
{
  public:
    PrimMST() {};
    PrimMST(Graph g):graph(g) {};
    ~PrimMST() {};
    list<Edge> tree();
    int cost();

  private:
    Graph graph;	// Graph used by Prim's Algorithm
    bool isSrcInSet (Edge e, list<Node> L);
    bool isDstInSet (Edge e, list<Node> L);
    Node getNodeInSet (int num, list<Node> L);
};

// Return true if edges 'e' src is in node list 'L' and false otherwise
bool PrimMST::isSrcInSet (Edge e, list<Node> L) {
  list<Node>::iterator found = find(L.begin(), L.end(), e.getSrc());
  if (found != L.end()) return true;
  return false;
}

// Return true if edges 'e' dst is in node list 'L' and false otherwise
bool PrimMST::isDstInSet (Edge e, list<Node> L) {
  list<Node>::iterator found = find(L.begin(), L.end(), e.getDst());
  if (found != L.end()) return true;
  return false;
}

// Return node with nodeNumber 'num' in node list 'L'
Node PrimMST::getNodeInSet (int num, list<Node> L) {
  Node n(num);
  list<Node>::iterator found = find(L.begin(), L.end(), n);
  if (found != L.end()) return (*found);
}

// Return a list<Edge> containing the list of edges in the minimum spanning tree found by Prim's Algorithm
list<Edge> PrimMST::tree() {
  list<Edge> mst;
  PriorityQueue<Edge> p;
  Node selected;
  Edge selectedEdge;
  
  // Create selectedNodes empty set
  // Create candidateNodes set containing all nodes
  list<Node> selectedNodes;
  list<Node> candidateNodes = graph.vertices();

  while (!candidateNodes.empty()) {
    // Add first node from candidateNodes set to selectedNodes set and remove it from candidateNodes
    selected = candidateNodes.front();
    selectedNodes.push_back(selected);
    candidateNodes.remove(selected);
    
    // Insert all edges from first node into priority queue
    list<Node> sEdges = graph.neighbors(selected.getNodeNumber());
    for(list<Node>::iterator j=sEdges.begin(); j != sEdges.end(); ++j) {
      Edge e(selected.getNodeNumber(),(*j).getNodeNumber(),(*j).getEdgeWeight());
      p.insert(e);
    }

    while(p.size() != 0) {
      selectedEdge = p.top();	// Remove the lower edge from the priority queue
      if ((isSrcInSet(selectedEdge,selectedNodes)) && (isDstInSet(selectedEdge,candidateNodes))) {
	mst.push_back(selectedEdge);
	// Add selected node to selectedNodes set and remove it from candidateNodes set
	selected = getNodeInSet(selectedEdge.getDst(),candidateNodes);
	selectedNodes.push_back(selected);
	candidateNodes.remove(selected);
	// Insert all edges from dst node into priority queue
	sEdges = graph.neighbors(selected.getNodeNumber());
	for(list<Node>::iterator j=sEdges.begin(); j != sEdges.end(); ++j) {
	  Edge e(selected.getNodeNumber(),(*j).getNodeNumber(),(*j).getEdgeWeight());
	  p.insert(e);
	}
      } 
    }
  }
     
  return mst;
}

// Return the cost of the minimum spanning tree found by Prim's Algorithm
int PrimMST::cost() {
  int mstCost=0;
  list<Edge> mst;
  mst = tree();
  for(list<Edge>::iterator i=mst.begin(); i != mst.end(); ++i)
    mstCost += (*i).getWeight();
  return mstCost;
}

//==============================================================================
// KruskalMST Class
// Implements Kruskal's Algorithm to find a minimum spanning tree in a graph
//==============================================================================
class KruskalMST
{
  public:
    KruskalMST() {};
    KruskalMST(Graph g):graph(g) {};
    ~KruskalMST() {};
    list<Edge> tree();
    int cost();

  private:
    Graph graph;	// Graph used by Kruskal's Algorithm
    map<int, int> set;	// Map used to store nodes set
    void makeSet(Node n);
    int findSet(Node n);
    void combineSet(Node u, Node v);
};

// Assign node 'n' to a new set 
void KruskalMST::makeSet(Node n) {
  static int setNumber = 0;
  set[n.getNodeNumber()] = setNumber++;
}

// Return the set of node 'n' 
int KruskalMST::findSet(Node n) {
  return set[n.getNodeNumber()]; 
}

// Move all nodes in the same set of node 'v' to the set of node 'u' 
void KruskalMST::combineSet(Node u, Node v) {
  int uSet = set[u.getNodeNumber()], vSet = set[v.getNodeNumber()];
  for (map<int,int>::iterator it=set.begin(); it != set.end(); ++it)
    if (it->second == vSet)
      set[it->first] = uSet;
}

// Return a list<Edge> containing the list of edges in the minimum spanning tree found by Kruskal's Algorithm
list<Edge> KruskalMST::tree() {
  list<Edge> mst;
  PriorityQueue<Edge> p;
  Node selected;
  Edge selectedEdge;
  
  // Add each node of the graph to its own set
  list<Node> allNodes = graph.vertices();  
  for(list<Node>::iterator i=allNodes.begin(); i != allNodes.end(); ++i)
    makeSet(*i);
  
  // Insert all edges from the graph into priority queue
  for(list<Node>::iterator i=allNodes.begin(); i != allNodes.end(); ++i) {
    list<Node> iEdges = graph.neighbors((*i).getNodeNumber());
    for(list<Node>::iterator j=iEdges.begin(); j != iEdges.end(); ++j) {
      Edge e((*i).getNodeNumber(),(*j).getNodeNumber(),(*j).getEdgeWeight());
      p.insert(e);
    }
  }
 
  while(p.size() != 0) {
    selectedEdge = p.top();	// Remove the lower edge from the priority queue
    if (findSet(selectedEdge.getSrc()) != findSet(selectedEdge.getDst())) {
      mst.push_back(selectedEdge);	// If src and dst are in different sets, then add this edge to mst
      combineSet(selectedEdge.getSrc(),selectedEdge.getDst());	// Combine both sets into one
    } 
  }
     
  return mst;
}

// Return the cost of the minimum spanning tree found by Kruskal's Algorithm
int KruskalMST::cost() {
  int mstCost=0;
  list<Edge> mst;
  mst = tree();
  for(list<Edge>::iterator i=mst.begin(); i != mst.end(); ++i)
    mstCost += (*i).getWeight();
  return mstCost;
}

//==============================================================================
// Monte Carlo Class
// Used to generate random graphs and run simulations 
//==============================================================================
class MonteCarlo
{
  public:
    MonteCarlo();
    ~MonteCarlo() {};
    Graph randomGraph(int vert, double density, int minDistEdge, int maxDistEdge);
    void run(Graph g);
  
  private:
};

// Constructor of MonteCarlo Class
// Initializes the seed of random number generator
MonteCarlo::MonteCarlo()
{
  srand(time(NULL));
}

// Return a random Graph generated with number of nodes, density and edge weight range informed
Graph MonteCarlo::randomGraph(int numVert, double density, int minDistEdge, int maxDistEdge)
{
  int randDistEdge;
  char srcVert, dstVert;
  
  Graph g(numVert);

  for (int i=0; i<g.V(); ++i)
    for (int j=i+1; j<g.V(); ++j)
    {
      double p = ((static_cast<double>(rand())) / RAND_MAX);	// Generate random probability
      if (p < density)	// If random probability is less than density, edge (i,j) will be set
      {
	randDistEdge = rand() % (maxDistEdge - minDistEdge) + minDistEdge; // Generate random edge weight 
	g.set_edge_value(i,j,randDistEdge);
      }
    }
  
  return g;
}

// Run a simulation finding the Prim's and Kruskal's MST in a given graph 
void MonteCarlo::run(Graph g)
{
  static int turn=0;
  
  cout << endl << "=== RUNNING SIMULATION No. " << ++turn << " ===" << endl;
  
  // Print out graph information
  double d = static_cast<double>(g.E())/((static_cast<double>(g.V())*static_cast<double>(g.V())-1)/2)*100;	// Calculate real density reached
  cout << "Vertices: " << g.V() << endl;
  cout << "Edges: " << g.E() << " (density: " << d << "%)" << endl;
  cout << "Nodes: " << g.vertices() << endl;
  cout << "Graph: " << g << endl;
  
  PrimMST prim(g);
  cout << "Prim MST: " << prim.tree() << endl;
  cout << "Prim MST cost: " << prim.cost() << endl;
  cout << endl;

  KruskalMST kruskal(g);
  cout << "Kruskal MST: " << kruskal.tree() << endl;
  cout << "Kruskal MST cost: " << kruskal.cost() << endl;
}

//==============================================================================
// Main Function
//==============================================================================
int main()
{
  MonteCarlo simulation;

  // Reading the example graph from file and calculating MST cost
  Graph g("sampletestdata");
  simulation.run(g);

  // Generating a rendom graph (10 nodes and 20% density) and calculating MST cost
  // Graph g2 = simulation.randomGraph(10,0.2,1,10);
  // simulation.run(g2);
    
  // Generating a rendom graph (10 nodes and 40% density) and calculating MST cost
  // Graph g3 = simulation.randomGraph(10,0.4,1,10);
  // simulation.run(g3);

  return 0;  
}