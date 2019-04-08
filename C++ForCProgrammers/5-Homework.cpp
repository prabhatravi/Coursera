// Hex game (with Monte Carlo moves)
// 
// Compile with: g++ 5-Homework.cpp -o 5-Homework -lboost_regex -std=c++11
//
// Brief explanation:
// The board is represented by a graph.
// Each board position corresponds to one node identified by a number and marked with one color
// (WHITE=available, BLUE=player1 and RED=player2).
// There is an edge connecting two nodes if those nodes correspond to neighbor positions in the board. 
// There are 4 extra virtual nodes (WEST, EAST, NORTH and SOUTH) connected to nodes of the start and end positions in the board.
// WEST and EAST start BLUE. NORTH and SOUTH start RED.
// Each node is initialized in its own set.
// At each move, the selected position is marked with player's color and the node set is combined with sets of all neighbors
// with the same color.
// If start virtual node and end virtual node are in the same set, so there is a path connecting them, and there is a winner.
//
// Computer plays using Monte Carlo moves at each round.
// To implement a Monte Carlo move, every available position in the board is evaluated to check what is the best next move. 
// A Monte Carlo probability (probability of win in a set of at least 1000 random simulations) is calculated and the position with
// higher probability is chosen.
// To optimize the probability calculation, a simulation is interrupted before its end if that position can't beat the best current
// probability anymore. 
// Many internal structures of the board and graph were modified to optimize Monte Carlo evaluation, resulting in significant
// time reductions (more than 50% of total time).
// As a result, in a 11x11 board, computer is able to choose the best position in less than 60s. 
//
#include <iostream>
#include <iomanip>
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

// SIMUL set the default number of simulations for each Monte Carlo move evaluation  
const int SIMUL = 1000;

// INFINIT is used to represent no edge/path between two nodes 
const int INFINIT = numeric_limits<int>::max();

// Node color definition 
// WHITE means available position
// BLACK means invalid color
// BLUE and RED represent player's colors
enum class Color {BLUE, RED, WHITE, BLACK};

// Overload operator << to print a Color
ostream &operator<<(ostream& output, const Color &c) {
  switch (c) {
    case Color::BLUE:	output << "X"; break;
    case Color::RED:	output << "O"; break;
    case Color::WHITE:	output << "."; break;
    default: 		output << "~"; break;
  }  
  return output;
}

// Game player definition 
// HUMAN means user
// COMPUTER means Monte Carlo moves
enum class gamePlayer {HUMAN, COMPUTER};

// Overload operator << to print gamePlayer
ostream &operator<<(ostream &output, const gamePlayer &p) {
  switch (p) {
    case gamePlayer::HUMAN:	output << "Human"; break;
    case gamePlayer::COMPUTER:	output << "Computer"; break;
    default:			output << "Unknown";
  }
  return output;
}

//==============================================================================
// Class Node
// Used to store information (number and color) about nodes of a graph 
//==============================================================================
class Node {
  public:
    Node():number(-1),color(Color::WHITE) {};
    Node(int n, Color c=Color::WHITE):number(n),color(c) {};
    ~Node() {};
    Node& operator=(const Node &n);
    bool operator==(const Node &n) const;
    friend ostream &operator<<(ostream &output,const Node &n);
    void setNumber(const int &n) { number = n; };	// Set node number
    int getNumber() const { return number; };		// Return node number
    void setColor(const Color &c) { color = c; };	// Set node color
    Color getColor() const { return color; };		// Return node color
    
  private:
    int number;
    Color color;
};

// Assignment operator definition for Node
Node& Node::operator=(const Node &n) {
  this->number = n.number;
  this->color = n.color;
  return *this;
}

// Return true if two nodes have the same node number and false otherwise
bool Node::operator==(const Node &n) const {
  return (this->number == n.number);
}

// Overload operator << to print Node
ostream &operator<<(ostream &output,const Node &n) {
    output << "(" << n.number << "," << n.color << ")";
  return output;
}

// Overload operator << to print list<Node>
ostream &operator<<(ostream &output, const list<Node> &L) {
  for (auto i:L)
    output << i << " ";
  return output;
}

//==============================================================================
// Graph Class
// Represents a Graph through an adjacency matrix
//==============================================================================
class Graph
{
  public:
    Graph() {};
    Graph(int numVertices, int initialValue);
    void setNodeColor(const int &x, const Color &c) { nodes[x].setColor(c); };			// Set the color of a given node
    Color getNodeColor(const int &x) const { return nodes[x].getColor(); };			// Return the color of a given node
    void setEdgeWeight(const int &x, const int &y, const int &weight);				// Set the edge weight between two given nodes
    int getEdgeWeight(const int &x, const int &y) const { return adjMatrix[x*numV+y]; };	// Return the edge weight between two given nodes
    bool adjacent(const int &x, const int &y) const { return (adjMatrix[x*numV+y]!=INFINIT); };	// Return true if two nodes are neighbors and false otherwise
    list<Node> neighbors(const int &x) const;
    int V() const { return numV; };
    int E() const { return numE; };
    list<Node> vertices() const;
    void show() const;

  private:
    int numV;			// Number of nodes of the Graph
    int numE;			// Number of edges of the Graph
    vector<Node> nodes;		// Vector with all nodes of the Graph		
    vector<int> adjMatrix;	// Adjacency matrix representing the Graph
};

// Initialize the graph setting all node's color to WHITE (default) and all edge's weight to INFINIT
Graph::Graph(int numVertices, int initialValue=INFINIT) {
  numV = numVertices;
  numE = 0;
  nodes.resize(numV);
  for (int x=0; x<numV; ++x)
    nodes[x].setNumber(x);
  adjMatrix.resize(numV*numV, initialValue);
}

// Set edge weight between two nodes
void Graph::setEdgeWeight(const int &x, const int &y, const int &weight) {
  if (adjMatrix[x*numV+y] == INFINIT)
    ++numE;
  adjMatrix[x*numV+y] = adjMatrix[y*numV+x] = weight;
}

// Return the neighbors list of a node
list<Node> Graph::neighbors(const int &x) const {
  list<Node> adjNodes;
  for (int y=0; y<numV;++y)
    if (adjMatrix[x*numV+y] != INFINIT)
      adjNodes.push_back(nodes[y]);
  return adjNodes;
}

// Return a list of all nodes in the Graph
list<Node> Graph::vertices() const {
  list<Node> allNodes;
  for (int x=0; x<numV;++x)
    allNodes.push_back(nodes[x]);
  return allNodes;
}

// Print out adjacency matrix representing the Graph (used for debug purposes)
void Graph::show() const {
  cout << setw(3) << " ";
  for (int y=0; y<numV;++y)
    cout << setw(3) << nodes[y].getNumber();
  cout << endl;
  for (int x=0; x<numV;++x) {
    cout << setw(3) << nodes[x].getNumber();
    for (int y=0; y<numV;++y)
      if (adjMatrix[x*numV + y] != INFINIT)
	cout << setw(3) << adjMatrix[x*numV + y];
      else
	cout << setw(3) << "-";
    cout << endl;
  }
}

//==============================================================================
// Coordinate Class
// Used to represent a coordinate (i,j) 
//==============================================================================
class Coord {
  public:
    Coord():i(0),j(0) {};
    Coord(const int &i, const int &j):i(i),j(j) {};
    Coord operator=(const Coord &c) { i=c.i; j=c.j; return *this; };
    bool operator==(const Coord &c) { return ((i==c.i) && (j==c.j)); };
    friend ostream &operator<<(ostream &output, const Coord &c);
    int geti() const { return i; };	// Return i component of a coordinate (i,j)
    int getj() const { return j; };	// Return j component of a coordinate (i,j)
    
  private:
    int i, j;
};

// Overload operator << to print a Coord element
ostream &operator<<(ostream &output, const Coord &c) {
  output << "(" << c.i << "," << c.j << ")";
  return output;
}

// Overload operator << to print a vector<Coord>
ostream &operator<<(ostream &output, const vector<Coord> &vec) {
  for (auto i:vec) 
    output << i << " ";
  return output;
}

//==============================================================================
// Board Class
// Used to represent the board of the hex game 
//==============================================================================
class Board {
  public:
    Board() {};
    Board(const int &d);	// Initialize the board with dimension (d x d)
    ~Board() {};
    void init(const int &d);
    friend ostream &operator<<(ostream &output, const Board &b);
    int getDimension() const { return dimension; };	// Return the dimension of the board
    int coordToPos(const Coord &c) const { return (c.geti()*dimension+c.getj()); };	// Transform coordinates (i,j) into positions in the board
    Color getPosColor(const Coord &c) const;
    void setPosColor(const Coord &coordinate, const Color &color);
    bool isCoordValid(const Coord &c) const;
    bool isCoordOccupied(const Coord &c) const;
    bool hasConnectedPath(const int &player) const;
    vector<Coord> availableCoord() const;
    
  private:
    int dimension;		// Dimension of the board (dimension x dimension)
    Graph g;			// Graph that represents the board
    int startNode[2];		// Store start node for player1 (0) and player2 (1)
    int endNode[2];		// Store end node for player1 (0) and player2 (1)
    vector<int> set;		// Vector used to store nodes set

    void combineSet(const int &u, const int &v);
    void combineNeighborSet(const int &u);
    void initSets();
    bool sameSet(const int &u, const int &v) const;
};

// Create a graph with d^2 nodes and 4 virtual nodes
// Each node of the graph represents one position of the board
// Each edge of the graph connects adjacent positions in the board
// Each virtual node represents a special node (WEST, EAST, NORTH and SOUTH)
Board::Board(const int &d) {
  dimension = d;
  const int edgeweight = 1;
  g = Graph(d*d+4);	// Initialize the graph with dimension^2 nodes + 4 virtual nodes
  for (int i=0; i<dimension; ++i)
    for (int j=0; j<dimension; ++j) {
      // neighbors is a vector containing coordinates of all positions adjacent to board position (i,j)
      vector<Coord> neighbors = {{i-1,j},{i-1,j+1},{i,j-1},{i+1, j-1}, {i+1, j},{i,j+1}};
      // For each valid neighbor position, create an edge between it and node (i,j)
      for (auto k:neighbors)
	if (isCoordValid(k))
	  g.setEdgeWeight(coordToPos({i,j}),coordToPos(k),edgeweight);
    }
  int vnodeWest, vnodeEast, vnodeNorth, vnodeSouth;
  vnodeWest = d*d;	// Node (d^2) represents WEST virtual node
  g.setNodeColor(vnodeWest,Color::BLUE);	// Assign BLUE color to WEST virtual node
  for (int i=0; i<d; ++i)			// Connect all nodes in the column 0 to WEST virtual node
    g.setEdgeWeight(vnodeWest,coordToPos({i,0}),edgeweight);
  vnodeEast = d*d+1;	// Node (d^2+1) represents EAST virtual node
  g.setNodeColor(vnodeEast,Color::BLUE);	// Assign BLUE color to EAST virtual node
  for (int i=0; i<d; ++i)			// Connect all nodes in the column 'd-1' to EAST virtual node
    g.setEdgeWeight(vnodeEast,coordToPos({i,d-1}),edgeweight);
  vnodeNorth = d*d+2;	// Node (d^2+2) represents NORTH virtual node
  g.setNodeColor(vnodeNorth,Color::RED);	// Assign RED color to NORTH virtual node
  for (int j=0; j<d; ++j)			// Connect all nodes in the row 0 to NORTH virtual node
    g.setEdgeWeight(vnodeNorth,coordToPos({0,j}),edgeweight);
  vnodeSouth = d*d+3;	// Node (d^2+3) represents SOUTH virtual node
  g.setNodeColor(vnodeSouth,Color::RED);	// Assign RED color to SOUTH virtual node
  for (int j=0; j<d; ++j)			// Connect all nodes in the row 'd-1' to SOUTH virtual node
    g.setEdgeWeight(vnodeSouth,coordToPos({d-1,j}),edgeweight);
  startNode[0]=vnodeWest;	// Assign WEST to start node of player1
  endNode[0]=vnodeEast;		// Assign EAST to end node of player1
  startNode[1]=vnodeNorth;	// Assign NORTH to start node of player2
  endNode[1]=vnodeSouth;	// Assign SOUTH to end node of player2  
  initSets();		// Assign each node of the graph to its own set
}

// Move all nodes in the same set of node 'v' to the set of node 'u' 
void Board::combineSet(const int &u, const int &v) {
  int uSet = set[u], vSet = set[v];
  for (auto &i:set)
    if (i == vSet)
      i = uSet;
}

// Add each node of the graph to its own set
void Board::initSets() {
  list<Node> allNodes = g.vertices();  
  set.resize(allNodes.size());
  for (auto i:allNodes)
    set[i.getNumber()] = i.getNumber();
}

// Move n and all neighbors of 'n' with n's color to the same set 
void Board::combineNeighborSet(const int &u) {
  for (auto i:g.neighbors(u))
    if (g.getNodeColor(u)==g.getNodeColor(i.getNumber()))
      combineSet(u,i.getNumber());
}

// Return true if two nodes are in the same set and false otherwise 
bool Board::sameSet(const int &u, const int &v) const {
  return (set[u] == set[v]);
}

// Return true if there is a path from startNode to endNode of that player 
bool Board::hasConnectedPath(const int &player) const {
  return sameSet(startNode[player],endNode[player]);
}

// Return the color of a board position 
Color Board::getPosColor(const Coord &c) const {
  return g.getNodeColor(coordToPos(c));
}

// Set the color of a board position 
void Board::setPosColor(const Coord &coordinate, const Color &color) {
  g.setNodeColor(coordToPos(coordinate), color);
  combineNeighborSet(coordToPos(coordinate));
}

// Return true if 'c' is a valid coordinate in the board an false otherwise 
bool Board::isCoordValid(const Coord &c) const {
  return (((c.geti()>=0) && (c.geti()<dimension)) && ((c.getj()>=0) && (c.getj()<dimension)));
}

// Return true if 'c' is an occupied position in the board and false otherwise 
bool Board::isCoordOccupied(const Coord &c) const {
  return ((getPosColor(c)==Color::RED) || ((getPosColor(c))==Color::BLUE));
}

// Return a vector with available (free) positions in the board 
vector<Coord> Board::availableCoord() const {
  vector<Coord> available;
  for (int i=0; i<dimension; ++i)
    for (int j=0; j<dimension; ++j)
      if (getPosColor({i,j})==Color::WHITE)
	available.push_back({i,j});
  return available;
}

// Overload operator << to print a board
ostream &operator<<(ostream &output, const Board &b) {
  // Show NORTH label
  output << setw(4*b.dimension/2+10) << "NORTH" << endl << endl;
  // Show column names (letters)
  output << setw(4) << " " << setw(2) << " ";
  for (int j=0; j<b.dimension; ++j)
    output << setw(4) << static_cast<char>(j+'A');
  output << endl << endl;
  // Show line numbers, board positions and line numbers again
  for (int i=0; i<b.dimension; ++i) {
    if (i==(b.dimension/2))
      output << setw(i*2+4) << "WEST";
    else
      output << setw(i*2+4) << " ";
    output << setw(4) << i+1 ;      
    for (int j=0; j<b.dimension; ++j)
      output << setw(4) << b.g.getNodeColor(b.coordToPos({i,j}));
    output << setw(4) << i+1 ;
    if (i==(b.dimension/2))
      output << setw(8) << "EAST";
    output << endl << endl;
  }
  output << setw(4) << " ";
  // Show column names (letters) again
  output << setw((b.dimension)*2+4) << " ";
  for (int j=0; j<b.dimension; ++j)
    output << setw(4) << static_cast<char>(j+'A');
  output << endl;
  // Show "SOUTH" label
  output << endl;
  output << setw(4*(b.dimension+3)) << "SOUTH" << endl;
  return output;
}

//==============================================================================
// ProgressBar Class
// Used to print out a progress bar
//==============================================================================
class ProgressBar {
  public:
    ProgressBar():barProgress(0.0),barWidth(70) {};
    void init(const string &text, const int &width, const double &progress=0.0);
    void draw() const;
    void update(const double &progress);
   
  private:
    double barProgress;			// Store the completed percentage
    int barWidth;			// Store the bar width
    string barText;			// Store the text showed before the progress bar (left side)
    double starttime, currenttime;	// Store the initial and elapsed time to show after the bar (right side)
};

// Initialize the progress bar attributes and start counting the time
void ProgressBar::init(const string &text, const int &width, const double &progress) {
  barWidth = width;
  barProgress = progress;
  barText = text;
  starttime = static_cast<double>(clock());
}
 
// Draw the progress bar according to its attributes
void ProgressBar::draw() const {
  cout << "\r" << barText << " [";
  int pos = barWidth * barProgress;
  for (int i = 0; i < barWidth; ++i)
    if (i < pos)
      cout << "#";
    else
      cout << " ";
  cout << "] " << setw(3) << int(barProgress * 100.0) << "% (" << fixed << setprecision(2) << currenttime << "s)";
  cout.flush();
}

// Update completed percentage and elapsed time of progress bar and draw it again
void ProgressBar::update(const double &progress) {
  if (progress<0.0)
    barProgress = 0.0;
  else
    if (progress>1.0)
      barProgress = 1.0;
    else
      barProgress = progress;
  currenttime = static_cast<double>(clock() - starttime) / CLOCKS_PER_SEC;
  draw();
}

//==============================================================================
// Game Class
// Used to control the game flow 
//==============================================================================
class Game {
  public:
    Game();
    ~Game() {};
    void init();
    
  private:
    Board b;			// Store the board used in the game
    int round;			// Store the round number
    int winner;			// Store the winner of the game (-1: No winner, 0: Player1 and 1: Player2) 
    gamePlayer player[2];	// Store what kind of player (HUMAN or COMPUTER) is Player1 (position 0) and Player2 (position 1) 

    int currentPlayer() const { return (round%2); };			// Return the current player
    Color currentColor() const { return static_cast<Color>(round%2); };	// Return the color of current player
    int playRound(const int &player, const Color &color, const Coord &pos);
    int currentMove() const { return (round/2)+1; };			// Return the current move number
    Coord readLocation() const;
    Coord mcLocation() const;
    Coord nextMove(const Game &g) const;
    double probMonteCarlo(Game gamecopy, const Coord &c, const double &bestprob, const int &numsim=SIMUL) const;
};

// Constructor
// Initialize some internal attributes to default value and seed random number generator
Game::Game() {
  player[0] = player[1] = gamePlayer::COMPUTER;
  round = 0;
  winner = -1;
  srand(time(NULL));
}

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
  cout << "* " << setw(8) << player[0] << "       *    " << static_cast<Color>(0) << "    * connects West-East            *" << endl;
  cout << "* " << setw(8) << player[1] << "       *    " << static_cast<Color>(1) << "    * connects North-South          *" << endl;
  cout << "************************************************************" << endl;
  cout << "-----------------------------------------------------------------------" << endl;
  // Create the board with users's specified dimension
  b = Board(dimension);
  // Play round by round until we have a winner
  while (winner==-1) {
    cout << "Board: " << endl;
    cout << b << endl;
    cout << "-----------------------------------------------------------------------" << endl;
    cout << "Move " << currentMove() << ": " << player[currentPlayer()] << "(" << currentColor() << ")" << endl;
    cout << "-----------------------------------------------------------------------" << endl;
    Coord currentPos;
    // Select if location will be entered by user input (HUMAN) or Monte Carlo choice (COMPUTER)
    if (player[currentPlayer()] == gamePlayer::HUMAN)  
      currentPos = readLocation();
    else
      currentPos = mcLocation();
    winner = playRound(currentPlayer(),currentColor(),currentPos);
  }
  // Show winner information
  cout << "-----------------------------------------------------------------------" << endl;
  cout << "Winner: " << player[winner] << "(" << static_cast<Color>(winner) << ")" << endl;
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
// Return the winner of the game (-1: No winner, 0: Player1 and 1: Player2)
int Game::playRound(const int &player, const Color &color, const Coord &pos) {
  b.setPosColor(pos,color);		// Set entered position with player's color
  round++;				// Increment round number
  if (b.hasConnectedPath(player))	// Check if start virtual node and end virtual node are connected
    return player;			// If there is a winner, return the winner number (0: Player1 and 1: Player2)
  return -1;				// If there is no winner, return -1
}

// Validate user location input and return the coordinate in the board chosen by player
Coord Game::readLocation() const {
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
	location = {lin-1,(col-'A')};
	if (!b.isCoordValid(location))	// Check if entered coordinate is valid in the board
	  valid=false;
      }
      if (!valid)	// If coordinate is not in a valid format, shows appropriate error message
	cout << "ERROR - Invalid position!" << endl;
    } while (!valid);
    if (b.isCoordOccupied(location)) {	// If coordinate is occupied in the board, shows appropriate error message
      cout << "ERROR - Location already occupied!" << endl;
      valid = false;
    }
  } while (!valid);
  return location;
}

// Return the coordinate in the board generated by MonteCarlo evaluation
Coord Game::mcLocation() const {
  Coord mcCoord = nextMove(*this);	// Get the best next move 
  cout << "Enter location: ";
  cout.flush();
  sleep(1);
  cout << static_cast<char>(mcCoord.getj()+'A');
  cout.flush();
  sleep(1);
  cout << mcCoord.geti()+1 << endl;
  cout.flush();
  sleep(1);  
  return mcCoord;
}

// Return the best next move based on MonteCarlo evaluation   
Coord Game::nextMove(const Game &g) const {
  double bestprob = -1.0, probMC = 0.0;
  Coord bestMove;
  ProgressBar bar;
  double visited=0.0, total=0.0;
  vector<Coord> candidates = g.b.availableCoord();	// Set all available positions in the board as candidates
  random_shuffle(candidates.begin(), candidates.end());	// Shuffle candidates (to be fair)  
  bar.init("Thinking",55,0.0);				// Create the progress bar and start counting time
  bar.draw();						// Draw an empty progress bar
  total=static_cast<double>(candidates.size());		// Set the total number of candidates to be evaluated
  for (auto i:candidates) {	
    probMC = probMonteCarlo(g, i, bestprob); 		// For each candidate position, evaluate its Monte Carlo probability
    if (bestprob < probMC) {				// If its Monte Carlo probability is the best known, set this candidate as the best move
      bestprob = probMC;
      bestMove = i;
    }
    visited++;						// Update the visited candidates counter
    bar.update(visited/total);				// Draw an updated progress bar
  }
  cout << endl;
  return bestMove;
}

// Return the MonteCarlo evaluation (probability of win in a set of random simulations)   
double Game::probMonteCarlo(Game gamecopy, const Coord &c, const double &bestprob, const int &numsim) const {
  int numwins=0, mcplayer=gamecopy.currentPlayer();
  Game g;
  vector<Coord> available, availablecopy;
  vector<Coord>::iterator posptr;
  
  gamecopy.winner = gamecopy.playRound(gamecopy.currentPlayer(),gamecopy.currentColor(),c);
  availablecopy=gamecopy.b.availableCoord();
//  cout << "Evaluating position: " << c;
//  clock_t start_time = clock();
  int i = 0;
  while ((i < numsim) && (((numsim-i)+numwins)>(bestprob*numsim))) {	// If this position can't beat bestprob, interrupt simulation
    g = gamecopy;
    available = availablecopy;							// Get a copy of all avalable positions in the board
    random_shuffle(available.begin(), available.end());				// Shuffle available positions
    posptr = available.begin();							// Point to the first random move
    while (g.winner==-1) {							// Play the game until there is a winner
      g.winner = g.playRound(g.currentPlayer(),g.currentColor(),(*posptr));	// Play in a random available position
      posptr++;									// Go to the next random position
    }
    if (g.winner == mcplayer)	// If the winner is the function caller, increments the number of wins 
      numwins++;
    ++i;
  }
//  cout << " (" << i << " simulations) - elapsed time: " << (static_cast<double>(clock() - start_time)) / CLOCKS_PER_SEC << "s - numwins: " << numwins <<  endl;
  return static_cast<double>(numwins)/static_cast<double>(numsim);
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