#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>

template <typename vertex>
class weighted_graph {

	public:

	//You will need to add some data members here
	//to actually represent the graph internally,
	//and keep track of whatever you need to.
	std::vector<std::vector<int> > adj_matrix;
	std::vector<vertex> vertices;
	size_t size;
	int numVertices;
	int numEdges;
	int totalWeight;
	
	//The graph_iterator class provides an iterator
	//over the vertices of the graph.
	//This is one of the harder parts, so if you're
	//not too comfortable with C++ leave this for last.
	//If you are, there are many ways of doing this,
	//as long as it passes the tests, it's okay.
	class graph_iterator {

		private:

		//You may need data members here.
		const weighted_graph<vertex> * owner;
		int position;

		public:
			graph_iterator(const weighted_graph &);
			graph_iterator(const weighted_graph &, size_t);
			~graph_iterator();
			graph_iterator operator=(const graph_iterator&);
			bool operator==(const graph_iterator&) const;
			bool operator!=(const graph_iterator&) const;
			graph_iterator operator++();
			graph_iterator operator++(int);
			const vertex operator*();
			const vertex* operator->();
	};

	//The neighbour_iterator class provides an iterator
	//over the neighbours of a given vertex. This is
	//probably harder (conceptually) than the graph_iterator.
	//Unless you know how iterators work.
	class neighbour_iterator {

		private:

		//You may need data members here.

		public:
			neighbour_iterator(const neighbour_iterator&);
			neighbour_iterator(const weighted_graph &, const vertex&);
			neighbour_iterator(const weighted_graph &, const vertex&, size_t);
			~neighbour_iterator();
			neighbour_iterator operator=(const neighbour_iterator& it);
			bool operator==(const neighbour_iterator&) const;
			bool operator!=(const neighbour_iterator&) const;
			neighbour_iterator operator++();
			neighbour_iterator operator++(int);			
			const std::pair<vertex, int> operator*();
			const std::pair<const vertex, int>* operator->();
	};

	public:


	weighted_graph();
	//weighted_graph(const size_t&); //A constructor for weighted_graph with initial size as a parameter.
	~weighted_graph(); //A destructor. Depending on how you do things, this may
					   //not be necessary.

	bool are_adjacent(const vertex&, const vertex&) const; //Returns true if the two vertices are
														   //adjacent, false otherwise.
	bool has_vertex(const vertex&) const; //Returns true if the passed in vertex is 
										  //a vertex of the graph, false otherwise.

	void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
	void add_edge(const vertex&, const vertex&, const int&); //Adds an edge between the two vertices
															 //with the given weight (as an int).

	void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
	void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.
	void set_edge_weight(const vertex&, const vertex&, const int&); //Changes the edge weight between the two
																	//vertices to the new weight (the int).
	int get_index(const vertex&); //Returns the index of a vertex in the graph.
	int get_edge_weight(const vertex&, const vertex&) const; //Returns the weight on the edge between the two vertices.
	int degree(const vertex&) const; //Returns the degree of the vertex.
	int weighted_degree(const vertex&); //Returns the sum of the weights on all the edges incident to the vertex.
	int num_vertices() const; //Returns the total number of vertices in the graph.
	int num_edges() const; //Returns the total number of edges in the graph (just the count, not the weight).
	int total_weight(); //Returns the sum of all the edge weights in the graph.

	std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
	std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex.

	graph_iterator begin(); //Returns a graph_iterator pointing to the start of the vertex set.
	graph_iterator end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

	neighbour_iterator neighbours_begin(const vertex&); //Returns a neighbour_iterator pointing to the start
														//of the neighbour set for the given vertex.
	neighbour_iterator neighbours_end(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end
													  //of the neighbour set for the given vertex.

	std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they
													//are visited in by a depth-first traversal starting at
													//the given vertex.
	std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they
													  //are visisted in by a breadth-first traversal starting
													  //at the given vertex.

	weighted_graph<vertex> mst(); //Returns a minimum spanning tree of the graph.

};

//Define all your methods down here (or move them up into the header, but be careful you don't double up).
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

// graph iterator class to keep track of position by accessing the underlying pointers
template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g){ 
    owner = &g;
    position = 0;
}
// copy constructor
template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g, size_t start_pos){
    owner = &g;
    position = start_pos;
}
template <typename vertex> weighted_graph<vertex>::graph_iterator::~graph_iterator(){}
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator=(const graph_iterator& it){ 
    auto g = graph_iterator( *owner );
    return g; 
}
template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator==(const graph_iterator& it) const { 
    return this->position == it.position && this->owner == it.owner; 
}
template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator!=(const graph_iterator& it) const { 
    return (this->position != it.position || this->owner != it.owner); 
}
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(){
    this->position++;
    return *this;
}
template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(int){
    this->position++;
    return *this; 
}
template <typename vertex> const vertex weighted_graph<vertex>::graph_iterator::operator*(){ 
	auto v = *(owner->vertices.begin() + position );
    return v; 
}
template <typename vertex> const vertex* weighted_graph<vertex>::graph_iterator::operator->(){ 
    auto v = (owner->vertices.begin() + position ); 
    return v; 
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u) {}
template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u, size_t start_pos) {}
template <typename vertex> weighted_graph<vertex>::neighbour_iterator::~neighbour_iterator() {}
template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator=(const neighbour_iterator& it) { auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }
template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator==(const neighbour_iterator& it) const { return false; }
template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator!=(const neighbour_iterator& it) const { return false; }
template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++() { auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }
template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++(int){ auto n = neighbour_iterator(weighted_graph<vertex>(), vertex()); return n; }			
template <typename vertex> const std::pair<vertex, int> weighted_graph<vertex>::neighbour_iterator::operator*(){ auto p = std::pair<vertex,int>(); return p; }
template <typename vertex> const std::pair<const vertex, int>* weighted_graph<vertex>::neighbour_iterator::operator->(){ return nullptr; }

template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::begin() {
	return graph_iterator(*this);
}

template <typename vertex>	typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::end() {
	auto it = graph_iterator(*this, this->vertices.size());
    return it;
}

template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_begin(const vertex& u) {
	return neighbour_iterator(*this, u);
}

template <typename vertex>	typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_end(const vertex& u) {
	auto n_it = neighbour_iterator(*this, u, this->vertices.size()-1 );
        return n_it;
}

template <typename vertex> weighted_graph<vertex>::weighted_graph(){
	//initialise all global variables and the graph's data structures
	this->size = 0;
    this->numVertices = 0;
	this->numEdges = 0;
	this->totalWeight = 0;    
    
    vertices = std::vector<vertex>();
    
    adj_matrix = std::vector<std::vector<int> >(size, std::vector<int>(size, 0));
    adj_matrix.resize(0);
}
// possible copy constructor for future use
/*template <typename vertex> weighted_graph<vertex>::weighted_graph(const size_t& size){
	this->size = size;
	this->numVertices = 0;
	this->numEdges = 0;
	this->totalWeight = 0;
	
    vertices = std::vector<vertex>();
	adj_matrix = std::vector<std::vector<int>>(size, std::vector<int>(size, 0));
	ad_matrix.resize(0);
}*/

template <typename vertex> weighted_graph<vertex>::~weighted_graph(){
	// Nothing needed here, no dynamic allocation used
}

template <typename vertex> bool weighted_graph<vertex>::has_vertex(const vertex& u) const {
    // loops through vertices vector, if i == passed in vertex return true
	for (int i = 0; i < vertices.size(); i++) {
        if (vertices[i] == u) {
            return true;
        }
    }
	return false;

}

template <typename vertex> bool weighted_graph<vertex>::are_adjacent(const vertex& u, const vertex& v) const {
	// compare the two passed in vertices to see that they're adjacent, else return false
	if ((u >= 0) && (u < size) && (v >= 0) && (v < size)) {
		return adj_matrix[u][v];
	}
	return false;
}

template <typename vertex> void weighted_graph<vertex>::add_vertex(const vertex& v) {
	// add vertex to vertices vector
	this->vertices.push_back(v);
  	// resize matrix
	adj_matrix.resize(numVertices+1);
    for( int i = 0; i < adj_matrix.size() ; i++){
		adj_matrix[i].resize(numVertices + 1, 0);
	}
    // increment vertex count    
    numVertices++;
    size++;
}

template <typename vertex> void weighted_graph<vertex>::add_edge(const vertex& u, const vertex& v, const int& weight) {
	// add edge by accessing the indexes of the vertices and assigning the weight
	if ((u != v) && (has_vertex(u)) && (has_vertex(v))){
		adj_matrix[get_index(u)][get_index(v)] = weight;
    	adj_matrix[get_index(v)][get_index(u)] = weight;
		// increment total weight counter
		totalWeight += weight;
	}
	//increment the edge counter
	numEdges++;
}


template <typename vertex> void weighted_graph<vertex>::remove_vertex(const vertex& u) {	
	// store the position so that the index can be used to remove vertex
	int position = get_index(u);
	// delete from vertices vector
	this->vertices.erase(this->vertices.begin() + position);
	// decrement counters
	numVertices--;
	size--;
	// resize adjacency matrix
	adj_matrix.erase( adj_matrix.begin() + position  );
        for( int i = 0; i < adj_matrix.size() ; i++){
            adj_matrix[i].erase(adj_matrix[i].begin() + position);
        }
}


template <typename vertex> void weighted_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
	// remove the edge by accessing the vertices indexes and assigning them the value of 0
	adj_matrix[get_index(u)][get_index(v)] = 0;
	adj_matrix[get_index(v)][get_index(u)] = 0;
	// decrement the edge counter
	numEdges--;
}

template <typename vertex> void weighted_graph<vertex>::set_edge_weight(const vertex& u, const vertex& v, const int& weight) {
	// set given vertices to passed in weight parameter
	adj_matrix[u][v] = weight;
    adj_matrix[v][u] = weight;
}

template <typename vertex> int weighted_graph<vertex>::get_index(const vertex& u) {
	// helper function to loop through vertices vector and return the index
	// else return -1 if not found
	for ( int i = 0; i < vertices.size(); i++  ) {
        if ( vertices[i] == u ) {
            return i;
		}
     }
	 return -1;
}

template <typename vertex> int weighted_graph<vertex>::get_edge_weight(const vertex& u, const vertex& v) const {
	// check that they are neighbours then return the edge weight stored in the matrix
	if ((u != v) && (has_vertex(u)) && (has_vertex(v))) {
		return adj_matrix[u][v];
	}
}

template <typename vertex> int weighted_graph<vertex>::degree(const vertex& u) const {
	//Returns the degree of the vertex
	//loop through every vertex but u
	//if not == 0, add one to degree count;
	//return degree
	int degree_count = 0;
	for (int i=0; i < adj_matrix.size(); i++) {
		if (adj_matrix[u][i] != 0) {
			degree_count++;
		}
	}
	return degree_count;
}

template <typename vertex> int weighted_graph<vertex>::weighted_degree(const vertex& u) {
	// same as above but tally the weight stored at the index of u
	int degree_weight = 0;
	for (int i=0; i < adj_matrix.size(); i++) {
			degree_weight += adj_matrix[get_index(u)][i];
	}
	return degree_weight;
}

template <typename vertex> int weighted_graph<vertex>::num_vertices() const {
	// return the vertex counter is incremented/decremented as needed
	return this->numVertices;
}

template <typename vertex> int weighted_graph<vertex>::num_edges() const {
	// return the edge counter that is incremented/decremented as needed
	return this->numEdges;
}

template <typename vertex> int weighted_graph<vertex>::total_weight() {
	// return the weight counter that is incremented/decremented as needed
	return this->totalWeight;
}

template <typename vertex>	std::vector<vertex> weighted_graph<vertex>::get_vertices() {
	// return the vector that vertices are added and removed from
	return this->vertices;
}

template <typename vertex>	std::vector<vertex> weighted_graph<vertex>::get_neighbours(const vertex& u) {
	// Create neighbours vector, loop through vertices list, 
	// check if other vertices are adjacent to u, if so add to neighbours and return
	std::vector<vertex> neighbours;
	for (int i = 0; i < vertices.size(); i++) {
		if (vertices[i] != u && are_adjacent(vertices[i], u)) {
			neighbours.push_back(vertices[i]);
		}
	}
	return neighbours;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::depth_first(const vertex& start_vertex){
	// traverse vertically (dft)
	bool visited[size];
	for (unsigned i = 0; i < size; i++){
		visited[i] = false;
	}
	// use stack (LIFO) for dft
	std::stack<vertex> unprocessed;
	unprocessed.push(start_vertex);
	
	std::vector<vertex> ordered;
	// as long as stack isn't empty, remove top element and check if visited and process
	// to ordered vector
	while (!unprocessed.empty()){
		int n = unprocessed.top();
		unprocessed.pop();
		if (!visited[n]){
			visited[n] = true;
			ordered.push_back(n);
			for (unsigned i = size; i != 0; i--){
				if (adj_matrix[n][i-1]){
					unprocessed.push(i-1);
				}
			}
		}
	}
	return ordered;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::breadth_first(const vertex& start_vertex){
	// traverse horizontally (bft)
	bool visited[size];
	for (unsigned i = 0; i < size; i++){
		visited[i] = false;
	}
	// use queue (FIFO) for bft
	std::queue<vertex> unprocessed;
	unprocessed.push(start_vertex);
	
	std::vector<vertex> ordered;
	// as long as queue isn't empty, remove first element and check if visited and process
	// to ordered vector
	while (!unprocessed.empty()){
		int n = unprocessed.front();
		unprocessed.pop();
		if (!visited[n]){
			visited[n] = true;
			ordered.push_back(n);
			for (unsigned i = 0; i < size; i++){
				if (adj_matrix[n][i]){
					unprocessed.push(i);
				}
			}
		}
	}
	return ordered;
}

template <typename vertex>	weighted_graph<vertex> weighted_graph<vertex>::mst() {
	// create tree
	weighted_graph<int> tree;
	//int minEdge;
	//add starting vertex to tree
	tree.add_vertex(vertices[0]);
	
	bool visited[size];
	for (unsigned i = 0; i < size; i++){
		visited[i] = false;
	}
	
	std::queue<std::pair<int,int> > unprocessed;
	//unprocessed.push({start, start});
	//while tree is less than vertices(numVertices) of graph
	/*while (tree.size < numVertices) {
		
	}*/
		// find the smallest edge between a tree vertex and a non-tree vertex
		// 
		// add the new vertex and the edge to the tree.
	return tree;
	//return weighted_graph<vertex>();
}


#endif
