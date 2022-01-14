#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>
// #include <boost/graph/graph_traits.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/graph/edge_list.hpp>
#include <boost/graph/page_rank.hpp>
#include <fstream>
#include <map>
#include <cassert>
#include <chrono>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/bc_clustering.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <string>
#include <math.h>
#include <boost/pending/disjoint_sets.hpp>
using namespace boost;

typedef std::pair<int,int> Edge;
typedef std::vector<Edge> EdgeList;
typedef disjoint_sets_with_storage< identity_property_map, identity_property_map, find_with_full_path_compression > DisjointSet;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;
typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;
namespace qi = boost::spirit::qi;
double max_centrality = boost::lexical_cast<double>(10);

#define NFIELD 4
#define GRAPH_SIZE 4915350
#define MAXW 128

Graph g;

std::map<int, int> hashMap;
std::map<int, int> revMap;
std::vector<int> originalPeelValues;
std::map<Edge, int> edgePeelValues;
int countOfVertices = 0;
int vertexCount = 1;
int firstTurn = 0;


std::vector<Edge> readEdges(){

    int a, b,c;
    char comma = ',';

    std::vector<Edge> edges;

    // Get starting timepoint
    auto start1 = std::chrono::high_resolution_clock::now();

    std::fstream f("./layer-5.csv");
    if (!f.is_open()) {
        std::cerr << "error: file open failed " << ".\n";
        return edges;
    }

    int count = 0;
    int verA = -1, verB = -1;
    for (;;) {          /* loop continually */
        f >> a >> comma >> b >> comma >> c;
        if (f.fail() || f.eof())   
            break;
        count++;
        if(hashMap.find(a) == hashMap.end()){
            hashMap[a] = vertexCount;
            verA = vertexCount;
            if(revMap.find(vertexCount) == revMap.end()){
                revMap[vertexCount] = a;
            }
            vertexCount++;
            
        }
        else{
            verA = hashMap[a];
        }
        if(hashMap.find(b) == hashMap.end()){
            hashMap[b] = vertexCount;
            verB = vertexCount;
            if(revMap.find(vertexCount) == revMap.end()){
                revMap[vertexCount] = b;
            }

            vertexCount++;
        }
        else{
            verB = hashMap[b];
        }
        Edge e = Edge(verA,verB);
        edgePeelValues[e] = 0;
        edges.push_back(e);

        // add_edge(verA, verB, g);
        // break;
        if(count > GRAPH_SIZE){
            break;
        }
        f.ignore (MAXW, '\n');
    }

    f.close();

    // Get ending timepoint
    auto stop1 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::seconds>(stop1 - start1);
    std::cout << "Time taken to read csv file: "<< duration1.count() << " seconds" << std::endl;

    return edges;
}

std::vector<Edge> readEdgesFromBinary(){
    std::vector<Edge> edges;
    std::fstream file;
    int records = 0;
    std::cout << "Loading from disk?" << std::endl;
    file.open("sample.raw", std::ios::in | std::ios::binary);
    file.seekg(0, std::ios::end);
    int nodes = (int)file.tellg() / sizeof(int);

    std::cout << "total Nodes: "<<nodes << " Edges: "<< nodes/2 <<std::endl;
    file.seekg(0, std::ios::beg);
    int verV,verU;

    int count = 0;
    int verA = -1, verB = -1;
    for(int i = 0; i < nodes/2; i = i + 1){
      file.read( (char*)(&verV), 4 );
      file.read( (char*)(&verU), 4 );
    //   if(records < 10){
    //     std::cout << "First 10 edges" << records << "Edge :" <<verV << " - " << verU << std::endl;
    //   }
    //   if(records > (nodes/2) - 10){
    //     std::cout << "Last 10 edges" << records << "Edge :" <<verV << " - " << verU << std::endl;
    //   }
        count++;
        if(hashMap.find(verV) == hashMap.end()){
            hashMap[verV] = vertexCount;
            verA = vertexCount;
            if(revMap.find(vertexCount) == revMap.end()){
                revMap[vertexCount] = verV;
            }
            vertexCount++;
            
        }
        else{
            verA = hashMap[verV];
        }
        if(hashMap.find(verU) == hashMap.end()){
            hashMap[verU] = vertexCount;
            verB = vertexCount;
            if(revMap.find(vertexCount) == revMap.end()){
                revMap[vertexCount] = verU;
            }

            vertexCount++;
        }
        else{
            verB = hashMap[verU];
        }




      Edge e = Edge(verA,verB);
      edges.push_back(e);
      records++;
    }
    file.close();

    return edges;
}

std::vector<int> getPeelValues(std::vector<int> &peelValues){
    
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;

    countOfVertices = num_vertices(g);

    std::vector<int> bin(countOfVertices+1),position(countOfVertices+1),vert(countOfVertices+1), degree(countOfVertices+1);
    int maxDegree = 0;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first){
        int ver = int(index[*vp.first]);
        int deg = boost::degree(*vp.first, g)/2;
        degree[ver] = deg;
        
        if(maxDegree < deg){
            maxDegree = deg;
        }
    }

    for(int d=0;d<=maxDegree;d++){
        bin[d] = 0;
    }

    for(int i = 1; i<=countOfVertices; i++){
        bin[degree[i]]++;
    }

    int start = 1,num;
    for(int d = 0; d<= maxDegree;d++){
        num = bin[d];
        bin[d] = start;
        start+=num;
    }

    for(int v = 1; v <=countOfVertices;v++){
        position[v] = bin[degree[v]];
        vert[position[v]] = v;
        bin[degree[v]]++;
    }
    for(int d = maxDegree; d >0 ; d--){
        bin[d] = bin[d-1];
    }
    bin[0] = 1;


    // n disjoint sets
    DisjointSet ds(num_vertices(g));
    for(int i =1 ; i< countOfVertices; i++){
        ds.make_set(vert[i]);
    }

    for(int i =1; i<=countOfVertices;i++){
        int v = vert[i];
        int counter = 0;
        auto neighbours = boost::adjacent_vertices(v, g);

        for (int u : make_iterator_range(neighbours)){

            // union u and v - whole graph connected components
            ds.union_set(u,v);
            
            if(degree[u] > degree[v]){
                int du = degree[u],pu = position[u];
                int pw = bin[du], w = vert[pw];
                if(u != w){
                    position[u] = pw;
                    vert[pu] = w;
                    position[w] = pu;
                    vert[pw] = u;
                }
                bin[du]++;
                degree[u]--;
            }
        }
    }

    peelValues =  degree;

    std::vector<int> vertexCC(countOfVertices + 1);
    for(int i =1;i<=countOfVertices;i++){
        vertexCC[vert[i]] = ds.find_set(i);
    }

    // std::cout<<"Completed"<<std::endl;
    return vertexCC;
}


void writeToAsciiFile(std::vector<int>peelValues, std::vector<int>vertexCC){

    // writing to ascii file
    std::ofstream myfile;
    std::string name( "v_p_c_in_ascii_file.csv" );
    std::cout<< "Creating a Peel value: "<<name<<std::endl;
    myfile.open (name);
    for(int i =1;i<peelValues.size();i++){
        
        myfile <<revMap[i] <<","<<peelValues[i]<<","<<vertexCC[i]<<"\n";
    }
    std::cout<< "Finished creating a Peel value: "<<name<<std::endl;
}


void writeToBinaryFile(std::vector<int> peelValues,std::vector<int> vertexCC){

    std::fstream my_binary_file;;
    my_binary_file.open("v_p_c.raw", std::ios::out | std::ios::binary);
    std::cout<< "Writing v,p,c to binary file: "<<std::endl;
    int count = 0;
    int arr [3];
    for(int i =1;i<=peelValues.size();i++){
        int a = revMap[i], b = peelValues[i], c = vertexCC[i];
        my_binary_file.write( (char*)(&a), 4);
        my_binary_file.write( (char*)(&b), 4);
        my_binary_file.write( (char*)(&c), 4);
    }
    my_binary_file.close();
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Edge> edges = readEdgesFromBinary();
    int turn = 0;
    if(edges.size() == 0){
        std::cout<< "Only one edge in the file"<<std::endl;
        return 1;
    }
    
    for (Edge & edge : edges) {
        add_edge(edge.first, edge.second, g);
    }
    std::cout << "Graph is prepared"<< std::endl;
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by read from binary and prepare graph: "<< duration.count() << " seconds" << std::endl;
    
    
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<int> peelValues, vertexCC;
    //pass vector by reference and also compute meta edges

        
    vertexCC = getPeelValues(peelValues);

    writeToAsciiFile(peelValues,vertexCC);
    writeToBinaryFile(peelValues,vertexCC);
    // Get ending timepoint
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time taken by peel values: "<< duration.count() << " seconds" << std::endl;

    
}

// this tool computes the peel value of a graph