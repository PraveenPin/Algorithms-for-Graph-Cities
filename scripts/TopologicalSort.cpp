#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/graph/edge_list.hpp>
#include <fstream>
#include <map>
#include <cassert>
#include <chrono>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/graph/bc_clustering.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <string>
#include <math.h>
using namespace boost;

typedef std::pair<int,int> Edge;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> Graph;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
typedef graph_traits<Graph>::vertex_descriptor Vertex;

int countOfVertices = 0;
int vertexCount = 1;
int firstTurn = 0;
void split(const std::string &s, char delim, std::vector<int> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(std::stoi(item));
    }
}
void getTopologicalSort(){
    char comma = ',';
    std::ifstream infile ("../Testcases/fpmeta.csv");
// input whole file, if edge of layer-5, add an directed edge from smaller vertex, to larger vertex
    std::string line;
    std::vector<Edge> edges;

    std::map<int, int> hashMap;
    std::map<int, int> revMap;
    std::map<Edge, int> weightOfEdge;
    int vertexCount = 0, verA = -1, verB = -1;
    while (std::getline(infile, line))
    {
        std::vector<int> row_values;
        split(line, comma, row_values);

        int a = row_values[0],b = row_values[1], weight = row_values[2];
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
        weightOfEdge[e] = weight;

        edges.push_back(e);

     }

    Graph metaGraph;
    
    for (Edge & edge : edges) {
        add_edge(edge.first, edge.second, metaGraph);
    }


    auto start1 = std::chrono::high_resolution_clock::now();

    std::cout << "Num vertices:" << num_vertices(metaGraph)<< std::endl;

    std::map<int,int> in_degree_map;
    std::map<int,int> zero_degree_map;
    std::vector<int> zero_degree_list;
    std::vector<std::vector<int>> orderedComponents;

    typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;

    for (vp = vertices(metaGraph); vp.first != vp.second; ++vp.first){
        int ver = *vp.first;
        int degree = in_degree(*vp.first, metaGraph);
        if(degree == 0){
            zero_degree_map[ver] = degree;  
            zero_degree_list.push_back(ver);
        }  
        else if(degree > 0){
            in_degree_map[ver] = degree;
        } 
    }

    while(!zero_degree_list.empty()){
        std::vector<int>  temp;
        orderedComponents.push_back(zero_degree_list);
        for (auto ver : zero_degree_list)
        {
            auto neighbours = boost::adjacent_vertices(ver, metaGraph);
            for (auto verU : make_iterator_range(neighbours)){
                in_degree_map[verU] -= 1;
                if(in_degree_map[verU] == 0){
                    temp.push_back(verU);
                }
            }
        }
        zero_degree_list = temp;
    }

    std::map<int,int> idx2level;
    std::vector<int> sizeOfEachLevel;
    for(int level = 0; level < orderedComponents.size(); level++){
        sizeOfEachLevel.push_back(orderedComponents[level].size());
        for(auto idx: orderedComponents[level]){
            idx2level[idx] = level;
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start1);
    std::cout << "Time taken by Toplogical sort and other preprocessing: "<< duration.count() << " seconds" << std::endl;
    

    std::ofstream myfile;
    std::string name( "../Ascii_files/fpmeta.span" );
    name+=".csv";
    std::cout<< "Creating a cc-fpmeta: "<<name<<std::endl;
    myfile.open (name);

    for (vp = vertices(metaGraph); vp.first != vp.second; ++vp.first){

        int verV = *vp.first;
        int vLevel = idx2level[verV];
        auto neighbours = boost::adjacent_vertices(*vp.first, metaGraph);
        for (int verU : make_iterator_range(neighbours)){
            int uLevel = idx2level[verU];
            if(uLevel != vLevel + 1){
                continue;
            }
            else{
                Edge e = Edge(verV, verU);
                int weight = weightOfEdge[e];
                myfile <<revMap[verV]<<","<<revMap[verU]<<","<<weight<<"\n";
            }
            
        }
    }
    std::cout<< "Finished creating a cc-fpmeta: "<<name<<std::endl;
    myfile.close();

}


int main()
{
    getTopologicalSort();
    return 1;

}

//this tool uses union find data structe to compute peel connected components with Meta Edges, meta nodes in single pass