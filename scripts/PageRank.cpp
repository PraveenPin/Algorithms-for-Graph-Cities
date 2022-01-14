#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/graph/edge_list.hpp>
#include <boost/graph/page_rank.hpp>
#include <fstream>
#include <map>
#include <cassert>
#include <chrono>
#include <boost/config.hpp>
#include <boost/throw_exception.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;

typedef std::pair<int,int> Edge;
typedef std::vector<Edge> EdgeList;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> Graph;
typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iter;

#define NFIELD 4
#define GRAPH_SIZE 4915350
#define MAXW 128

int main()
{
    int a, b,c;
    char comma = ',';

    std::map<int, int> hashMap;
    std::map<int, int> revMap;


    // Get starting timepoint
    auto start1 = std::chrono::high_resolution_clock::now();

    std::fstream f("../Testcases/layer-5.csv");
    if (!f.is_open()) {
        std::cerr << "error: file open failed " << ".\n";
        return 1;
    }

    std::vector<Edge> edges;
    int count = 0;
    int vertexCount = 0;
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
        // std::cout <<"Old edges: "<<a<<" "<<b<<" New Edges: "<<verA<<" "<<verB<<std::endl;
        edges.push_back(Edge(verA,verB));
        // break;
        if(count > GRAPH_SIZE){
            break;
        }
        f.ignore (MAXW, '\n');
    }

    f.close();


    Graph g(GRAPH_SIZE);

    for(int i = 0; i< GRAPH_SIZE;i++){
        add_edge(edges[i].first, edges[i].second, g);
    }

    // Get ending timepoint
    auto stop1 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::seconds>(stop1 - start1);
    std::cout << "Time taken to read csv file: "<< duration1.count() << " seconds" <<" Vertex count: "<<vertexCount<< std::endl;

    




    // Get starting timepoint
    auto start2 = std::chrono::high_resolution_clock::now();
    
    //Page Rank stuff
    std::vector<double> ranks(num_vertices(g));
//try iterations = 10, compare with 20, estimate mean square error and check
    page_rank(g, make_iterator_property_map(ranks.begin(), get(boost::vertex_index, g)));

    // Get ending timepoint
    auto stop2 = std::chrono::high_resolution_clock::now();

    auto duration2 = std::chrono::duration_cast<std::chrono::seconds>(stop2 - start2);
    std::cout << "Time taken by pagerank: "<< duration2.count() << " seconds" << std::endl;

    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
    std::pair<vertex_iter, vertex_iter> vp;
    

    std::fstream my_binary_file;;
    my_binary_file.open("../binaries/vertex_page_rank.raw", std::ios::out | std::ios::binary);
    std::cout<< "Writing v,p to binary file: "<<std::endl;
    
    for (vp = vertices(g); vp.first != vp.second; ++vp.first){
        int ver = int(index[*vp.first]);
        
        if(revMap.find(ver)!= revMap.end()){
            int a = revMap[ver], b = ranks[*vp.first];
            my_binary_file.write( (char*)(&a), 4);
            my_binary_file.write( (char*)(&b), 4);
            // std::cout << "For old vertex : "<< revMap[ver] << " and new vertex: " << ver<< " Rank : "<< ranks[*vp.first] << std::endl;
        }
    }
    my_binary_file.close();
    std::cout<< "Completed writing page rank values to binary file"<<std::endl;

}