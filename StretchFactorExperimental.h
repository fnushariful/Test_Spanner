#ifndef STRETCHFACTOR_STRETCHFACTOREXPERIMENTAL_H
#define STRETCHFACTOR_STRETCHFACTOREXPERIMENTAL_H

#include <algorithm> // swap
#include <deque>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>


#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/squared_distance_2.h> //for 2D functions

#include <omp.h>

#include "Utilities.h"
#include "DelaunayL2.h"

namespace spanners {

    using namespace std;

    template<typename PointIterator, typename EdgeIterator>
    number_t StretchFactorExpDijk(PointIterator pointsBegin,
                                  PointIterator pointsEnd,
                                  EdgeIterator edgesBegin,
                                  EdgeIterator edgesEnd,
                                  const size_t numberOfThreads = 4) {
        vector<Point> P(pointsBegin, pointsEnd);
        const vector<Edge> E(edgesBegin, edgesEnd);

        const index_t n = P.size();
        cout<<"Point Size "<<n<<endl;
        vector<unordered_set<index_t>> G(n),
                DelG(n);


        for(const Edge& e : E) {

            G[e.first].insert(e.second);
            G[e.second].insert(e.first);
        }

        vector<Edge> edgesOfDT;
        vector<index_t> index;
        spatialSort<K>(P, index);

        //Step 1: Construct Delaunay triangulation
        DelaunayL2 DT; //DelaunayTriangulationSFH DT(P, edgesOfDT);
        /*Add IDs to the vertex handle. IDs are the number associated to the vertex, also maped as an index in handles.
          (i.e. Vertex with the ID of 10 will be in location [10] of handles.)*/
        FaceHandle hint;
        //cout<<"del:";
        for (size_t entry : index) {
            auto vh = DT.insert(P[entry], hint);
            hint = vh->face();
            vh->info() = entry;
        }
        VertexHandle v_inf = DT.infinite_vertex();

        // Convert DT to adjacency list
        for(auto vit=DT.finite_vertices_begin(); vit != DT.finite_vertices_end(); ++vit) {
            VertexCirculator N = DT.incident_vertices(vit),
                    done(N);
            do {
                if( N != v_inf )
                    DelG[vit->info()].insert(N->info());
            } while( ++N != done );
        }

        struct bfsElement {
            typedef index_t level_t;

            index_t vertex;
            level_t level;
            bfsElement(const index_t& vertex, const level_t& level )
             : vertex(vertex), level(level) {}
        };
        typedef std::deque<bfsElement> BfsQueue;

        typedef std::multimap<number_t, index_t> ShortestPathsQueue;
        typedef ShortestPathsQueue::iterator SPQHandle;


        vector<number_t> stretchFactorOfG(numberOfThreads,0.0);
        vector<pair<index_t, index_t>> worstPairOfG(numberOfThreads);
        vector<unordered_map<index_t,number_t>> tracker(n);

//#pragma omp parallel for num_threads(numberOfThreads)
        for( index_t u = 0; u < n; u++ ) {
            cout<<" Starting U is "<<u<<endl;

            // BFS variables /////////////////////////

            size_t lvl = 0;
            BfsQueue bfs;
            bfs.emplace_back(u,lvl);
            unordered_set<index_t> frontier;
            vector<bool> known(n,false);
            known[u] = true;

            // Dijkstra variables ////////////////////

            unordered_map<index_t, number_t> shortestPathLength(n);
            shortestPathLength[u] = 0.0;

            ShortestPathsQueue open;
            unordered_map<index_t, SPQHandle> openHandle(n);
            openHandle[u] = open.emplace(shortestPathLength[u], u);

            // Iterative Deepening Loop
            number_t t_u = 0.0, t_lvl = 0.0;

            do {
                ++lvl;
                t_u = t_lvl;
                cout<<"Level value "<<t_u<<endl;
                t_lvl = 0.0;

                while(!bfs.empty() && bfs.begin()->level < lvl) {
                    cout<<"BFS "<<bfs.empty()<<" BFS level begin "<<bfs.begin()->level<<" BFS lvl value "<<lvl<<endl;
                    auto v = bfs.begin()->vertex;
                    bfs.pop_front();
                    cout<<"current node V "<<v<<endl;
                    for( auto w : DelG[v] ) {
                        cout<<"W "<<w<<endl;
                        if(!known[w]) {
                            bfs.emplace_back(w,lvl);
                            known[w] = true;
                            if(!contains(shortestPathLength, w)) {
                                frontier.insert(w);
                                cout<<"frontier "<<endl;
                                copy(frontier.begin(),
                                     frontier.end(),
                                     std::ostream_iterator<int>(std::cout, " "));
                                cout<<"\n frontier "<<endl;
                            } else {
                                auto t_w = shortestPathLength[w] / getDistance(P[u], P[w]);

                                t_lvl = CGAL::max(t_lvl,t_w);
                                cout<<"t_w "<<t_w<<endl;
                            }
                        }
                    }

                }
                cout<<"BFS endl"<<endl;

                while(!frontier.empty()) {
                    assert(!open.empty());
                    auto nextShortestPath = open.begin();
                    cout<<nextShortestPath->first<<" nextShortestPath "<<nextShortestPath->second<<endl;
                    index_t v = nextShortestPath->second;
                    shortestPathLength[v] = nextShortestPath->first;
                    open.erase(nextShortestPath);
                    openHandle.erase(v);

                    if( contains(frontier,v) ) {
                        cout<<"frontiner erasing "<<v<<endl;
                        frontier.erase(v);
                        auto t_v = shortestPathLength[v] / getDistance(P[u], P[v]);
                        cout<<"U "<<u<<"frontiner "<<v<<"t_v "<<t_v<<endl;
                        t_lvl = CGAL::max(t_lvl, t_v);
                        cout<<"t_lvl "<<t_lvl<<endl;
                    }

                    for( auto w : G[v] ) {
                        cout<<v<<" adj In Disktra  "<<w<<endl;
                        if( !contains(shortestPathLength, w) ) {
                            number_t newDist = shortestPathLength[v] + getDistance(P[w], P[v]);
                            cout<<v<<" newDist "<<w<<" "<<newDist<<endl;

                            if(!contains(openHandle, w) || newDist < openHandle[w]->first ) {
                                if (contains(openHandle, w)) {
                                    cout<<"erasing "<<w<<endl;
                                    open.erase(openHandle[w]);
                                }
                                openHandle[w] = open.emplace(newDist, w);
                                cout<<"OpenHandle "<<openHandle[w]->first<<" "<<openHandle[w]->second<<endl;
                            }
                        }
                    }
                }
                cout<<"All end "<<endl;
            } while(!bfs.empty() && !(t_u > t_lvl)); ////

            // Done searching neighborhood of u

            cout<<"U level SF "<<u<<" "<<t_u<<endl;
            if (t_u > stretchFactorOfG[omp_get_thread_num()]) {
//                worstPairOfG[omp_get_thread_num()].first  = uWorstPair.first;
//                worstPairOfG[omp_get_thread_num()].second = uWorstPair.second;
                stretchFactorOfG[omp_get_thread_num()]    = t_u;
            }
        }

        //index_t worstPairU = worstPairOfG[0].first, worstPairV = worstPairOfG[0].second;
//        number_t finalSf = stretchFactorOfG[0];
//        for(index_t i = 1; i < numberOfThreads; i++)
//            if( stretchFactorOfG[i] >  finalSf ) {
//                finalSf = stretchFactorOfG[i];
//                worstPairU = worstPairOfG[i].first;
//                worstPairV = worstPairOfG[i].second;
//            }

        //cout << "Heuristic worst pair: " <<  worstPairU << ", " << worstPairV << endl;
        return *std::max_element(stretchFactorOfG.begin(),stretchFactorOfG.end());
    }

} // spanners

#endif //STRETCHFACTOR_STRETCHFACTOREXPERIMENTAL_H
