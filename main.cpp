#include "GraphPrinter.h"
#include "PointGenerators.h"
#include "StretchFactorExact.h"
#include "StretchFactorExperimental.h"
#include "WSPD.h"
#include "point_partition.h"
#include "neighbor.h"

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <CGAL/Point_set_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/assertions.h>

using namespace std;
using namespace neighbor;

#define debug(x)         cerr<<__LINE__<<" "<<#x<<" "<<x<<endl;

// Traits
typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;

typedef K::FT                                                FT;
typedef K::Point_2                                           Point;
typedef K::Segment_2                                         Segment;
//typedef CGAL::Alpha_shape_vertex_base_2<K>                   Vb;
//typedef CGAL::Alpha_shape_face_base_2<K>                     Fb;
//typedef CGAL::Triangulation_data_structure_2<Vb,Fb>          Tds;
//typedef CGAL::Delaunay_triangulation_2<K,Tds>                Triangulation_2;
//typedef CGAL::Alpha_shape_2<Triangulation_2>                 Alpha_shape_2;
//typedef Alpha_shape_2::Alpha_shape_edges_iterator            Alpha_shape_edges_iterator;
//
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef pair<double, int> pPair;

typedef CGAL::Point_set_2<K>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K>::Vertex_handle  Vertex_handle;
typedef K::Point_2                           Point_2;
CGAL::Point_set_2<K> PSet;

vector<spanners::Edge> totalEdges;

typedef pair<string,string> Option;
typedef vector<Option> OptionsList;

vector<Point> adjacentList[1000];
map<pair<double, double>, int> cordinateID;
map<int,pair<double,double>> cordinateIDtoPoint;

//template <class OutputIterator>
//void alpha_edges( const Alpha_shape_2& A, OutputIterator out)
//{
//    Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
//            end = A.alpha_shape_edges_end();
//    for( ; it!=end; ++it)
//        *out++ = A.segment(*it);
//}

struct Vertex{
    Point_2 p;
    int id;
    vector<Vertex>* adj;
    bool operator== (Vertex const& b) const{
        return this->p.x() == b.p.x() && this->p.y() == b.p.y();
    }
};

struct Pair{
    Vertex u;
    Vertex v;
    double weight;
};



double eucDist(const Vertex& u, const Vertex& v){
    return sqrt(pow(u.p.x()-v.p.x(),2) + pow(u.p.y()-v.p.y(),2));
}

double eucdistanceCalculate(double x1, double y1, double x2, double y2)
{
    double x = x1 - x2; //calculating number to square in next step
    double y = y1 - y2;
    double dist;

    dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
    dist = sqrt(dist);

    return dist;
}

class comparatorPair{
public:
    bool operator() (const Pair &a, const Pair &b){
        return eucDist(a.u, a.v) > eucDist(b.u, b.v);
    }
};

void stretchFactor(const vector<Vertex>& points, vector<Pair>& s){
    Pair ins {points[0], points[1], eucDist(points[0],points[1])};
    s.emplace_back(ins);
    for(int i = 0; i < points.size(); ++i){
        for(int j = i + 1; j < points.size(); ++j){
            if(eucDist(points[i], points[j]) > s[0].weight){
                s.clear();
                ins = {points[i], points[j], eucDist(points[0],points[1])};
                s.emplace_back(ins);
            }else if(eucDist(points[i], points[j]) == s[0].weight){
                ins = {points[i], points[j], eucDist(points[0],points[1])};
                s.emplace_back(ins);
            }
        }
    }
}

void print(const vector<Vertex>& points, const vector<spanners::Edge> edge2){
    vector<Point_2> points2;
    //unordered_set<Point_2> val;

    for(const auto& x : points){
        spanners::Point ins{x.p.x(), x.p.y()};
        points2.emplace_back(ins);
    }
     double x= spanners::StretchFactorDijkstraReduction(points2.begin(),points2.end(),edge2.begin(),edge2.end());
    cout<<" SF "<<x<<endl;
    /*
    transform(points.begin(),points.end(),back_inserter(points2), [](const auto &value){
        return value.p;
    });
    */

    /*
    vector<spanners::Edge> edge2;
    for(const auto& x : edges){
        spanners::Edge ins{x.u.id, x.v.id};
        edge2.emplace_back(ins);
    }
     */


    /*
    transform(edges.begin(), edges.end(), edge2, [](const auto &value){
        return make_pair(value.u.p, value.v.p);
    });
    */

    spanners::GraphPrinter printer("./", "output","article");
    printer.autoscale(points2.begin(), points2.end());
    printer.drawEdges(edge2.begin(), edge2.end(), points2);

    //printer.drawVertices(points2.begin(), points2.end());
    printer.display();

    vector<Pair> sFactorPoints;
    stretchFactor(points, sFactorPoints);
}

void initPoints(vector<Point> &points, int numOfPts){
    double lowerBound = 0;
    double upperBound = 40;
    uniform_real_distribution<double> unif(lowerBound, upperBound);
    default_random_engine re;
    srand(time(NULL));
    int r = rand() % 1000;
    for(int i = 0; i < r; ++i){
        double toss = unif(re);
    }
    for(int i = 0; i < numOfPts; ++i){
        double x = unif(re);
        double y = unif(re);
        int a = rand() % 2;
        if(a == 0){
            x *= -1;
        }
        a = rand() % 2;
        if(a == 0){
            y *= -1;
        }
        bool exists = false;
        Point ins {x,y};
        for(auto z : points){
            if(z == ins) exists = true;
        }
        if(!exists){
            points.emplace_back(ins);
        }else{
            --i;
        }
    }
}

void initPoints(vector<Vertex> &verts, int numOfPoints){
    srand(time(0));
    for(int i = 0; i < numOfPoints; i++){
        double x = rand() % 20;
        if(x != 0){
            if(rand() % 2) x *= -1;
        }

        //x += static_cast<double>(rand())/100;

        double y = rand() % 20;
        if(y != 0){
            if(rand() % 2) y *= -1;
        }

        //y += static_cast<double>(rand())/100;

        Vertex a = {Point_2{x,y}};
        bool exists = false;
        for(const auto& x : verts){
            exists = (a==x);
        }
        if(!exists){
            Vertex ins {Point_2{x,y},i, new vector<Vertex>};
            verts.emplace_back(ins);
        }else{
            --i;
        }
    }
}

vector<int> adjList[5000];
int adjListSize;

double BFS(Vertex start, Vertex end, vector<Pair> edges)
{
    //cout<<"start "<<start.p.x()<<" "<<start.p.y()<<" end "<<end.p.x()<<" "<<end.p.y()<<endl;
//    for( auto x: edges )
//    {
//        cout<<cordinateID[make_pair(x.u.p.x(),x.u.p.y())]<<"-->";
//        for( auto y: *x.u.adj )
//        {
//            cout<<" "<<cordinateID[make_pair(y.p.x(),y.p.y())];
//        }
//        cout<<endl;
//    }
    //cout<<"start "<<start.p<<" "<<cordinateID[make_pair(start.p.x(),start.p.y())]<<" end "<<end.p<<" "<<cordinateID[make_pair(end.p.x(),end.p.y())]<<endl;
//        cout<<"adjList "<<endl;
//    for( int i=0; i<16; i++ )
//    {
//        cout<<i<<" ";
//        for( int j=0; j<adjList[i].size(); j++ )
//        {
//           cout<<" "<<adjList[i][j];
//        }
//        cout<<endl;
//    }
    double distance[5000];
    double parent[5000];
    int check[5000];

    memset(check,0,sizeof check);
    //memset(distance,0,sizeof distance);
    for( int i=0; i<5000; i++ ) distance[i] = 5000;

    queue<Vertex> bfsQueue;
    bfsQueue.push(start);

    distance[cordinateID[make_pair(start.p.x(),start.p.y())]] = 0;
    check[cordinateID[make_pair(start.p.x(),start.p.y())]] = 1;
    parent[cordinateID[make_pair(start.p.x(),start.p.y())]] = -1;

    while(!bfsQueue.empty())
    {
        const Vertex current = bfsQueue.front();
        bfsQueue.pop();
        //if(!check[current.id] ) {
        Vertex tmp;
            //cout<<"Current id "<<cordinateID[make_pair(current.p.x(),current.p.y())]<<" "<<current.p.x()<<" "<<current.p.y()<<endl;
                //cout<<" edges info "<<tmp.p<<" "<<tmp.id<<" "<<tmp.adj->size()<<endl;

                    for(int i=0; i<adjList[cordinateID[make_pair(current.p.x(),current.p.y())]].size(); i++ ) {
                        int next = adjList[cordinateID[make_pair(current.p.x(),current.p.y())]][i];
                        Vertex tmpPoint = {Point {cordinateIDtoPoint[next].first,cordinateIDtoPoint[next].second}};
                        //cout<<"Next "<<next<<" check "<<check[next]<<endl;
                        //if( !check[next] ) {
                        if(distance[next] > distance[cordinateID[make_pair(current.p.x(),current.p.y())]] + eucDist(current,tmpPoint)){
                        //cout<<"Next "<<next<<" "<<cordinateIDtoPoint[next].first<<" "<<cordinateIDtoPoint[next].second<<endl;
                        //Vertex tmpPoint = {Point {cordinateIDtoPoint[next].first,cordinateIDtoPoint[next].second}};
                        //cout<<"tmpPoint "<<tmpPoint.p.x()<<" "<<tmpPoint.p.y()<<endl;
                        //cout<<"Distance "<<eucDist(tmpPoint,current)<<endl;
                        //distance[cordinateID[make_pair(tmpPoint.p.x(),tmpPoint.p.y())]] = distance[cordinateID[make_pair(current.p.x(),current.p.y())]] + eucDist(tmpPoint,current);
                        distance[next] = distance[cordinateID[make_pair(current.p.x(),current.p.y())]] + eucDist(current,tmpPoint);
                        //cout<<"Next distance "<<distance[cordinateID[make_pair(tmpPoint.p.x(),tmpPoint.p.y())]]<<endl;
                        bfsQueue.emplace(tmpPoint);
                        check[next] = 1;
                        parent[next] = cordinateID[make_pair(current.p.x(),current.p.y())];
                        //cout<<"Check "<<check[next]<<endl;
                        }
                    }

//        if(cordinateID[make_pair(current.p.x(),current.p.y())] == cordinateID[make_pair(end.p.x(),end.p.y())])
//        {
//            //cout<<" start "<<cordinateID[make_pair(current.p.x(),current.p.y())]<<" end "<<cordinateID[make_pair(end.p.x(),end.p.y())]<<endl;
////            for( int i=0; i<14; i++ )
////                cout<<i<<" "<<parent[i]<<endl;
//            return distance[cordinateID[make_pair(end.p.x(),end.p.y())]];
//        }
        }
    return distance[cordinateID[make_pair(end.p.x(),end.p.y())]];
}

double shortestPath_FG(const vector<Vertex> &points, const Vertex &u, const Vertex &v, vector<vector<double>>& m){
    priority_queue<pPair, vector<pPair>, greater<>> p;
    vector<double> distances(points.size(), INT16_MAX);
    distances[u.id] = 0;
    p.push(make_pair(0, u.id));
    while(!p.empty()){
        int cur = p.top().second;
        p.pop();
        //cout<<"current "<<cur<<endl;

        for(const auto& x : *points[cur].adj){
            //cout<<"x id value "<<x.id<<" "<<points[cur].adj->size()<<endl;
            if(distances[x.id] > distances[cur] + eucDist(points[cur], points[x.id])){
                distances[x.id] = distances[cur] + eucDist(points[cur], points[x.id]);
                m[u.id][x.id] = distances[x.id];
                m[x.id][u.id] = distances[x.id];
                p.push(make_pair(distances[x.id], x.id));
            }
        }
    }
    return distances[v.id];
}


double shortestPath(const vector<Vertex> &points, vector<Pair> &edges, const Vertex &u, const Vertex &v){
    priority_queue<pPair, vector<pPair>, greater<>> p;
    vector<double> distances(points.size(), INT16_MAX);
    distances[u.id] = 0;
    p.push(make_pair(0, u.id));
   // cout<<" Points "<<points.size()<<endl;
   // cout<<" inside shortestPath "<<u.id<<endl;
    while(!p.empty()){
        int cur = p.top().second;
       // cout<<" current "<<cur<<endl;
        p.pop();
        for(const auto& x : *points[cur].adj){
            //cout<<"x id value "<<x.id<<" "<<points[cur].adj->size()<<endl;
            if(distances[x.id] > distances[cur] + eucDist(points[cur], points[x.id])){
                distances[x.id] = distances[cur] + eucDist(points[cur], points[x.id]);
                p.push(make_pair(distances[x.id], x.id));
            }
        }
    }
    return distances[v.id];
}


void FG_Greedy(const vector<Vertex>& points, const vector<Pair>& pairs, vector<Pair>& edges, vector<vector<double>>& m, double t){

    for(const auto& x : pairs){
        //cout<<"Pairs "<<pairs.size()<<" "<<x.v.id<<" "<<x.u.id<<endl;
        if(m[x.u.id][x.v.id] > t * eucDist(x.u, x.v)){
            shortestPath_FG(points, x.u, x.v, m);
        }
        if(m[x.u.id][x.v.id] > t * eucDist(x.u, x.v)){
            edges.emplace_back(x);
            x.u.adj->emplace_back(x.v);
            x.v.adj->emplace_back(x.u);
        }
    }
}

bool compareFunction(pair<double,double> a, pair<double,double> b)
{
    return a.second<b.second;
}

void BFSHeuristic(int totalPoint)
{
   for( int i=0; i<totalPoint; i++ )
   {
       vector<pair<double,double>> tmp;
       //cout<<"old adj "<<endl;
       for( int j=0; j<adjList[i].size(); j++ )
       {
           double sourcetoAdjPointDist = eucdistanceCalculate(cordinateIDtoPoint[i].first,cordinateIDtoPoint[i].second,cordinateIDtoPoint[adjList[i][j]].first,cordinateIDtoPoint[adjList[i][j]].second);
           tmp.emplace_back(make_pair(adjList[i][j],sourcetoAdjPointDist));
           //cout<<adjList[i][j]<<" ";
       }
       //cout<<endl;
       sort(tmp.begin(),tmp.end(),compareFunction);
       //cout<<"BFSHeuristic "<<endl;
//       for( auto x: tmp )
//       {
//           cout<<x.first<<" "<<x.second<<endl;
//       }
//       cout<<"NewAdj "<<endl;
//       for( int j=0; j<adjList[i].size(); j++ )
//       {
//           cout<<adjList[i][j]<<" ";
//          adjList[i][j] = tmp[j].first;
//          cout<<adjList[i][j]<<endl;
//       }
   }
}

void Org_Greedy(const vector<Vertex>& points, const vector<Pair>& pairs, vector<Pair>& edges, double t){
    //int cnt =0;
    //cout<<"Inside Org_Greedy "<<pairs.size()<<endl;
    for(const auto& x : pairs){
        //cout<<"Inside Org_Greedy shortestPath "<<x.u.p<<" "<<x.v.p<<" "<<x.u.id<<" "<<x.v.id<<endl;
        //cout<<shortestPath(points, edges, x.u, x.v)<<" shortestPath value"<<endl;
        if(shortestPath(points, edges, x.u, x.v) > t * eucDist(x.u,x.v)) {
            //cout<<"Inside Org_Greedy shortestPath inside "<<endl;
            edges.emplace_back(x);
            x.u.adj->emplace_back(x.v);
            x.v.adj->emplace_back(x.u);
            adjList[cordinateID[make_pair(x.u.p.x(),x.u.p.y())]].emplace_back(cordinateID[make_pair(x.v.p.x(),x.v.p.y())]);
            adjList[cordinateID[make_pair(x.v.p.x(),x.v.p.y())]].emplace_back(cordinateID[make_pair(x.u.p.x(),x.u.p.y())]);
            //adjListSize++;
        }
        //sort(adjList[cordinateID[make_pair(x.u.p.x(),x.u.p.y())]].begin(),adjList[cordinateID[make_pair(x.u.p.x(),x.u.p.y())]].end());
        //sort(adjList[cordinateID[make_pair(x.v.p.x(),x.v.p.y())]].begin(),adjList[cordinateID[make_pair(x.v.p.x(),x.v.p.y())]].end());
    }
//    cout<<"adjList "<<endl;
//    for( int i=0; i<20; i++ )
//    {
//        cout<<i<<" ";
//        for( int j=0; j<adjList[i].size(); j++ )
//        {
//           cout<<" "<<adjList[i][j];
//        }
//        cout<<endl;
//    }
}


void sortPairs(const vector<Vertex> &points, vector<Pair> &pairs){
    priority_queue<Pair, vector<Pair>, comparatorPair> p;
    for(int i = 0; i < points.size(); i++){
        for(int j = i + 1; j < points.size(); j++){
            //cout<<i<<" "<<j<<" "<<points[i].p<<" "<<points[j].p<<" "<<points[i].id<<" "<<points[j].id<<endl;
            Pair ins = {points[i], points[j]};
            ins.weight = eucDist(points[i], points[j]);
            //cout<<"weight "<<ins.weight<<endl;
            p.push(ins);
        }
    }
    //cout<<"pairs size "<<p.size()<<endl;
    int prioprityQueueSize = p.size();
    for(int i = 0; i <prioprityQueueSize; i++){
        pairs.emplace_back(p.top());
        p.pop();
    }
    //cout<<"Pairs size after pop "<<pairs.size()<<endl;
}

vector<spanners::Edge> shfitEdgeID(vector<spanners::Edge> edgeA, size_t cellSize)
{
    vector<spanners::Edge> shifitedEdge;

    for( auto x: edgeA )
    {
        x.first = x.first + cellSize;
        x.second = x.second + cellSize;
        spanners::Edge ins{x.first, x.second};
        shifitedEdge.emplace_back(ins);
    }
    cout<<"cellB size "<<edgeA.size()<<" cellA size "<<cellSize<<endl;
    for(auto x: shifitedEdge)
    {
        cout<<"After shifted "<<x.first<<" "<<x.second<<endl;
    }
    return shifitedEdge;
}


double StretchFactorCalculation(vector<Point> points, vector<spanners::Edge> edges)
{
    double sF = spanners::StretchFactorDijkstraReduction(points.begin(),points.end(),edges.begin(),edges.end());
    //cout<<sF<<endl;
    return sF;
}


int main() {

    int n = 1000;

    vector<Point> points;

    spanners::RandomPointGenerator_2 pg;
    //pg.generatePointsInsideASquare(n, 200, points);
    pg.generatePointsInsideADisc(n,5,points);

    spanners::CellToPointsMap cells;
    //spanners::CellToLeadersMap leadersMap;
    spanners::point_partition(points, cells);

    cout<<"AFter point Partition"<<endl;
    long rows, cols;
    rows = cols = -1;
    for(const auto &x : cells){
        if(x.first.first > rows) rows = x.first.first;
        if(x.first.second > cols) cols = x.first.second;
    }

    cout<<rows<<" "<<cols<<endl;

    vector<pair<spanners::GridCell,spanners::GridCell>> neighborPair = neighbor::findingNeighbourPair(cells,rows,cols);

    for( auto x: neighborPair )
    {
        cout<<x.first.pointsInsideTheCell.size()<<" "<<x.second.pointsInsideTheCell.size()<<endl;
    }


    clock_t st, ed;
    st = clock();
    //vector<spanners::Edge> totalEdges;

    spanners::GraphPrinter printer("./", "output", "article");
    Point_2 p;

    vector <Vertex> pointsnA;
    vector <Pair> edges;
    vector<Pair> alledgesAsPair;

    vector <Vertex> pointsB;
    double t = 1.25;
    //vector <Vertex> points, firstPointsSetV;
    vector <Pair> pairs;

    unsigned int pointSize = 200;
    //double s = 400;
    vector <Point> squarePoints, gridA, gridB, sPoints, firstPointsSet, testPoints;
    vector <spanners::Edge> edge2;
    spanners::RandomPointGenerator_2 gP;


    int i = 0;

    edges.clear();
    edge2.clear();
    pairs.clear();

    //gP.generatePointsInsideASquare(1000,1000,testPoints);


    //printer.drawVertices(testPoints.begin(),testPoints.end());
    //printer.display();

    gP.generatePointsInsideASquareNormal(pointSize, 1, squarePoints);
    gP.generatePointsInsideASquareNormal(pointSize, 1, sPoints);

    transform(sPoints.begin(), sPoints.end(), back_inserter(squarePoints), [](const Point &p) {
        return Point(p.x() - 10, p.y() - 10);
    });
    printer.autoscale(squarePoints.begin(), squarePoints.end());
    //cout<<"Hello "<<endl;
    i = 0;
    int cnt = 1;
    vector <Vertex> cellA, cellB, allpoints, totalpoint, totPoints;
    int id = 0;

    //WSPD::WellSepPairDec wellSepPairDec(3,squarePoints);

    for (auto x: squarePoints) {
        Vertex ins{Point_2{x.x(), x.y()}, i++, new vector <Vertex>};
        //cout<<"Points "<<x.x()<<" "<<x.y()<<endl;
        totPoints.emplace_back(ins);
    }
    edges.clear();
    edge2.clear();
    pairs.clear();
    //cout<<"size of all points "<<totPoints.size()<<endl;
    //vector<vector<double>> mm(totPoints.size(), std::vector<double>(totPoints.size(), INT16_MAX));
    //sortPairs(totPoints, pairs);
    //Org_Greedy(totPoints,pairs,edges,1.25);

    for (const auto &x: edges) {
        spanners::Edge ins{x.u.id, x.v.id};
        edge2.emplace_back(ins);
    }
    //print(totPoints,edge2);

    cnt = 0;
    i = 0;
    for (auto x: squarePoints) {
        if (cordinateID.find(make_pair(x.x(), x.y())) == cordinateID.end()) {
            //cout<<"id "<<id<<endl;
            cordinateID[make_pair(x.x(), x.y())] = id;
            cordinateIDtoPoint[id] = make_pair(x.x(), x.y());
            //cout<<"id "<<id<<endl;
            //cout<<"Cordinate to Point "<<cordinateIDtoPoint[id].first<<" "<<cordinateIDtoPoint[id].second<<endl;
            //cout<<"Id "<<cordinateID[make_pair(cordinateIDtoPoint[id].first, cordinateIDtoPoint[id].second)]<<endl;
            ++id;
        }
        if (cnt == pointSize) i = 0;
        Vertex ins{Point_2{x.x(), x.y()}, i++, new vector <Vertex>};
        //cout<<x.x()<<" "<<x.y()<<" "<<i<<endl;
        if (cnt >= 0 and cnt < pointSize) {
            cellA.emplace_back(ins);
            gridA.emplace_back(x);
            allpoints.emplace_back(ins);
            cnt++;
            //cout<<"Points A "<<x.x()<<" "<<x.y()<<" "<<i<<endl;
        } else if (cnt >= pointSize and cnt < pointSize+pointSize ) {
            cellB.emplace_back(ins);
            gridB.emplace_back(x);
            allpoints.emplace_back(ins);
            cnt++;
            //cout<<x.x()<<" setB "<<x.y()<<" "<<i<<endl;
            //cout<<"Points B "<<x.x()<<" "<<x.y()<<" "<<i<<endl;
        }
    }
//    for (auto x: cordinateID) {
//        cout << "[ " << x.first.first << " " << x.first.second << "] " << x.second << endl;
//    }
    i = 0;
    //vector<spanners::Edge> edge2;

    edges.clear();
    pairs.clear();
    edge2.clear();
    //cout<<"size of all points "<<allpoints.size()<<endl;
    //vector<vector<double>> nmm(allpoints.size(), std::vector<double>(allpoints.size(), INT16_MAX));
    //sortPairs(allpoints, pairs);
    //Org_Greedy(allpoints,pairs,edges,1.25);

    for (const auto &x: edges) {
        spanners::Edge ins{x.u.id, x.v.id};
        edge2.emplace_back(ins);
    }

    //double sF = spanners::StretchFactorDijkstraReduction(squarePoints.begin(),squarePoints.end(),edges.begin(),edges.end());
    //cout<<"Stretch Factor "<<sF<<endl;

    //print(allpoints,edge2);


    //cout<<"cell B size is "<<cellB.size()<<endl;

    PSet.insert(gridA.begin(), gridA.end());

    set <Point_2> borderpointA;

    for (auto x: gridB) {
        Vertex_handle v1 = PSet.nearest_neighbor({x.x(), x.y()});
        borderpointA.insert({v1->point().x(), v1->point().y()});
    }
    //cout<<" borderpointA "<<borderpointA.size()<<endl;
    PSet.clear();

    set <Point_2> borderpointB;

    PSet.insert(gridB.begin(), gridB.end());

    for (auto x: gridA) {
        Vertex_handle v1 = PSet.nearest_neighbor({x.x(), x.y()});
        //cout<<v1->point()<<endl;
        borderpointB.insert({v1->point().x(), v1->point().y()});
    }

//    for( auto x: borderpointA) { cout<<x.x()<<" "<<x.y()<<endl;
//    printer.drawVertexWithLabel(x.x(),x.y(),"A");
//    }
//
//    for( auto x: borderpointB ) { cout<<x.x()<<" "<<x.y()<<endl;
//        printer.drawVertexWithLabel(x.x(),x.y(),"B");
//    }

    //cout<<" boderpointB "<<borderpointB.size()<<endl;

    vector <Vertex> mixedpoint;

    i = 0;
    vector <Point_2> allborderpoint;

    for (auto x: borderpointA) {
        Vertex ins{Point_2{x.x(), x.y()}, i++, new vector <Vertex>};
        mixedpoint.emplace_back(ins);
        allborderpoint.emplace_back(Point_2{x.x(), x.y()});
    }

    for (auto x: borderpointB) {
        Vertex ins{Point_2{x.x(), x.y()}, i++, new vector <Vertex>};
        mixedpoint.emplace_back(ins);
        allborderpoint.emplace_back(Point_2{x.x(), x.y()});
    }
    i = 0;

    // vector<spanners::Edge> edge2;
    edge2.clear();

    edges.clear();
    pairs.clear();

    vector <vector<double>> mmmm(cellA.size(), std::vector<double>(cellA.size(), INT16_MAX));
    sortPairs(cellA, pairs);
    //cout<<"Before stitching "<<cellA.size()<<endl;
    //cout<<"Pairs size "<<pairs.size()<<endl;
    Org_Greedy(cellA, pairs, edges, t);
    //cout<<" CellA after FG_Greedy"<<edges.size()<<endl;

    for (const auto &x: edges) {
        spanners::Edge ins{x.u.id, x.v.id};
        spanners::Edge inss(cordinateID[make_pair(x.u.p.x(), x.u.p.y())], cordinateID[make_pair(x.v.p.x(), x.v.p.y())]);
        //cout<<" A edge map ID "<<cordinateID[make_pair(x.u.p.x(),x.u.p.y())]<<" "<<cordinateID[make_pair(x.v.p.x(),x.v.p.y())]<<endl;
        edge2.emplace_back(ins);
        totalEdges.emplace_back(inss);
        alledgesAsPair.emplace_back(x);
        // cout<<"edges "<<x.u.id<<" "<<x.v.id<<endl;
        //cout<<"Unique id "<<cordinateID[make_pair(x.u.p.x(),x.u.p.y())]<<" "<<cordinateID[make_pair(x.v.p.x(),x.v.p.y())]<<endl;
    }

    //printer.autoscale(gridA.begin(),gridA.end());
    printer.drawEdges(edge2.begin(), edge2.end(), gridA);



    //printer.drawVertices(allborderpoint.begin(),allborderpoint.end());
    //printer.display();

    edge2.clear();
    edges.clear();
    pairs.clear();
    vector <vector<double>> mmmmmm(cellB.size(), std::vector<double>(cellB.size(), INT16_MAX));
    sortPairs(cellB, pairs);

    Org_Greedy(cellB, pairs, edges, t);
    //cout<<" cellB edges size is "<<edges.size()<<endl;
    for (const auto &x: edges) {
        spanners::Edge ins{x.u.id, x.v.id};

        edge2.emplace_back(ins);
        spanners::Edge inss(cordinateID[make_pair(x.u.p.x(), x.u.p.y())], cordinateID[make_pair(x.v.p.x(), x.v.p.y())]);
        //cout<<" B edge map ID "<<cordinateID[make_pair(x.u.p.x(),x.u.p.y())]<<" "<<cordinateID[make_pair(x.v.p.x(),x.v.p.y())]<<endl;
        // cout<<"edges "<<x.u.id<<" "<<x.v.id<<endl;
        //cout<<"Unique ID "<<cordinateID[make_pair(x.u.p.x(),x.u.p.y())]<<" "<<cordinateID[make_pair(x.v.p.x(),x.v.p.y())]<<endl;
        totalEdges.emplace_back(inss);
        alledgesAsPair.emplace_back(x);
    }
    //print(cellB,edge2);
    //printer.autoscale(gridB.begin(),gridB.end());
    printer.drawEdges(edge2.begin(), edge2.end(), gridB);
    //printer.drawVertices(gridB.begin(),gridB.end());
    //printer.display();


    edges.clear();
    pairs.clear();
    edge2.clear();
    vector <vector<double>> m(mixedpoint.size(), std::vector<double>(mixedpoint.size(), INT16_MAX));
    sortPairs(mixedpoint, pairs);
    Org_Greedy(mixedpoint, pairs, edges, t);

    vector <spanners::Edge> edge22;
    for (const auto &x: edges) {
        spanners::Edge ins{x.u.id, x.v.id};
        //printer.drawVertexWithLabel(x.u.id,x.v.id,"A");
        spanners::Edge insss(cordinateID[make_pair(x.u.p.x(), x.u.p.y())],
                             cordinateID[make_pair(x.v.p.x(), x.v.p.y())]);
        //cout<<" border edge map ID "<<cordinateID[make_pair(x.u.p.x(),x.u.p.y())]<<" "<<cordinateID[make_pair(x.v.p.x(),x.v.p.y())]<<endl;

        spanners::Edge inss(cordinateID[make_pair(x.u.p.x(), x.u.p.y())], cordinateID[make_pair(x.v.p.x(), x.v.p.y())]);
        edge2.emplace_back(ins);
        edge22.emplace_back(insss);
        //cout<<"mixed point id edges "<<x.u.id<<" "<<x.v.id<<endl;
        //cout<<"mixed point id Unique ID "<<cordinateID[make_pair(x.u.p.x(),x.u.p.y())]<<" "<<cordinateID[make_pair(x.v.p.x(),x.v.p.y())]<<endl;
        totalEdges.emplace_back(insss);
        alledgesAsPair.emplace_back(x);
    }

    string activeEdgeColor =  "ff0000";
    double activeEdgeWidth = 0.5;

    OptionsList activeEdgeOptionsSF = { // active edge options
            {"color",      activeEdgeColor},
            {"line width", to_string(activeEdgeWidth)}
    };

    printer.drawEdges(edge2.begin(),edge2.end(),allborderpoint,activeEdgeOptionsSF);
    //printer.drawVerticesWithInfo(allborderpoint.begin(),allborderpoint.end());


    cout<<"Given SF "<<t<<endl;
    cout<<"stretch factor before "<<StretchFactorCalculation(squarePoints,totalEdges)<<endl;


    edge2.clear();
    vector <Point_2> extraStretchFactorEdges;
    vector<Vertex> StretchFactorVertex;
    BFSHeuristic(squarePoints.size());
    for( auto cellAPoint: cellA )
    {
        for( auto cellBPoint: cellB )
        {
            double shortestDistance = BFS(cellAPoint,cellBPoint,edges);
            double ecludian = eucDist(cellAPoint,cellBPoint);
            //cout<<"shortestDistance "<<shortestDistance<<" ecludian "<<ecludian<<endl;
            if( shortestDistance > t * ecludian )
            {
                edge2.emplace_back(cordinateID[make_pair(cellAPoint.p.x(), cellAPoint.p.y())],
                                   cordinateID[make_pair(cellBPoint.p.x(), cellBPoint.p.y())]);
                cout<<"edges id main "<<cordinateID[make_pair(cellAPoint.p.x(), cellAPoint.p.y())]<<" "<<
                        cordinateID[make_pair(cellBPoint.p.x(), cellBPoint.p.y())]<<endl;
                extraStretchFactorEdges.emplace_back(Point_2{cellAPoint.p.x(), cellAPoint.p.y()});
                extraStretchFactorEdges.emplace_back(Point_2{cellBPoint.p.x(), cellBPoint.p.y()});
                StretchFactorVertex.emplace_back(cellAPoint);
                StretchFactorVertex.emplace_back(cellBPoint);
                //cout<<"CellAPoint "<<cellAPoint.p.x()<<" "<<cellAPoint.p.y()<<endl;
                //cout<<"CellBPoint "<<cellBPoint.p.x()<<" "<<cellBPoint.p.y()<<endl;
                spanners::Edge ins{cordinateID[make_pair(cellAPoint.p.x(), cellAPoint.p.y())],
                        cordinateID[make_pair(cellBPoint.p.x(), cellBPoint.p.y())]};
                totalEdges.emplace_back(ins);

                adjList[cordinateID[make_pair(cellAPoint.p.x(),cellAPoint.p.y())]].emplace_back(cordinateID[make_pair(cellBPoint.p.x(),cellBPoint.p.y())]);
                adjList[cordinateID[make_pair(cellBPoint.p.x(),cellBPoint.p.y())]].emplace_back(cordinateID[make_pair(cellAPoint.p.x(),cellAPoint.p.y())]);
                //cout<<"Shortest Distance "<<shortestDistance<<" "<<ecludian<<endl;
            }
            //break;
        }
        //break;
    }

    activeEdgeColor =  "0000ff";
    activeEdgeWidth = 0.5;

    activeEdgeOptionsSF = { // active edge options
            {"color",      activeEdgeColor},
            {"line width", to_string(activeEdgeWidth)}
    };
    printer.setCordinateIDtoPointforGraphPrinter(cordinateIDtoPoint);
    printer.setCordinateIDforGraphPrinter(cordinateID);
    printer.drawEdgesforStretchFactor(edge2.begin(),edge2.end(),extraStretchFactorEdges,activeEdgeOptionsSF);
    //printer.drawVerticesWithID(extraStretchFactorEdges.begin(),extraStretchFactorEdges.end());
    //printer.drawVerticesWithID(squarePoints.begin(),squarePoints.end());

//    cout<<"adjList "<<endl;
//    for( int i=0; i<14; i++ )
//    {
//        cout<<i<<" ";
//        for( int j=0; j<adjList[i].size(); j++ )
//        {
//           cout<<" "<<adjList[i][j];
//        }
//        cout<<endl;
//    }

    cout<<"E/P "<<totalEdges.size()/squarePoints.size()<<endl;
    cout<<"stretch factor after BFS "<<StretchFactorCalculation(squarePoints,totalEdges)<<endl;

    printer.display();


    //double SFE = spanners::StretchFactorExpDijk(squarePoints.begin(), squarePoints.end(),totalEdges.begin(),totalEdges.end(),4);
    //cout<<SFE<<endl;

//    for( auto x: allborderpoint)
//    {
//        cout<<x.x()<<" "<<x.y()<<endl;
//    }
    //cout<<"stretch factor "<<StretchFactorCalculation(allborderpoint,edge22)<<endl;

   // cout << "stretch factor " << StretchFactorCalculation(squarePoints, totalEdges) << endl;
//    cout << "Stretch factor " << endl;
//    vector <spanners::Edge> tmpEdges;
//    tmpEdges = totalEdges;
//
//    int tmpEdgesLength = tmpEdges.size();
//    double mi = INT_MAX;
//    vector <spanners::Edge> largeStretchFactorEdges;
//    vector <Point> largeStretchFactorPoints;
//    double xTmpFirst, yTmpFirst, xTmpSecond, yTmpSecond;
//    for (int i = 0; i < tmpEdgesLength; i++) {
//        spanners::Edge tmpOne = tmpEdges.at(0);
//        cout << "Deleting index 0 edge " << tmpOne.first << " " << tmpOne.second << endl;
//        tmpEdges.erase(tmpEdges.begin() + 0);
//        double sF = StretchFactorCalculation(squarePoints, tmpEdges);
//        cout << "stretch factor " << sF << endl;
//        mi = min(mi, sF);
//        //if( sF <= mi )
//        tmpEdges.emplace_back(tmpOne);
//        if (sF > mi) {
//            largeStretchFactorEdges.emplace_back(tmpOne);
//            for (auto x: cordinateID) {
//                if (x.second == tmpOne.first) {
//                    xTmpFirst = x.first.first;
//                    yTmpFirst = x.first.second;
//                    largeStretchFactorPoints.emplace_back(Point_2{xTmpFirst, yTmpFirst});
//                } else if (x.second == tmpOne.second) {
//                    xTmpSecond = x.first.first;
//                    yTmpSecond = x.first.second;
//                    largeStretchFactorPoints.emplace_back(Point_2{xTmpSecond, yTmpSecond});
//                }
//            }
//        }
//    }
//    cout<<mi<<" minimum "<<endl;
//
//
//    //printer.display();
//    //cout<<"MixedPoint sizes "<<edges.size()<<endl;
//    //printer.autoscale(allborderpoint.begin(),allborderpoint.end());
//
//    typedef pair<string,string> Option;
//    typedef vector<Option> OptionsList;
//
//
//    string activeEdgeColor =     "0000ff";
//    double activeEdgeWidth = 0.5;
//
//    OptionsList activeEdgeOptions = { // active edge options
//            {"color",      activeEdgeColor},
//            {"line width", to_string(activeEdgeWidth)}
//    };
//
//    printer.drawEdges(edge2.begin(),edge2.end(),allborderpoint,activeEdgeOptions);
//    //printer.drawVertices(allborderpoint.begin(),allborderpoint.end());
//
//    activeEdgeColor =     "ff0000";
//    activeEdgeWidth = 0.5;
//
//    OptionsList activeEdgeOptionsSF = { // active edge options
//            {"color",      activeEdgeColor},
//            {"line width", to_string(activeEdgeWidth)}
//    };
//
//    printer.drawEdges(largeStretchFactorEdges.begin(),largeStretchFactorEdges.end(),largeStretchFactorPoints,activeEdgeOptionsSF);
//
//    printer.display();

//       string activeEdgeColor =  "ff0000";
//       double activeEdgeWidth = 0.5;
//
//       OptionsList activeEdgeOptionsSF = { // active edge options
//               {"color",      activeEdgeColor},
//               {"line width", to_string(activeEdgeWidth)}
//       };
//
//    printer.drawEdges(edge2.begin(),edge2.end(),allborderpoint,activeEdgeOptionsSF);
    //printer.display();
    //BFS for finding the bad pairs
//    for( auto x: cellA)
//    {
//        for( auto y: cellB)
//        {
//            cout<<x.p<<" "<<y.p<<endl;
//            double xShortDistanceY = BFS(x,y,alledgesAsPair);
//            double xEuclideanY = eucDist(x,y);
//            cout<<"short "<<xShortDistanceY<<" Euclide "<<xEuclideanY<<endl;
//        }
//    }

//     for( int i=0; i<squarePoints.size(); i++ )
//     {
//         cout<<i<<" ";
//        for( int j=0; j<adjList[i].size(); j++ )
//        {
//           cout<<adjList[i][j]<<" ";
//        }
//        cout<<endl;
//     }

    ed = clock();
    //print(allpoints,totalEdges);

    //printer.display();
    cout<<"time duration "<<(double) (ed-st)/CLOCKS_PER_SEC<<endl;

    return 0;

}
