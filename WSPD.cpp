//
// Created by justin on 12/28/21.
//

#include "WSPD.h"
#include <CGAL/Bbox_2.h>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <queue>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

using namespace WSPD;
using namespace std;

bool WellSepPairDec::CalcS(double s, Pair in) {
    double radiusV, radiusW;
    double orientation;
    double distance;
    CGAL::Bbox_2 boxV = CGAL::bbox_2(in.Sv.begin(), in.Sv.end());
    CGAL::Bbox_2 boxW = CGAL::bbox_2(in.Sw.begin(), in.Sw.end());
    double vx, vy, wx, wy;
    vx = boxV.xmax() - boxV.xmin();
    vy = boxV.ymax() - boxV.ymin();
    wx = boxW.xmax() - boxW.xmin();
    wy = boxW.ymax() - boxW.ymin();
    Point centerV{vx/2, vy/2};
    Point centerW{wx/2, wy/2};
}

WellSepPairDec::WellSepPairDec(double s, const vector<Point> &points){
    this->root = new SplitTree(points);

    this->s = s;
    ExtractPairs(this->root, this->pairs);
}

void WellSepPairDec::InsertInternalPairs(SplitTree *root, queue<Pair> &q) {
    SplitTree *iter = root;
    if(iter->s1 != nullptr && iter->s2 != nullptr){
        q.push(Pair{iter->s1->pSet, iter->s2->pSet});
        cout<<"Queue "<<q.size()<<endl;
        InsertInternalPairs(iter->s1, q);
        InsertInternalPairs(iter->s2, q);
    }else if(iter->s1 != nullptr && iter->s2 == nullptr){
        InsertInternalPairs(iter->s1, q);
    }else if(iter->s1 == nullptr && iter->s2 != nullptr){
        InsertInternalPairs(iter->s2, q);
    }else return;
}

void WellSepPairDec::ExtractPairs(SplitTree *root, vector<Pair> &pairs){
    queue<Pair> q;
    InsertInternalPairs(this->root, q);
    while(!q.empty()){
        Pair cur = q.front();
        q.pop();
        if(CalcS(this->s, cur)){
            this->pairs.emplace_back(cur);
        }else{
            InsertNonWellSeparatedPairs(cur, q);
        }
    }
}

void WellSepPairDec::InsertNonWellSeparatedPairs(Pair cur, queue<Pair> &q) {
    CGAL::Bbox_2 boxSv = CGAL::bbox_2(cur.Sv.begin(), cur.Sv.end());
    CGAL::Bbox_2 boxSw = CGAL::bbox_2(cur.Sw.begin(), cur.Sw.end());
    double dxSv, dySv, dxSw, dySw;
    dxSv = boxSv.xmax() - boxSv.xmin();
    dySv = boxSv.ymax() - boxSv.ymin();
    dxSw = boxSw.xmax() - boxSw.xmin();
    dySw = boxSw.ymax() - boxSw.ymin();
    vector<Point> *larger, *smaller;
    double *XL, *YL;
    double largestV, largestW;
    if(dxSv >= dySv){
        largestV = dxSv;
    }else{
        largestV = dySv;
    }
    if(dxSw >= dySw){
        largestW = dxSw;
    }else{
        largestW = dySw;
    }
    if(largestV >= largestW){
        larger = &cur.Sv;
        smaller = &cur.Sw;
        XL = &dxSv;
        YL = &dySv;
    }else{
        larger = &cur.Sw;
        smaller = &cur.Sv;
        XL = &dxSw;
        YL = &dySw;
    }
    double div;
    vector<Point> w1, w2;
    if(*XL >= *YL){
        div = *XL / 2;
        for(auto p : *larger){
            if(p.x() <= div){
                w1.emplace_back(p);
            }else{
                w2.emplace_back(p);
            }
        }
    }else{
        div = *YL / 2;
        for(auto p : *larger){
            if(p.y() <= div){
                w1.emplace_back(p);
            }else{
                w2.emplace_back(p);
            }
        }
    }
    q.push(Pair{*smaller, w1});
    q.push(Pair{*smaller, w2});
}

SplitTree::SplitTree(const vector<Point> &points){
    cout<<"Checked"<<endl;
    this->pSet = points;
    if(this->pSet.size() == 1){
        this->s1 = this->s2 = nullptr;
    }else{
        boundingBox = CGAL::bbox_2(points.begin(), points.end());
        boundaries[0] = boundingBox.xmin();
        boundaries[1] = boundingBox.xmax();
        boundaries[2] = boundingBox.ymin();
        boundaries[3] = boundingBox.ymax();
        double dy = boundingBox.ymax() - boundingBox.ymin();
        double dx = boundingBox.xmax() - boundingBox.xmin();
        vector<Point> s1Points, s2Points;
        if(dy >= dx){
            double mid = dy/2;
            for(auto p : this->pSet){
                if(p.y() >= mid){
                    s1Points.emplace_back(p);
                }else{
                    s2Points.emplace_back(p);
                }
            }
        }else{
            double mid = dx/2;
            for(auto p : this->pSet){
                if(p.x() >= mid){
                    s1Points.emplace_back(p);
                }else{
                    s2Points.emplace_back(p);
                }
            }
        }
        if(s1Points.empty()){
            this->s1 = nullptr;
        }else{
            this->s1 = new SplitTree(s1Points);
        }
        if(s2Points.empty()){
            this->s2 = nullptr;
        }else{
            this->s2 = new SplitTree(s2Points);
        }
    }
}

