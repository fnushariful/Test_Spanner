//
// Created by justin on 12/28/21.
//

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Bbox_2.h>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <queue>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

#ifndef MERGESPANNER_WSPD_H
#define MERGESPANNER_WSPD_H
namespace WSPD {
    using namespace std;

    struct Pair{
        vector<Point> Sv, Sw;
    };

    class SplitTree {
    public:
        double boundaries[4];
        SplitTree *s1, *s2;
        vector<Point> pSet;
        CGAL::Bbox_2 boundingBox;
        SplitTree(const vector<Point> &points);
    };

    class WellSepPairDec{
    public:
        double s;
        SplitTree* root;
        vector<Pair> pairs;
        WellSepPairDec(double s, const vector<Point> &points);
        void InsertInternalPairs(SplitTree* root, queue<Pair> &q);
        void ExtractPairs(SplitTree* root, vector<Pair> &pairs);
        bool CalcS(double s, Pair in);
        void InsertNonWellSeparatedPairs(Pair cur, queue<Pair> &q);
    };

}
#endif //MERGESPANNER_WSPD_H
