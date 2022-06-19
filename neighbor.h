//
// Created by shariful on 1/6/22.
//

#include <iostream>
#include <vector>
#include <chrono>
#include <queue>
#include <limits>
#include <cmath>
#include <random>
#include <unordered_set>
#include <boost/heap/fibonacci_heap.hpp> // ordering
#include <CGAL/Bbox_2.h>

#include "Utilities.h"

#ifndef SPANNEREXP_NEIGHBOR_H
#define SPANNEREXP_NEIGHBOR_H

#define validGrid(nx,ny,row,col) nx>=0 && nx<row && ny>=0 && ny<col
using namespace std;
using namespace spanners;

namespace neighbor{

    typedef std::pair<long, long> longIntPair;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2 Point;

//    class GridCell {
//    public:
//        vector<spanners::Point> pointsInsideTheCell;
//        longIntPair neighbors[8] = {};
//        Point leader;
//
//        void addNeighbor(int index, longIntPair neighbor) {
//            this->neighbors[index] = neighbor;
//        }
//    };

    //typedef unordered_map<longIntPair, GridCell, boost::hash<longIntPair>> CellToPointsMap;

    int currentX, currentY;
    int row = 0, column = 0;
    vector<pair<spanners::GridCell,spanners::GridCell>> findingNeighbourPair(spanners::CellToPointsMap partitionedPointSets, long row, long cols) {
        cout<<"In findingNeighbourPair "<<endl;
        vector<pair<spanners::GridCell,spanners::GridCell>> neighbourPair;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < cols; j++) {
                currentX = i, currentY = j;

                int upX = currentX + 1;
                int upY = currentY ;

                int rightX = currentX;
                int rightY = currentY + 1 ;

                int leftX = currentX;
                int leftY = currentY - 1;

                int rightUpDiagonalX = currentX + 1;
                int rightUpDiagonalY = currentY + 1;

                int leftUpDiagonalX = currentX - 1;
                int leftUpDiagonalY = currentY - 1;


                cout << "merging "<<currentX<<","<<currentY<<" "<<endl;

//                if (validGrid(upX, upY, row, column) and validGrid(rightX, rightY, row, column) and
//                    validGrid(leftX, leftY, row, column) and
//                    validGrid(leftUpDiagonalX, leftUpDiagonalY, row, column) and
//                    validGrid(rightUpDiagonalX, rightUpDiagonalY, row, column))
                {

                    if ( validGrid(upX, upY, row, cols) and validGrid(rightX, rightY, row, cols) and partitionedPointSets[make_pair(upX, upY)].pointsInsideTheCell.size() > 0 and
                        partitionedPointSets[make_pair(rightX, rightY)].pointsInsideTheCell.size() > 0) {
                        //cout<<"neighborPair "<<endl;
                        neighbourPair.emplace_back(make_pair(partitionedPointSets[make_pair(currentX, currentY)],
                                                             partitionedPointSets[make_pair(upX, upY)]));
                        neighbourPair.emplace_back(make_pair(partitionedPointSets[make_pair(currentX, currentY)],
                                                             partitionedPointSets[make_pair(rightX, rightY)]));
                    } else if ( validGrid(upX, upY, row, cols) and validGrid(rightX, rightY, row, cols) and validGrid(rightUpDiagonalX, rightUpDiagonalY, row, cols)
                    and partitionedPointSets[make_pair(upX, upY)].pointsInsideTheCell.size() == 0
                               and partitionedPointSets[make_pair(rightX, rightY)].pointsInsideTheCell.size() == 0
                               and partitionedPointSets[make_pair(rightUpDiagonalX,
                                                                  rightUpDiagonalY)].pointsInsideTheCell.size() > 0) {
                        neighbourPair.emplace_back(make_pair(partitionedPointSets[make_pair(currentX, currentY)],
                                                             partitionedPointSets[make_pair(rightUpDiagonalX,
                                                                                            rightUpDiagonalY)]));
                    } else if (validGrid(upX, upY, row, cols) and validGrid(leftUpDiagonalX, leftUpDiagonalY, row, cols) and validGrid(leftX, leftY, row, cols) and
                    partitionedPointSets[make_pair(upX, upY)].pointsInsideTheCell.size() == 0
                               and partitionedPointSets[make_pair(leftX, leftY)].pointsInsideTheCell.size() == 0
                               and partitionedPointSets[make_pair(leftUpDiagonalX,
                                                                  leftUpDiagonalY)].pointsInsideTheCell.size() > 0) {
                        neighbourPair.emplace_back(make_pair(partitionedPointSets[make_pair(currentX, currentY)],
                                                             partitionedPointSets[make_pair(leftUpDiagonalX,
                                                                                            leftUpDiagonalY)]));
                    }
                }
            }
        }
        return neighbourPair;
    }
}



#endif //SPANNEREXP_NEIGHBOR_H
