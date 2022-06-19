//
// Created by lucas.futurist on 10/9/21.
//

#ifndef NNSPANNER_H
#define NNSPANNER_H

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

#endif //CLION_PROJECT_POINT_PARTITION_H

namespace spanners {

    using namespace std;

    typedef std::pair<long, long> longIntPair;

    class GridCell {
    public:
        vector<spanners::Point> pointsInsideTheCell;
        longIntPair neighbors[8] = {};
        Point leader;

        void addNeighbor(int index, longIntPair neighbor) {
            this->neighbors[index] = neighbor;
        }
    };

    const longIntPair standard = std::make_pair(-1,-1);

    typedef unordered_map<longIntPair, GridCell, boost::hash<longIntPair>> CellToPointsMap;
    typedef vector<vector<spanners::Point>> partitioned_points;
    typedef CGAL::Bbox_2 box;

    // takes a point set and partitions into 2D array of
    void point_partition(vector<Point> &P, CellToPointsMap &P_k, int size = 500) {

        // get the bounding box
        const size_t N = P.size();
        box bounding_box = CGAL::bbox_2(P.begin(), P.end());

        auto min_x = bounding_box.xmin();
        auto min_y = bounding_box.ymin();
        auto max_x = bounding_box.xmax();
        auto max_y = bounding_box.ymax();

        size_t numPartitions = ceil(N / size);
        bool squareFound = false;
        long sqrt_k = 0;
        long boxes = 1;

        while (!squareFound) {
            ++sqrt_k;
            boxes = sqrt_k * sqrt_k;
            squareFound = boxes >= numPartitions;
        }

        //  cout << "number of boxes: " << boxes << endl;
        //  cout << "number of desired partitions (lower bound): " << numPartitions << endl;
        cout << "height: " << sqrt_k << endl << endl;
        cout << "min_x: " << min_x << endl;
        cout << "max_x: " << max_x << endl << endl;
        cout << "min_y: " << min_y << endl;
        cout << "max_y: " << max_y << endl << endl;

        spanners::number_t total_y = (max_y - min_y);
        spanners::number_t total_x = (max_x - min_x);

        spanners::number_t length_y = total_y / sqrt_k;
        spanners::number_t length_x = total_x / sqrt_k;

        // adjustment translation constants
        const spanners::number_t adj_x = abs(min_x);
        const spanners::number_t adj_y = abs(min_y);

        bool partitionComplete = false;

        long dimX = 0;
        long dimY = 0;

        while (!partitionComplete) {

            P_k.clear();

            // cout << "test " << endl;

            size_t x_boxes = ceil(total_x / length_x);
            size_t y_boxes  = ceil(total_y / length_y);;
            boxes = x_boxes * y_boxes;

            for (auto v: P) {

                auto v_x = v.x() + adj_x;
                auto v_y = v.y() + adj_y;

                // determine x/y box coords
                long box_x = floor(v_x / length_x);
                long box_y = floor(v_y / length_y);

                dimX = (box_x > dimX) ? box_x : dimX;
                dimY = (box_y > dimY) ? box_y : dimY;

                bool localMax_x = abs(v_x - (box_x * length_x)) < EPSILON && v_x > 0;
                bool localMax_y = abs(v_y - (box_y * length_y)) < EPSILON && v_y > 0;

                box_x -= localMax_x;
                box_y -= localMax_y;

              /*  if (abs(v_x - (box_x * length_x)) < EPSILON)
                    --box_x;

                if (abs(v_y - (box_y * length_y)) < EPSILON)
                    --box_y;

                if (abs(v_x - min_x) < EPSILON)
                    ++box_x;

                if (abs(v_y - min_y) < EPSILON)
                    ++box_y; */

                // <x, y> --> <column, row> coordinate pair
                longIntPair position = std::make_pair(box_x, box_y);
                //Vertex ins{Point_2{v_x, v_y}, i++, new vector<Vertex>};
                P_k[position].pointsInsideTheCell.push_back(v);


                // account for boundary case
           /*     if (box_y == sqrt_k)
                    --box_y;
                if (box_x == sqrt_k)
                    --box_x; */
            }

            partitionComplete = true;

            for (auto cell: P_k) {

                if (cell.second.pointsInsideTheCell.size() > size) {

                    partitionComplete = false;
                    break;

                }
            }

            if (!partitionComplete) {

                length_x = 0.95 * length_x;
                length_y = 0.95 * length_y;

            }

        }

        bool occupiedTable[dimY+1][dimX+1];

        for (auto& cell : P_k) {
            //cout << "(" << cell.first.first << ", " << cell.first.second << ") " << endl;

            for (auto point : cell.second.pointsInsideTheCell) {
                //cout << " (" << point.x() << ", " << point.y() << ") ";
            }

            cout << endl;


            longIntPair position = cell.first;
            long x_pos = position.first;
            long y_pos = position.second;

            occupiedTable[y_pos][x_pos] = 1;

            auto points = cell.second.pointsInsideTheCell;
            double sumX = 0, sumY = 0;

            for (auto point : points) {
                sumX += point.x();
                sumY += point.y();
            }

            double meanX = sumX / points.size();
            double meanY = sumY / points.size();
            Point meanPoint(meanX, meanY);

            double minDistance = spanners::INF;

            for (auto point : points) {

                double distance = sqrt(CGAL::squared_distance(point, meanPoint));

                if (distance < minDistance) {
                    minDistance = distance;
                    cell.second.leader = point;
                }

            }

        }

        cout << endl << "occupiedTable size" << (dimY+1)*(dimX+1) << endl;
        for (int y = dimY; y >= 0; y--) {
            for (int x = 0; x < dimX+1; x++)
                cout << occupiedTable[y][x] << " ";
            cout << endl;
        }

        // identify neighbors
        for (auto& cell : P_k) {

            // identify position of cell
            long cellX = cell.first.first;
            long cellY = cell.first.second;

            longIntPair key = cell.first;

            int i[8] = {0, 1, 1, 1, 0, -1, -1, -1};
            int j[8] = {1, 1, 0, -1, -1, -1, 0, 1};

            // identify neighbors, beware the arrayIndexOutOfBound exception

            for (size_t cone = 0; cone < 8; cone++) {

                bool boundaryCase = (cellX == dimX && (cone == 0 || cone == 1 || cone == 7)) ||
                                    (cellX == 0 && (cone == 3 || cone == 4 || cone == 5)) ||
                                    (cellY == 0 && (cone == 5 || cone == 6 || cone == 7)) ||
                                    (cellY == dimY && (cone == 1 || cone == 2 || cone == 3));

                if (boundaryCase) {
                    //P_k[key].addNeighbor(cone, standard);
                    cell.second.addNeighbor(cone, standard);
                    continue;
                }

                long y_pos = cellY + i[cone];
                long x_pos = cellX + j[cone];
                longIntPair neighbor = (occupiedTable[y_pos][x_pos]) ? std::make_pair(x_pos, y_pos) : standard;

                cell.second.addNeighbor(cone, neighbor);
            }

        }
        cout << endl;

        for (auto cell : P_k) {

            cout << "(" << cell.first.first << ", " << cell.first.second << ")" << " --> neighbors" << endl << "    ";
            for (auto neighbor : cell.second.neighbors) {
                cout << "(" << neighbor.first << ", " << neighbor.second << ") ";
            }
            cout << endl;

        }

        cout << endl << "cardinality of each grid cell" << endl;

        long total = 0;
        long minCardinality = 500;
        long maxCardinality = 0;
        for (auto cell: P_k) {

            auto position = cell.first;
            auto cardinality = cell.second.pointsInsideTheCell.size();

            maxCardinality = (cardinality > maxCardinality) ? cardinality : maxCardinality;
            minCardinality = (cardinality < minCardinality) ? cardinality : minCardinality;

            cout << "(" << position.first << ", " << position.second << ") --> " << cardinality << endl;

            total += cardinality;

        }

        cout << endl << "occupied grids: " << P_k.size() << endl;
        cout << "unoccupied grids: " << (dimX+1)*(dimY+1) - P_k.size()  << endl;
        cout << "total grids: " << pow((sqrt(boxes)), 2) << " (" << sqrt(boxes) << " x " << sqrt(boxes)  << ") " << endl;
        cout << "% occupied: " << (double((P_k.size())) / double(pow((sqrt(boxes)), 2))) * 100 << "% " << endl;
        cout << "points: " << total << endl;
        cout << "minimum: " << minCardinality << endl;
        cout << "maximum: " << maxCardinality << endl;

    }

}