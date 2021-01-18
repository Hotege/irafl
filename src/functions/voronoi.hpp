#if !defined(_FUNCTIONS_VORONOI_HPP_)
#define _FUNCTIONS_VORONOI_HPP_

#include <irafldef.h>
#include "functions.hpp"
#include <cstdint>
#include <ctime>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace IRAFLFunc
{
    namespace name
    {
        static std::string VoronoiName = "voronoi";
    }
    namespace voronoi
    {
        typedef std::pair<double, double> coord_t;
        typedef std::vector<coord_t> vtvec_t;
        typedef struct _VP
        {
            double r;
            uint32_t k;
        } VP;
        double RandDouble()
        {
            union
            {
                double_t d;
                uint64_t u64;
                uint16_t u16[sizeof(uint64_t) / sizeof(uint16_t)];
            } x;
            for (uint32_t i = 0; i < sizeof(uint64_t) / sizeof(uint16_t); ++i)
                x.u16[i] = (rand() << 1) | (rand() & 1);
            x.u64 = (x.u64 >> 12) | 0x3FF0000000000000ULL;
            return x.d - 1;
        }
        void IRAFL_LOCAL GetPoisson(vtvec_t& coords, const double& r, const uint32_t& k, const long& width, const long& height)
        {
            // according to https://blog.csdn.net/gouki04/article/details/100143747
            // step 0
            uint32_t n = 2;
            double cellSize = r / glm::sqrt(n);
            uint32_t cols = (uint32_t)glm::ceil(width / cellSize);
            uint32_t rows = (uint32_t)glm::ceil(height / cellSize);
            std::vector<glm::dvec2> cells;
            std::vector<std::vector<int>> grids(rows, std::vector<int>(cols, -1));
            // step 1
            auto v0 = glm::dvec2(RandDouble() * width, RandDouble() * height);
            uint32_t col = (uint32_t)glm::floor(v0.x / cellSize);
            uint32_t row = (uint32_t)glm::floor(v0.y / cellSize);
            int32_t v0Idx = cells.size();
            cells.push_back(v0);
            grids[row][col] = (int32_t)v0Idx;
            std::vector<int32_t> activeList;
            activeList.push_back(v0Idx);
            // step 2
            while (activeList.size() > 0)
            {
                auto rid = rand() % activeList.size();
                auto viIdx = activeList[rid];
                auto& vi = cells[viIdx];
                bool found = false;
                for (uint32_t i = 0; i < k; i++)
                {
TryFindNewVertex:
                    auto dir = [&]()->glm::dvec2 {
                        glm::dvec2 rd;
                        double r1 = r, r2 = 2 * r;
                        double a = 1 / (r2 * r2 - r1 * r1);
                        double rr = glm::sqrt(RandDouble() / a + r1 * r1);
                        double rt = RandDouble() * glm::two_pi<double>();
                        rd.x = rr * glm::cos(rt);
                        rd.y = rr * glm::sin(rt);
                        return rd;
                    }();
                    auto vk = vi + dir;
                    if (vk.x < 0 || vk.x >= width || vk.y < 0 || vk.y >= height)
                        goto TryFindNewVertex;
                    col = (uint32_t)glm::floor(vk.x / cellSize);
                    row = (uint32_t)glm::floor(vk.y / cellSize);
                    if (grids[row][col] != -1)
                        continue;
                    bool ok = true;
                    int32_t minR = (int32_t)glm::floor((vk.y - r) / cellSize);
                    int32_t maxR = (int32_t)glm::floor((vk.y + r) / cellSize);
                    int32_t minC = (int32_t)glm::floor((vk.x - r) / cellSize);
                    int32_t maxC = (int32_t)glm::floor((vk.x + r) / cellSize);
                    for (int32_t oR = minR; oR <= maxR; oR++)
                    {
                        if (oR < 0 || oR >= (int32_t)rows)
                            continue;
                        for (int32_t oC = minC; oC <= maxC; oC++)
                        {
                            if (oC < 0 || oC >= (int32_t)cols)
                                continue;
                            auto vjIdx = grids[oR][oC];
                            if (vjIdx != -1)
                            {
                                auto& vj = cells[vjIdx];
                                auto dist = glm::length(vj - vk);
                                if (dist < r)
                                {
                                    ok = false;
                                    goto EndOfDistanceCheck;
                                }
                            }
                        }
                    }
EndOfDistanceCheck:
                    if (ok)
                    {
                        auto vkIdx = cells.size();
                        cells.push_back(vk);
                        grids[row][col] = vkIdx;
                        activeList.push_back(vkIdx);
                        found = true;
                        break;
                    }
                }
                if (!found)
                    for (auto it = activeList.begin(); it != activeList.end();)
                        if ((*it) == viIdx)
                        {
                            it = activeList.erase(it);
                            break;
                        }
                        else
                            it++;
            }
            for (auto p : cells)
                coords.push_back({ p.x, p.y });
        }
        namespace Delaunay
        {
            typedef struct _cell_t
            {
                coord_t A;
                coord_t B;
                coord_t C;

                void Sort()
                {
                    if (A < B && B < C)
                        ;
                    else if (A < C && C < B)
                    {
                        coord_t t;
                        t.first = B.first;
                        t.second = B.second;
                        B.first = C.first;
                        B.second = C.second;
                        C.first = t.first;
                        C.second = t.second;
                    }
                    else if (B < A && A < C)
                    {
                        coord_t t;
                        t.first = A.first;
                        t.second = A.second;
                        A.first = B.first;
                        A.second = B.second;
                        B.first = t.first;
                        B.second = t.second;
                    }
                    else if (B < C && C < A)
                    {
                        coord_t t;
                        t.first = A.first;
                        t.second = A.second;
                        A.first = B.first;
                        A.second = B.second;
                        B.first = C.first;
                        B.second = C.second;
                        C.first = t.first;
                        C.second = t.second;
                    }
                    else if (C < A && A < B)
                    {
                        coord_t t;
                        t.first = C.first;
                        t.second = C.second;
                        C.first = B.first;
                        C.second = B.second;
                        B.first = A.first;
                        B.second = A.second;
                        A.first = t.first;
                        A.second = t.second;
                    }
                    else if (C < B && B < A)
                    {
                        coord_t t;
                        t.first = A.first;
                        t.second = A.second;
                        A.first = C.first;
                        A.second = C.second;
                        C.first = t.first;
                        C.second = t.second;
                    }
                }
                bool operator==(const _cell_t& c)
                {
                    return A == c.A && B == c.B && C == c.C;
                }
            } cell_t;
            typedef std::vector<cell_t> cvec_t;
            void IRAFL_LOCAL ConvexHull(vtvec_t*& res, const vtvec_t* _vts)
            {
                res = new vtvec_t;
                double minY = double(0x7FFFFFFF);
                auto minIt = _vts->end();
                for (auto it = _vts->begin(); it != _vts->end(); it++)
                    if ((*it).second < minY)
                    {
                        minY = (*it).second;
                        minIt = it;
                    }
                auto baseIt = minIt;
                auto searchIt = _vts->begin();
                if (minIt == _vts->end())
                    searchIt++;
                do
                {
                    double dot = 1.1;
                    for (auto it = _vts->begin(); it != _vts->end(); it++)
                        if (it != searchIt && it != baseIt)
                        {
                            glm::dvec3 base((*searchIt).first - (*baseIt).first, (*searchIt).second - (*baseIt).second, 0);
                            glm::dvec3 n((*it).first - (*baseIt).first, (*it).second - (*baseIt).second, 0);
                            auto c = glm::cross(base, n);
                            if (c.z <= 0)
                            {
                                searchIt = it;
                                dot = c.z;
                            }
                        }
                    baseIt = searchIt;
                    res->push_back(*searchIt);
                    for (auto it = _vts->begin(); it != _vts->end(); it++)
                    {
                        bool found = true;
                        for (auto p : *res)
                            if (p == *it)
                            {
                                found = false;
                                break;
                            }
                        if (found)
                        {
                            searchIt = it;
                            break;
                        }
                    }
                } while (baseIt != minIt);
            }
            bool IRAFL_LOCAL FindCommon(coord_t& va, coord_t& vb, coord_t& opposite1, coord_t& opposite2, const cell_t& t1, const cell_t& t2)
            {
                // AB - AB
                if ((t1.A == t2.A && t1.B == t2.B) || (t1.A == t2.B && t1.B == t2.A))
                {
                    va = t1.A;
                    vb = t1.B;
                    opposite1 = t1.C;
                    opposite2 = t2.C;
                    return true;
                }
                // AB - BC
                if ((t1.A == t2.B && t1.B == t2.C) || (t1.A == t2.C && t1.B == t2.B))
                {
                    va = t1.A;
                    vb = t1.B;
                    opposite1 = t1.C;
                    opposite2 = t2.A;
                    return true;
                }
                // AB - CA
                if ((t1.A == t2.C && t1.B == t2.A) || (t1.A == t2.A && t1.B == t2.C))
                {
                    va = t1.A;
                    vb = t1.B;
                    opposite1 = t1.C;
                    opposite2 = t2.B;
                    return true;
                }
                // BC - AB
                if ((t1.B == t2.A && t1.C == t2.B) || (t1.B == t2.B && t1.C == t2.A))
                {
                    va = t1.B;
                    vb = t1.C;
                    opposite1 = t1.A;
                    opposite2 = t2.C;
                    return true;
                }
                // BC - BC
                if ((t1.B == t2.B && t1.C == t2.C) || (t1.B == t2.C && t1.C == t2.B))
                {
                    va = t1.B;
                    vb = t1.C;
                    opposite1 = t1.A;
                    opposite2 = t2.A;
                    return true;
                }
                // BC - CA
                if ((t1.B == t2.C && t1.C == t2.A) || (t1.B == t2.A && t1.C == t2.C))
                {
                    va = t1.B;
                    vb = t1.C;
                    opposite1 = t1.A;
                    opposite2 = t2.B;
                    return true;
                }
                // CA - AB
                if ((t1.C == t2.A && t1.A == t2.B) || (t1.C == t2.B && t1.A == t2.A))
                {
                    va = t1.C;
                    vb = t1.A;
                    opposite1 = t1.B;
                    opposite2 = t2.C;
                    return true;
                }
                // CA - BC
                if ((t1.C == t2.B && t1.A == t2.C) || (t1.C == t2.C && t1.A == t2.B))
                {
                    va = t1.C;
                    vb = t1.A;
                    opposite1 = t1.B;
                    opposite2 = t2.A;
                    return true;
                }
                // CA - CA
                if ((t1.C == t2.C && t1.A == t2.A) || (t1.C == t2.A && t1.A == t2.C))
                {
                    va = t1.C;
                    vb = t1.A;
                    opposite1 = t1.B;
                    opposite2 = t2.B;
                    return true;
                }
                return false;
            }
            bool IRAFL_LOCAL IsVertexInCircumcircle(const coord_t& A, const coord_t& B, const coord_t& C, const coord_t& P)
            {
                double x1 = A.first, y1 = A.second;
                double x2 = B.first, y2 = B.second;
                double x3 = C.first, y3 = C.second;
                glm::dvec2 mid2 = { (x2 + x1) / 2, (y2 + y1) / 2 };
                glm::dvec2 mid3 = { (x3 + x1) / 2, (y3 + y1) / 2 };
                glm::dmat2 mm2 = {
                    { x2 - x1, -(y2 - y1) },
                    { mid2.y, mid2.x }
                };
                glm::dmat2 mm3 = {
                    { x3 - x1, -(y3 - y1) },
                    { mid3.y, mid3.x }
                };
                double detm2 = glm::determinant(mm2);
                double detm3 = glm::determinant(mm3);
                glm::dvec2 eqVec = { detm2, detm3 };
                glm::dmat2 eqMat = {
                    { x2 - x1, y2 - y1 },
                    { x3 - x1, y3 - y1 }
                };
                auto eqMatInv = glm::inverse(eqMat);
                auto center = glm::transpose(eqMatInv) * eqVec;
                auto radius = glm::distance(center, { x1, y1 });
                auto rd = glm::distance(center, { P.first, P.second });
                return rd < radius;
            }
            bool IRAFL_LOCAL IsVertexInCell(const cell_t& c, const coord_t& vt)
            {
                double x1 = c.A.first, y1 = c.A.second;
                double x2 = c.B.first, y2 = c.B.second;
                double x3 = c.C.first, y3 = c.C.second;
                glm::dvec2 t = { vt.first - x1, vt.second - y1 };
                glm::dvec2 b2 = { x2 - x1, y2 - y1 };
                glm::dvec2 b3 = { x3 - x1, y3 - y1 };
                glm::dmat2 eqMatAB = {
                    { b2.x, b3.x },
                    { b2.y, b3.y }
                };
                glm::dvec2 eqVecAB = { t.x, t.y };
                auto eqMatABInv = glm::inverse(eqMatAB);
                auto ab = glm::transpose(eqMatABInv) * eqVecAB;
                glm::dmat2 eqMatLG = {
                    { t.x, b2.x - b3.x },
                    { t.y, b2.y - b3.y }
                };
                glm::dvec2 eqVecLG = { b2.x, b2.y };
                auto eqMatLGInv = glm::inverse(eqMatLG);
                auto lg = glm::transpose(eqMatLGInv) * eqVecLG;
                return ab.x >= 0 && ab.y >= 0 && lg.x >= 1 && ab.x * ab.y > 0;
            }
            coord_t IRAFL_LOCAL FindCircumCircle(double& radius, const coord_t& A, const coord_t& B, const coord_t& C)
            {
                double x1 = A.first, y1 = A.second;
                double x2 = B.first, y2 = B.second;
                double x3 = C.first, y3 = C.second;
                glm::dvec2 mid2 = { (x2 + x1) / 2, (y2 + y1) / 2 };
                glm::dvec2 mid3 = { (x3 + x1) / 2, (y3 + y1) / 2 };
                glm::dmat2 mm2 = {
                    { x2 - x1, -(y2 - y1) },
                    { mid2.y, mid2.x }
                };
                glm::dmat2 mm3 = {
                    { x3 - x1, -(y3 - y1) },
                    { mid3.y, mid3.x }
                };
                double detm2 = glm::determinant(mm2);
                double detm3 = glm::determinant(mm3);
                glm::dvec2 eqVec = { detm2, detm3 };
                glm::dmat2 eqMat = {
                    { x2 - x1, y2 - y1 },
                    { x3 - x1, y3 - y1 }
                };
                auto eqMatInv = glm::inverse(eqMat);
                auto center = glm::transpose(eqMatInv) * eqVec;
                radius = glm::distance(center, { x1, y1 });
                return { center.x, center.y };
            }
            void IRAFL_LOCAL Triangulate(cvec_t& _cells, const vtvec_t& _coords, const long& width, const long& height)
            {
                const vtvec_t* coords = &_coords;
                cvec_t* cells = new cvec_t;
                vtvec_t* vtsConvexHull = nullptr;
                ConvexHull(vtsConvexHull, coords);
                vtvec_t rests;
                cvec_t standby;
                for (auto p : *coords)
                {
                    bool found = false;
                    for (auto q : *vtsConvexHull)
                        if (p == q)
                        {
                            found = true;
                            break;
                        }
                    if (!found)
                        rests.push_back(p);
                }
                bool noInnerVertices = false;
TAG_SPLIT:
                if (rests.size() > 0)
                {
                    auto itICHFinal = vtsConvexHull->end();
                    itICHFinal--;
                    for (auto it = vtsConvexHull->begin(); it != vtsConvexHull->end(); it++)
                    {
                        auto itNext = it;
                        itNext++;
                        coord_t v1 = *it;
                        coord_t v2;
                        if (it != itICHFinal)
                            v2 = *itNext;
                        else
                        {
                            if (noInnerVertices)
                                break;
                            else
                                v2 = *(vtsConvexHull->begin());
                        }
                        cvec_t acpt;
                        cell_t single = { *(rests.begin()), v1, v2 };
                        single.Sort();
                        standby.push_back(single);
                        while (!standby.empty())
                        {
                            auto ins = *(standby.begin());
                            standby.erase(standby.begin());
                            bool engage = true;
                            for (auto itc = cells->begin(); itc != cells->end(); itc++)
                            {
                                auto old = *itc;
                                if (ins == old)
                                {
                                    engage = false;
                                    break;
                                }
                                coord_t va, vb, o1, o2;
                                bool common = FindCommon(va, vb, o1, o2, ins, old);
                                if (common)
                                {
                                    bool inCircumcircle = IsVertexInCircumcircle(old.A, old.B, old.C, o1);
                                    if (inCircumcircle)
                                    {
                                        cell_t c1 = { o1, o2, va };
                                        c1.Sort();
                                        cell_t c2 = { o1, o2, vb };
                                        c2.Sort();
                                        standby.push_back(c1);
                                        standby.push_back(c2);
                                        cells->erase(itc);
                                        engage = false;
                                        break;
                                    }
                                }
                            }
                            if (engage)
                                acpt.push_back(ins);
                        }
                        while (!acpt.empty())
                        {
                            auto ins = *(acpt.begin());
                            acpt.erase(acpt.begin());
                            cells->push_back(ins);
                        }
                    }
                    rests.erase(rests.begin());
                    while (!rests.empty())
                    {
                        auto vt = *(rests.begin());
                        rests.erase(rests.begin());
                        for (auto itc = cells->begin(); itc != cells->end(); itc++)
                        {
                            auto c = *itc;
                            if (IsVertexInCell(c, vt))
                            {
                                cell_t cAB = { vt, c.A, c.B };
                                cAB.Sort();
                                cell_t cBC = { vt, c.B, c.C };
                                cBC.Sort();
                                cell_t cCA = { vt, c.C, c.A };
                                cCA.Sort();
                                standby.push_back(cAB);
                                standby.push_back(cBC);
                                standby.push_back(cCA);
                                cells->erase(itc);
                                break;
                            }
                        }
                        cvec_t acpt;
                        while (!standby.empty())
                        {
                            auto ins = *(standby.begin());
                            standby.erase(standby.begin());
                            bool engage = true;
                            for (auto itc = cells->begin(); itc != cells->end(); itc++)
                            {
                                auto old = *itc;
                                if (ins == old)
                                {
                                    engage = false;
                                    break;
                                }
                                coord_t va, vb, o1, o2;
                                bool common = FindCommon(va, vb, o1, o2, ins, old);
                                if (common)
                                {
                                    bool inCircumcircle = IsVertexInCircumcircle(old.A, old.B, old.C, o1);
                                    if (inCircumcircle)
                                    {
                                        cell_t c1 = { o1, o2, va };
                                        c1.Sort();
                                        cell_t c2 = { o1, o2, vb };
                                        c2.Sort();
                                        standby.push_back(c1);
                                        standby.push_back(c2);
                                        cells->erase(itc);
                                        engage = false;
                                        break;
                                    }
                                }
                            }
                            if (engage)
                                acpt.push_back(ins);
                        }
                        while (!acpt.empty())
                        {
                            auto ins = *(acpt.begin());
                            acpt.erase(acpt.begin());
                            cells->push_back(ins);
                        }
                    }
                }
                else
                {
                    rests.push_back(vtsConvexHull->back());
                    vtsConvexHull->pop_back();
                    noInnerVertices = true;
                    goto TAG_SPLIT;
                }
                for (auto it = cells->begin(); it != cells->end();)
                {
                    double radius = 0;
                    auto cnt = FindCircumCircle(radius, (*it).A, (*it).B, (*it).C);
                    if (cnt.first < 0 || cnt.first >= width || cnt.second < 0 || cnt.second >= height || radius <= 0)
                        it = cells->erase(it);
                    else
                        it++;
                }
                /*for (auto p : *cells)
                {
                    MKTriangle* tri = MKTriangle::Create();
                    for (MKVertex* q : *vertices)
                        if (q->x == p.A.first && q->y == p.A.second)
                        {
                            tri->A = q;
                            break;
                        }
                    for (MKVertex* q : *vertices)
                        if (q->x == p.B.first && q->y == p.B.second)
                        {
                            tri->B = q;
                            break;
                        }
                    for (MKVertex* q : *vertices)
                        if (q->x == p.C.first && q->y == p.C.second)
                        {
                            tri->C = q;
                            break;
                        }
                    tri->Sort();
                    bool found = false;
                    for (MKTriangle* q : *triangles)
                        if (*q == *tri)
                        {
                            found = true;
                            break;
                        }
                    if (!found)
                        triangles->push_back(tri);
                }*/
                _cells = *cells;
                delete vtsConvexHull;
                vtsConvexHull = nullptr;
                delete cells;
                cells = nullptr;
            }
        }
        void IRAFL_LOCAL voronoi(unsigned char*& dst, const unsigned char* src, const long& w, const long& h, const void* params)
        {
            VP* vp = (VP*)params;
            vtvec_t coords;
            GetPoisson(coords, vp->r, vp->k, w, h);
            Delaunay::cvec_t triangles;
            Delaunay::Triangulate(triangles, coords, w, h);
            memset(dst, -1, w * h * 4);
            bool* flags = new bool[w * h];
            memset(flags, 0, w * h);
            long* belongs = new long[w * h];
            memset(belongs, -1, sizeof(long) * w * h);
            std::vector<double> sums[3] = { std::vector<double>(coords.size(), 0), std::vector<double>(coords.size(), 0), std::vector<double>(coords.size(), 0) };
            long vpos = 0;
            for (auto p : coords)
            {
                vtvec_t neighbors;
                auto contain = [](const vtvec_t& vts, const coord_t& vt)->bool
                {
                    for (auto p : vts)
                        if (p == vt)
                            return true;
                    return false;
                };
                double xl = p.first, xr = p.first, yt = p.second, yb = p.second;
                for (auto q : triangles)
                    if (q.A == p)
                    {
                        xl = std::min(std::min(q.B.first, xl), q.C.first);
                        xr = std::max(std::max(q.B.first, xr), q.C.first);
                        yt = std::min(std::min(q.B.second, yt), q.C.second);
                        yb = std::max(std::max(q.B.second, yb), q.C.second);
                        if (!contain(neighbors, q.B))
                            neighbors.push_back(q.B);
                        if (!contain(neighbors, q.C))
                            neighbors.push_back(q.C);
                    }
                    else if (q.B == p)
                    {
                        xl = std::min(std::min(q.C.first, xl), q.A.first);
                        xr = std::max(std::max(q.C.first, xr), q.A.first);
                        yt = std::min(std::min(q.C.second, yt), q.A.second);
                        yb = std::max(std::max(q.C.second, yb), q.A.second);
                        if (!contain(neighbors, q.C))
                            neighbors.push_back(q.C);
                        if (!contain(neighbors, q.A))
                            neighbors.push_back(q.A);
                    }
                    else if (q.C == p)
                    {
                        xl = std::min(std::min(q.A.first, xl), q.B.first);
                        xr = std::max(std::max(q.A.first, xr), q.B.first);
                        yt = std::min(std::min(q.A.second, yt), q.B.second);
                        yb = std::max(std::max(q.A.second, yb), q.B.second);
                        if (!contain(neighbors, q.A))
                            neighbors.push_back(q.A);
                        if (!contain(neighbors, q.B))
                            neighbors.push_back(q.B);
                    }
                for (long y = ceil(yt); y <= floor(yb); y++)
                    for (long x = ceil(xl); x <= floor(xr); x++)
                    {
                        if (x == 554 && y == 107)
                        {
                            int njdkn = 0;
                        }
                        double d2p = (x - p.first) * (x - p.first) + (y - p.second) * (y - p.second);
                        bool match = true;
                        for (auto n : neighbors)
                            if ((x - n.first) * (x - n.first) + (y - n.second) * (y - n.second) < d2p)
                            {
                                match = false;
                                break;
                            }
                        if (match)
                        {
                            flags[y * w + x] = true;
                            belongs[y * w + x] = vpos;
                            sums[0][vpos] += src[(y * w + x) * 4 + 0];
                            sums[1][vpos] += src[(y * w + x) * 4 + 1];
                            sums[2][vpos] += src[(y * w + x) * 4 + 2];
                        }
                    }
                vpos++;
            }
            for (long y = 0; y < h; y++)
                for (long x = 0; x < w; x++)
                    if (!flags[y * w + x])
                    {
                        double d2 = 1e20;
                        long idx = 0, midx = 0;
                        for (auto p : coords)
                        {
                            double d2p = (x - p.first) * (x - p.first) + (y - p.second) * (y - p.second);
                            if (d2p < d2)
                            {
                                d2 = d2p;
                                midx = idx;
                            }
                            idx++;
                        }
                        belongs[y * w + x] = midx;
                        sums[0][midx] += src[(y * w + x) * 4 + 0];
                        sums[1][midx] += src[(y * w + x) * 4 + 1];
                        sums[2][midx] += src[(y * w + x) * 4 + 2];
                    }
            std::vector<int> count(coords.size(), 0);
            for (long y = 0; y < h; y++)
                for (long x = 0; x < w; x++)
                    count[belongs[y * w + x]]++;
            for (long y = 0; y < h; y++)
                for (long x = 0; x < w; x++)
                {
                    dst[(y * w + x) * 4 + 0] = sums[0][belongs[y * w + x]] / count[belongs[y * w + x]];
                    dst[(y * w + x) * 4 + 1] = sums[1][belongs[y * w + x]] / count[belongs[y * w + x]];
                    dst[(y * w + x) * 4 + 2] = sums[2][belongs[y * w + x]] / count[belongs[y * w + x]];
                }
            delete[]flags;
            flags = nullptr;
            delete[]belongs;
            belongs = nullptr;
        }
    }
}

#endif
