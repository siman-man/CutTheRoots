#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits.h>
#include <time.h>
#include <string>
#include <cassert>
#include <string.h>
#include <sstream>
#include <set>
#include <cstdio>
#include <cfloat>
#include <cstdlib>
#include <cmath>
#include <stack>
#include <queue>

using namespace std;

typedef long long ll;

const int MAX_NP = 105;
const int MAX_W = 1024;
const int MAX_H = 1024;

const int COUNTER_CLOCKWISE =     1;
const int CLOCKWISE         =    -1;
const int ONLINE_BACK       =     2;
const int ONLINE_FRONT      =    -2;
const int ON_SEGMENT        =     0;
const double EPS            = 1e-10;

unsigned long long xor128(){
  static unsigned long long rx=123456789, ry=362436069, rz=521288629, rw=88675123;
  unsigned long long rt = (rx ^ (rx<<11));
  rx=ry; ry=rz; rz=rw;
  return (rw=(rw^(rw>>19))^(rt^(rt>>8)));
}

class Vector {
  public:
    int id;
    int value;
    int dist;
    int depth;
    double x, y;
    vector<int> neighbor;

    Vector(double y = 0.0, double x = 0.0) {
      this->value = 0;
      this->dist = 0;
      this->depth = 0;
      this->y = y;
      this->x = x;
    }

    bool operator==(const Vector& v) const {
      return (x == v.x && y == v.y);
    }

    bool operator<(const Vector& v) const {
      return (x != v.x)? x < v.x : y < v.y;
    }

    Vector operator+(Vector p) {
      return Vector(y + p.y, x + p.x);
    }

    Vector operator-(Vector p) {
      return Vector(y - p.y, x - p.x);
    }

    Vector operator*(double k) {
      return Vector(k * y, k * x);
    }

    Vector operator/(double k) {
      return Vector(y / k, x / k);
    }

    double norm() {
      return x*x + y*y;
    }

    double abs() {
      return sqrt(norm());
    }
};

struct Circle {
  Vector center;
  double radius;
};

struct Point {
  Vector *p1;
  Vector *p2;
  Vector *p3;
  double dist;

  Point(Vector *p1, Vector *p2, Vector *p3, double dist){
    this->p1 = p1;
    this->p2 = p2;
    this->p3 = p3;
    this->dist = dist;
  }

  bool operator >(const Point &e) const{
    return dist > e.dist;
  }    
};

class Triangle{
  public:
    const Vector *p1, *p2, *p3;

    bool operator==(const Triangle& t) const{
      return (*p1 == *t.p1 && *p2 == *t.p2 && *p3 == *t.p3) ||
        (*p1 == *t.p1 && *p2 == *t.p3 && *p3 == *t.p2) ||
        (*p1 == *t.p2 && *p2 == *t.p3 && *p3 == *t.p1) ||
        (*p1 == *t.p2 && *p2 == *t.p1 && *p3 == *t.p3) ||
        (*p1 == *t.p3 && *p2 == *t.p1 && *p3 == *t.p2) ||
        (*p1 == *t.p3 && *p2 == *t.p2 && *p3 == *t.p1);
    }

    bool operator<(const Triangle& t) const{
      return !(*getMinVertex() == *t.getMinVertex())? 
        *getMinVertex() < *t.getMinVertex() :
        !(*getMidVertex() == *t.getMidVertex())?
        *getMidVertex() < *t.getMidVertex() :
        *getMaxVertex() < *t.getMaxVertex();
    }

    bool hasCommonPoints(const Triangle& t) const{
      return *p1 == *t.p1 || *p1 == *t.p2 || *p1 == *t.p3 ||
        *p2 == *t.p1 || *p2 == *t.p2 || *p2 == *t.p3 ||
        *p3 == *t.p1 || *p3 == *t.p2 || *p3 == *t.p3;
    }

  private:
    inline const Vector* getMinVertex() const{
      return *p1 < *p2 && *p1 < *p3 ? p1 : (*p2 < *p3)? p2 : p3;
    }

    inline const Vector* getMidVertex() const{
      return ((*p1 < *p2 && *p2 < *p3) || (*p3 < *p2 && *p2 < *p1))? p2 :
        ((*p2 < *p3 && *p3 < *p1) || (*p1 < *p3 && *p3 < *p2))? p3 : p1;
    }

    inline const Vector* getMaxVertex() const{
      return (*p2 < *p1 && *p3 < *p1)? p1 : (*p3 < *p2)? p2 : p3;
    }
};

Vector vectorList[MAX_NP];
vector<Vector> vertexList;

class Delaunay2d{
  public:
    typedef const set<Vector>               ConstVertexSet;
    typedef ConstVertexSet::const_iterator  ConstVertexIter;

    typedef set<Triangle>                   TriangleSet;
    typedef set<Triangle>::iterator         TriangleIter;

    typedef map<Triangle, bool>             TriangleMap;

    static void getDelaunayTriangles(ConstVertexSet &pVertexSet, TriangleSet *triangleSet){
      Triangle hugeTriangle;{
        double maxX, maxY; maxX = maxY = DBL_MIN;
        double minX, minY; minX = minY = DBL_MAX;

        for(ConstVertexIter it = pVertexSet.begin(); it != pVertexSet.end(); it++){
          double y = it->y;
          double x = it->x;

          maxX = max(maxX, x);
          minX = min(minX, x);

          maxY = max(maxY, y);
          minY = min(minY, y);
        }

        double centerX = (maxX - minX) * 0.5;
        double centerY = (maxY - minY) * 0.5;

        double dx = maxX - centerX;
        double dy = maxY - centerY;
        double radius = sqrt(dx*dx + dy*dy) + 1.0;

        Vector *p1 = &vectorList[MAX_NP-1];
        p1->x = centerX - sqrt(3.0) * radius;
        p1->y = centerY - radius;

        Vector *p2 = &vectorList[MAX_NP-2];
        p2->x = centerX + sqrt(3.0) * radius;
        p2->y = centerY - radius;

        Vector *p3 = &vectorList[MAX_NP-3];
        p3->x = centerX;
        p3->y = centerY + 2.0 * radius;

        hugeTriangle.p1 = p1;
        hugeTriangle.p2 = p2;
        hugeTriangle.p3 = p3;
      }

      triangleSet->insert(hugeTriangle);

      for(ConstVertexIter vIter = pVertexSet.begin(); vIter != pVertexSet.end(); vIter++){
        const Vector *p = &*vIter;

        TriangleMap rddcMap;

        for(TriangleIter tIter = triangleSet->begin(); tIter != triangleSet->end();){
          Triangle t = *tIter;

          Circle c;{
            double x1 = t.p1->x; double y1 = t.p1->y;
            double x2 = t.p2->x; double y2 = t.p2->y;
            double x3 = t.p3->x; double y3 = t.p3->y;

            double m = 2.0 * ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
            double x = ((y3-y1)*(x2*x2-x1*x1+y2*y2-y1*y1)
                +(y1-y2)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;
            double y = ((x1-x3)*(x2*x2-x1*x1+y2*y2-y1*y1)
                +(x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1))/m;

            c.center.x = x;
            c.center.y = y;

            double dx = t.p1->x - x;
            double dy = t.p1->y - y;
            double radius = sqrt(dx*dx + dy*dy);

            c.radius = radius;
          }

          double dx = c.center.x - p->x;
          double dy = c.center.y - p->y;
          double dist = sqrt(dx*dx + dy*dy);

          if(dist < c.radius){
            Triangle t1;
            t1.p1 = p; t1.p2 = t.p1; t1.p3 = t.p2;
            addElementToRedundanciesMap(&rddcMap, t1);

            Triangle t2;
            t2.p1 = p; t2.p2 = t.p2; t2.p3 = t.p3;
            addElementToRedundanciesMap(&rddcMap, t2);

            Triangle t3;
            t3.p1 = p; t3.p2 = t.p3; t3.p3 = t.p1;
            addElementToRedundanciesMap(&rddcMap, t3);

            triangleSet->erase(tIter++);
          }else{
            tIter++;
          }
        }

        for(TriangleMap::iterator iter = rddcMap.begin(); iter != rddcMap.end(); iter++){
          if(iter->second){
            triangleSet->insert(iter->first);
          }
        }
      }

      for(TriangleIter tIter = triangleSet->begin(); tIter != triangleSet->end();){
        if(hugeTriangle.hasCommonPoints(*tIter)){
          triangleSet->erase(tIter++);
        } else {
          tIter++;
        }
      }
    };

  private:
    static inline void addElementToRedundanciesMap(TriangleMap *pRddcMap, Triangle &t) {
      TriangleMap::iterator it = pRddcMap->find(t);

      if (it != pRddcMap->end() && it->second) {
        pRddcMap->erase(it);
        pRddcMap->insert(TriangleMap::value_type(t, false));
      } else {
        pRddcMap->insert(TriangleMap::value_type(t, true));
      }
    }
};

inline double dot(Vector &a, Vector &b){
  return a.x * b.x + a.y * b.y;
}

inline double cross(Vector &a, Vector &b){
  return a.x * b.y - a.y * b.x;
}

inline int ccw(Vector &p0, Vector &p1, Vector &p2){
  Vector a = p1 - p0;
  Vector b = p2 - p0;

  if(cross(a, b) > EPS) return COUNTER_CLOCKWISE;
  if(cross(a, b) < -EPS) return CLOCKWISE;
  if(dot(a, b) < -EPS) return ONLINE_BACK;
  if(a.norm() < b.norm()) return ONLINE_FRONT;

  return ON_SEGMENT;
}

inline bool intersect(Vector &p1, Vector &p2, Vector &p3, Vector &p4){
  return ((ccw(p1, p2, p3) * ccw(p1, p2, p4) <= 0) && (ccw(p3, p4, p1) * ccw(p3, p4, p2) < 0));
}

int randomNum[MAX_H+10];

struct Line {
  int fromY;
  int fromX;
  int toY;
  int toX;

  Line(int fromY = -1, int fromX = -1, int toY = -1, int toX = -1) {
    this->fromY = fromY;
    this->fromX = fromX;
    this->toY = toY;
    this->toX = toX;
  }
};

struct Edge {
  int fromY;
  int fromX;
  int toY;
  int toX;
  int removed;
  ll hashCode;

  Edge(int fromY = -1, int fromX = -1, int toY = -1, int toX = -1) {
    this->fromY = fromY;
    this->fromX = fromX;
    this->toY = toY;
    this->toX = toX;

    this->removed = 0;

    this->hashCode = randomNum[fromY] ^ randomNum[fromX] ^ randomNum[toY] ^ randomNum[toX];
  }
};

struct Root {
  int id;
  int from;
  int to;
  int depth;
  int value;
  int aid;
  int removed;

  Root(int from = -1, int to = -1) {
    this->from = from;
    this->to = to;

    this->depth = 0;
    this->aid = -1;
    this->value = 0;
    this->removed = 0;
  }
};

int g_NP;
vector<Edge> edgeList;
vector<Root> aliveRootList;
Root rootList[105000];
int rootListSize;
int g_aliveRootSize;
int g_edgeListSize;

class CutTheRoots {
  public:
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
      fprintf(stderr, "NP = %d\n", NP);
      fprintf(stderr, "points = %ld\n", points.size());
      fprintf(stderr, "roots = %ld\n", roots.size());

      vector<int> ret;
      init();

      g_NP = NP;
      int PS = points.size()/2;
      rootListSize = roots.size()/2;
      set<Vector> vertices;
      set<Triangle> triangles;

      for (int i = 0; i < PS; ++i) {
        Vector v;
        v.id = i;
        v.x = points[2*i];
        v.y = points[2*i+1];

        //fprintf(stderr,"fromY = %d, fromX = %d\n", (int)v.y, (int)v.x);
        vertexList.push_back(v);

        if(i < NP) {
          vertices.insert(v);
        }
      }

      int depthLimit = 3;

      if (NP <= 50) {
        depthLimit++;
      }
      if (NP <= 30) {
        depthLimit++;
      }
      if (NP <= 15) {
        depthLimit++;
      }
      if (NP <= 10) {
        depthLimit++;
      }

      for (int i = NP; i < PS; i++) {
        int j = roots[2*(i-NP)];
        int fromX = points[2*j];
        int fromY = points[2*j+1];

        int k = roots[2*(i-NP)+1];
        int toX = points[2*k];
        int toY = points[2*k+1];

        Vector *from = getVertex(j);
        Vector *to = getVertex(k);

        Root root(j, k);
        root.id = i;
        int dist = calcDistDetail(fromY, fromX, toY, toX);

        from->neighbor.push_back(i);

        to->dist = dist;
        to->depth = from->depth + 1;

        root.depth = to->depth;
        rootList[i] = root;

        if (to->depth <= depthLimit && dist > 5) {
          rootList[i].aid = g_aliveRootSize;
          g_aliveRootSize++;
          aliveRootList.push_back(root);
        }
      }

      for(int i = 0; i < NP; i++){
        searchRoot(i);
      }

      // ドロネー三角分割
      Delaunay2d::getDelaunayTriangles(vertices, &triangles);

      fprintf(stderr,"triangle num = %lu\n", triangles.size());

      set<Triangle>::iterator it = triangles.begin();

      map<ll, bool> checkList;

      while(it != triangles.end()){
        Triangle t = (*it);

        const Vector *v1 = t.p1;
        const Vector *v2 = t.p2;
        const Vector *v3 = t.p3;

        Edge edge1((int)v1->y, (int)v1->x, (int)v2->y, (int)v2->x);
        Edge edge2((int)v1->y, (int)v1->x, (int)v3->y, (int)v3->x);
        Edge edge3((int)v2->y, (int)v2->x, (int)v3->y, (int)v3->x);

        if (!checkList[edge1.hashCode]) {
          edgeList.push_back(edge1);
          checkList[edge1.hashCode] = true;
        }
        if (!checkList[edge2.hashCode]) {
          edgeList.push_back(edge2);
          checkList[edge2.hashCode] = true;
        }
        if (!checkList[edge3.hashCode]) {
          edgeList.push_back(edge3);
          checkList[edge3.hashCode] = true;
        }

        it++;
      }

      int esize = edgeList.size();
      int crossLimit = 3;

      for(int i = 0; i < NP-1; i++) {
        Vector v1 = vertexList[i];

        for(int j = i+1; j < NP; j++) {
          Vector v2 = vertexList[j];
          Edge newEdge(v1.y, v1.x, v2.y, v2.x);
          Vector p1(newEdge.fromY, newEdge.fromX);
          Vector p2(newEdge.toY, newEdge.toX);

          int crossCount = 0;

          for(int k = 0; k < esize && crossCount <= crossLimit; k++) {
            Edge *edge = getEdge(k);
            Vector p3(edge->fromY, edge->fromX);
            Vector p4(edge->toY, edge->toX);

            if (intersect(p1, p2, p3, p4)) {
              crossCount++;
            }
          }

          if (crossCount <= crossLimit) {
            edgeList.push_back(newEdge);
          }
        }
      }

      g_edgeListSize = edgeList.size();

      int rsize = g_aliveRootSize;
      fprintf(stderr, "edge size = %d, root size = %d\n", g_edgeListSize, rsize);
      vector<Line> lines;

      for(int i = 0; i < NP; i++) {
        Line line = getBestLine();

        if (line.fromY == -1) {
          continue;
        }

        //fprintf(stderr,"fromY = %d, fromX = %d, toY = %d, toX = %d\n", line.fromY, line.fromX, line.toY, line.toX);

        removeRoot(line);
        removeEdge(line);

        lines.push_back(line);
        lines = cleanLine(lines);
      }

      lines = cleanLine(lines);
      int lsize = lines.size();

      fprintf(stderr,"line size = %d\n", lsize);

      for(int i = 0; i < lsize; i++) {
        Line line = lines[i];

        ret.push_back(line.fromX);
        ret.push_back(line.fromY);
        ret.push_back(line.toX);
        ret.push_back(line.toY);
      }

      return ret;
    }

    Line getBestLine() {
      Line bestLine;
      int maxValue = INT_MIN;
      int limit = 100;

      if (g_NP <= 75) {
        limit = 200;
      }
      if (g_NP <= 60) {
        limit = 300;
      }
      if (g_NP <= 40) {
        limit = 400;
      }
      if (g_NP <= 25) {
        limit = 800;
      }
      if (g_NP <= 20) {
        limit = 1000;
      }

      for(int i = 0; i < limit; i++) {
        /*
        int idA = xor128() % g_edgeListSize;
        int idB = xor128() % g_edgeListSize;

        while(idA == idB) {
          idA = xor128() % g_edgeListSize;
          idB = xor128() % g_edgeListSize;
        }

        Edge *edgeA = getEdge(idA);
        Edge *edgeB = getEdge(idB);

        int cyA = (edgeA->fromY + edgeA->toY) / 2 + xor128()%300 - 150;
        int cxA = (edgeA->fromX + edgeA->toX) / 2 + xor128()%300 - 150;

        int cyB = (edgeB->fromY + edgeB->toY) / 2 + xor128()%300 - 150;
        int cxB = (edgeB->fromX + edgeB->toX) / 2 + xor128()%300 - 150;
        */

        int cyA = xor128()%MAX_H;
        int cxA = xor128()%MAX_W;

        int cyB = xor128()%MAX_H;
        int cxB = xor128()%MAX_W;

        Edge edgeC(cyA, cxA, cyB, cxB);
        Line line = edge2line(edgeC);

        int removeValue = 0;
        removeValue = removeRoot(line, true);
        int removeCount = removeEdge(line, true);
        int eval = 200 * removeCount - removeValue/10;

        if(removeCount > 0 && maxValue < eval) {
          maxValue = eval;
          bestLine = line;
        }
      }

      return bestLine;
    }

    vector<Line> cleanLine(vector<Line> lines) {
      vector<Line> result;
      int lsize = lines.size();
      int cleanCount = 0;

      for(int i = 0; i < lsize; i++) {
        Line line = lines[i];

        if (cleanEdge(line, true)) {
          cleanEdge(line);
          cleanCount++;
        } else {
          result.push_back(line);
        }
      }

      //fprintf(stderr,"clean line count = %d\n", cleanCount);

      return result;
    }

    int removeEdge(Line line, bool evalMode = false) {
      int removeCount = 0;
      Vector p1(line.fromY, line.fromX);
      Vector p2(line.toY, line.toX);

      for(int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (edge->removed > 0 && evalMode) {
          continue;
        }

        Vector p3(edge->fromY, edge->fromX);
        Vector p4(edge->toY, edge->toX);

        if (intersect(p1, p2, p3, p4)) {
          removeCount++;

          if (!evalMode) {
            edge->removed++;
          }
        }
      }

      return removeCount;
    }

    bool cleanEdge(Line line, bool evalMode = false) {
      Vector p1(line.fromY, line.fromX);
      Vector p2(line.toY, line.toX);

      for(int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (edge->removed >= 2 && evalMode) {
          continue;
        }

        Vector p3(edge->fromY, edge->fromX);
        Vector p4(edge->toY, edge->toX);

        if (intersect(p1, p2, p3, p4)) {
          if (edge->removed == 1 && evalMode) {
            return false;
          }

          if (!evalMode) {
            edge->removed--;
          }
        }
      }

      return true;
    }

    int removeRoot(Line line, bool evalMode = false) {
      int removeValue = 0;
      Vector p1(line.fromY, line.fromX);
      Vector p2(line.toY, line.toX);

      int rsize = g_aliveRootSize;

      for(int i = 0; i < rsize; i++) {
        Root *root = getAliveRoot(i);

        if (root->removed > 0 && evalMode) {
          continue;
        }

        Vector p3 = vertexList[root->from];
        Vector p4 = vertexList[root->to];

        if (intersect(p1, p2, p3, p4)) {
          removeValue += max(p3.value, p4.value);

          if (!evalMode) {
            root->removed++;
            cleanRoot(root->to);
          }
        }
      }

      return removeValue;
    }

    int searchRoot(int rootId) {
      Vector *v = getVertex(rootId);
      int nsize = v->neighbor.size();
      int value = 0;

      for(int i = 0; i < nsize; i++) {
        int rid = v->neighbor[i];
        Root *root = getRoot(rid);

        value += searchRoot(root->to);
      }

      value += v->dist;
      v->value = value;

      return value;
    }

    void cleanRoot(int rootId) {
      Vector *v = getVertex(rootId);
      int nsize = v->neighbor.size();

      for(int i = 0; i < nsize; i++) {
        int rid = v->neighbor[i];

        if (g_aliveRootSize <= rid) continue;
        Root *r = getRoot(rid);
        if (r->aid < 0) continue;
        Root *root = getAliveRoot(r->aid);
        root->removed++;

        cleanRoot(root->to);
      }
    }

    void init() {
      g_aliveRootSize = 0;

      for(int i = 0; i < MAX_H+10; i++) {
        randomNum[i] = xor128();
      }
    }

    Line edge2line(Edge edge) {
      Line line;
      int fromY, fromX, toY, toX;

      if (edge.fromX > edge.toX) {
        fromY = edge.toY;
        fromX = edge.toX;
        toY = edge.fromY;
        toX = edge.fromX;
      } else {
        fromY = edge.fromY;
        fromX = edge.fromX;
        toY = edge.toY;
        toX = edge.toX;
      }

      int dy = toY - fromY;
      int dx = toX - fromX;

      if (dx == 0) {
        line.fromX = fromX;
        line.fromY = 0;
        line.toX = toX;
        line.toY = MAX_H;
      } else if (dy == 0) {
        line.fromY = fromY;
        line.fromX = 0;
        line.toY = toY;
        line.toX = MAX_W;
      } else {
        double slope = dy / (double)dx;
        int b = (int)(fromY - slope * fromX);
        int c = (int)(slope * MAX_W + b);

        //fprintf(stderr,"b = %d, c = %d, slope = %4.2f\n", b, c, slope);

        if (b < 0) {
          line.fromY = 0;
          line.fromX = -1 * (int)(b / slope);

          if(c <= MAX_H) {
            line.toY = c;
            line.toX = MAX_W;
          } else {
            line.toY = MAX_H;
            line.toX = (int)((MAX_W - b) / slope);
          }
        } else if (0 <= b && b <= MAX_H) {
          line.fromY = b;
          line.fromX = 0;

          if (c < 0) {
            line.toY = 0;
            line.toX = -1 * (int)(b / slope);
          } else if (0 <= c && c <= MAX_H) {
            line.toY = c;
            line.toX = MAX_W;
          } else {
            line.toY = MAX_H;
            line.toX = (int)((MAX_W - b) / slope);
          }
        } else {
          line.fromY = MAX_H;
          line.fromX = (int)((MAX_W - b) / slope);

          if (c < 0) {
            line.toY = 0;
            line.toX = -1 * (int)(b / slope);
          } else {
            line.toY = c;
            line.toX = MAX_W;
          }
        }
      }

      /*
      assert(0 <= line.fromX && line.fromX <= MAX_W);
      assert(0 <= line.fromY && line.fromY <= MAX_H);
      assert(0 <= line.toX && line.toX <= MAX_W);
      assert(0 <= line.toY && line.toY <= MAX_H);
      */

      return line;
    }

    inline int calcDist(int y1, int x1, int y2, int x2) {
      int dy = y1 - y2;
      int dx = x1 - x2;
      return dy*dy + dx*dx;
    }

    inline int calcDistDetail(int y1, int x1, int y2, int x2) {
      int dy = y1 - y2;
      int dx = x1 - x2;
      return (int)sqrt(dy*dy + dx*dx);
    }

    Edge* getEdge(int id) {
      return &edgeList[id];
    }

    Root* getRoot(int id) {
      return &rootList[id];
    }

    Root* getAliveRoot(int id) {
      return &aliveRootList[id];
    }

    Vector* getVertex(int id) {
      return &vertexList[id];
    }
};

int main() {
  int NP; cin >> NP;
  int Npoints; cin >> Npoints;
  vector<int> points(Npoints);
  for(int i = 0; i < Npoints; i++){ cin >> points[i]; }
  int Nroots; cin >> Nroots;
  vector<int> roots(Nroots);
  for(int i = 0; i < Nroots; i++){ cin >> roots[i]; }
  CutTheRoots cr;
  vector<int> ret = cr.makeCuts(NP, points, roots);
  cout << ret.size() << endl;
  for (int i = 0; i < ret.size(); ++i) { cout << ret[i] << endl; }
  cout.flush();
}
