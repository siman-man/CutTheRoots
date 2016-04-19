#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <queue>

using namespace std;

typedef long long ll;

const int MAX_W = 1024;
const int MAX_H = 1024;
const int MAX_PS = 105110;

unsigned long long xor128() {
  static unsigned long long rx=123456789, ry=362436069, rz=521288629, rw=88675123;
  unsigned long long rt = (rx ^ (rx<<11));
  rx=ry; ry=rz; rz=rw;
  return (rw=(rw^(rw>>19))^(rt^(rt>>8)));
}

class Vector {
  public:
    int id;
    int x, y;
    double value;
    double dist;
    int depth;
    vector<int> roots;

    Vector(int y = -1, int x = -1) {
      this->value = 0.0;
      this->dist = 0.0;
      this->depth = 0;
      this->y = y;
      this->x = x;
    }
};

vector<Vector> g_vertexList;

inline bool intersect(ll y1, ll x1, ll y2, ll x2, ll y3, ll x3, ll y4, ll x4) {
  return (((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*((x2-x1)*(y4-y1)-(y2-y1)*(x4-x1)) <= 0);
}

inline bool ontheline(int ay, int ax, int by, int bx, int cy, int cx) {
  return ((bx-ax)*(cy-ay)-(by-ay)*(cx-ax) == 0);
}

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
  double length;
  int removed;

  Edge(int fromY = -1, int fromX = -1, int toY = -1, int toX = -1) {
    this->fromY = fromY;
    this->fromX = fromX;
    this->toY = toY;
    this->toX = toX;
    this->removed = 0;
  }
};

struct Root {
  int id;
  int from;
  int to;
  int depth;
  double value;
  double length;
  int removed;

  Root(int from = -1, int to = -1) {
    this->from = from;
    this->to = to;

    this->depth = 0;
    this->length = 0.0;
    this->value = 0.0;
    this->removed = 0;
  }
};

struct Node {
  Edge edge;
  double eval;

  Node(Edge &edge, double value){ 
    this->edge = edge;
    this->eval = value;
  }
    
  bool operator >(const Node &e) const {
    return eval < e.eval;
  }
};

int g_NP;
int g_PS;
int g_RC;
int g_depthLimit;
int g_branchBonus;
double g_refineRate;
double g_minValue;
int g_tryLimit;
int g_cutLimit;
int g_rootListSize;
int g_ARC;
int g_edgeListSize;
int g_updated[MAX_PS];
int g_time;
Line g_line;

vector<Edge> g_edgeList;
vector<int> g_activeRootList;
Root rootList[MAX_PS];

priority_queue<Node, vector<Node>, greater<Node> > pque;

class CutTheRoots {
  public:
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
      vector<int> ret;
      g_PS = points.size()/2;
      g_RC = roots.size()/2;

      fprintf(stderr, "NP = %d\n", NP);
      fprintf(stderr, "points = %d\n", g_PS);
      fprintf(stderr, "roots = %d\n", g_RC);

      init(NP);

      for (int i = 0; i < g_PS; i++) {
        Vector v;
        v.id = i;
        v.x = points[2*i];
        v.y = points[2*i+1];

        g_vertexList.push_back(v);
      }

      for (int i = NP; i < g_PS; i++) {
        int j = roots[2*(i-NP)];
        int fromX = points[2*j];
        int fromY = points[2*j+1];

        int k = roots[2*(i-NP)+1];
        int toX = points[2*k];
        int toY = points[2*k+1];

        Vector *from = getVertex(j);
        Vector *to = getVertex(k);

        double dist = calcDistDetail(fromY, fromX, toY, toX);

        Root root(j, k);
        root.id = i;
        root.length = dist;

        from->roots.push_back(root.id);

        to->dist = dist;
        to->depth = from->depth + 1;

        root.depth = to->depth;
        rootList[g_rootListSize++] = root;

        if (to->depth <= g_depthLimit && dist > g_cutLimit) {
          g_activeRootList.push_back(root.id);
        }
      }

      g_ARC = g_activeRootList.size();

      directTryLimit();

      for (int i = 0; i < NP-1; i++) {
        Vector *v1 = getVertex(i);

        for (int j = i+1; j < NP; j++) {
          Vector *v2 = getVertex(j);
          double length = calcDistDetail(v1->y, v1->x, v2->y, v2->x);

          if (length >= 1000) continue;

          Edge newEdge(v1->y, v1->x, v2->y, v2->x);
          newEdge.length = length;

          g_edgeList.push_back(newEdge);
        }
      }

      g_edgeListSize = g_edgeList.size();

      int rsize = g_ARC;
      fprintf(stderr, "edge size = %d, root size = %d\n", g_edgeListSize, rsize);

      while (!finishCheck()) {
        refreshJar();
        Edge edge = getBestEdge();
        edge = refineEdge(edge);

        updateLine(edge.fromY, edge.fromX, edge.toY, edge.toX);
        double removeValue = removeRoot(g_line);
        int removeCount = removeEdge(g_line);

        pque.push(Node(edge, removeValue / (double)removeCount));
      }

      vector<Edge> edges = cleanEdges();

      for (int i = 0; i < edges.size(); i++) {
        Edge edge = edges[i];

        ret.push_back(edge.fromX);
        ret.push_back(edge.fromY);
        ret.push_back(edge.toX);
        ret.push_back(edge.toY);
      }

      return ret;
    }

    bool finishCheck() {
      for (int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (edge->removed == 0) return false;
      }

      return true;
    }

    Edge getBestEdge() {
      Edge bestEdge;
      g_minValue = DBL_MAX;
      int y1, x1, y2, x2;

      for (int i = 0; i < g_tryLimit; i++) {
        y1 = xor128()%MAX_H;
        x1 = xor128()%MAX_W;

        y2 = xor128()%MAX_H;
        x2 = xor128()%MAX_W;

        updateLine(y1, x1, y2, x2);
        if (lineOnThePlant()) continue;

        int removeEdgeCount = removeEdgeEval(g_line);

        if (removeEdgeCount > 0) {
          double removeValue = 0.0;
          removeValue = removeRootEval(g_line, (double) removeEdgeCount);
          double eval = removeValue / (double) removeEdgeCount;

          if (g_minValue > eval) {
            g_minValue = eval;
            bestEdge.fromY = y1;
            bestEdge.fromX = x1;
            bestEdge.toY = y2;
            bestEdge.toX = x2;
          }
        }
      }

      return bestEdge;
    }

    Edge refineEdge(Edge &currentEdge) {
      Edge bestEdge = currentEdge;
      g_minValue = DBL_MAX;
      int y1, x1, y2, x2;
      int cy = (currentEdge.fromY + currentEdge.toY) / 2;
      int cx = (currentEdge.fromX + currentEdge.toX) / 2;
      int refineLimit = g_refineRate * g_tryLimit;
      int firstLimit = 0.33 * refineLimit;
      int secondLimit = 0.66 * refineLimit;

      for (int i = 0; i < refineLimit; i++) {
        if (i == 0) {
          y1 = currentEdge.fromY;
          x1 = currentEdge.fromX;
          y2 = currentEdge.toY;
          x2 = currentEdge.toX;
        } else if (firstLimit < i) {
          y1 = currentEdge.fromY;
          x1 = currentEdge.fromX;
          y2 = xor128()%MAX_H;
          x2 = xor128()%MAX_W;
        } else if (secondLimit < i) {
          y1 = xor128()%MAX_H;
          x1 = xor128()%MAX_W;
          y2 = currentEdge.toY;
          x2 = currentEdge.toX;
        } else {
          y1 = cy;
          x1 = cx;
          y2 = xor128()%MAX_H;
          x2 = xor128()%MAX_W;
        }

        updateLine(y1, x1, y2, x2);
        if (lineOnThePlant()) continue;

        int removeEdgeCount = removeEdgeEval(g_line);

        if (removeEdgeCount > 0) {
          double removeValue = 0.0;
          removeValue = removeRootEval(g_line, (double) removeEdgeCount);
          double eval = removeValue / (double) removeEdgeCount;

          if (g_minValue > eval) {
            g_minValue = eval;
            bestEdge.fromY = y1;
            bestEdge.fromX = x1;
            bestEdge.toY = y2;
            bestEdge.toX = x2;
          }
        }
      }

      return bestEdge;
    }

    inline bool lineOnThePlant() {
      for (int i = 0; i < g_NP; i++) {
        Vector *v = getVertex(i);

        if(ontheline(g_line.fromY, g_line.fromX, g_line.toY, g_line.toX, v->y, v->x)) return true;
      }

      return false;
    }

    vector<Edge> cleanEdges() {
      vector<Edge> result;
      priority_queue<Node, vector<Node>, greater<Node> > pqueCopy;
      int qsize = pque.size();

      fprintf(stderr,"pque size = %d\n", qsize);

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();
        Edge edge = node.edge;
        updateLine(edge.fromY, edge.fromX, edge.toY, edge.toX);

        if (cleanEdge(g_line, true)) {
          cleanEdge(g_line);
        } else {
          pqueCopy.push(node);
          result.push_back(edge);
        }
      }

      pque = pqueCopy;
      return result;
    }

    int removeEdge(const Line &line) {
      int removeCount = 0;

      for (int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (intersect(line.fromY, line.fromX, line.toY, line.toX, edge->fromY, edge->fromX, edge->toY, edge->toX)) {
          removeCount += 1001 - edge->length;

          edge->removed++;
        }
      }

      return removeCount;
    }

    int removeEdgeEval(const Line &line) {
      int removeCount = 0;

      for (int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (edge->removed > 0) continue;

        if (intersect(line.fromY, line.fromX, line.toY, line.toX, edge->fromY, edge->fromX, edge->toY, edge->toX)) {
          removeCount += 1001 - edge->length;
        }
      }

      return removeCount;
    }

    bool cleanEdge(const Line &line, bool evalMode = false) {
      for(int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (edge->removed >= 2 && evalMode) continue;

        if (intersect(line.fromY, line.fromX, line.toY, line.toX, edge->fromY, edge->fromX, edge->toY, edge->toX)) {
          if (edge->removed == 1 && evalMode) return false;

          if (!evalMode) {
            edge->removed--;
          }
        }
      }

      return true;
    }

    double removeRoot(const Line &line) {
      double removeValue = 0.0;
      g_time++;
      vector<int> arl;

      for (int i = 0; i < g_ARC; i++) {
        int rid = getActiveRoot(i);
        Root *root = getRoot(rid);

        Vector *p3 = getVertex(root->from);
        Vector *p4 = getVertex(root->to);

        if (g_updated[root->id] != g_time && intersect(line.fromY, line.fromX, line.toY, line.toX, p3->y, p3->x, p4->y, p4->x)) {
          removeValue += root->value + root->length;

          g_updated[root->id] = g_time;
          root->removed++;
          cleanRoot(root->to);
        }

        if (p4->depth <= g_depthLimit && root->removed == 0 && root->length > g_cutLimit) {
          arl.push_back(rid);
        }
      }

      g_activeRootList = arl;
      g_ARC = g_activeRootList.size();

      return removeValue;
    }

    double removeRootEval(const Line &line, double removeEdgeCount) {
      double removeValue = 0.0;

      for (int i = 0; i < g_ARC; i++) {
        int rid = getActiveRoot(i);
        Root *root = getRoot(rid);

        Vector *p3 = getVertex(root->from);
        Vector *p4 = getVertex(root->to);

        if (intersect(line.fromY, line.fromX, line.toY, line.toX, p3->y, p3->x, p4->y, p4->x)) {
          removeValue += root->value + root->length;
        }

        if (g_minValue < removeValue / removeEdgeCount) break;
      }

      return removeValue;
    }

    void refreshJar() {
      for (int i = 0; i < g_rootListSize; i++) {
        Root *root = getRoot(i);
        root->value = 0.0;
      }

      for (int i = 0; i < g_NP; i++) {
        searchRoot(i);
      }

      for (int i = g_NP; i < g_PS; i++) {
        Root *root = getRoot(i);
        Vector *v1 = getVertex(root->from);
        Vector *v2 = getVertex(root->to);

        root->value = min(v1->value, v2->value);
      }
    }

    double searchRoot(int rootId) {
      Vector *v = getVertex(rootId);
      v->value = 0.0;

      int nsize = v->roots.size();

      for (int i = 0; i < nsize; i++) {
        int rid = v->roots[i];
        Root *root = getRoot(rid);

        if (root->removed > 0) continue;

        v->value += searchRoot(root->to) + g_branchBonus;
      }

      v->value += v->dist;

      return v->value;
    }

    void cleanRoot(int rootId) {
      Vector *v = getVertex(rootId);

      int nsize = v->roots.size();

      for (int i = 0; i < nsize; i++) {
        int rid = v->roots[i];

        Root *root = getRoot(rid);
        root->removed++;
        g_updated[root->id] = g_time;

        cleanRoot(root->to);
      }
    }

    void init(int NP) {
      g_NP = NP;
      g_ARC = 0;
      g_rootListSize = g_NP;

      g_time = 1;
      memset(g_updated, 0, sizeof(g_updated));

      g_branchBonus = 5;
      g_cutLimit = 2;
      g_refineRate = 0.1;

      if (NP >= 75) {
        g_depthLimit = 6;
      } else if (NP >= 50) {
        g_depthLimit = 7;
      } else if (NP >= 30) {
        g_depthLimit = 7;
      } else {
        g_depthLimit = 9;
      }
    }

    void directTryLimit() {
      if (g_NP >= 100) {
        g_tryLimit = 1100;
      } else if (g_NP >= 95) {
        g_tryLimit = 1200;
      } else if (g_NP >= 90) {
        g_tryLimit = 1500;
      } else if (g_NP >= 85) {
        g_tryLimit = 1600;
      } else if (g_NP >= 80) {
        g_tryLimit = 1900;
      } else if (g_NP >= 75) {
        g_tryLimit = 2200;
      } else if (g_NP >= 70) {
        g_tryLimit = 2500;
      } else if (g_NP >= 65) {
        g_tryLimit = 3000;
      } else if (g_NP >= 60) {
        g_tryLimit = 3500;
      } else if (g_NP >= 55) {
        g_tryLimit = 4000;
      } else if (g_NP >= 50) {
        g_tryLimit = 4500;
      } else if (g_NP >= 45) {
        g_tryLimit = 6000;
      } else if (g_NP >= 40) {
        g_tryLimit = 8000;
      } else if (g_NP >= 35) {
        g_tryLimit = 12000;
      } else if (g_NP >= 30) {
        g_tryLimit = 18000;
      } else if (g_NP >= 25) {
        g_tryLimit = 20000;
      } else if (g_NP >= 20) {
        g_tryLimit = 20000;
      } else if (g_NP >= 15) {
        g_tryLimit = 25000;
      } else {
        g_tryLimit = 40000;
      }

      g_tryLimit *= 6;

      g_tryLimit *= (1.0 - g_refineRate);
    }

    void updateLine(int fromY, int fromX, int toY, int toX) {
      if (fromX > toX) {
        fromX ^= toX;
        toX ^= fromX;
        fromX ^= toX;

        fromY ^= toY;
        toY ^= fromY;
        fromY ^= toY;
      }

      int dy = toY - fromY;
      int dx = toX - fromX;

      if (dx == 0) {
        g_line.fromY = 0;
        g_line.fromX = fromX;
        g_line.toY = MAX_H;
        g_line.toX = toX;
      } else if (dy == 0) {
        g_line.fromY = fromY;
        g_line.fromX = 0;
        g_line.toY = toY;
        g_line.toX = MAX_W;
      } else {
        int extendL = 1+max(fromX/dx, fromY/dy);
        int extendR = 1+max((MAX_W-fromX)/dx, (MAX_H-fromY)/dy);
        g_line.fromX = fromX - dx*extendL;
        g_line.fromY = fromY - dy*extendL;
        g_line.toX = fromX + dx*extendR;
        g_line.toY = fromY + dy*extendR;
      }
    }

    inline double calcDistDetail(int y1, int x1, int y2, int x2) {
      int dy = y2 - y1;
      int dx = x2 - x1;
      return sqrt(dy*dy + dx*dx);
    }

    inline Edge* getEdge(int id) {
      return &g_edgeList[id];
    }

    inline Root* getRoot(int id) {
      return &rootList[id];
    }

    inline int getActiveRoot(int id) {
      return g_activeRootList[id];
    }

    inline Vector* getVertex(int id) {
      return &g_vertexList[id];
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
