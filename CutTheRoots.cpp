#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <string.h>
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

class UnionFind {
  public:

  int parent[MAX_PS];
  int rank[MAX_PS];

  void init(int n) {
    for (int i = 0; i < n; i++) {
      parent[i] = i;
      rank[i] = 0;
    }
  }

  int find(int x) {
    if (parent[x] == x) {
      return x;
    } else {
      return parent[x] = find(parent[x]);
    }
  }

  void unite(int x, int y) {
    x = find(x);
    y = find(y);
    if (x == y) return;

    if (rank[x] < rank[y]) {
      parent[x] = y;
    } else {
      parent[y] = x;
      if (rank[x] == rank[y]) rank[x]++;
    }
  }

  bool same(int y, int x) {
    return find(x) == find(y);
  }
};

class Vector {
  public:
    int id;
    int x, y;
    double value;
    double dist;
    int depth;
    unordered_set<int> roots;

    Vector(int y = -1, int x = -1) {
      this->value = 0.0;
      this->dist = 0.0;
      this->depth = 0;
      this->y = y;
      this->x = x;
    }

    bool operator<(const Vector& v) const {
      return (x != v.x)? x < v.x : y < v.y;
    }

    Vector operator-(Vector p) {
      return Vector(y - p.y, x - p.x);
    }

    double norm() {
      return x*x + y*y;
    }
};

vector<Vector> g_vertexList;
typedef vector<Vector> Polygon;

inline int dot(Vector &a, Vector &b){
  return a.x * b.x + a.y * b.y;
}

inline int cross(Vector &a, Vector &b){
  return a.x * b.y - a.y * b.x;
}

inline double calcArea(Polygon s) {
  double area = 0.0;
  int size = s.size();

  for (int i = 0; i < size; i++) {
    area += cross(s[i], s[(i+1)%size]);
  }

  return 0.5 * fabs(area);
}

inline bool ccw(int ay, int ax, int by, int bx, int cy, int cx) {
  return ((bx-ax) * (cy-ay) - (by-ay) * (cx-ax) >= 0);
}

inline bool intersect(ll y1, ll x1, ll y2, ll x2, ll y3, ll x3, ll y4, ll x4) {
  return (((x2-x1) * (y3-y1) - (y2-y1) * (x3-x1)) * ((x2-x1) * (y4-y1) - (y2-y1) * (x4-x1)) <= 0);
}

inline bool ontheline(int ay, int ax, int by, int bx, int cy, int cx) {
  return ((bx-ax) * (cy-ay) - (by-ay) * (cx-ax) == 0);
}

Polygon andrewScan(Polygon s) {
  Polygon u, l;

  if(s.size() < 3) return s;
  sort(s.begin(), s.end());

  u.push_back(s[0]);
  u.push_back(s[1]);

  l.push_back(s[s.size() - 1]);
  l.push_back(s[s.size() - 2]);

  for (int i = 2; i < s.size(); i++) {
    for (int n = u.size(); n >= 2 && ccw(u[n-2].y, u[n-2].x, u[n-1].y, u[n-1].x, s[i].y, s[i].x); n--) {
      u.pop_back();
    }
    u.push_back(s[i]);
  }

  for (int i = s.size() -3; i >= 0; i--) {
    for (int n = l.size(); n >= 2 && ccw(l[n-2].y, l[n-2].x, l[n-1].y, l[n-1].x, s[i].y, s[i].x); n--) {
      l.pop_back();
    }
    l.push_back(s[i]);
  }

  reverse(l.begin(), l.end());

  for (int i = u.size() - 2; i >= 1; i--) {
    l.push_back(u[i]);
  }

  return l;
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
  double removeValue;

  Node(Edge &edge, double value){ 
    this->edge = edge;
    this->removeValue = value;
  }
    
  bool operator >(const Node &e) const {
    return removeValue < e.removeValue;
  }
};

int g_NP;
int g_PS;
int g_RC;
int g_depthLimit;
int g_branchBonus;
double g_randomRate;
int g_tryLimit;
int g_cutLimit;
int g_rootListSize;
int g_ARC;
int g_edgeListSize;
int g_vsize;
int g_updated[MAX_PS];
int g_time;
Line g_line;

vector<Edge> g_edgeList;
vector<int> g_activeRootList;
vector<int> g_convexHullVertex;
Root rootList[MAX_PS];
UnionFind uf;

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

        from->roots.insert(root.id);

        to->dist = dist;
        to->depth = from->depth + 1;

        if (to->depth > 0) {
          uf.unite(j, k);
        }

        root.depth = to->depth;
        rootList[g_rootListSize++] = root;

        if (to->depth <= g_depthLimit && dist > g_cutLimit) {
          g_activeRootList.push_back(root.id);
        }
      }

      g_ARC = g_activeRootList.size();
      unordered_map<int, Polygon> polygons;

      directTryLimit();

      for (int i = NP; i < g_PS; i++) {
        int id = uf.find(i);

        Vector v = g_vertexList[i];

        polygons[id].push_back(v);
      }

      unordered_map<int, Polygon>::iterator pit = polygons.begin();
      unordered_map<int, bool> vcheck;

      while (pit != polygons.end()) {
        Polygon polygon = (*pit).second;
        Polygon totsu = andrewScan(polygon);
        pit++;

        if (totsu.size() <= 2) continue;

        vector<Root> eds = polygon2roots(totsu);
        int ees = eds.size();

        for (int i = 0; i < ees; i++) {
          Root rt = eds[i];

          if (!vcheck[rt.from]) {
            g_convexHullVertex.push_back(rt.from);
            vcheck[rt.from] = true;
          }

          if (!vcheck[rt.to]) {
            g_convexHullVertex.push_back(rt.to);
            vcheck[rt.to] = true;
          }
        }
      }

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
      g_vsize = g_convexHullVertex.size();

      int rsize = g_ARC;
      fprintf(stderr, "edge size = %d, root size = %d\n", g_edgeListSize, rsize);

      while (!finishCheck()) {
        refreshJar();
        Edge edge = getBestEdge();

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
      double maxValue = -DBL_MAX;
      int y1, x1, y2, x2;
      int checkpoint = g_randomRate * g_tryLimit;

      for (int i = 0; i < g_tryLimit; i++) {
        if (i <= checkpoint) {
          y1 = xor128()%MAX_H;
          x1 = xor128()%MAX_W;

          y2 = xor128()%MAX_H;
          x2 = xor128()%MAX_W;
        } else {
          int a = g_convexHullVertex[xor128()%g_vsize];
          int b = g_convexHullVertex[xor128()%g_vsize];

          while (a == b) {
            a = g_convexHullVertex[xor128()%g_vsize];
            b = g_convexHullVertex[xor128()%g_vsize];
          }

          Vector *v1 = getVertex(a);
          Vector *v2 = getVertex(b);

          y1 = v1->y;
          x1 = v1->x;
          y2 = v2->y;
          x2 = v2->x;
        }

        updateLine(y1, x1, y2, x2);
        if (lineOnThePlant()) continue;

        int removeEdgeCount = removeEdgeEval(g_line);

        if (removeEdgeCount > 0) {
          double removeValue = 0.0;
          removeValue = removeRootEval(g_line);
          double eval = -1 * removeValue / (double) removeEdgeCount;

          if (maxValue < eval) {
            maxValue = eval;
            bestEdge.fromY = y1;
            bestEdge.fromX = x1;
            bestEdge.toY = y2;
            bestEdge.toX = x2;
          }
        }
      }

      return bestEdge;
    }

    bool lineOnThePlant() {
      for (int i = 0; i < g_NP; i++) {
        Vector *v = getVertex(i);

        if(ontheline(g_line.fromY, g_line.fromX, g_line.toY, g_line.toX, v->y, v->x)) {
          return true;
        }
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
          p3->value -= p4->value;
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

    double removeRootEval(const Line &line) {
      double removeValue = 0.0;

      for (int i = 0; i < g_ARC; i++) {
        int rid = getActiveRoot(i);
        Root *root = getRoot(rid);

        Vector *p3 = getVertex(root->from);
        Vector *p4 = getVertex(root->to);

        if (intersect(line.fromY, line.fromX, line.toY, line.toX, p3->y, p3->x, p4->y, p4->x)) {
          removeValue += root->value + root->length;
        }
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

      vector<int> convexList;

      for (int i = 0; i < g_vsize; i++) {
        int vid = g_convexHullVertex[i];
        Vector *v = getVertex(vid);

        if (v->value > 0) {
          convexList.push_back(vid);
        }
      }

      g_convexHullVertex = convexList;
      g_vsize = g_convexHullVertex.size();
    }

    double searchRoot(int rootId) {
      Vector *v = getVertex(rootId);
      v->value = 0.0;
      double value = 0.0;

      unordered_set<int>::iterator it = v->roots.begin();

      while(it != v->roots.end()) {
        int rid = (*it);
        it++;
        Root *root = getRoot(rid);

        if (root->removed > 0) continue;

        value += searchRoot(root->to) + g_branchBonus;
      }

      value += v->dist;
      v->value = value;

      return value;
    }

    void cleanRoot(int rootId) {
      Vector *v = getVertex(rootId);

      unordered_set<int>::iterator it = v->roots.begin();

      while (it != v->roots.end()) {
        int rid = (*it);
        it++;

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

      for(int i = 0; i < MAX_H+10; i++) {
        randomNum[i] = xor128();
      }

      uf.init(g_PS);
      g_branchBonus = 5;
      g_cutLimit = 2;
      g_randomRate = 0.1;

      if (NP >= 75) {
        g_depthLimit = 4;
      } else if (NP >= 50) {
        g_depthLimit = 4;
      } else if (NP >= 30) {
        g_depthLimit = 5;
      } else if (NP >= 15) {
        g_depthLimit = 6;
      } else {
        g_depthLimit = 7;
      }
    }

    void directTryLimit() {
      if (g_NP >= 90) {
        g_tryLimit = 1000;
      } else if (g_NP >= 85) {
        g_tryLimit = 1300;
      } else if (g_NP >= 75) {
        g_tryLimit = 1600;
      } else if (g_NP >= 70) {
        g_tryLimit = 1800;
      } else if (g_NP >= 60) {
        g_tryLimit = 1900;
      } else if (g_NP >= 50) {
        g_tryLimit = 2000;
      } else if (g_NP >= 40) {
        g_tryLimit = 3200;
      } else if (g_NP >= 30) {
        g_tryLimit = 3500;
      } else if (g_NP >= 25) {
        g_tryLimit = 3800;
      } else if (g_NP >= 20) {
        g_tryLimit = 5000;
      } else if (g_NP >= 15) {
        g_tryLimit = 6500;
      } else {
        g_tryLimit = 20000;
      }
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

    vector<Root> polygon2roots(Polygon &pol) {
      vector<Root> roots;
      int psize = pol.size();

      for(int i = 0; i < psize; i++) {
        Vector v1 = pol[i%psize];
        Vector v2 = pol[(i+1)%psize];
        Root root;

        if (v1.depth <= v2.depth) {
          root.from = v1.id;
          root.to = v2.id;
        } else {
          root.from = v2.id;
          root.to = v1.id;
        }

        roots.push_back(root);
      }

      return roots;
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
