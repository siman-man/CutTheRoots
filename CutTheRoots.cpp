#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits.h>
#include <string>
#include <cassert>
#include <string.h>
#include <set>
#include <cstdio>
#include <cfloat>
#include <cstdlib>
#include <cmath>
#include <queue>

using namespace std;

typedef long long ll;

const int MAX_NP = 105;
const int MAX_W = 1024;
const int MAX_H = 1024;
const int MAX_PS = 105000;

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

class UnionFind {
  int parent[MAX_PS];
  int rank[MAX_PS];

  public:

  void init(int n){
    for(int i = 0; i < n; i++){
      parent[i] = i;
      rank[i] = 0;
    }
  }

  int find(int x){
    if (parent[x] == x) {
      return x;
    } else {
      return parent[x] = find(parent[x]);
    }
  }

  void unite(int x, int y){
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
    double x, y;
    int id;
    int value;
    int dist;
    int depth;
    int branchCount;
    set<int> roots;

    Vector(double y = 0.0, double x = 0.0) {
      this->value = 0;
      this->dist = 0;
      this->depth = 0;
      this->branchCount = 0;
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

int norm2(int y, int x) {
  return y*y + x*x;
}

vector<Vector> vertexList;

typedef vector<Vector> Polygon;

inline int dot(Vector &a, Vector &b){
  return a.x * b.x + a.y * b.y;
}

inline int dot2(int ay, int ax, int by, int bx) {
  return ax * bx + ay * by;
}

inline int cross(Vector &a, Vector &b){
  return a.x * b.y - a.y * b.x;
}

inline int cross2(int ay, int ax, int by, int bx) {
  return ax * by - ay * bx;
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

inline int ccw2(int ay, int ax, int by, int bx, int cy, int cx) {
  int y1 = by - ay;
  int x1 = bx - ax;
  int y2 = cy - ay;
  int x2 = cx - ax;

  if(cross2(y1, x1, y2, x2) > EPS) return COUNTER_CLOCKWISE;
  if(cross2(y1, x1, y2, x2) < -EPS) return CLOCKWISE;
  //if(dot2(y1, x1, y2, x2) < -EPS) return ONLINE_BACK;
  //if(norm2(y1, x1) < norm2(y2, x2)) return ONLINE_FRONT;

  return ON_SEGMENT;
}

inline bool intersect(Vector &p1, Vector &p2, Vector &p3, Vector &p4){
  return ((ccw(p1, p2, p3) * ccw(p1, p2, p4) <= 0) && (ccw(p3, p4, p1) * ccw(p3, p4, p2) < 0));
}

inline bool intersect_fast(int y1, int x1, int y2, int x2, int y3, int x3, int y4, int x4) {
  return ((ccw2(y1, x1, y2, x2, y3, x3) * ccw2(y1, x1, y2, x2, y4, x4) <= 0) && (ccw2(y3, x3, y4, x4, y1, x1) * ccw2(y3, x3, y4, x4, y2, x2) < 0));
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
    for (int n = u.size(); n >= 2 && ccw(u[n-2], u[n-1], s[i]) != CLOCKWISE; n--) {
      u.pop_back();
    }
    u.push_back(s[i]);
  }

  for (int i = s.size() -3; i >= 0; i--) {
    for (int n = l.size(); n >= 2 && ccw(l[n-2], l[n-1], s[i]) != CLOCKWISE; n--) {
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

/**
 * IN 2
 * ON 1
 * OUT 0
 */
int contains(Polygon g, Vector v) {
  int n = g.size();
  bool x = false;

  for (int i = 0; i < n; i++) {
    Vector a = g[i] - v;
    Vector b = g[(i+1)%n] - v;

    if(abs(cross(a, b)) < EPS && dot(a, b) < EPS) return 1;
    if(a.y > b.y) swap(a, b);
    if(a.y < EPS && EPS < b.y && cross(a, b) > EPS) x = !x;
  }
  return (x ? 2 : 0);
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
  int from;
  int to;
  int removed;
  ll hashCode;

  Edge(int fromY = -1, int fromX = -1, int toY = -1, int toX = -1) {
    this->fromY = fromY;
    this->fromX = fromX;
    this->toY = toY;
    this->toX = toX;

    this->from = -1;
    this->to = -1;
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

struct Node {
  Edge edge;
  int removeValue;

  Node(Edge &edge, int value){ 
    this->edge = edge;
    this->removeValue = value;
  }
    
  bool operator >(const Node &e) const{
    return removeValue < e.removeValue;
  }    
};

int g_NP;
int g_PS;
int g_depthLimit;
vector<Edge> edgeList;
vector<Root> activeRootList;
vector<int> convexHullVertex;
Root rootList[105000];
int g_rootListSize;
int g_activeRootSize;
int g_edgeListSize;
Line g_line;
UnionFind uf;

priority_queue<Node, vector<Node>, greater<Node> > pque;

class CutTheRoots {
  public:
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
      fprintf(stderr, "NP = %d\n", NP);
      fprintf(stderr, "points = %ld\n", points.size());
      fprintf(stderr, "roots = %ld\n", roots.size());

      vector<int> ret;

      g_NP = NP;
      g_PS = points.size()/2;

      init();

      for (int i = 0; i < g_PS; ++i) {
        Vector v;
        v.id = i;
        v.x = points[2*i];
        v.y = points[2*i+1];

        vertexList.push_back(v);
      }

      g_depthLimit = 4;

      if (NP <= 50) {
        g_depthLimit++;
      }
      if (NP <= 30) {
        g_depthLimit++;
      }
      if (NP <= 15) {
        g_depthLimit++;
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

        Root root(j, k);
        root.id = i;
        int dist = calcDistDetail(fromY, fromX, toY, toX);

        from->roots.insert(root.id);

        to->dist = dist;
        to->depth = from->depth + 1;

        if (to->depth > 1) {
          uf.unite(j, k);
        }

        root.depth = to->depth;
        rootList[g_rootListSize++] = root;

        if (to->depth <= g_depthLimit && dist > 3) {
          rootList[i].aid = g_activeRootSize;
          g_activeRootSize++;
          activeRootList.push_back(root);
        }
      }

      for (int i = 0; i < NP; i++) {
        searchRoot(i);
      }

      map<int, bool> checkListV;
      map<int, Polygon> polygons;

      for (int i = NP; i < g_PS; i++) {
        int id = uf.find(i);

        Vector v = vertexList[i];

        polygons[id].push_back(v);
      }

      map<int, Polygon>::iterator pit = polygons.begin();
      map<int, bool> vcheck;

      while (pit != polygons.end()) {
        Polygon pol = (*pit).second;
        Polygon totsu = andrewScan(pol);
        vector<Root> eds = polygon2roots(totsu);
        int ees = eds.size();

        if (ees >= 3) {
          for (int i = 0; i < ees; i++) {
            Root rt = eds[i];

            if (!vcheck[rt.from]) {
              convexHullVertex.push_back(rt.from);
              vcheck[rt.from] = true;
            }
            if (!vcheck[rt.to]) {
              convexHullVertex.push_back(rt.to);
              vcheck[rt.to] = true;
            }
          }
        }

        pit++;
      }

      for (int i = NP; i < g_PS; i++) {
        Root *root = getRoot(i);
        Vector *v1 = getVertex(root->from);
        Vector *v2 = getVertex(root->to);

        root->value = min(v1->value, v2->value);

        if (root->aid != -1) {
          Root *aroot = getActiveRoot(root->aid);
          aroot->value = root->value;
        }
      }

      for (int i = 0; i < NP-1; i++) {
        Vector *v1 = getVertex(i);

        for (int j = i+1; j < NP; j++) {
          Vector *v2 = getVertex(j);
          Edge newEdge(v1->y, v1->x, v2->y, v2->x);
          newEdge.from = i;
          newEdge.to = j;

          edgeList.push_back(newEdge);
        }
      }

      g_edgeListSize = edgeList.size();

      int rsize = g_activeRootSize;
      fprintf(stderr, "edge size = %d, root size = %d\n", g_edgeListSize, rsize);

      for (int i = 0; i < NP; i++) {
        Edge edge = getBestEdge();

        if (edge.fromY == -1) {
          continue;
        }

        edge2line(edge);

        int removeValue = removeRoot(g_line);
        removeEdge(g_line);

        pque.push(Node(edge, removeValue));
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

    Edge getBestEdge() {
      Edge bestEdge;
      int maxValue = INT_MIN;
      int limit = 100;

      if (g_NP >= 90) {
        limit = 400;
      } else if (g_NP >= 75) {
        limit = 600;
      } else if (g_NP >= 60) {
        limit = 800;
      } else if (g_NP >= 40) {
        limit = 1000;
      } else if (g_NP >= 30) {
        limit = 1200;
      } else if (g_NP >= 25) {
        limit = 1200;
      } else {
        limit = 1500;
      }

      int vsize = convexHullVertex.size();

      for (int i = 0; i < limit; i++) {
        int cyA, cxA, cyB, cxB;

        if (xor128()%100 <= 10) {
          cyA = xor128()%MAX_H;
          cxA = xor128()%MAX_W;

          cyB = xor128()%MAX_H;
          cxB = xor128()%MAX_W;
        } else {
          int a = convexHullVertex[xor128()%vsize];
          int b = convexHullVertex[xor128()%vsize];

          while (a == b) {
            a = convexHullVertex[xor128()%vsize];
            b = convexHullVertex[xor128()%vsize];
          }

          Vector *v1 = getVertex(a);
          Vector *v2 = getVertex(b);

          cyA = v1->y;
          cxA = v1->x;
          cyB = v2->y;
          cxB = v2->x;
        }

        Edge edge(cyA, cxA, cyB, cxB);
        edge2line(edge);

        int removeValue = 0;
        removeValue = removeRootEval(g_line);
        int removeCount = removeEdge(g_line, true);
        int eval = 200 * removeCount - removeValue;

        if(removeCount > 0 && maxValue < eval) {
          maxValue = eval;
          bestEdge = edge;
        }
      }

      return bestEdge;
    }

    vector<Edge> cleanEdges() {
      vector<Edge> result;
      int cleanCount = 0;
      priority_queue<Node, vector<Node>, greater<Node> > pqueCopy;

      while (!pque.empty()) {
        Node node = pque.top(); pque.pop();
        Edge edge = node.edge;
        edge2line(edge);

        if (cleanEdge(g_line, true)) {
          cleanEdge(g_line);
          cleanCount++;
        } else {
          pqueCopy.push(node);
          result.push_back(edge);
        }
      }

      pque = pqueCopy;
      return result;
    }

    int removeEdge(Line &line, bool evalMode = false) {
      int removeCount = 0;

      for (int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (edge->removed > 0 && evalMode) {
          continue;
        }

        if (intersect_fast(line.fromY, line.fromX, line.toY, line.toX, edge->fromY, edge->fromX, edge->toY, edge->toX)) {
          removeCount++;

          if (!evalMode) {
            edge->removed++;
          }
        }
      }

      return removeCount;
    }

    bool cleanEdge(Line &line, bool evalMode = false) {
      for(int i = 0; i < g_edgeListSize; i++) {
        Edge *edge = getEdge(i);

        if (edge->removed >= 2 && evalMode) {
          continue;
        }

        if (intersect_fast(line.fromY, line.fromX, line.toY, line.toX,
              edge->fromY, edge->fromX, edge->toY, edge->toX)) {
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

    int removeRoot(Line &line) {
      int removeValue = 0;

      for(int i = 0; i < g_activeRootSize; i++) {
        Root *root = getActiveRoot(i);

        Vector *p3 = getVertex(root->from);
        Vector *p4 = getVertex(root->to);

        if (intersect_fast(line.fromY, line.fromX, line.toY, line.toX, p3->y, p3->x, p4->y, p4->x)) {
          removeValue += root->value;

          root->removed++;
          cleanRoot(root->to);
        }
      }

      return removeValue;
    }

    int removeRootEval(Line &line) {
      int removeValue = 0;

      for(int i = 0; i < g_activeRootSize; i++) {
        Root *root = getActiveRoot(i);

        if (root->removed > 0) continue;

        Vector *p3 = getVertex(root->from);
        Vector *p4 = getVertex(root->to);

        if (intersect_fast(line.fromY, line.fromX, line.toY, line.toX, p3->y, p3->x, p4->y, p4->x)) {
          removeValue += root->value;
        }
      }

      return removeValue;
    }

    int searchRoot(int rootId) {
      Vector *v = getVertex(rootId);
      int value = 0;

      set<int>::iterator it = v->roots.begin();

      while(it != v->roots.end()) {
        int rid = (*it);
        Root *root = getRoot(rid);

        value += searchRoot(root->to);
        it++;
      }

      value += v->dist + 5;
      v->value = value;

      return value;
    }

    void cleanRoot(int rootId) {
      Vector *v = getVertex(rootId);

      set<int>::iterator it = v->roots.begin();

      while(it != v->roots.end()) {
        int rid = (*it);
        it++;

        if (g_activeRootSize <= rid) continue;
        Root *r = getRoot(rid);

        if (r->aid < 0) continue;
        Root *root = getActiveRoot(r->aid);
        root->removed++;
        root->value = 0.0;

        cleanRoot(root->to);
      }
    }

    void init() {
      g_activeRootSize = 0;
      g_rootListSize = g_NP;

      for(int i = 0; i < MAX_H+10; i++) {
        randomNum[i] = xor128();
      }

      uf.init(g_PS);
    }

    void edge2line(Edge &edge) {
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
        int ox = fromX;
        int oy = fromY;
        int extendL = 1+max(ox/dx, oy/dy);
        int extendR = 1+max((MAX_W-ox)/dx, (MAX_H-oy)/dy);
        g_line.fromX = ox - dx*extendL;
        g_line.fromY = oy - dy*extendL;
        g_line.toX = ox + dx*extendR;
        g_line.toY = oy + dy*extendR;
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

    inline int calcDistDetail(int y1, int x1, int y2, int x2) {
      int dy = y1 - y2;
      int dx = x1 - x2;
      return (int)sqrt(dy*dy + dx*dx);
    }

    inline Edge* getEdge(int id) {
      return &edgeList[id];
    }

    inline Root* getRoot(int id) {
      return &rootList[id];
    }

    inline Root* getActiveRoot(int id) {
      return &activeRootList[id];
    }

    inline Vector* getVertex(int id) {
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
