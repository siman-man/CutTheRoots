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
#include <cstdlib>
#include <cmath>
#include <stack>
#include <queue>

using namespace std;

typedef long long ll;

const int MAX_NP = 105;
const int MAX_W = 1024;
const int MAX_H = 1024;

unsigned long long xor128(){
  static unsigned long long rx=123456789, ry=362436069, rz=521288629, rw=88675123;
  unsigned long long rt = (rx ^ (rx<<11));
  rx=ry; ry=rz; rz=rw;
  return (rw=(rw^(rw>>19))^(rt^(rt>>8)));
}

struct Point {
  int id;
  int y;
  int x;

  Point(int id = 0, int y = -1, int x = -1){
    this->id = id;
    this->y = y;
    this->x = x;
  }
};

struct Plant {
  int id;
  int cid;
  int y;
  int x;

  Plant(int id = 0, int y = -1, int x = -1){
    this->id = id;
    this->cid = 0;
    this->y = y;
    this->x = x;
  }
};

Plant plantList[MAX_NP];

Point p1;
Point p2;

int g_NP;

int g_clusterList[MAX_NP];
int g_clusterId = 1;
int g_targetClusterId = 0;
int g_clusterCount[MAX_NP];

vector<int> g_plantIdList;

class CutTheRoots {
  public:
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
      fprintf(stderr, "NP = %d\n", NP);
      fprintf(stderr, "points = %ld\n", points.size());
      fprintf(stderr, "roots = %ld\n", roots.size());

      vector<int> ret;
      init();

      g_NP = NP;

      for (int i = 0; i < NP; ++i) {
        int x = points[2*i];
        int y = points[2*i+1];

        plantList[i] = Plant(i, y, x);
      }

      for(int i = 0; i < g_NP-1; i++){
        setup();

        if(g_plantIdList.size() <= 1) {
          break;
        }

        k_means();

        int cy = (p1.y + p2.y) / 2;
        int cx = (p1.x + p2.x) / 2;

        int dy = p2.y - p1.y;
        int dx = p2.x - p1.x;

        if (dy == 0) {
          ret.push_back(cx);
          ret.push_back(cy);
          ret.push_back(cx);
          ret.push_back(cy+1);
        } else if (dx == 0) {
          ret.push_back(cx);
          ret.push_back(cy);
          ret.push_back(cx+1);
          ret.push_back(cy);
        } else {
          double slope = -1 * dy / (double)dx;
          double rslope = round(dx / (double)dy);

          fprintf(stderr,"dy = %d, dx = %d, slope = %4.2f, rslope = %4.2f\n", dy, dx, slope, rslope);
          fprintf(stderr,"slope = %d\n", (100*dy / dx));

          int nx = cx+50;
          int ny = cy-(50*dx)/dy;

          if (nx > MAX_W || ny > MAX_H) {
            int diff = max(nx, ny) - MAX_H;
            ret.push_back(max(0, cx-diff));
            ret.push_back(max(0, cy-diff));
            ret.push_back(max(0, cx+50-diff));
            ret.push_back(max(0, cy-(50*dx)/dy-diff));
          } else if (nx < 0 || ny < 0 ) {
            int diff = -1 * min(nx, ny);
            ret.push_back(min(1024, cx+diff));
            ret.push_back(min(1024, cy+diff));
            ret.push_back(min(1024, cx+50+diff));
            ret.push_back(max(0, cy-(50*dx)/dy+diff));
          } else {
            ret.push_back(cx);
            ret.push_back(cy);
            ret.push_back(cx+50);
            ret.push_back(max(0, cy-(50*dx)/dy));
          }
        }
      }

      setup();
      return ret;
    }

    void init() {
      memset(g_clusterList, 0, sizeof(g_clusterList));
    }

    void setup() {
      g_plantIdList.clear();
      memset(g_clusterCount, 0, sizeof(g_clusterCount));
      int maxCount = 1;
      g_targetClusterId = -1;

      for(int i = 0; i < g_NP; i++) {
        Plant *p = getPlant(i);
        g_clusterCount[p->cid]++;

        if(g_clusterCount[p->cid] > maxCount) {
          maxCount = g_clusterCount[p->cid];
          g_targetClusterId = p->cid;
        }
      }

      for(int i = 0; i < g_NP; i++) {
        Plant *p = getPlant(i);

        fprintf(stderr,"%d ", g_clusterCount[i]);

        if(g_targetClusterId == p->cid) {
          g_plantIdList.push_back(p->id);
        }
      }

      fprintf(stderr,"\n");
    }

    /**
     * set random point
     */
    void setRandomPoint() {
      fprintf(stderr,"setRandomPoint =>\n");
      int id1 = 0;
      int id2 = 0;
      int size = g_plantIdList.size();

      while(id1 == id2) {
        id1 = g_plantIdList[xor128()%size];
        id2 = g_plantIdList[xor128()%size];
      }

      Plant *plantA = getPlant(id1);
      Plant *plantB = getPlant(id2);

      p1.y = plantA->y;
      p1.x = plantA->x;

      p2.y = plantB->y;
      p2.x = plantB->x;
    }

    void k_means() {
      setRandomPoint();

      fprintf(stderr,"cluster id = %d, cluster count = %ld\n", g_targetClusterId, g_plantIdList.size());

      while(true) {
        if(!separeteNP()) {
          break;
        }
        replacePoint();
      }

      g_clusterId++;
      fprintf(stderr,"p1.y = %d, p1.x = %d\n", p1.y, p1.x);
      fprintf(stderr,"p2.y = %d, p2.x = %d\n", p2.y, p2.x);
    }

    /**
     * k-means(2)
     */
    bool separeteNP(){
      bool update = false;
      int size = g_plantIdList.size();

      fprintf(stderr,"size: %d, separeteNP =>\n", size);
      assert(g_targetClusterId != g_clusterId);

      for(int i = 0; i < size; i++) {
        int id = g_plantIdList[i];
        Plant *p = getPlant(id);
        int d1 = calcDist(p->y, p->x, p1.y, p1.x);
        int d2 = calcDist(p->y, p->x, p2.y, p2.x);
        int temp = p->cid;

        if(d1 <= d2) {
          g_clusterList[id] = g_targetClusterId;
          p->cid = g_targetClusterId;
        } else {
          g_clusterList[id] = g_clusterId;
          p->cid = g_clusterId;
        }

        if (temp != g_clusterList[id]) {
          update = true;
        }
      }

      return update;
    }

    void replacePoint() {
      fprintf(stderr,"replacePoint =>\n");
      int c1Count = 0;
      int c2Count = 0;
      int sumY1 = 0;
      int sumX1 = 0;
      int sumY2 = 0;
      int sumX2 = 0;
      int size = g_plantIdList.size();

      for(int i = 0; i < size; i++){
        int id = g_plantIdList[i];
        Plant *p = getPlant(id);

        if(g_clusterList[id] == g_targetClusterId) {
          sumY1 += p->y;
          sumX1 += p->x;
          c1Count++;
        } else {
          sumY2 += p->y;
          sumX2 += p->x;
          c2Count++;
        }
      }

      assert(c1Count > 0 && c2Count > 0);

      p1.y = sumY1 / c1Count;
      p1.x = sumX1 / c1Count;

      p2.y = sumY2 / c2Count;
      p2.x = sumX2 / c2Count;
    }

    int calcDist(int y1, int x1, int y2, int x2) {
      int dy = y1 - y2;
      int dx = x1 - x2;
      return dy*dy + dx*dx;
    }

    Plant* getPlant(int id) {
      return &plantList[id];
    }
};

int main() {
  int NP;
  cin >> NP;

  int Npoints;
  cin >> Npoints;
  vector<int> points(Npoints);
  for(int i = 0; i < Npoints; i++){
    cin >> points[i];
  }

  int Nroots;
  cin >> Nroots;
  vector<int> roots(Nroots);
  for(int i = 0; i < Nroots; i++){
    cin >> roots[i];
  }

  CutTheRoots cr;
  vector<int> ret = cr.makeCuts(NP, points, roots);

  cout << ret.size() << endl;
  for (int i = 0; i < ret.size(); ++i) {
    cout << ret[i] << endl;
  }
  cout.flush();
}

