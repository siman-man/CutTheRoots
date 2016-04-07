#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits.h>
#include <time.h>
#include <string>
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

class CutTheRoots {
  public:
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
      fprintf(stderr, "NP = %d\n", NP);
      fprintf(stderr, "points = %ld\n", points.size());
      fprintf(stderr, "roots = %ld\n", roots.size());

      vector<int> xs(NP);
      vector<int> ys(NP);
      for (int i = 0; i < NP; ++i) {
        xs[i] = points[2 * i];
        ys[i] = points[2 * i + 1];
      }

      sort(xs.begin(), xs.end());
      sort(ys.begin(), ys.end());
      vector<int> ret(4 * (NP - 1));
      for (int i = 0; i < NP - 1; ++i) {
        int x = (xs[i] + xs[i + 1]) / 2;
        int y = (ys[i] + ys[i + 1]) / 2;
        ret[4 * i] = x;
        ret[4 * i + 1] = 0;
        ret[4 * i + 2] = x;
        ret[4 * i + 3] = 1024;
        /*
        ret[8 * i] = x;
        ret[8 * i + 1] = 0;
        ret[8 * i + 2] = x;
        ret[8 * i + 3] = 1024;

        ret[8 * i + 4] = 0;
        ret[8 * i + 5] = y;
        ret[8 * i + 6] = 1024;
        ret[8 * i + 7] = y;
        */
      }
      return ret;
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

