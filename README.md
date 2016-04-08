Problem: CutTheRoots
問題: CutTheRoots

Problem Statement
問題文
Several seeds were planted in one jar. The grown-up plants have to be separated and moved to other jars, one plant per jar.
いくつかの植物の種が1つのツボの中に植えれらています。成長した植物は1つのツボにつき1つの植物となるように他のツボに移さなければいけません
While growing, the plants developed a massive intertwining root system which has to be separated as well with as little damage as possible. 
成長していく中で、巨大に絡み合った根系を形成した植物をなるべく最小のダメージで切り離す必要があります。

The jar with plants is represented as a circle on a plane. Each plant is represented as a 2D root system, 
ツボの中の植物は平面上の円形で表現されます. 各植物は平面の根系を持ちます

since we are not interested in its stem (it doesn't need to be cut when the plants are separated). 
幹については考慮しなくて良いです

The root system of one plant starts with a plant base point - the point where root system connects to the stem. 
根系は植物の基盤部分(幹に接続している部分）から始まり、
Each root grows either from the base point of a plant or from the end of another root. 
各根っこは他の根っこから伸びて成長していきます

To separate the plants, you need to cut the soil in the jar. 
植物を分けるためには、あなたはツボから土ごと切り分ける必要があります。
Each cut is modeled as an infinite straight line which goes through points (X1, Y1) and (X2, Y2). 
切り取る際には(X1, Y1)と(X2, Y2)の2点を通る無限の長さの直線が使用されます
Cut lines cut every root they intersect. 
その線で横切られた根っこは切り取られます
Two plants are considered to be separated as long as there exists at least one line which separates their base points. 
二つの植物は、その基点を分離する少なくとも1つのラインが存在するように分離することが考えられています

Your code must implement a single method int[] makeCuts(int NP, int[] points, int[] roots). 
あなたは次のメソッドmakeCutsを実装します。
The parameters give you the following information:
パラメータについては以下のとおりです

NP - the number of plants.
NP - 植物の数
points - the coordinates of endpoints of roots. Endpoint i (0-based) has coordinates (point[2*i], point[2*i+1]). 
各植物の根っこの座標
Note that there will be more than NP points in this array; the first NP points will be base points of plants.
注意: pointsの数はNPより多くなる可能性があります、NPは植物の幹の数です

roots - the roots, represented with indices of their endpoints. 
roots - 根っこ、どの植物に対応しているかを示す
Root j (0-based) connects endpoints roots[2*j] and roots[2*j+1]. 

It is guaranteed that each root endpoint belongs to exactly one plant 
必ずどの植物に対応しているかは決まります
(i.e. no chain of roots exists which would connect base points of two different plants). 

Roots of one plant or different plants can intersect each other, but not share an endpoint.


You must return a vector <int> representing a set of at most 4*NP lines which separate all plants. If your 
最大で4*NPの数の分割直線を指定できます。
return describes L lines, it must contain exactly 4*L elements, line k (0-based) is defined as going through points 
もしLを指定した場合は4*Lの数の要素を含む必要があります。
(ret[4*k], ret[4*k+1]) and (ret[4*k+2], ret[4*k+3]). 


The score for an individual test case is the total length of roots which stay attached to the base point of their plants,
各テストケースのスコアは、各植物に残っている根っこの長さの総長になります。
divided by total length of all roots before cutting. If there exist two plants which are not separated with the cuts, or
カットする前の総長、                                     もし2つ以上の植物が残っていた場合にはスコアは0点です
the return value was invalid in any other way, the score for that test case is 0. 
また、不正な値が入っていた場合にも0点です。
Your overall score will be the sum of YOUR/MAX over all test cases, where YOUR is your raw score on a test case, 
最終的なスコアは（あなたのスコア/そのテストケースでの最大点)となります
and MAX is the maximal score achieved by anyone on that test case (test cases on which you scored 0 don't 
contribute to the sum).Tools and Information
An offline tester is available here. You can use it to test/debug your solution locally. 
You can also check its source code for exact implementation of test case generation and score calculation. 
That page also contains links to useful information and sample solutions in several languages.

 
Definition
      
Class:  CutTheRoots
Method: makeCuts
Parameters: int, vector <int>, vector <int>
Returns:  vector <int>
Method signature: vector <int> makeCuts(int NP, vector <int> points, vector <int> roots)
(be sure your method is public)
    
 
Notes
- This match is rated.
この試合はレートがつきます
- The time limit is 10 seconds and the memory limit is 1024MB for a single test case.
制限時間は10秒です、メモリの制限は1024MBです
- The number of plants (NP) will be less than 105 and there will be at least 5 plants. There will be less than 105000 roots and at least 50.
NPの値は最大で105、最小で5の値をとります。根っこは最大で105000、最小で50の値を取ります
- X1,Y1,X2 and Y2 must be between 0 and 1024. All plant and root coordinates will also fall within this range.
X1,Y1,X2,Y2の値は最大で1024、最小で0の値をとります、この範囲外の値は無効となります。


## 考察

### 未知のもの

* 植物を分割する際に、根っこの切断を最小にするための切り分け方

### 与えられているもの

* 各植物の座標
* 根っこの座標
* 根っこの長さ
* 根っこがどの植物に属しているか

### 条件

* 必ず1つずつ植物を切り分けなければいけない

### 類似問題

* 去年のSmallPolygonがそれに近い
* グラフの最大カット

### 問題の置き換え


### アプローチ

植物の幹同士を接続したグラフを作成、これを全てバラバラになるように切断する。
これで最低限の回答が得られる。
