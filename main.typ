#import "@preview/showybox:2.0.1": showybox

#import "@preview/codly:0.2.0": *
#show: codly-init.with()

#import "@preview/lovelace:0.2.0": *
#show: setup-lovelace
#let comment(body) = {
  let centered-asterisk = "*" // "\u{2217}"
  let s = "/" + centered-asterisk + " " + body + " " + centered-asterisk + "/"
  let t = text(fill: gray.darken(30%),
  font: "Source Code Pro",
  s)
  t
}
#set page(
  paper: "a4",
  numbering: "1",
  // flipped: true,
  // columns: 2,
  // margin: 3%,
  header: [TCR The minute his head is in view, hit it with the rock!]
)
#counter(page).update(1)

#codly(
  languages: (
    cpp: (
      name: "C++",
      icon: "", //text(font: "FiraCode Nerd Font", "\u{e61d}"),
      color: rgb("#2995df")
    ),
  )
)
#set heading(numbering: "1.")
#show outline.entry.where(
  level: 1
): it => {
  v(12pt, weak: true)
  strong(it)
}
#outline(indent: auto)


Compilation Command

`g++ -std=c++11 -O2 -Wall test.cpp -o test`

= Template

```cpp
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using ld = long double;
using vi = vector<ll>;
typedef pair<ll, ll> pii;
#define FOR(a, c) for (int(a) = 0; (a) < (c); (a)++) 
#define all(x) begin(x), end(x)
#define allRev(x) rbegin(x), rend(x)
#define FAST_IO                     \
  ios_base::sync_with_stdio(false); \
  cin.tie(0);                       \
  cout.tie(0);

int main() {
  // freopen("input.txt", "r", stdin); // redirect file to STDIN
}
```
128bit
```cpp
using lll = __int128;
ostream& operator<<( ostream& o, __int128 n ) {
	auto t = n<0 ? -n : n; char b[128], *d = end(b);
	do *--d = '0'+t%10, t /= 10; while (t);
	if(n<0) *--d = '-';
	o.rdbuf()->sputn(d,end(b)-d);
	return o;
}
```

== split str
```cpp
vector<string> split_str(string s, const char delim) {
  vector<string> out;
  istringstream ss(s);
  string e;
  while (getline(ss, e, delim)){
    out.push_back(e);
  }
  return out;
}
```

== Vectors
```cpp
// insert
v.insert(v.begin() + index, value); // head or index | O(n)
v.push_back(value);	                // tail | O(1) access
auto v_head = v.front();            // head | O(1)
auto v_index = v.at(index);         // index (or v[i]) | O(1)
auto v_tail = v.back();	            // tail | O(1)
// remove
v.erase(v.begin() + index);         // head or index | O(n)
v.pop_back();                       // tail | O(1)
// clear - common to all data structures
v.clear();                          // O(n)
```

== Lists
```cpp
list<int> l;
int value = 1;
int index = 1;
// insert
l.push_front(value);                // head | O(1)
l.insert(l.begin(), index, value);  // index | O(n)
l.push_back(value);                 // tail | O(1)
// access
auto v_head = l.front();            // head | O(1)
auto v_value = *next(l.begin(), 1); // index | O(n)
auto v_tail = l.back();             // tail | O(1)
// remove
auto pos = next(l.begin(),index);
l.pop_front();                      // head | O(1)
l.erase(pos);                       // index | O(1) + O(n)
l.pop_back();                       // tail | O(1)

list<int> l = {1,2,3};
list<int> l2;
l2.splice(l2.begin(), l);           // transfer elements between list
l.remove(1);                        // remove an element by value
l.unique();                         // remove duplicates
l.merge(l2);                        // merge two sorted lists
```

== Pairs and tuples
```cpp
pair<int,int> p = {1,2};
cout << p.first << ' ' << p.second << '\n';
pair<int,int> p2 = {3,4};
p.swap(p2);
tuple<int,int,int> t = {1,2,3};
cout << get<0>(t) << get<1>(t) << get<2>(t) << '\n';
```

== Map
Map: $O log n$

Unordered Map: $O 1$

```cpp
using pss = pair<string,string>;
int main() {
  map<string,string> m;
  // unordered_map<string,string> m;
  m.insert(pss("key", "value")); // insert
  string value = m.at("key"); // access (or m["key"])
  m.erase("key"); // remove
  bool exists = (m.find("key") != m.end()); // find
					    // iterate
  for(auto &it: m) {
    cout << it.first << ' ' << it.second;
  }
}
```

== Set
Set: $O log n$

Unordered Set: $O 1$
```cpp
int main() {
  set<int> s;
  // unordered_set<int> s;
  int value = 1;
  s.insert(value); // insert
  s.erase(value); // remove
  bool exists = (s.find(value) != s.end()); // find
					    // iterate
  for(auto &i: s) {
    cout << i;
  }
  cout << '\n';
}
```


= Math

$
  x = frac(- b plus.minus sqrt(b^2 - 4 a c) , 2 a)\
  sin(v + w) = sin(v) cos(w) + cos(v)sin(w)\
  cos(v + w) = cos(v)cos(w) + sin(v)sin(w)
$
== Law of sines
$
  sin(A) / a = sin(B)/b = sin(C)/c
$

== Law of cosines
$
  a^2 = b^2 + c^2 - 2 b c cos(alpha)
$

== Complex numbers
$
  r = sqrt(a^2 + b^2)\
  theta = arctan(b/a)\
  a = r cos theta\
  b = r sin theta\
  z^n = r^n (cos(n theta) + i sin(n theta) = r e^(i n theta))\
  root(n, r) e^(i (frac(theta + 2 pi k, n)))
$

== Permutations
$
  n "elements of length" r\
  "with repititions" n^r\
  "without repititions" frac(n!, (n - r)!)\
  "The number of distinguishable permutations that can be formed from a collection of"\ n "objects,"
  " where object" o_i "appears" k_i "times is:"\
  frac(n!,k_1! k_2! k_3! ...)
$

e.g. MISSISSIPPI: 1M 4I 4S 2P
$ n = 1 + 4 + 4 + 2 = 11\ 11!/( 1!4!4!2! ) = 34 650 $

== Combinations (order doesn t matter)

$ ""_n C _r = frac(n!, r!(n - r)!) = vec(n, r) $

with repetitions:

$ ""_(n + r - 1)C_r $

== Sequences
$
  A_n = n/2(2a_1 + (n-1)d)
$
$
  G_n = a_1(1-r^n)/(1-r)
$

= Number Theory
== GCD
```cpp
template <typename T>
T gcd(T a, T b) {
  if (b == 0) {
    return a;
  }
  return gcd(b, a % b);
}
```
Iter
```cpp
template <typename T>
T gcd(T a, T b) {
  T result = min(a, b);
  while (result > 0) {
    if (a % result == 0 && b % result == 0) {
        break;
    }
    result--;
  }
  return result;
}
```

== LCM

$ a * b = gcd(a, b) * "lcm"(a, b) $
```cpp
T lcm(T a, Tb) {
  return (a*b) / gcd(a, b);
}
```

== Primes
```cpp
bool isPrime(int n) {
  if (n < 2) {
    return false;
  }
  for (int x = 2; x * x <= n; ++x) {
    if (n % x == 0) {
      return false;
    }
  }
  return true;
}

vector<int> primeFactorisation(int n) {
  vector<int> facts;
  for (int x = 2; x * x <= n; ++x) {
    while (n % x == 0) {
      facts.push_back(x);
      n /= x;
    }
  }
  if (n > 1) {
    facts.push_back(n);
  }
  return facts;
}
```

== Complex Numbers

```cpp
constexpr double pi = 3.1415926535897932384626433; // or std::acos(-1)
struct Complex { using T = Complex; double u,v;
  Complex(double u=0, double v=0) : u{u}, v{v} {}
  T operator+(T r) const { return {u+r.u, v+r.v}; }
  T operator-(T r) const { return {u-r.u, v-r.v}; }
  T operator*(T r) const { return {u*r.u - v*r.v, u*r.v + v*r.u}; }
  T operator/(T r) const {
    auto norm = r.u*r.u+r.v*r.v;
    return {(u*r.u + v*r.v)/norm, (v*r.u - u*r.v)/norm};
  }
  T operator*(double r) const { return T{u*r, v*r}; }
  T operator/(double r) const { return T{u/r, v/r}; }
  T inv() const { return T{1,0}/ *this; }
  T conj() const { return T{u, -v}; }
  static T root(ll k){ return {cos(2*pi/k), sin(2*pi/k)}; }
  bool zero() const { return max(abs(u), abs(v)) < 1e-6; }
};
```
== Matrix

```cpp
#define REP(i, n) for(auto i = decltype(n)(0); i < (n); i++)
using T = double;
constexpr T EPS = 1e-8;
template<int R, int C>
using M = array<array<T,C>,R>;	// matrix
template<int R, int C>
T ReducedRowEchelonForm(M<R,C> &m, int rows) {	// return the determinant
  int r = 0; T det = 1;							// MODIFIES the input
  for(int c = 0; c < rows && r < rows; c++) {
    int p = r;
    for(int i=r+1; i<rows; i++) if(abs(m[i][c]) > abs(m[p][c])) p=i;
    if(abs(m[p][c]) < EPS){	det = 0; continue; }
    swap(m[p], m[r]);		det = -det;
    T s = 1.0 / m[r][c], t;	det *= m[r][c];
    REP(j,C) m[r][j] *= s;				// make leading term in row 1
    REP(i,rows) if (i!=r){ t = m[i][c]; REP(j,C) m[i][j] -= t*m[r][j]; }
    ++r;
  }
  return det;
}
bool error, inconst;	// error => multiple or inconsistent
template<int R,int C>	// Mx = a; M:R*R, v:R*C => x:R*C
M<R,C> solve(const M<R,R> &m, const M<R,C> &a, int rows){
  M<R,R+C> q;
  REP(r,rows){
    REP(c,rows) q[r][c] = m[r][c];
    REP(c,C) q[r][R+c] = a[r][c];
  }
  ReducedRowEchelonForm<R,R+C>(q,rows);
  M<R,C> sol; error = false, inconst = false; 
  REP(c,C) for(auto j = rows-1; j >= 0; --j){
    T t=0; bool allzero=true;
    for(auto k = j+1; k < rows; ++k)
      t += q[j][k]*sol[k][c], allzero &= abs(q[j][k]) < EPS;
    if(abs(q[j][j]) < EPS)
      error = true, inconst |= allzero && abs(q[j][R+c]) > EPS;
    else sol[j][c] = (q[j][R+c] - t) / q[j][j]; // usually q[j][j]=1
  }
  return sol;
}
```

= Sorting and Searching

```cpp
// Binary predicate
int compare(const void* ap, const void* bp) {
  // Typecasting
  const int* a = (int*)ap;
  const int* b = (int*)bp;
  if (*a < *b)
    return -1;
  else if (*a > *b)
    return 1;
  else
    return 0;
}
int main() {
  vector<int> vec;
  auto first = vec.begin();
  auto last = vec.end();
  sort(first, last);
  // reverse: sort(vec.rbegin(), vec.rend())
  // sorting arrays: a[]; sort(a, a+n)
  int key = 1;
  bool isPresent = binary_search(first, last, key);
  int* p1 = (int*)bsearch(&key, vec.data(), vec.size(), sizeof(int), compare)
}
```

= Graph
== Dijkstra

```cpp
// Time Complexity: O(E log V)
// Space Complexity: O(V + E)
#define ll long long
#define fi first
#define se second
typedef pair<ll, ll> ii;
int main() {
  ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
  ll v, e; cin >> v >> e;
  vector<ii> adj[v];
  ll dist[v];
  for (ll i = 0; i < v; i++) dist[i] = INT64_MAX;
  for (ll i = 0; i < e; i++) {
    ll u, v, wi; cin >> u >> v >> wi; u--, v--;
    adj[u].push_back(ii(v, wi));
    adj[v].push_back(ii(u, wi));
  }
  priority_queue<ii> q; // (dist, node)
  q.push(ii(-0, 0));
  while(!q.empty()) {
    ii cur = q.top(); q.pop();
    ll curDist = -cur.fi, curNode = cur.se;
    if (curDist > dist[curNode]) continue;
    for (ll i = 0; i < adj[curNode].size(); i++) {
      ii nx = adj[curNode][i];
      ll nxDist = curDist + nx.se;
      if (nxDist < dist[nx.fi]) {
        dist[nx.fi] = nxDist;
        q.push(ii(-nxDist, nx.fi));
      }
    }   
  }
  cout << dist[v-1] << '\n';
}
```

== Floyd Warshall

```cpp
int main() {
  ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
  ll t; cin >> t;
  while (t--) {
    string s; cin >> s;
    ll n = s.size();
    ll adj[n][n]; adj[0][0] = 0;
    for (ll j=1; j<n; j++) adj[0][j] = (s[j]=='Y'?1:1e9);
    for (ll i=1; i<n; i++) {
      cin >> s;
      for (ll j=0; j<n; j++) {
        if (i == j) adj[i][j] = 0;
        else adj[i][j] = (s[j]=='Y'?1:1e9);
      }
    }
    for (ll k=0; k<n; k++) {
      for (ll i=0; i<n; i++) {
        for (ll j=0; j<n; j++) {
          adj[i][j] = min(adj[i][j], adj[i][k] + adj[k][j]);
        }
      }
    }
    ll mxVal = 0, mxInd = 0;
    for (ll i=0; i<n; i++) {
      ll ans = 0;
      for (ll j=0; j<n; j++) {
        if (adj[i][j] == 2) ans++;
      }
      if (ans > mxVal) {
        mxVal = ans;
        mxInd = i;
      }
    }
    cout << mxInd << ' ' << mxVal << '\n';
  }
}
```

== Breadth First Search

```cpp
void bfs(vector<vector<int>>& adj, int s, vector<bool>& visited) {
  queue<int> q;
  visited[s] = true;
  q.push(s);
  while (!q.empty()) {
    int curr = q.front();
    q.pop();
    for (int x : adj[curr]) {
      if (!visited[x]) {
        visited[x] = true;
        q.push(x);
      }
    }
  }
}
void addEdge(vector<vector<int>>& adj, int u, int v) {
  adj[u].push_back(v);
  adj[v].push_back(u); // Undirected Graph
}

int main() {
  // Number of vertices in the graph
  int V = 5;
  // Adjacency list representation of the graph
  vector<vector<int>> adj(V);
  addEdge(adj, 0, 1);
  addEdge(adj, 0, 2);
  addEdge(adj, 1, 3);
  addEdge(adj, 1, 4);
  addEdge(adj, 2, 4);
  // Mark all the vertices as not visited
  vector<bool> visited(V, false);
  // Perform BFS traversal starting from vertex 0
  bfs(adj, 0, visited);
  return 0;
}
```
== Depth First Search

```cpp
// Function to add an edge to the adjacency list
void addEdge(vector<vector<int>> &adj, int s, int t){
  // Add edge from vertex s to t
  adj[s].push_back(t); 
  // Due to undirected Graph
  adj[t].push_back(s);
}

// Recursive function for DFS traversal
void DFSRec(vector<vector<int>> &adj, vector<bool> &visited, int s){
  // Mark the current vertex as visited
  visited[s] = true;
  // Print the current vertex
  cout << s << " ";
  // Recursively visit all adjacent vertices that are not visited yet
  for (int i : adj[s])
    if (visited[i] == false)
      DFSRec(adj, visited, i);
}

// Main DFS function that initializes the visited array and call DFSRec
void DFS(vector<vector<int>> &adj, int s){
  vector<bool> visited(adj.size(), false);
  // Call the recursive DFS function
  DFSRec(adj, visited, s);
}

int main(){
  int V = 5; // Number of vertices in the graph
  // Create an adjacency list for the graph
  vector<vector<int>> adj(V);
  // Define the edges of the graph
  vector<vector<int>> edges={{1, 2},{1, 0},{2, 0},{2, 3},{2, 4}};
  // Populate the adjacency list with edges
  for (auto &e : edges)
      addEdge(adj, e[0], e[1]);
  int source = 1; // Starting vertex for DFS
  cout << "DFS from source: " << source << endl;
  DFS(adj, source); // Perform DFS starting from the source vertex
  return 0;
}
```
= Data Srtuctures
== Trie

```cpp
const int ALPHABET_SIZE = 26;
inline int mp(char c) { return c - 'a'; }
struct Node {
	Node* ch[ALPHABET_SIZE];
	bool isleaf = false;
	Node() {
		for(int i = 0; i < ALPHABET_SIZE; ++i) ch[i] = nullptr;
	}
	void insert(string &s, int i = 0) {
		if (i == s.length()) isleaf = true;
		else {
			int v = mp(s[i]);
			if (ch[v] == nullptr)
				ch[v] = new Node();
			ch[v]->insert(s, i + 1);
		}
	}
	bool contains(string &s, int i = 0) {
		if (i == s.length()) return isleaf;
		else {
			int v = mp(s[i]);
			if (ch[v] == nullptr) return false;
			else return ch[v]->contains(s, i + 1);
		}
	}
	void cleanup() {
		for (int i = 0; i < ALPHABET_SIZE; ++i)
			if (ch[i] != nullptr) {
				ch[i]->cleanup();
				delete ch[i];
			}
	}
};
```

== Fenwick
```cpp
template <class T>
struct FenwickTree {		// use 1 based indices!!!
	int n; vector<T> tree;
	FenwickTree(int n) : n(n) { tree.assign(n + 1, 0); }
	T query(int l, int r) { return query(r) - query(l - 1); }
	T query(int r) {
		T s = 0;
		for(; r > 0; r -= (r & (-r))) s += tree[r];
		return s;
	}
	void update(int i, T v) {
		for(; i <= n; i += (i & (-i))) tree[i] += v;
	}
};
```

== Heap
```cpp
vector<int> v1 = { 20, 30, 40, 25, 15 };
// Converting vector into a heap
make_heap(v1.begin(), v1.end());
// Displaying the maximum element of heap using front()
cout << v1.front() << endl;
// add element
v1.push_back(1);
push_heap(v1.begin(), v1.end());
// remove element
// using pop_heap() function to move the largest element to the end
pop_heap(v1.begin(), v1.end());
// actually removing the element from the heap using pop_back()
v1.pop_back();
```

== Linked List
```cpp
list<int> l = {1, 2, 3};
l.push_front(0);
l.push_back(4);
l.unique();
```


== Search Tree
```cpp
template <typename T> struct Node {
   T data;
   Node *left;
   Node *right;
};
template <typename T> Node<T> *newNode(T data) {
   Node<T> *n = new Node<T>;
   n->data = data;
   n->left = n->right = nullptr;
   return n;
}
template <typename T> Node<T> *insertNode(Node<T> *root, int data) {
   if (root == nullptr) {
      return newNode(data);
   }
   if (data < root->data) {
      root->left = insertNode(root->left, data);
   } else if (data > root->data) {
      root->right = insertNode(root->right, data);
   }
   return root;
}
template <typename T> Node<T> *searchNode(Node<T> *root, int key) {
   if (root == nullptr || root->data == key) {
      return root;
   }
   if (root->data < key) {
      return searchNode(root->right, key);
   }
   return searchNode(root->left, key);
}
template <typename T> Node<T> *minValueNode(Node<T> *node) {
   Node<T> *current = node;
   while (current && current->left != nullptr) {
      current = current->left;
   }
   return current;
}
template <typename T> Node<T> *deleteNode(Node<T> *root, int data) {
   if (root == nullptr)
      return root;
   if (data < root->data) {
      root->left = deleteNode(root->left, data);
   } else if (data > root->data) {
      root->right = deleteNode(root->right, data);
   } else {
      if (root->left == nullptr) {
         Node<T> *temp = root->right;
         delete root;
         return temp;
      } else if (root->right == nullptr) {
         Node<T> *temp = root->left;
         delete root;
         return temp;
      }
      Node<T> *temp = minValueNode(root->right);
      root->data = temp->data;
      root->right = deleteNode(root->right, temp->data);
   }
   return root;
}
```


= Libraries
== BigInt

```cpp
class BigInt{
 string digits;
 public:
    //Constructors:
    BigInt(unsigned long long n = 0);
    BigInt(string &);
    BigInt(const char *);
    BigInt(BigInt &);
    //Helper Functions:
    friend void divide_by_2(BigInt &a);
    friend bool Null(const BigInt &);
    friend int Length(const BigInt &);
    int operator[](const int)const;
    /* * * * Operator Overloading * * * */
    //Direct assignment
    BigInt &operator=(const BigInt &);
    //Post/Pre - Incrementation
    BigInt &operator++();
    BigInt operator++(int temp);
    BigInt &operator--();
    BigInt operator--(int temp);
    //Addition and Subtraction
    friend BigInt &operator+=(BigInt &, const BigInt &);
    friend BigInt operator+(const BigInt &, const BigInt &);
    friend BigInt operator-(const BigInt &, const BigInt &);
    friend BigInt &operator-=(BigInt &, const BigInt &);
    //Comparison operators
    friend bool operator==(const BigInt &, const BigInt &);
    friend bool operator!=(const BigInt &, const BigInt &);
    friend bool operator>(const BigInt &, const BigInt &);
    friend bool operator>=(const BigInt &, const BigInt &);
    friend bool operator<(const BigInt &, const BigInt &);
    friend bool operator<=(const BigInt &, const BigInt &);
    //Multiplication and Division
    friend BigInt &operator*=(BigInt &, const BigInt &);
    friend BigInt operator*(const BigInt &, const BigInt &);
    friend BigInt &operator/=(BigInt &, const BigInt &);
    friend BigInt operator/(const BigInt &, const BigInt &);
    //Modulo
    friend BigInt operator%(const BigInt &, const BigInt &);
    friend BigInt &operator%=(BigInt &, const BigInt &);
    //Power Function
    friend BigInt &operator^=(BigInt &,const BigInt &);
    friend BigInt operator^(BigInt &, const BigInt &);
    //Square Root Function
    friend BigInt sqrt(BigInt &a);
    //Read and Write
    friend ostream &operator<<(ostream &,const BigInt &);
    friend istream &operator>>(istream &, BigInt &);
};

BigInt::BigInt(string & s){
  digits = "";
  int n = s.size();
  for (int i = n - 1; i >= 0;i--){
    if(!isdigit(s[i]))
      throw("ERROR");
    digits.push_back(s[i] - '0');
  }
}
BigInt::BigInt(unsigned long long nr){
  do{
    digits.push_back(nr % 10);
    nr /= 10;
  } while (nr);
}
BigInt::BigInt(const char *s){
  digits = "";
  for (int i = strlen(s) - 1; i >= 0;i--)
  {
    if(!isdigit(s[i]))
      throw("ERROR");
    digits.push_back(s[i] - '0');
  }
}
BigInt::BigInt(BigInt & a){
  digits = a.digits;
}
 
bool Null(const BigInt& a){
  if(a.digits.size() == 1 && a.digits[0] == 0)
    return true;
  return false;
}
int Length(const BigInt & a){
  return a.digits.size();
}
int BigInt::operator[](const int index)const{
  if(digits.size() <= index || index < 0)
    throw("ERROR");
  return digits[index];
}
bool operator==(const BigInt &a, const BigInt &b){
  return a.digits == b.digits;
}
bool operator!=(const BigInt & a,const BigInt &b){
  return !(a == b);
}
bool operator<(const BigInt&a,const BigInt&b){
  int n = Length(a), m = Length(b);
  if(n != m) return n < m;
  while(n--)
    if(a.digits[n] != b.digits[n])
      return a.digits[n] < b.digits[n];
  return false;
}
bool operator>(const BigInt&a,const BigInt&b){
  return b < a;
}
bool operator>=(const BigInt&a,const BigInt&b){
  return !(a < b);
}
bool operator<=(const BigInt&a,const BigInt&b){
  return !(a > b);
}
 
BigInt &BigInt::operator=(const BigInt &a){
  digits = a.digits;
  return *this;
}
 
BigInt &BigInt::operator++(){
  int i, n = digits.size();
  for (i = 0; i < n && digits[i] == 9;i++)
    digits[i] = 0;
  if(i == n)
    digits.push_back(1);
  else
     digits[i]++;
  return *this;
}
BigInt BigInt::operator++(int temp){
  BigInt aux;
  aux = *this;
  ++(*this);
  return aux;
}
 
BigInt &BigInt::operator--(){
  if(digits[0] == 0 && digits.size() == 1)
    throw("UNDERFLOW");
  int i, n = digits.size();
  for (i = 0; digits[i] == 0 && i < n;i++)
    digits[i] = 9;
  digits[i]--;
  if(n > 1 && digits[n - 1] == 0)
    digits.pop_back();
  return *this;
}
BigInt BigInt::operator--(int temp){
  BigInt aux;
  aux = *this;
  --(*this);
  return aux;
}
 
BigInt &operator+=(BigInt &a,const BigInt& b){
  int t = 0, s, i;
  int n = Length(a), m = Length(b);
  if(m > n)
    a.digits.append(m - n, 0);
  n = Length(a);
  for (i = 0; i < n;i++){
    if(i < m)
      s = (a.digits[i] + b.digits[i]) + t;
    else
      s = a.digits[i] + t;
    t = s / 10;
    a.digits[i] = (s % 10);
  }
  if(t)
    a.digits.push_back(t);
  return a;
}
BigInt operator+(const BigInt &a, const BigInt &b){
  BigInt temp;
  temp = a;
  temp += b;
  return temp;
}
 
BigInt &operator-=(BigInt&a,const BigInt &b){
   if(a < b)
     throw("UNDERFLOW");
   int n = Length(a), m = Length(b);
   int i, t = 0, s;
   for (i = 0; i < n;i++){
     if(i < m)
       s = a.digits[i] - b.digits[i]+ t;
     else
       s = a.digits[i]+ t;
     if(s < 0)
       s += 10,
       t = -1;
     else
       t = 0;
     a.digits[i] = s;
   }
   while(n > 1 && a.digits[n - 1] == 0)
     a.digits.pop_back(),
     n--;
   return a;
}
BigInt operator-(const BigInt& a,const BigInt&b){
  BigInt temp;
  temp = a;
  temp -= b;
  return temp;
}
 
BigInt &operator*=(BigInt &a, const BigInt &b)
{
  if(Null(a) || Null(b)){
     a = BigInt();
     return a;
  }
  int n = a.digits.size(), m = b.digits.size();
  vector<int> v(n + m, 0);
  for (int i = 0; i < n;i++)
    for (int j = 0; j < m;j++){
      v[i + j] += (a.digits[i] ) * (b.digits[j]);
    }
  n += m;
  a.digits.resize(v.size());
  for (int s, i = 0, t = 0; i < n; i++) {
    s = t + v[i];
    v[i] = s % 10;
    t = s / 10;
    a.digits[i] = v[i] ;
  }
  for (int i = n - 1; i >= 1 && !v[i];i--)
     a.digits.pop_back();
  return a;
}
BigInt operator*(const BigInt&a,const BigInt&b){
  BigInt temp;
  temp = a;
  temp *= b;
  return temp;
}
 
BigInt &operator/=(BigInt& a,const BigInt &b){
  if(Null(b))
    throw("Arithmetic Error: Division By 0");
  if(a < b){
    a = BigInt();
    return a;
  }
  if(a == b){
    a = BigInt(1);
    return a;
  }
  int i, lgcat = 0, cc;
  int n = Length(a), m = Length(b);
  vector<int> cat(n, 0);
  BigInt t;
  for (i = n - 1; t * 10 + a.digits[i]  < b;i--){
    t *= 10;
    t += a.digits[i] ;
  }
  for (; i >= 0; i--){
    t = t * 10 + a.digits[i];
    for (cc = 9; cc * b > t;cc--);
    t -= cc * b;
    cat[lgcat++] = cc;
  }
  a.digits.resize(cat.size());
  for (i = 0; i < lgcat;i++)
    a.digits[i] = cat[lgcat - i - 1];
  a.digits.resize(lgcat);
  return a;
}
BigInt operator/(const BigInt &a,const BigInt &b){
  BigInt temp;
  temp = a;
  temp /= b;
  return temp;
}
 
BigInt &operator%=(BigInt& a,const BigInt &b){
  if(Null(b))
    throw("Arithmetic Error: Division By 0");
  if(a < b){
    return a;
  }
  if(a == b){
    a = BigInt();
    return a;
  }
  int i, lgcat = 0, cc;
  int n = Length(a), m = Length(b);
  vector<int> cat(n, 0);
  BigInt t;
  for (i = n - 1; t * 10 + a.digits[i] < b;i--){
    t *= 10;
    t += a.digits[i];
  }
  for (; i >= 0; i--){
    t = t * 10 + a.digits[i];
    for (cc = 9; cc * b > t;cc--);
    t -= cc * b;
    cat[lgcat++] = cc;
  }
  a = t;
  return a;
}
BigInt operator%(const BigInt &a,const BigInt &b){
  BigInt temp;
  temp = a;
  temp %= b;
  return temp;
}
 
BigInt &operator^=(BigInt & a,const BigInt & b){
  BigInt Exponent, Base(a);
  Exponent = b;
  a = 1;
  while(!Null(Exponent)){
    if(Exponent[0] & 1)
      a *= Base;
    Base *= Base;
    divide_by_2(Exponent);
  }
  return a;
}
BigInt operator^(BigInt & a,BigInt & b){
  BigInt temp(a);
  temp ^= b;
  return temp;
}
 
void divide_by_2(BigInt & a){
  int add = 0;
  for (int i = a.digits.size() - 1; i >= 0;i--){
    int digit = (a.digits[i] >> 1) + add;
    add = ((a.digits[i] & 1) * 5);
    a.digits[i] = digit;
  }
  while(a.digits.size() > 1 && !a.digits.back())
    a.digits.pop_back();
}
 
BigInt sqrt(BigInt & a){
  BigInt left(1), right(a), v(1), mid, prod;
  divide_by_2(right);
  while(left <= right){
    mid += left;
    mid += right;
    divide_by_2(mid);
    prod = (mid * mid);
    if(prod <= a){
      v = mid;
      ++mid;
      left = mid;
    }
    else{
      --mid;
      right = mid;
    }
    mid = BigInt();
  }
  return v;
}
 
istream &operator>>(istream &in,BigInt&a){
  string s;
  in >> s;
  int n = s.size();
  for (int i = n - 1; i >= 0;i--){
    if(!isdigit(s[i]))
      throw("INVALID NUMBER");
    a.digits[n - i - 1] = s[i];
  }
  return in;
}
 
ostream &operator<<(ostream &out,const BigInt &a){
  for (int i = a.digits.size() - 1; i >= 0;i--)
    cout << (short)a.digits[i];
  return cout;
}
```

= Pseudo Code
#pseudocode-list(
   indentation-guide-stroke: .5pt + gray.darken(20%),
)[
  + *Dijkstra (G, v)*
  - *input* connected weighted graph $G$ with node $v$
  - *output* function $d$ yielding for every node the length of a shortest path to $v$
  -  $S <- "nodes"(G)$ #comment[ initialize ToDo list ]
  + *For All* {$u in "nodes"(G)$}
    + $d[u] <-$ *if* $u = v$ *then* 0 *else* $infinity$
    + #comment[initialize $d$]
  + *While* $S$ is not empty
    + $u <-$ node in $S$ such that $d[u]$ is minimal
    + remove $u$ from $S$
    + *For All* $z \in S$ with $(u, z) in "edges"(G)$
      + $d[u] <- min(d[z], d[u] + "weight"[u][z])$
  + *Return* d
]


