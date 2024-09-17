#import "@preview/showybox:2.0.1": showybox

#import "@preview/codly:0.2.0": *
#show: codly-init.with()

#import "@preview/lovelace:0.2.0": *
#show: setup-lovelace

#set page(
  paper: "a4"
)
#set page(numbering: "1")
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
#set page(
  header: [
    #grid(columns: (50%, 50%),
        grid.cell( "Linear Algebra for Computer Science", ),
        grid.cell( "Matteo Bongiovanni (S5560349)",align: right ),
      )
  Homework 4
  ]
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
#define FOR(a, c) for (int(a) = 0; (a) < (c); (a)++) 
#define FAST_IO                     \
  ios_base::sync_with_stdio(false); \
  cin.tie(0);                       \
  cout.tie(0);

int main() {
}
```

= Number Theory
== GCD
```cpp
template <typename T>
inline T gcd(T a, T b) {
  if (b == 0) {
    return a;
  }
  return gcd(b, a % b);
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

vector<int> getFactors(int n) {
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

= Graph
== Dijkstra

```cpp
// Time Complexity: O(E log V)
// Space Complexity: O(V + E)
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define fi first
#define se second

typedef pair<int, int> ii;
signed main() {
    ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);
    int v, e; cin >> v >> e;

    vector<ii> adj[v];
    int dist[v];

    for (int i = 0; i < v; i++) dist[i] = INT64_MAX;
    for (int i = 0; i < e; i++) {
        int u, v, wi; cin >> u >> v >> wi; u--, v--;
        adj[u].push_back(ii(v, wi));
        adj[v].push_back(ii(u, wi));
    }

    priority_queue<ii> q; // (dist, node)
    q.push(ii(-0, 0));
    while(!q.empty()) {
        ii cur = q.top(); q.pop();
        int curDist = -cur.fi, curNode = cur.se;

        if (curDist > dist[curNode]) continue;

        for (int i = 0; i < adj[curNode].size(); i++) {
            ii nx = adj[curNode][i];
            int nxDist = curDist + nx.se;
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

=======
== Trie

```cpp
string lower(string s) {
   string out = "";
   for (char c : s) {
      out += tolower(c);
   }
   return out;
}
class TrieNode {
 public:
   TrieNode *children[26];
   bool end_of_word;
   char letter;
   TrieNode() {
      end_of_word = false;
      for (int i = 0; i < 26; i++) {
         children[i] = NULL;
      }
      letter = '\0';
   }
};

class Trie {
 public:
   TrieNode root;
   void Insert(string str) {
     str = lower(str);
      TrieNode *current = &root;
      for (size_t i = 0; i < str.size(); i++) {
         if (current->children[str.at(i) - 'a'] == NULL) {
            current->children[str.at(i) - 'a'] = new TrieNode;
            current->children[str.at(i) - 'a']->letter = str.at(i);
         }
         current = current->children[str.at(i) - 'a'];
      }
      current->end_of_word = true;
   }
   TrieNode *Search(string str) {
     str = lower(str);
      TrieNode *current = &root;
      for (size_t i = 0; i < str.size(); i++) {
         if (current->children[str.at(i) - 'a']) {
            current = current->children[str.at(i) - 'a'];
         } else {
            current = NULL;
            break;
         }
      }
      return current;
   }
};

```

= Libraries
== BigInt

```cpp
#include <bits/stdc++.h>
 
using namespace std;
 
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
    if(n != m)
        return n < m;
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
    for (int s, i = 0, t = 0; i < n; i++)
    {
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

