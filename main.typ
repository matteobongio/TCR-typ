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


= Template

```cpp
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using ld = long double;

int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0); cout.tie(0);
}
```

