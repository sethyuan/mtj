# cav/mtj

A [Matrix Toolkits Java](https://github.com/fommil/matrix-toolkits-java) clojure wrapper and core.matrix implementation.

The goal is to implement as much core.matrix protocols as possible and, provide extra API for MTJ features that core.matrix does not have (yet).

## Installation

Not released to Clojars yet, because much implementation is still missing. It'll take some time before I can finish all necessary implementations.

## Example

Just use it as a core.matrix implementation, like this:

```clojure
(require '[clojure.core.matrix :as m])

(m/set-current-implementation :mtj)

(m/matrix [1 2 3])
...
```

## API

Extra features API that core.matrix does not have will go here.

## License

(The MIT License)

Copyright (c) 2014 Seth Yuan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
