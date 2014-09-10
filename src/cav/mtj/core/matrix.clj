(ns cav.mtj.core.matrix
  "Wrapping MTJ for core.matrix."
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.protocols :refer :all]
            [clojure.core.matrix.implementations :refer [register-implementation]]
            [clojure.core.matrix.utils :refer [error]])
  (:import [no.uib.cipr.matrix
            Matrix Vector DenseVector DenseMatrix Matrices
            BandMatrix LowerSymmBandMatrix
            PermutationMatrix]))

(def imp
  (reify
    PImplementation
    (implementation-key [m] :mtj)
    (meta-info [m] {:doc "MTJ as an implementation of core.matrix."})
    (supports-dimensionality? [m dimensions]
      (or (= dimensions 1) (= dimensions 2)))
    (new-vector [m length] (DenseVector. (int length)))
    (new-matrix [m rows columns] (DenseMatrix. (int rows) (int columns)))
    (new-matrix-nd [m shape]
      (case (count shape)
        1 (new-vector m (first shape))
        2 (new-matrix m (first shape) (second shape))
        nil))
    (construct-matrix [m data]
      (cond
        (number? data) (double data)
        (instance? Vector data) data
        (instance? Matrix data) data   
        (instance? (Class/forName "[D") data) (DenseVector. ^doubles data)
        (instance? (Class/forName "[[D") data) (DenseMatrix. ^"[[D" data)
        (is-vector? data)
          (let [^long size (dimension-count data 0)
                v (DenseVector. size)]
            (dotimes [i size]
              (.set v i ^double (get-1d data i)))
            v)
        (sequential? data)
          (let [mat (DenseMatrix. (count data) (count (first data)))]
            (loop [i 0, rows data]
              (when-let [r (first rows)]
                (loop [j 0, vs r]
                  (when-let [v (first vs)]
                    (.set mat i j v)
                    (recur (inc j) (rest vs))))
                (recur (inc i) (rest rows))))
            mat)))

    PRowColMatrix
    (column-matrix [m data]
      (cond
        (instance? Vector data) (DenseMatrix. ^Vector data)
        (is-vector? data)
          (let [^long size (dimension-count data 0)
                mat (DenseMatrix. size 1)]
            (dotimes [i size]
              (.set mat i 0 (get-1d data i)))
            mat)
        :else (error "data has to be a 1D vector.")))
    (row-matrix [m data]
      (if (is-vector? data)
        (let [^long size (dimension-count data 0)
              mat (DenseMatrix. 1 size)]
          (dotimes [i size]
            (.set mat 0 i (get-1d data i)))
          mat)
        (error "data has to be a 1D vector.")))

    PSparse
      (sparse-coerce [m data] nil)

    PDense
      (dense-coerce [m data] data)

    PCoercion
      (coerce-param [m data]
        (cond
          (instance? Vector data) data
          (instance? Matrix data) data
          :else nil))

    PSpecialisedConstructors
      (identity-matrix [m dims] (Matrices/identity (int dims)))
      (diagonal-matrix [m diagonal-values]
        (let [^long size (dimension-count diagonal-values 0)
              mat (LowerSymmBandMatrix. size 0)]
          (dotimes [i size]
            (.set mat i i (get-1d diagonal-values i)))
          mat))

    PPermutationMatrix
      (permutation-matrix [m permutation]
        (PermutationMatrix. (int-array permutation)))))

(extend-protocol PDimensionInfo
  Vector
  (dimensionality [v] 1)
  (get-shape [v] [(.size v)])
  (is-scalar? [v] false)
  (is-vector? [v] true)
  (dimension-count [v ^long dimension-number]
    (case dimension-number 
      0 (.size v)
      (error "Illegal dimension number.")))

  Matrix
  (dimensionality [m] 2)
  (get-shape [m] [(.numRows m) (.numColumns m)])
  (is-scalar? [m] false)
  (is-vector? [m] false)
  (dimension-count [m ^long dimension-number]
    (case dimension-number 
      0 (.numRows m)
      1 (.numColumns m)
      (error "Illegal dimension number."))))

(extend-protocol PIndexedAccess
  Vector
  (get-1d [m row] (.get m (int row)))
  (get-2d [m row column] (error "1D vector must be accessed by only 1 index."))
  (get-nd [m indexes]
    (if (= 1 (count indexes))
      (get-1d m (first indexes))
      (error "1D vector must be accessed by only 1 index.")))

  Matrix
  (get-1d [m row] (error "2D matrix must be accessed by row and column."))
  (get-2d [m row column] (.get m (int row) (int column)))
  (get-nd [m indexes]
    (if (= 2 (count indexes))
      (get-2d m (first indexes) (second indexes))
      (error "2D matrix must be accessed by row and column."))))

(extend-protocol PIndexedSetting
  Vector
  (set-1d [m row v]
    (doto (.copy m)
      (.set (int row) (double v))))
  (set-2d [m row column v] (error "1D vector must be set with only 1 index."))
  (set-nd [m indexes v]
    (if (= 1 (count indexes))
      (set-1d m (first indexes) v)
      (error "1D vector must be set with only 1 index.")))
  (is-mutable? [m] true)

  Matrix
  (set-1d [m row v] (error "2D matrix must be set with 2 indexes."))
  (set-2d [m row column v]
    (doto (.copy m)
      (.set (int row) (int column) (double v))))
  (set-nd [m indexes v]
    (if (= 2 (count indexes))
      (set-2d m (first indexes) (second indexes) v)
      (error "2D matrix must be set with 2 indexes.")))
  (is-mutable? [m] true))

(extend-protocol PIndexedSettingMutable
  Vector
  (set-1d! [m row v] (.set m (int row) (double v)))
  (set-2d! [m row column v] (error "1D vector must be set with only 1 index."))
  (set-nd! [m indexes v]
    (if (= 1 (count indexes))
      (set-1d m (first indexes) v)
      (error "1D vector must be set with only 1 index.")))
  (is-mutable? [m] true)

  Matrix
  (set-1d! [m row v] (error "2D matrix must be set with 2 indexes."))
  (set-2d! [m row column v] (.set m (int row) (int column) (double v)))
  (set-nd! [m indexes v]
    (if (= 2 (count indexes))
      (set-2d m (first indexes) (second indexes) v)
      (error "2D matrix must be set with 2 indexes.")))
  (is-mutable? [m] true))

(extend-protocol PMatrixCloning
  Vector
  (clone [m] (.copy m))
  
  Matrix
  (clone [m] (.copy m)))

(extend-protocol PTypeInfo
  Vector
  (element-type [m] Double/TYPE)

  Matrix
  (element-type [m] Double/TYPE))

(extend-protocol PArrayMetrics
  Vector
  (nonzero-count [m] (Matrices/cardinality m))

  Matrix
  (nonzero-count [m] (Matrices/cardinality m)))

(extend-protocol PMutableMatrixConstruction
  Vector
  (mutable-matrix [m] m)

  Matrix
  (mutable-matrix [m] m))

(extend-protocol PSameShape
  Vector
  (same-shape? [a b]
    (and (instance? Vector b) (= (.size a) (.size ^Vector b))))

  Matrix
  (same-shape? [a b]
    (and (instance? Matrix b)
         (let [b ^Matrix b]
           (and (= (.numRows a) (.numRows b))
                (= (.numColumns a) (.numColumns b)))))))

(extend-protocol PMatrixSlices
  Vector
  (get-row [m i] (error "Cannot access row of 1D vector."))
  (get-column [m i] (error "Cannot access column of 1D vector."))
  (get-major-slice [m i] (.get m (int i)))
  (get-slice [m dimension i]
    (if (= dimension 0)
      (.get m (int i))
      (error "Only dimension 0 is supported for 1D vector.")))

  Matrix
  (get-row [m i]
    ;; Alternative that returns a new DenseVector each time.
    ;; (let [row (DenseVector. (.numColumns m))]
    ;;   (dotimes [j (.size row)]
    ;;     (.set row j (.get m i j)))
    ;;   row)
    (Matrices/getSubMatrix m (int-array [i]) (int-array (range (.numColumns m)))))
  (get-column [m i]
    ;; Alternative that returns a new DenseVector each time.
    ;; (Matrices/getColumn m i)
    (Matrices/getSubMatrix m (int-array (range (.numRows m))) (int-array [i])))
  (get-major-slice [m i] (get-row m i))
  (get-slice [m dimension i]
    (case (long dimension)
      0 (get-row m i)
      1 (get-column m i)
      (error "Only dimensions 0 and 1 are supported for a 2D matrix."))))

(extend-protocol PSliceJoin
  Vector
  (join [m ^Vector a]
    (let [m-size (.size m)
          a-size (.size a)
          v (DenseVector. (+ m-size a-size))]
      (dotimes [i m-size] (.set v i (.get m i)))
      (dotimes [i a-size] (.set v (+ i m-size) (.get a i)))
      v))

  Matrix
  (join [m ^Matrix a]
    (let [m-col-size (.numColumns m)
          a-col-size (.numColumns a)]
      (if (= m-col-size a-col-size)
        (let [m-row-size (.numRows m)
              a-row-size (.numRows a)
              mat (DenseMatrix. (+ m-row-size a-row-size) m-col-size)]
          (dotimes [i m-row-size]
            (dotimes [j m-col-size]
              (.set mat i j (.get m i j))))
          (dotimes [i a-row-size]
            (dotimes [j a-col-size]
              (.set mat (+ i m-row-size) j (.get a i j))))
          mat)
        (error "Column size not compatible.")))))

(extend-protocol PSliceJoinAlong
  Vector
  (join-along [m a dim]
    (if (= dim 0)
      (join m a)
      (error "1D vector can only be joined on dimension 0.")))

  Matrix
  (join-along [m ^Matrix a dim]
    (case (long dim)
      0 (join m a)
      1 (let [m-row-size (.numRows m)
              a-row-size (.numRows a)]
          (if (= m-row-size a-row-size)
            (let [m-col-size (.numColumns m)
                  a-col-size (.numColumns a)
                  mat (DenseMatrix. m-row-size (+ m-col-size a-col-size))]
              (dotimes [j m-col-size]
                (dotimes [i m-row-size]
                  (.set mat i j (.get m i j))))
              (dotimes [j a-col-size]
                (dotimes [i a-row-size]
                  (.set mat i (+ j m-col-size) (.get a i j))))
              mat)
            (error "Row size not compatible.")))
      (error "2D matrix can only be joined on dimensions 0 and 1."))))

(extend-protocol PSubVector
  Vector
  (subvector [m start length]
    (Matrices/getSubVector m (int-array (range start (+ start length))))))

(extend-protocol PMatrixSubComponents
  Matrix
  (main-diagonal [m]
    (let [col-size (.numColumns m)
          v (DenseVector. col-size)]
      (dotimes [i col-size]
        (.set v i (.get m i i)))
      v)))

(extend-protocol PZeroCount
  Vector
  (zero-count [v]
    (- (.size v) (nonzero-count v)))

  Matrix
  (zero-count [m]
    (- (* (.numRows m) (.numColumns m)) (nonzero-count m))))

(extend-protocol PAssignment
  Vector
  (assign! [v source]
    (if (number? source)
      (fill! v (double source))
      (error "Only a scalar can be assigned to a 1D vector.")))

  Matrix
  (assign! [m source]
    (cond
      (number? source) (fill! m (double source))
      (and (is-vector? source)
           (= (dimension-count source 0) (dimension-count m 1)))
        (do
          (dotimes [i (.numRows m)] (set-row! m i source))
          m)
      :else (error "source has to have dimensionality 0 or 1 for a 2D matrix."))))

(extend-protocol PImmutableAssignment
  Vector
  (assign [v ^double source]
    (fill! (.copy v) source))

  Matrix
  (assign [m source]
    (assign! (.copy m) source)))

(extend-protocol PMutableFill
  Vector
  (fill! [v ^double value]
    (dotimes [i (.size v)]
      (.set v i value))
    v)

  Matrix
  (fill! [m ^double value]
    (dotimes [i (.numRows m)]
      (dotimes [j (.numColumns m)]
        (.set m i j value)))
    m))

(extend-protocol PDoubleArrayOutput
  Vector
  (to-double-array [v] (Matrices/getArray v))
  (as-double-array [m] nil)

  Matrix
  (to-double-array [m] (Matrices/getArray m))
  (as-double-array [m] nil))

(extend-protocol PObjectArrayOutput
  Vector
  (to-object-array [v] (Matrices/getArray v))
  (as-object-array [m] nil)

  Matrix
  (to-object-array [m] (Matrices/getArray m))
  (as-object-array [m] nil))

(extend-protocol PMatrixMultiply
  Vector
  (matrix-multiply [v a]
    (cond
      (number? a) (element-multiply v a)
      (instance? Vector a) (.dot v ^Vector a)
      (instance? Matrix a)
        (let [^Matrix a a
              len (.size v)
              res-size (.numColumns a) 
              row-matrix (DenseMatrix. 1 len)
              res-matrix (DenseMatrix. 1 res-size)
              res-vector (DenseVector. res-size)]
          (dotimes [i len]
            (.set row-matrix 0 i (.get v i)))
          (.mult row-matrix a res-matrix)
          (dotimes [i res-size]
            (.set res-vector i (.get res-matrix 0 i)))
          res-vector)
      :else (error "Can only multiply scalar, vector and matrix.")))
  (element-multiply [v a]
    (cond
      (number? a) (.scale (.copy v) (double a))
      (instance? Vector a)
        (let [^Vector a a
              size (.size v)]
          (if (= (.size a) size)
            (let [res (DenseVector. size)]
              (dotimes [i size]
                (.set res i (* (.get v i) (.get a i))))
              res)
            (error "shape of a is not compatible.")))
      :else (error "shape of a is not compatible.")))

  Matrix
  (matrix-multiply [m a]
    (cond
      (number? a) (element-multiply m a)
      (instance? Vector a) (.mult m ^Vector a (DenseVector. (.numRows m)))
      (instance? Matrix a)
        (let [^Matrix a a]
          (.mult m a (DenseMatrix. (.numRows m) (.numColumns a))))
      :else (error "Can only multiply scalar, vector and matrix.")))
  (element-multiply [m a]
    (cond
      (number? a) (.scale (.copy m) (double a))
      (instance? Vector a)
        (let [^Vector a a
              row-size (.numRows m)
              col-size (.numColumns m)]
          (if (= (.size a) col-size)
            (let [res-matrix (DenseMatrix. row-size col-size)]
              (dotimes [i row-size]
                (dotimes [j col-size]
                  (.set res-matrix i j (* (.get m i j) (.get a j)))))
              res-matrix)
            (error "shape of a is not compatible.")))
      (instance? Matrix a)
        (let [^Matrix a a
              row-size (.numRows m)
              col-size (.numColumns m)]
          (if (and (= (.numRows a) row-size) (= (.numColumns a) col-size))
            (let [res-matrix (DenseMatrix. row-size col-size)]
              (dotimes [i row-size]
                (dotimes [j col-size]
                  (.set res-matrix i j (* (.get m i j) (.get a i j)))))
              res-matrix)
            (error "shape of a is not compatible.")))
      :else (error "shape of a is not compatible."))))

(extend-protocol PMatrixProducts
  Vector
  (inner-product [m ^Vector a] (matrix-multiply m a))
  (outer-product [m a]
    (let [c (column-matrix imp m)
          r (row-matrix imp a)]
      (matrix-multiply c r))))

(extend-protocol PRowSetting
  Matrix
  (set-row [m i row] (set-row! (.copy m) i row))
  (set-row! [m i row]
    (if (is-vector? row)
      (if (= (dimension-count row 0) (.numColumns m))
        (do
          (dotimes [j (.numColumns m)]
            (.set m i j (get-1d row j)))
          m)
        (error "row has a incompatible shape."))
      (error "row must be a 1D vector."))))

(extend-protocol PColumnSetting
  Matrix
  (set-column [m i column] (set-column! (.copy m) i column))
  (set-column! [m i column]
    (if (is-vector? column)
      (if (= (dimension-count column 0) (.numRows m))
        (do
          (dotimes [j (.numRows m)]
            (.set m j i (get-1d column j)))
          m)
        (error "column has a incompatible shape."))
      (error "column must be a 1D vector."))))

(register-implementation imp)

(comment
  
(m/set-current-implementation :mtj)
(m/equals (m/matrix :mtj [[1 2] [3 4]]) (m/matrix :mtj [[1 2] [3 4]]))
(m/mmul (m/diagonal-matrix :mtj [1 2]) 3)
(m/mmul (m/matrix :mtj [[1 2] [3 4]]) 2)
(m/mmul (m/matrix :mtj [[1 2] [3 4]]) (m/matrix :mtj [2 3]))
(m/mmul (m/matrix :mtj [[1 2] [3 4]]) (m/matrix :mtj [[2 3] [2 1]]))
(m/mmul (m/matrix :vectorz [[1 2] [3 4]]) (m/matrix :vectorz [[2 3] [2 1]]))
(m/mmul (m/matrix :mtj [1 2 3]) 3)
(m/mmul (m/matrix :mtj [1 2 3]) (m/matrix :mtj [3 2 1]))
(m/mmul (m/matrix :mtj [2 3]) (m/matrix :mtj [[3 2 1]
                                              [4 5 6]]))
(m/mmul (m/matrix :vectorz [2 3]) (m/matrix :vectorz [[3 2 1]
                                              [4 5 6]]))
(m/mmul (m/row-matrix :mtj [1 2 3]) (m/matrix :mtj [1 2 3]))
(m/inner-product (m/matrix :mtj [1 2 3]) (m/matrix :mtj [4 3 2]))
(m/emul (m/matrix :mtj [1 2 3]) (m/matrix :mtj [4 5 6]))
(m/emul (m/matrix :mtj [[1 2] [3 4]]) (m/matrix :mtj [[2 2] [3 3]]))
(m/emul (m/matrix :mtj [[1 2] [3 4]]) (m/matrix :mtj [2 2]))

  )
