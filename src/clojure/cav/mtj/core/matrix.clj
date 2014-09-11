(ns cav.mtj.core.matrix
  "Wrapping MTJ for core.matrix."
  (:require [clojure.core.matrix :as m]
            [clojure.core.matrix.protocols :refer :all]
            [clojure.core.matrix.implementations :refer [register-implementation
                                                         get-canonical-object]]
            [clojure.core.matrix.utils :refer [error]])
  (:import [no.uib.cipr.matrix
            Matrix Vector DenseVector DenseMatrix Matrices
            Vector$Norm Matrix$Norm
            BandMatrix LowerSymmBandMatrix
            PermutationMatrix]
           [cav.mtj RowVectorView ColumnVectorView]))

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
                res (DenseVector. size)]
            (dotimes [i size]
              (.set res i (double (get-1d data i))))
            res)
        (sequential? data)
          (let [mat (DenseMatrix. (count data) (count (first data)))]
            (loop [i 0, rows data]
              (when-let [row (first rows)]
                (loop [j 0, vs row]
                  (when-let [v (first vs)]
                    (.set mat i j (double v))
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

    PComputeMatrix
    (compute-matrix [m shape f]
      (case (count shape)
        1 (let [res (DenseVector. (int (first shape)))]
            (dotimes [i (first shape)]
              (.set res i (double (f i))))
            res)
        2 (let [res (DenseMatrix. (int (first shape)) (int (second shape)))]
            (dotimes [i (first shape)]
              (dotimes [j (second shape)]
                (.set res i j (double (f i j)))))
            res)
        (error "Shape not supported.")))

    PSparse
    (sparse-coerce [m data] nil)

    PDense
    (dense-coerce [m data] data)

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
  (set-1d [m row v] (set-1d! (.copy m) row v))
  (set-2d [m row column v] (error "1D vector must be set with only 1 index."))
  (set-nd [m indexes v]
    (if (= 1 (count indexes))
      (set-1d! (.copy m) (first indexes) v)
      (error "1D vector must be set with only 1 index.")))
  (is-mutable? [m] true)

  Matrix
  (set-1d [m row v] (error "2D matrix must be set with 2 indexes."))
  (set-2d [m row column v] (set-2d! (.copy m) row column v))
  (set-nd [m indexes v]
    (if (= 2 (count indexes))
      (set-2d! (.copy m) (first indexes) (second indexes) v)
      (error "2D matrix must be set with 2 indexes.")))
  (is-mutable? [m] true))

(extend-protocol PIndexedSettingMutable
  Vector
  (set-1d! [m ^long row ^double v] (.set m row v))
  (set-2d! [m row column v] (error "1D vector must be set with only 1 index."))
  (set-nd! [m indexes v]
    (if (= 1 (count indexes))
      (set-1d! m (first indexes) v)
      (error "1D vector must be set with only 1 index.")))

  Matrix
  (set-1d! [m row v] (error "2D matrix must be set with 2 indexes."))
  (set-2d! [m ^long row ^long column ^double v] (.set m row column v))
  (set-nd! [m indexes v]
    (if (= 2 (count indexes))
      (set-2d! m (first indexes) (second indexes) v)
      (error "2D matrix must be set with 2 indexes."))))

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
  (mutable-matrix [m] (.copy m))

  Matrix
  (mutable-matrix [m] (.copy m)))

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

(extend-protocol PValidateShape
  Vector
  (validate-shape [m] (get-shape m))

  Matrix
  (validate-shape [m] (get-shape m)))

(extend-protocol PMatrixSlices
  Vector
  (get-row [m i] (error "Cannot access row of 1D vector."))
  (get-column [m i] (error "Cannot access column of 1D vector."))
  (get-major-slice [m ^long i] (.get m i))
  (get-slice [m dimension ^long i]
    (if (= dimension 0)
      (.get m i)
      (error "Only dimension 0 is supported for 1D vector.")))

  Matrix
  (get-row [m i] (RowVectorView. m i))
  (get-column [m i] (ColumnVectorView. m i))
  (get-major-slice [m i] (get-row m i))
  (get-slice [m ^long dimension i]
    (case dimension
      0 (get-row m i)
      1 (get-column m i)
      (error "Only dimensions 0 and 1 are supported for a 2D matrix."))))

(extend-protocol PSliceJoin
  Vector
  (join [m ^Vector a]
    (let [m-size (.size m)
          a-size (.size a)
          res (DenseVector. (+ m-size a-size))]
      (dotimes [i m-size] (.set res i (.get m i)))
      (dotimes [i a-size] (.set res (+ i m-size) (.get a i)))
      res))

  Matrix
  (join [m ^Matrix a]
    (let [m-col-size (.numColumns m)
          a-col-size (.numColumns a)]
      (if (= m-col-size a-col-size)
        (let [m-row-size (.numRows m)
              a-row-size (.numRows a)
              res (DenseMatrix. (+ m-row-size a-row-size) m-col-size)]
          (dotimes [i m-row-size]
            (dotimes [j m-col-size]
              (.set res i j (.get m i j))))
          (dotimes [i a-row-size]
            (dotimes [j a-col-size]
              (.set res (+ i m-row-size) j (.get a i j))))
          res)
        (error "Column size not compatible.")))))

(extend-protocol PSliceJoinAlong
  Vector
  (join-along [m a dim]
    (if (= dim 0)
      (join m a)
      (error "1D vector can only be joined on dimension 0.")))

  Matrix
  (join-along [m ^Matrix a ^long dim]
    (case dim 
      0 (join m a)
      1 (let [m-row-size (.numRows m)
              a-row-size (.numRows a)]
          (if (= m-row-size a-row-size)
            (let [m-col-size (.numColumns m)
                  a-col-size (.numColumns a)
                  res (DenseMatrix. m-row-size (+ m-col-size a-col-size))]
              (dotimes [j m-col-size]
                (dotimes [i m-row-size]
                  (.set res i j (.get m i j))))
              (dotimes [j a-col-size]
                (dotimes [i a-row-size]
                  (.set res i (+ j m-col-size) (.get a i j))))
              res)
            (error "Row size is not compatible.")))
      (error "2D matrix can only be joined on dimensions 0 and 1."))))

(extend-protocol PSubVector
  Vector
  (subvector [m start length]
    (Matrices/getSubVector m (int-array (range start (+ start length))))))

(extend-protocol PMatrixSubComponents
  Matrix
  (main-diagonal [m]
    (let [col-size (.numColumns m)
          res (DenseVector. col-size)]
      (dotimes [i col-size]
        (.set res i (.get m i i)))
      res)))

(extend-protocol PZeroCount
  Vector
  (zero-count [m]
    (- (.size m) (nonzero-count m)))

  Matrix
  (zero-count [m]
    (- (* (.numRows m) (.numColumns m)) (nonzero-count m))))

(extend-protocol PAssignment
  Vector
  (assign! [m source]
    (cond
      (number? source) (fill! m source)
      (same-shape? m source) (do
                               (dotimes [i (.size m)]
                                 (.set m i (.get ^Vector source i)))
                               m) 
      :else (error "Incompatible shape.")))

  Matrix
  (assign! [m source]
    (cond
      (number? source) (fill! m source)
      (and (is-vector? source)
           (= (dimension-count source 0) (dimension-count m 1)))
        (do
          (dotimes [i (.numRows m)]
            (set-row! m i source))
          m)
      (same-shape? m source) (do
                               (dotimes [i (.numRows m)]
                                 (dotimes [j (.numColumns m)]
                                   (.set m i j (.get ^Matrix source i j))))
                               m)
      :else (error "source has to have dimensionality 0 or 1 for a 2D matrix."))))

(extend-protocol PImmutableAssignment
  Vector
  (assign [m source] (assign! (.copy m) source))

  Matrix
  (assign [m source] (assign! (.copy m) source)))

(extend-protocol PMutableFill
  Vector
  (fill! [m ^double value]
    (dotimes [i (.size m)]
      (.set m i value))
    m)

  Matrix
  (fill! [m ^double value]
    (dotimes [i (.numRows m)]
      (dotimes [j (.numColumns m)]
        (.set m i j value)))
    m))

(extend-protocol PDoubleArrayOutput
  Vector
  (to-double-array [m] (Matrices/getArray m))
  (as-double-array [m] nil)

  Matrix
  (to-double-array [m] (Matrices/getArray m))
  (as-double-array [m] nil))

(extend-protocol PObjectArrayOutput
  Vector
  (to-object-array [m] (Matrices/getArray m))
  (as-object-array [m] nil)

  Matrix
  (to-object-array [m] (Matrices/getArray m))
  (as-object-array [m] nil))

(extend-protocol PMatrixMultiply
  Vector
  (matrix-multiply [m a]
    (cond
      (number? a) (element-multiply m a)
      (instance? Vector a) (.dot m ^Vector a)
      (instance? Matrix a)
        (let [^Matrix a a
              len (.size m)
              res-size (.numColumns a) 
              row-matrix (DenseMatrix. 1 len)
              res-matrix (DenseMatrix. 1 res-size)
              res (DenseVector. res-size)]
          (dotimes [i len]
            (.set row-matrix 0 i (.get m i)))
          (.mult row-matrix a res-matrix)
          (dotimes [i res-size]
            (.set res i (.get res-matrix 0 i)))
          res)
      :else (error "Can only multiply scalar, vector and matrix.")))
  (element-multiply [m a] (element-multiply! (.copy m) a))

  Matrix
  (matrix-multiply [m a]
    (cond
      (number? a) (element-multiply m a)
      (instance? Vector a) (.mult m ^Vector a (DenseVector. (.numRows m)))
      (instance? Matrix a)
        (let [^Matrix a a]
          (.mult m a (DenseMatrix. (.numRows m) (.numColumns a))))
      :else (error "Can only multiply scalar, vector and matrix.")))
  (element-multiply [m a] (element-multiply! (.copy m) a)))

(extend-protocol PMatrixMultiplyMutable
  Vector
  (element-multiply! [m a]
    (cond
      (number? a) (.scale m (double a))
      (instance? Vector a)
        (let [^Vector a a
              size (.size m)]
          (if (= (.size a) size)
            (do
              (dotimes [i size]
                (.set m i (* (.get m i) (.get a i))))
              m)
            (error "shape of a is not compatible.")))
      :else (error "shape of a is not compatible.")))

  Matrix
  (element-multiply! [m a]
    (cond
      (number? a) (.scale m (double a))
      (instance? Vector a)
        (let [^Vector a a
              row-size (.numRows m)
              col-size (.numColumns m)]
          (if (= (.size a) col-size)
            (do
              (dotimes [i row-size]
                (dotimes [j col-size]
                  (.set m i j (* (.get m i j) (.get a j)))))
              m)
            (error "shape of a is not compatible.")))
      (instance? Matrix a)
        (let [^Matrix a a
              row-size (.numRows m)
              col-size (.numColumns m)]
          (if (and (= (.numRows a) row-size) (= (.numColumns a) col-size))
            (do
              (dotimes [i row-size]
                (dotimes [j col-size]
                  (.set m i j (* (.get m i j) (.get a i j)))))
              m)
            (error "shape of a is not compatible.")))
      :else (error "shape of a is not compatible."))))

(extend-protocol PMatrixProducts
  Vector
  (inner-product [m ^Vector a] (.dot m a))
  (outer-product [m ^Vector a]
    (let [c (column-matrix imp m)
          r (row-matrix imp a)]
      (matrix-multiply c r))))

(extend-protocol PValueEquality
  Vector
  (value-equals [m a] (matrix-equals m a))
  
  Matrix
  (value-equals [m a] (matrix-equals m a)))

(extend-protocol PMatrixEquality
  Vector
  (matrix-equals [a ^Vector b]
    (let [size (.size a)]
      (if (= size (.size b))
        (loop [i 0]
          (if (< i size)
            (if (= (.get a i) (.get b i))
              (recur (inc i))
              false)
            true))
        false)))

  Matrix
  (matrix-equals [a ^Matrix b]
    (let [row-size (.numRows a)
          col-size (.numColumns a)]
      (if (and (= (.numRows b) row-size)
               (= (.numColumns b) col-size))
        (loop [i 0]
          (if (< i row-size)
            (if (loop [j 0]
                  (if (< j col-size)
                    (if (= (.get a i j) (.get b i j))
                      (recur (inc j))
                      false)
                    true))
              (recur (inc i))
              false)
            true))
        false))))

(defn- === [a b eps]
  (and (>= a (- b eps)) (<= a (+ b eps))))

(extend-protocol PMatrixEqualityEpsilon
  Vector
  (matrix-equals-epsilon [a ^Vector b eps]
    (let [size (.size a)]
      (if (= size (.size b))
        (loop [i 0]
          (if (< i size)
            (if (=== (.get a i) (.get b i) eps)
              (recur (inc i))
              false)
            true))
        false)))

  Matrix
  (matrix-equals-epsilon [a ^Matrix b eps]
    (let [row-size (.numRows a)
          col-size (.numColumns a)]
      (if (and (= (.numRows b) row-size)
               (= (.numColumns b) col-size))
        (loop [i 0]
          (if (< i row-size)
            (if (loop [j 0]
                  (if (< j col-size)
                    (if (=== (.get a i j) (.get b i j) eps)
                      (recur (inc j))
                      false)
                    true))
              (recur (inc i))
              false)
            true))
        false))))

(extend-protocol PAddScaled
  Vector
  (add-scaled [m a factor] (add-scaled! (.copy m) a factor))
  
  Matrix
  (add-scaled [m a factor] (add-scaled! (.copy m) a factor)))

(extend-protocol PAddScaledMutable
  Vector
  (add-scaled! [m ^Vector a ^double factor] (.add m factor a))
  
  Matrix
  (add-scaled! [m ^Matrix a ^double factor] (.add m factor a)))

(extend-protocol PMatrixDivide
  Vector
  (element-divide
    ([m] (element-divide! (.copy m)))
    ([m a] (element-divide! (.copy m) a)))

  Matrix
  (element-divide
    ([m] (element-divide! (.copy m)))
    ([m a] (element-divide! (.copy m) a))))

(extend-protocol PMatrixDivideMutable
  Vector
  (element-divide!
    ([m]
     (dotimes [i (.size m)]
       (.set m i (/ 1.0 (.get m i))))
     m)
    ([m a]
     (cond
       (number? a) (element-multiply! m (/ 1.0 (double a)))
       (same-shape? m a) (do
                           (dotimes [i (.size m)]
                             (.set m i (/ (.get m i) (.get ^Vector a i))))
                           m)
       :else (error "Incompatible shape."))))

  Matrix
  (element-divide!
    ([m]
     (dotimes [i (.numRows m)]
       (dotimes [j (.numColumns m)]
         (.set m i j (/ 1.0 (.get m i j)))))
     m)
    ([m a]
     (cond
       (number? a) (element-multiply! m (/ 1.0 (double a)))
       (and (instance? Vector a) (= (.size ^Vector a) (.numColumns m)))
         (do
           (dotimes [i (.numRows m)]
             (dotimes [j (.numColumns m)]
               (.set m i j (/ (.get m i j) (.get ^Vector a j)))))
           m)
       (same-shape? m a) (do
                           (dotimes [i (.numRows m)]
                             (dotimes [j (.numColumns m)]
                               (.set m i j (/ (.get m i j) (.get ^Matrix a i j)))))
                           m)
       :else (error "Incompatible shape.")))))

(extend-protocol PMatrixScaling
  Vector
  (scale [m a] (scale! (.copy m) a))
  (pre-scale [m a] (scale! (.copy m) a))

  Matrix
  (scale [m a] (scale! (.copy m) a))
  (pre-scale [m a] (scale! (.copy m) a)))

(extend-protocol PMatrixMutableScaling
  Vector
  (scale! [m a] (.scale m a))
  (pre-scale! [m a] (scale! m a))

  Matrix
  (scale! [m a] (.scale m a))
  (pre-scale! [m a] (scale! m a)))

(extend-protocol PMatrixAdd
  Vector
  (matrix-add [m a] (matrix-add! (.copy m) a))
  (matrix-sub [m a] (matrix-sub! (.copy m) a))
  
  Matrix
  (matrix-add [m a] (matrix-add! (.copy m) a))
  (matrix-sub [m a] (matrix-add! (.copy m) a)))

(extend-protocol PMatrixAddMutable
  Vector
  (matrix-add! [m a]
    (cond
      (number? a)
        (let [a (double a)]
          (dotimes [i (.size m)]
            (.set m i (+ (.get m i) a)))
          m)
      (instance? Vector a) (.add m ^Vector a)
      :else (error "Incompatible shape.")))
  (matrix-sub! [m a]
    (cond
      (number? a) (matrix-add! m (- a))
      (instance? Vector a) (.add m -1.0 ^Vector a)
      :else (error "Incompatible shape.")))
  
  Matrix
  (matrix-add! [m a]
    (cond
      (number? a)
        (let [a (double a)]
          (dotimes [i (.numRows m)]
            (dotimes [j (.numColumns m)]
              (.set m i j (+ (.get m i j) a))))
          m)
      (instance? Vector a)
        (let [a ^Vector a]
          (dotimes [i (.numRows m)]
            (dotimes [j (.numColumns m)]
              (.set m i j (+ (.get m i j) (.get a j)))))
          m)
      (instance? Matrix a) (.add m ^Matrix a)
      :else (error "Incompatible shape.")))
  (matrix-sub! [m a]
    (cond
      (number? a) (matrix-add! m (- a))
      (instance? Vector a)
        (let [a ^Vector a]
          (dotimes [i (.numRows m)]
            (dotimes [j (.numColumns m)]
              (.set m i j (- (.get m i j) (.get a j)))))
          m)
      (instance? Matrix a) (.add m -1 ^Matrix a)
      :else (error "Incompatible shape."))))

(defn- compute-indices [pair size]
  (if (nil? pair)
    (range size)
    (let [[start len] pair] (range start (+ start len)))))

(extend-protocol PSubMatrix
  Matrix
  (submatrix [m dim-ranges]
    (if (= (count dim-ranges) 2)
      (Matrices/getSubMatrix
        m
        (int-array (compute-indices (first dim-ranges) (.numRows m)))
        (int-array (compute-indices (second dim-ranges) (.numColumns m))))
      (error "Only dimension 0 and 1 are supported."))))

(extend-protocol PTranspose
  Vector
  (transpose [m] m)

  Matrix
  (transpose [m] (.transpose m (DenseMatrix. (.numColumns m) (.numRows m)))))

(extend-protocol PTransposeInPlace
  Matrix
  (transpose! [m]
    (if (.isSquare m)
      (.transpose m)
      (error "Matrix has to be squared."))))

(extend-protocol PNumerical
  Vector
  (numerical? [m] true)

  Matrix
  (numerical? [m] true))

(extend-protocol PVectorOps
  Vector
  (vector-dot [a ^Vector b]
    (if (= (.size a) (.size b))
      (.dot a b)
      (error "Incompatible shape.")))
  (length [a] (.norm a Vector$Norm/Two))
  (length-squared [a] (let [l (.norm a Vector$Norm/Two)] (* l l)))
  (normalise [a] (normalise! (.copy a))))

(extend-protocol PMutableVectorOps
  Vector
  (normalise! [a]
    (element-divide! a (.norm a Vector$Norm/Two))))

(extend-protocol PMatrixOps
  Matrix
  (trace [m]
    (if (.isSquare m)
      (loop [i 0, sum 0.0]
        (if (< i (.numRows m))
          (recur (inc i) (+ sum (.get m i i)))
          sum))
      (error "Square matrix required.")))
  (determinant [m]
    (->> m
         (coerce-param (get-canonical-object :ndarray-double))
         (determinant)))
  (inverse [m]
    (let [iden (Matrices/identity (.numColumns m))
          res (.copy iden)]
      (.solve m iden res))))

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
(m/div 2 (m/matrix :mtj [1 2 3]))
(m/sqrt (m/matrix :mtj [[1 2] [3 4]]))

  )
