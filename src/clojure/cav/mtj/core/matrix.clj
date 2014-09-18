(ns cav.mtj.core.matrix
  "Wrapping MTJ for core.matrix."
  (:use [clojure.core.matrix.protocols])
  (:require [clojure.core.matrix]
            [clojure.core.matrix.implementations :as impl]
            [clojure.core.matrix.utils :refer [error same-shape-object?]])
  (:import [no.uib.cipr.matrix
            DenseVector DenseMatrix Matrices
            Vector$Norm Matrix$Norm
            LowerSymmBandMatrix
            PermutationMatrix
            SVD QR DenseLU DenseCholesky EVD SymmDenseEVD]
           [cav.mtj Mat Vec RowVectorView ColumnVectorView]))

(definline ^:private array? [a] `(not (is-scalar? ~a)))

(definline ^:private zero-dim? [a] `(= (dimensionality ~a) 0))

(def ^:private mtj-impl (Mat. (DenseMatrix. 1 1)))

(extend-type Mat
  PImplementation
  (implementation-key [a] :mtj)
  (meta-info [a] {:doc "MTJ as an implementation of core.matrix."})
  (supports-dimensionality? [a dimensions] (<= 1 dimensions 2))
  (new-vector [a length] (Vec. (DenseVector. (int length))))
  (new-matrix [a rows columns] (Mat. (DenseMatrix. (int rows) (int columns))))
  (new-matrix-nd [a shape]
    (case (count shape)
      1 (new-vector a (first shape))
      2 (new-matrix a (first shape) (second shape))
      nil))
  (construct-matrix [a data]
    (cond
      (zero-dim? data) (double (get-0d data))
      (instance? Vec data) data
      (instance? Mat data) data   
      (instance? (Class/forName "[D") data) (Vec. (DenseVector. ^doubles data))
      (instance? (Class/forName "[[D") data) (Mat. (DenseMatrix. ^"[[D" data))
      (and (sequential? data) (number? (first data))) (Vec. (DenseVector. (double-array data)))
      (sequential? data) (Mat. (DenseMatrix. ^"[[D" (into-array (map #(double-array %) data))))

      ;; Following implementation is more general, but known to be slow for
      ;; some matrix implementations like seqs, because of get-1d/get-2d usage.

      ;; (is-vector? data) (let [^int size (dimension-count data 0)
      ;;                         res (DenseVector. size)]
      ;;                     (dotimes [i size]
      ;;                       (.set res i (double (get-1d data i))))
      ;;                     (Vec. res))
      ;; (array? data) (if (= (dimensionality data) 2)
      ;;                 (let [^int row-size (dimension-count data 0)
      ;;                       ^int col-size (dimension-count data 1)
      ;;                       res (DenseMatrix. row-size col-size)]
      ;;                   (dotimes [i row-size]
      ;;                     (dotimes [j col-size]
      ;;                       (.set res i j (get-2d data i j))))
      ;;                   (Mat. res))
      ;;                 (error "Incompatible shape."))
      ))

  PRowColMatrix
  (column-matrix [a data]
    (cond
      (instance? Vec data) (Mat. (DenseMatrix. (.mtj ^Vec data)))
      (is-vector? data) (Mat. (DenseMatrix. ^"[[D" (into-array (map #(double-array [%]) data))))
      :else (error "data has to be a 1D vector.")))
  (row-matrix [a data]
    (if (is-vector? data)
      (Mat. (DenseMatrix. ^"[[D" (into-array [(double-array data)])))
      (error "data has to be a 1D vector.")))

  PComputeMatrix
  (compute-matrix [a shape f]
    (case (count shape)
      1 (let [res (DenseVector. (int (first shape)))]
          (dotimes [i (first shape)]
            (.set res i (double (f i))))
          (Vec. res))
      2 (let [res (DenseMatrix. (int (first shape)) (int (second shape)))]
          (dotimes [i (first shape)]
            (dotimes [j (second shape)]
              (.set res i j (double (f i j)))))
          (Mat. res))
      (error "Incompatible shape.")))

  PSparse
  (sparse-coerce [a data] nil)

  PDense
  (dense-coerce [a data] data)

  PCoercion
  (coerce-param [a param]
    (cond
      (instance? Vec param) param
      (instance? Mat param) param
      (is-vector? param)
        (do
          (println "WARNING: You are coercing a non-MTJ object, this hurts performance!")
          (construct-matrix a param))
      (array? param)
        (do
          (println "WARNING: You are coercing a non-MTJ object, this hurts performance!")
          (construct-matrix a param))
      :else nil))

  PSpecialisedConstructors
  (identity-matrix [a dims] (Mat. (Matrices/identity (int dims))))
  (diagonal-matrix [a diagonal-values]
    (let [^int size (dimension-count diagonal-values 0)
          res (LowerSymmBandMatrix. size 0)]
      (loop [i 0, vs diagonal-values]
        (when (< i size)
          (.set res i i (first vs))
          (recur (inc i) (rest vs))))
      (Mat. res)))

  PPermutationMatrix
  (permutation-matrix [a permutation]
    (Mat. (PermutationMatrix. (int-array permutation)))))

(extend-protocol PDimensionInfo
  Vec
  (dimensionality [a] 1)
  (get-shape [a] [(.size a)])
  (is-scalar? [a] false)
  (is-vector? [a] true)
  (dimension-count [a ^long dimension-number]
    (case dimension-number 
      0 (.size a)
      (error "Illegal dimension number.")))

  Mat
  (dimensionality [a] 2)
  (get-shape [a] [(.numRows a) (.numColumns a)])
  (is-scalar? [a] false)
  (is-vector? [a] false)
  (dimension-count [a ^long dimension-number]
    (case dimension-number 
      0 (.numRows a)
      1 (.numColumns a)
      (error "Illegal dimension number."))))

(extend-protocol PIndexedAccess
  Vec
  (get-1d [a row] (.get a (int row)))
  (get-2d [a row column] (error "1D vector must be accessed by only 1 index."))
  (get-nd [a indexes]
    (if (= 1 (count indexes))
      (get-1d a (first indexes))
      (error "1D vector must be accessed by only 1 index.")))

  Mat
  (get-1d [a row] (error "2D matrix must be accessed by row and column."))
  (get-2d [a row column] (.get a (int row) (int column)))
  (get-nd [a indexes]
    (if (= 2 (count indexes))
      (get-2d a (first indexes) (second indexes))
      (error "2D matrix must be accessed by row and column."))))

(extend-protocol PIndexedSetting
  Vec
  (set-1d [a row v] (set-1d! (.copy a) row v))
  (set-2d [a row column v] (error "1D vector must be set with only 1 index."))
  (set-nd [a indexes v]
    (if (= 1 (count indexes))
      (set-1d! (.copy a) (first indexes) v)
      (error "1D vector must be set with only 1 index.")))
  (is-mutable? [a] true)

  Mat
  (set-1d [a row v] (error "2D matrix must be set with 2 indexes."))
  (set-2d [a row column v] (set-2d! (.copy a) row column v))
  (set-nd [a indexes v]
    (if (= 2 (count indexes))
      (set-2d! (.copy a) (first indexes) (second indexes) v)
      (error "2D matrix must be set with 2 indexes.")))
  (is-mutable? [a] true))

(extend-protocol PIndexedSettingMutable
  Vec
  (set-1d! [a ^long row ^double v]
    (do
      (.set a row v)
      a))
  (set-2d! [a row column v] (error "1D vector must be set with only 1 index."))
  (set-nd! [a indexes v]
    (if (= 1 (count indexes))
      (set-1d! a (first indexes) v)
      (error "1D vector must be set with only 1 index.")))

  Mat
  (set-1d! [a row v] (error "2D matrix must be set with 2 indexes."))
  (set-2d! [a ^long row ^long column ^double v]
    (do
      (.set a row column v)
      a))
  (set-nd! [a indexes v]
    (if (= 2 (count indexes))
      (set-2d! a (first indexes) (second indexes) v)
      (error "2D matrix must be set with 2 indexes."))))

(extend-protocol PMatrixCloning
  Vec
  (clone [a] (.copy a))

  Mat
  (clone [a] (.copy a)))

(extend-protocol PTypeInfo
  Vec
  (element-type [a] Double/TYPE)

  Mat
  (element-type [a] Double/TYPE))

(extend-protocol PSparse
  Vec
  (sparse [a] a)

  Mat
  (sparse [a] a))

(extend-protocol PDense
  Vec
  (dense [a] a)

  Mat
  (dense [a] a))

(extend-protocol PArrayMetrics
  Vec
  (nonzero-count [a] (Matrices/cardinality (.mtj a)))

  Mat
  (nonzero-count [a] (Matrices/cardinality (.mtj a))))

(extend-protocol PMutableMatrixConstruction
  Vec
  (mutable-matrix [a] (.copy a))

  Mat
  (mutable-matrix [a] (.copy a)))

(extend-protocol PSameShape
  Vec
  (same-shape? [a b]
    (cond
      (instance? Vec b) (= (.size a) (.size ^Vec b))
      :else (same-shape-object? (get-shape a) (get-shape b))))

  Mat
  (same-shape? [a b]
    (cond
      (instance? Mat b) (let [^Mat b b]
                             (and (= (.numRows a) (.numRows b))
                                  (= (.numColumns a) (.numColumns b))))
      :else (same-shape-object? (get-shape a) (get-shape b)))))

(extend-protocol PValidateShape
  Vec
  (validate-shape [a] (get-shape a))

  Mat
  (validate-shape [a] (get-shape a)))

(extend-protocol PMatrixSlices
  Vec
  (get-row [a i] (error "Cannot access row of 1D vector."))
  (get-column [a i] (error "Cannot access column of 1D vector."))
  (get-major-slice [a ^long i] (.get a i))
  (get-slice [a dimension ^long i]
    (if (= dimension 0)
      (.get a i)
      (error "Only dimension 0 is supported for 1D vector.")))

  Mat
  (get-row [a i] (Vec. (RowVectorView. (.mtj a) i)))
  (get-column [a i] (Vec. (ColumnVectorView. (.mtj a) i)))
  (get-major-slice [a i] (get-row a i))
  (get-slice [a ^long dimension i]
    (case dimension
      0 (get-row a i)
      1 (get-column a i)
      (error "Only dimensions 0 and 1 are supported for a 2D matrix."))))

(extend-protocol PSliceJoin
  Vec
  (join [a x]
    (let [m-size (.size a)
          ^int size (+ m-size (dimension-count x 0))
          res (DenseVector. size)]
      (dotimes [i m-size] (.set res i (.get a i)))
      (loop [i m-size, vs x]
        (when (< i size)
          (.set res i (double (first vs)))
          (recur (inc i) (rest vs))))
      (Vec. res)))

  Mat
  (join [a x]
    (let [m-col-size (.numColumns a)
          ^int x-col-size (dimension-count x 1)]
      (if (= m-col-size x-col-size)
        (let [m-row-size (.numRows a)
              ^int x-row-size (dimension-count x 0)
              row-size (+ m-row-size x-row-size)
              res (DenseMatrix. row-size m-col-size)]
          (dotimes [i m-row-size]
            (dotimes [j m-col-size]
              (.set res i j (.get a i j))))
          (dotimes [i x-row-size]
            (dotimes [j m-col-size]
              (.set res (+ m-row-size i) j (get-2d x i j))))
          (Mat. res))
        (error "Column size not compatible.")))))

(extend-protocol PSliceJoinAlong
  Vec
  (join-along [a x dim]
    (if (= dim 0)
      (join a x)
      (error "1D vector can only be joined on dimension 0.")))

  Mat
  (join-along [a x ^long dim]
    (case dim 
      0 (join a x)
      1 (let [m-row-size (.numRows a)
              ^int a-row-size (dimension-count x 0)]
          (if (= m-row-size a-row-size)
            (let [m-col-size (.numColumns a)
                  ^int a-col-size (dimension-count x 1)
                  res (DenseMatrix. m-row-size (+ m-col-size a-col-size))]
              (dotimes [j m-col-size]
                (dotimes [i m-row-size]
                  (.set res i j (.get a i j))))
              (dotimes [j a-col-size]
                (dotimes [i a-row-size]
                  (.set res i (+ j m-col-size) (double (get-2d x i j)))))
              (Mat. res))
            (error "Row size is not compatible.")))
      (error "2D matrix can only be joined on dimensions 0 and 1."))))

(extend-protocol PSubVector
  Vec
  (subvector [a start length]
    (Vec. (Matrices/getSubVector (.mtj a) (int-array (range start (+ start length)))))))

(extend-protocol PMatrixSubComponents
  Mat
  (main-diagonal [a]
    (let [size (min (.numRows a) (.numColumns a))
          res (DenseVector. size)]
      (dotimes [i size]
        (.set res i (.get a i i)))
      (Vec. res))))

(extend-protocol PZeroCount
  Vec
  (zero-count [a]
    (- (.size a) (nonzero-count a)))

  Mat
  (zero-count [a]
    (- (* (.numRows a) (.numColumns a)) (nonzero-count a))))

(extend-protocol PAssignment
  Vec
  (assign! [a source]
    (cond
      (zero-dim? source) (fill! a (get-0d source))
      (same-shape? a source) (let [size (.size a)]
                               (loop [i 0, vs source]
                                 (when (< i size)
                                   (.set a i (double (first vs)))
                                   (recur (inc i) (rest vs))))
                               a) 
      :else (error "Incompatible shape.")))
  (assign-array! [a arr] (assign! a arr))

  Mat
  (assign! [a source]
    (cond
      (zero-dim? source) (fill! a (get-0d source))
      (and (is-vector? source)
           (= (dimension-count source 0) (dimension-count a 1)))
        (do
          (dotimes [i (.numRows a)]
            (set-row! a i source))
          a)
      (same-shape? a source) (do
                               (dotimes [i (.numRows a)]
                                 (dotimes [j (.numColumns a)]
                                   (.set a i j (double (get-2d source i j)))))
                               a)
      :else (error "source has to have dimensionality 0 or 1 for a 2D matrix.")))
  (assign-array! [a arr]
    (let [row-size (.numRows a)
          col-size (.numColumns a)]
      (loop [i 0, j 0, vs arr]
        (when (< i row-size)
          (if (< j col-size)
            (do
              (.set a i j (first vs))
              (recur i (inc j) (rest vs)))
            (recur (inc i) 0 vs))))
      a)))

(extend-protocol PImmutableAssignment
  Vec
  (assign [a source] (assign! (.copy a) source))

  Mat
  (assign [a source] (assign! (.copy a) source)))

(extend-protocol PMutableFill
  Vec
  (fill! [a ^double value]
    (dotimes [i (.size a)]
      (.set a i value))
    a)

  Mat
  (fill! [a ^double value]
    (dotimes [i (.numRows a)]
      (dotimes [j (.numColumns a)]
        (.set a i j value)))
    a))

(extend-protocol PDoubleArrayOutput
  Vec
  (to-double-array [a] (Matrices/getArray (.mtj a)))
  (as-double-array [a] nil)

  Mat
  (to-double-array [a]
    (let [row-size (.numRows a)
          col-size (.numColumns a)
          res (double-array (* row-size col-size))]
      (dotimes [i row-size]
        (dotimes [j col-size]
          (aset res (+ (* i col-size) j) (.get a i j))))
      res))
  (as-double-array [a] nil))

(extend-protocol PObjectArrayOutput
  Vec
  (to-object-array [a] (Matrices/getArray (.mtj a)))
  (as-object-array [a] nil)

  Mat
  (to-object-array [a] (to-double-array a))
  (as-object-array [a] nil))

(extend-protocol PMatrixMultiply
  Vec
  (matrix-multiply [a x]
    (cond
      (zero-dim? x) (element-multiply! (.copy a) (get-0d x))
      (instance? Vec x) (.dot a ^Vec x)
      (instance? Mat x)
        (let [^Mat x x
              len (.size a)
              res-size (.numColumns x) 
              row-matrix (DenseMatrix. 1 len)
              res-matrix (DenseMatrix. 1 res-size)
              res (DenseVector. res-size)]
          (dotimes [i len]
            (.set row-matrix 0 i (.get a i)))
          (.mult row-matrix (.mtj x) res-matrix)
          (dotimes [i res-size]
            (.set res i (.get res-matrix 0 i)))
          (Vec. res))
      (is-vector? x) (matrix-multiply a (coerce-param mtj-impl x))
      (array? x) (matrix-multiply a (coerce-param mtj-impl x))
      :else (error "Can only multiply scalar, vector and matrix.")))
  (element-multiply [a x] (element-multiply! (.copy a) x))

  Mat
  (matrix-multiply [a x]
    (cond
      (zero-dim? x) (element-multiply! (.copy a) (get-0d x))
      (instance? Vec x) (.mult a ^Vec x (Vec. (DenseVector. (.numRows a))))
      (instance? Mat x)
        (let [^Mat x x]
          (.mult a x (Mat. (DenseMatrix. (.numRows a) (.numColumns x)))))
      (is-vector? x) (matrix-multiply a (coerce-param mtj-impl x))
      (array? x) (matrix-multiply a (coerce-param mtj-impl x))
      :else (error "Can only multiply scalar, vector and matrix.")))
  (element-multiply [a x] (element-multiply! (.copy a) x)))

(extend-protocol PMatrixMultiplyMutable
  Vec
  (element-multiply! [a x]
    (cond
      (zero-dim? x) (.scale a (double (get-0d x)))
      (instance? Vec x)
        (let [^Vec x x
              size (.size a)]
          (if (= (.size x) size)
            (do
              (dotimes [i size]
                (.set a i (* (.get a i) (.get x i))))
              a)
            (error "shape of x is not compatible.")))
      (is-vector? x)
        (let [size (.size a)]
          (if (= (dimension-count x 0) size)
            (do
              (loop [i 0, vs x]
                (when (< i size)
                  (.set a i (* (.get a i) (first vs)))
                  (recur (inc i) (rest vs))))
              a)
            (error "shape of x is not compatible.")))
      :else (error "shape of x is not compatible.")))

  Mat
  (element-multiply! [a x]
    (cond
      (zero-dim? x) (.scale a (double (get-0d x)))
      (instance? Vec x)
        (let [^Vec x x
              row-size (.numRows a)
              col-size (.numColumns a)]
          (if (= (.size x) col-size)
            (do
              (dotimes [i row-size]
                (dotimes [j col-size]
                  (.set a i j (* (.get a i j) (.get x j)))))
              a)
            (error "shape of x is not compatible.")))
      (instance? Mat x)
        (let [^Mat x x
              row-size (.numRows a)
              col-size (.numColumns a)]
          (if (and (= (.numRows x) row-size) (= (.numColumns x) col-size))
            (do
              (dotimes [i row-size]
                (dotimes [j col-size]
                  (.set a i j (* (.get a i j) (.get x i j)))))
              a)
            (error "shape of x is not compatible.")))
      (is-vector? x)
        (let [row-size (.numRows a)
              col-size (.numColumns a)]
          (if (= (dimension-count x 0) col-size)
            (do
              (dotimes [i row-size]
                (loop [j 0, vs x]
                  (when (< j col-size)
                    (.set a i j (* (.get a i j) (first x)))
                    (recur (inc j) (rest vs)))))
              a)
            (error "shape of x is not compatible.")))
      (array? x)
        (let [row-size (.numRows a)
              col-size (.numColumns a)]
          (if (and (= (dimension-count x 0) row-size)
                   (= (dimension-count x 1) col-size))
            (do
              (dotimes [i row-size]
                (dotimes [j col-size]
                  (.set a i j (* (.get a i j) (get-2d x i j)))))
              a)
            (error "shape of x is not compatible.")))
      :else (error "shape of x is not compatible."))))

(extend-protocol PMatrixProducts
  Vec
  (inner-product [a x]
    (cond
      (instance? Vec x) (.dot a ^Vec x)
      (is-vector? x)
        (do
          (println "WARNING: Operations with mixed matrices/vectors may not yield the fastest result.")
          (matrix-multiply a x))
      :else (error "Need vectors for inner product.")))
  (outer-product [a x]
    (let [c (column-matrix mtj-impl a)
          r (row-matrix mtj-impl x)]
      (matrix-multiply c r))))

(extend-protocol PAddScaled
  Vec
  (add-scaled [a x factor] (add-scaled! (.copy a) x factor))

  Mat
  (add-scaled [a x factor] (add-scaled! (.copy a) x factor)))

(extend-protocol PAddScaledMutable
  Vec
  (add-scaled! [a x ^double factor]
    (cond
      (instance? Vec x) (.add a factor ^Vec x)
      :else (matrix-add! a (scale x factor))))

  Mat
  (add-scaled! [a x ^double factor]
    (cond
      (instance? Mat x) (.add a factor ^Mat x)
      :else (matrix-add! a (scale x factor)))))

(extend-protocol PMatrixDivide
  Vec
  (element-divide
    ([a] (element-divide! (.copy a)))
    ([a x] (element-divide! (.copy a) x)))

  Mat
  (element-divide
    ([a] (element-divide! (.copy a)))
    ([a x] (element-divide! (.copy a) x))))

(extend-protocol PMatrixDivideMutable
  Vec
  (element-divide!
    ([a]
     (dotimes [i (.size a)]
       (.set a i (/ 1.0 (.get a i))))
     a)
    ([a x]
     (cond
       (zero-dim? x) (element-multiply! a (/ 1.0 (double (get-0d x))))
       (same-shape? a x) (let [size (.size a)]
                           (loop [i 0, vs x]
                             (when (< i size)
                               (.set a i (/ (.get a i) (first vs)))
                               (recur (inc i) (rest vs))))
                           a)
       :else (error "Incompatible shape."))))

  Mat
  (element-divide!
    ([a]
     (dotimes [i (.numRows a)]
       (dotimes [j (.numColumns a)]
         (.set a i j (/ 1.0 (.get a i j)))))
     a)
    ([a x]
     (cond
       (zero-dim? x) (element-multiply! a (/ 1.0 (double (get-0d x))))
       (and (is-vector? x) (= (dimension-count x 0) (.numColumns a)))
         (let [col-size (.numColumns a)]
           (dotimes [i (.numRows a)]
             (loop [j 0, vs x]
               (when (< j col-size)
                 (.set a i j (/ (.get a i j) (first vs)))
                 (recur (inc j) (rest vs)))))
           a)
       (same-shape? a x) (do
                           (dotimes [i (.numRows a)]
                             (dotimes [j (.numColumns a)]
                               (.set a i j (/ (.get a i j) (get-2d x i j)))))
                           a)
       :else (error "Incompatible shape.")))))

(extend-protocol PMatrixScaling
  Vec
  (scale [a x] (scale! (.copy a) x))
  (pre-scale [a x] (scale! (.copy a) x))

  Mat
  (scale [a x] (scale! (.copy a) x))
  (pre-scale [a x] (scale! (.copy a) x)))

(extend-protocol PMatrixMutableScaling
  Vec
  (scale! [a x] (.scale a x))
  (pre-scale! [a x] (scale! a x))

  Mat
  (scale! [a x] (.scale a x))
  (pre-scale! [a x] (scale! a x)))

(extend-protocol PMatrixAdd
  Vec
  (matrix-add [a x] (matrix-add! (.copy a) x))
  (matrix-sub [a x] (matrix-sub! (.copy a) x))

  Mat
  (matrix-add [a x] (matrix-add! (.copy a) x))
  (matrix-sub [a x] (matrix-add! (.copy a) x)))

(extend-protocol PMatrixAddMutable
  Vec
  (matrix-add! [a x]
    (cond
      (instance? Vec x) (.add a ^Vec x)
      (zero-dim? x) (let [x (double (get-0d x))]
                      (dotimes [i (.size a)]
                        (.add a i x))
                      a)
      (is-vector? x)
        (let [size (.size a)]
          (println "WARNING: Operations with mixed matrices/vectors may not yield the fastest result.")
          (loop [i 0, vs x]
            (when (< i size)
              (.add a i (double (first vs)))
              (recur (inc i) (rest vs))))
          a)
      :else (error "Incompatible shape.")))
  (matrix-sub! [a x]
    (cond
      (instance? Vec x) (.add a -1.0 ^Vec x)
      (zero-dim? x) (matrix-add! a (- (get-0d x)))
      (is-vector? x)
        (let [size (.size a)]
          (println "WARNING: Operations with mixed matrices/vectors may not yield the fastest result.")
          (loop [i 0, vs x]
            (when (< i size)
              (.add a i (double (- (first vs))))
              (recur (inc i) (rest vs))))
          a)
      :else (error "Incompatible shape.")))

  Mat
  (matrix-add! [a x]
    (cond
      (instance? Mat x) (.add a ^Mat x)
      (zero-dim? x) (let [x (get-0d x)]
                      (dotimes [i (.numRows a)]
                        (dotimes [j (.numColumns a)]
                          (.add a i j x)))
                      a)
      (is-vector? x) (let [col-size (.numColumns a)]
                       (dotimes [i (.numRows a)]
                         (loop [j 0, vs x]
                           (when (< j col-size)
                             (.add a i j (first vs))
                             (recur (inc j) (rest vs)))))
                       a)
      (same-shape? a x) (do
                          (println "WARNING: Operations with mixed matrices/vectors may not yield the fastest result.")
                          (dotimes [i (.numRows a)]
                            (dotimes [j (.numColumns a)]
                              (.add a i j (get-2d x i j))))
                          a)
      :else (error "Incompatible shape.")))
  (matrix-sub! [a x]
    (cond
      (instance? Mat x) (.add a -1.0 ^Mat x)
      (zero-dim? x) (matrix-add! a (- (get-0d x)))
      (is-vector? x) (let [col-size (.numColumns a)]
                       (dotimes [i (.numRows a)]
                         (loop [j 0, vs x]
                           (when (< j col-size)
                             (.add a i j (- (first vs)))
                             (recur (inc j) (rest vs)))))
                       a)
      (same-shape? a x) (do
                          (println "WARNING: Operations with mixed matrices/vectors may not yield the fastest result.")
                          (dotimes [i (.numRows a)]
                            (dotimes [j (.numColumns a)]
                              (.add a i j (- (get-2d x i j)))))
                          a)
      :else (error "Incompatible shape."))))

(defn- compute-indices [pair size]
  (if (nil? pair)
    (range size)
    (let [[start len] pair] (range start (+ start len)))))

(extend-protocol PSubMatrix
  Mat
  (submatrix [a dim-ranges]
    (Mat. (Matrices/getSubMatrix
            (.mtj a)
            (int-array (compute-indices (first dim-ranges) (.numRows a)))
            (int-array (compute-indices (second dim-ranges) (.numColumns a)))))))

(extend-protocol PTranspose
  Vec
  (transpose [a] a)

  Mat
  (transpose [a] (.transpose a (Mat. (DenseMatrix. (.numColumns a) (.numRows a))))))

(extend-protocol PTransposeInPlace
  Mat
  (transpose! [a]
    (if (.isSquare a)
      (.transpose a)
      (error "Matrix has to be squared."))))

(extend-protocol PNumerical
  Vec
  (numerical? [a] true)

  Mat
  (numerical? [a] true))

(extend-protocol PVectorOps
  Vec
  (vector-dot [a b] (inner-product a b))
  (length [a] (.norm a Vector$Norm/Two))
  (length-squared [a] (let [l (.norm a Vector$Norm/Two)] (* l l)))
  (normalise [a] (normalise! (.copy a))))

(extend-protocol PMutableVectorOps
  Vec
  (normalise! [a]
    (element-divide! a (.norm a Vector$Norm/Two))))

(extend-protocol PMatrixOps
  Mat
  (trace [a]
    (if (.isSquare a)
      (let [size (.numRows a)]
        (loop [i 0, sum 0.0]
          (if (< i size)
            (recur (inc i) (+ sum (.get a i i)))
            sum)))
      (error "Square matrix required.")))
  (determinant [a]
    (->> a
         (coerce-param (impl/get-canonical-object :ndarray-double))
         (determinant)))
  (inverse [a]
    (let [^Mat iden (identity-matrix mtj-impl (.numColumns a))
          res (.copy iden)]
      (.solve a iden res))))

(extend-protocol PRowSetting
  Mat
  (set-row [a i row] (set-row! (.copy a) i row))
  (set-row! [a i row]
    (cond
      (zero-dim? row) (let [x (double (get-0d row))]
                        (dotimes [j (.numColumns a)]
                          (.set a j i x))
                        a)
      (is-vector? row) (let [col-size (.numColumns a)]
                         (loop [j 0, vs row]
                           (when (< j col-size)
                             (.set a i j (first vs))
                             (recur (inc j) (rest vs))))
                         a)
      :else (error "row has a incompatible shape."))))

(extend-protocol PColumnSetting
  Mat
  (set-column [a i column] (set-column! (.copy a) i column))
  (set-column! [a i column]
    (cond
      (zero-dim? column) (let [x (double (get-0d column))]
                           (dotimes [j (.numRows a)]
                             (.set a j i x))
                           a)
      (is-vector? column) (let [row-size (.numRows a)]
                            (loop [j 0, vs column]
                              (when (< j row-size)
                                (.set a j i (first vs))
                                (recur (inc j) (rest vs))))
                            a)
      :else (error "column has a incompatible shape."))))

(extend-protocol PElementCount
  Vec
  (element-count [a] (.size a))

  Mat
  (element-count [a] (* (.numRows a) (.numColumns a))))

(extend-protocol PFunctionalOperations
  Vec
  (element-seq [a] (map #(.get a %) (range (.size a))))
  (element-map
    ([a f] (element-map! (.copy a) f))
    ([a f x] (element-map! (.copy a) f x))
    ([a f x more] (element-map! (.copy a) f x more)))
  (element-map!
    ([a f]
     (dotimes [i (.size a)]
       (.set a i (double (f (.get a i)))))
     a)
    ([a f x]
     (let [size (.size a)]
       (loop [i 0, vs x]
         (when (< i size)
           (.set a i (double (f (.get a i) (first x))))
           (recur (inc i) (rest vs)))))
     a)
    ([a f x more]
     (let [xs (cons a (cons x more))
           size (.size a)]
       (loop [i 0, arrs xs]
         (when (< i size)
           (.set a i (double (apply f (mapv #(first %) arrs))))
           (recur (inc i) (mapv rest arrs))))
       a)))
  (element-reduce
    ([a f] (reduce f (element-seq a)))
    ([a f init] (reduce f init (element-seq a))))

  Mat
  (element-seq [a]
    (for [i (range (.numRows a))
          j (range (.numColumns a))]
      (.get a i j)))
  (element-map
    ([a f] (element-map! (.copy a) f))
    ([a f x] (element-map! (.copy a) f x))
    ([a f x more] (element-map! (.copy a) f x more)))
  (element-map!
    ([a f]
     (dotimes [i (.numRows a)]
       (dotimes [j (.numColumns a)]
         (.set a i j (double (f (.get a i j))))))
     a)
    ([a f x]
     (dotimes [i (.numRows a)]
       (dotimes [j (.numColumns a)]
         (.set a i j (double (f (.get a i j) (get-2d x i j))))))
     a)
    ([a f x more]
     (let [xs (cons a (cons x more))]
     (dotimes [i (.numRows a)]
       (dotimes [j (.numColumns a)]
         (.set a i j (double (apply f (mapv #(get-2d % i j) xs))))))
       a)))
  (element-reduce
    ([a f] (reduce f (element-seq a)))
    ([a f init] (reduce f init (element-seq a)))))

(extend-protocol PSelect
  Vec
  (select [a [index]]
    (cond
      (zero-dim? index) (.get a (int (get-0d index)))
      (array? index) (Vec. (Matrices/getSubVector (.mtj a) (int-array index)))
      (= index :all) a))

  Mat
  (select [a [row col]]
    (let [row-index (cond
                      (zero-dim? row) (int-array [(get-0d row)])
                      (sequential? row) (int-array row)
                      (= row :all) (int-array (range (.numRows a))))
          col-index (cond
                      (zero-dim? col) (int-array [(get-0d col)])
                      (sequential? col) (int-array col)
                      (= col :all) (int-array (range (.numColumns a))))]
      (Mat. (Matrices/getSubMatrix (.mtj a) row-index col-index)))))

(extend-protocol PSetSelection
  Vec
  (set-selection [a args values]
    (let [res (.copy a)]
      (assign! (select res args) values)
      res))

  Mat
  (set-selection [a args values]
    (let [res (.copy a)]
      (assign! (select res args) values)
      res)))

(extend-protocol PIndicesAccess
  Vec
  (get-indices [a indices]
    (let [size (count indices)
          res (DenseVector. size)]
      (loop [i 0, indices indices]
        (when-let [idx (first indices)]
          (if (sequential? idx)
            (if (= (count idx) 1)
              (.set res i (.get a (first idx)))
              (error "Incompatible shape of index."))
            (.set res i (.get a idx)))
          (recur (inc i) (rest indices))))
      (Vec. res)))

  Mat
  (get-indices [a indices]
    (let [size (count indices)
          res (DenseVector. size)]
      (loop [i 0, indices indices]
        (when-let [[row col] (first indices)]
          (.set res i (.get a row col))
          (recur (inc i) (rest indices))))
      (Vec. res))))

(extend-protocol PIndicesSetting
  Vec
  (set-indices [a indices values] (set-indices! (.copy a) indices values))
  (set-indices! [a indices values]
    (loop [indices indices, values values]
      (when-let [idx (first indices)]
        (if (sequential? idx)
          (if (= (count idx) 1)
            (.set a (int (first idx)) (double (first values)))
            (error "Incompatible shape of index."))
          (.set a (int idx) (double (first values))))
        (recur (rest indices) (rest values))))
    a)

  Mat
  (set-indices [a indices values] (set-indices! (.copy a) indices values))
  (set-indices! [a indices values]
    (loop [indices indices, values values]
      (when-let [[row col] (first indices)]
        (.set a row col (first values))
        (recur (rest indices) (rest values))))
    a))

(extend-protocol PNorm
  Vec
  (norm [a p]
    (condp = p
      1 (.norm a Vector$Norm/One)
      2 (.norm a Vector$Norm/Two)
      Double/POSITIVE_INFINITY (.norm a Vector$Norm/Infinity)
      (error "p-norm not supported.")))

  Mat
  (norm [a p]
    (condp = p
      1 (.norm a Matrix$Norm/One)
      2 (.norm a Matrix$Norm/Frobenius)
      Double/POSITIVE_INFINITY (.norm a Matrix$Norm/Infinity)
      Long/MAX_VALUE (.norm a Matrix$Norm/Maxvalue)
      (error "p-norm not supported."))))

(extend-protocol PSolveLinear
  Mat
  (solve [a b]
    (cond
      (is-vector? b) (let [^Vec b (coerce-param mtj-impl b)
                           res (.copy b)]
                       (.solve a b res)
                       res)
      (array? b) (let [^Mat b (coerce-param mtj-impl b)
                       res (.copy b)]
                   (.solve a b res)
                   res)
      :else (error "Can only solve for vectors and matrices."))))

(extend-protocol PQRDecomposition
  Mat
  (qr [a {:keys [return]}]
    (let [res (QR/factorize (.mtj a))]
      (into {} (mapv (fn [k]
                       [k (case k
                            :Q (Mat. (.getQ res))
                            :R (Mat. (.getR res)))])
                     return)))))

(extend-protocol PCholeskyDecomposition
  Mat
  (cholesky [a {:keys [return]}]
    (let [res (DenseCholesky/factorize (.mtj a))]
      (into {} (mapv (fn [k]
                       [k (case k
                            :L (Mat. (.getU res))
                            :L* (transpose (Mat. (.getU res))))])
                     return)))))

(extend-protocol PLUDecomposition
  Mat
  (lu [a {:keys [return]}]
    (let [res (DenseLU/factorize (.mtj a))]
      (into {} (mapv (fn [k]
                       [k (case k
                            :L (Mat. (.getL res))
                            :U (Mat. (.getU res))
                            :P (Mat. (.getP res)))])
                     return)))))

(extend-protocol PSVDDecomposition
  Mat
  (svd [a {:keys [return]}]
    (let [res (if (some #{:U :V*} return)
                (SVD/factorize (.mtj a))
                (-> (SVD. (.numRows a) (.numColumns a) false)
                    (.factor (.copy (.mtj a)))))]
      (into {} (mapv (fn [k]
                       [k (case k
                            :U (Mat. (.getU res))
                            :V* (Mat. (.getVt res))
                            :S (.getS res))])
                     return)))))

(extend-protocol PEigenDecomposition
  Mat
  (eigen [a {:keys [return symmetric]
             :or {return [:Q :A], symmetric false}}]
    (let [res (if symmetric
                (SymmDenseEVD/factorize (.mtj a))
                (EVD/factorize (.mtj a)))]
      (into {} (mapv (fn [k]
                       [k (case k
                            :Q (Mat. (if symmetric
                                       (.getEigenvectors ^SymmDenseEVD res)
                                       (.getRightEigenvectors ^EVD res)))
                            :A (diagonal-matrix mtj-impl
                                                (if symmetric
                                                  (.getEigenvalues ^SymmDenseEVD res)
                                                  (.getRealEigenvalues ^EVD res))))])
                     return)))))

(impl/register-implementation (Mat. (DenseMatrix. 1 1)))

(comment
  
(require '[clojure.core.matrix :as m])
(require '[clojure.core.matrix.linear :as l])
(require '[criterium.core :refer [quick-bench]])

(def w 1000)
(def h 1000)
(def numbers (for [_ (range w)]
               (for [_ (range h)]
                 (rand))))

(def clatrix (m/matrix :clatrix numbers))
(def mtj (m/matrix :mtj numbers))

(quick-bench (l/svd clatrix))
(quick-bench (l/svd mtj))

  )
