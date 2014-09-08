(ns cav.mtj.core.matrix
  "Wrapping MTJ for core.matrix."
  (:require [clojure.core.matrix :as m])
  (:use [clojure.core.matrix.protocols]
        [clojure.core.matrix.implementations :only [register-implementation]]
        [clojure.core.matrix.utils :only [error]])
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
          (is-scalar? data) (double data)
          (instance? Vector data) data
          (instance? Matrix data) data   
          (instance? (Class/forName "[D") data) (DenseVector. ^doubles data)
          (instance? (Class/forName "[[D") data) (DenseMatrix. ^"[[D" data)
          (is-vector? data)
            (let [^long size (dimension-count data 0)
                  v (DenseVector. size)]
              (loop [i 0]
                (when (< i size)
                  (.set v i ^double (get-1d data i))
                  (recur (inc i))))
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
              (loop [i 0]
                (when (< i size)
                  (.set mat i 0 (get-1d data i))
                  (recur (inc i))))
              mat)
          :else (error "data has to be a 1D vector.")))
      (row-matrix [m data]
        (if (is-vector? data)
          (let [^long size (dimension-count data 0)
                mat (DenseMatrix. 1 size)]
            (loop [i 0]
              (when (< i size)
                (.set mat 0 i (get-1d data i))
                (recur (inc i))))
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
      (identity-matrix  [m dims] (Matrices/identity (int dims)))
      (diagonal-matrix  [m diagonal-values]
        (let [size (dimension-count diagonal-values 0)
              mat (LowerSymmBandMatrix. (int size) 0)]
          (loop [i 0, vs diagonal-values]
            (when-let [v (first vs)]
              (.set mat i i v)
              (recur (inc i) (rest vs))))
          mat))

    PPermutationMatrix
      (permutation-matrix [m permutation]
        (PermutationMatrix. (int-array permutation)))
    ))

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
  (same-shape [a b]
    (and (instance? Vector b) (= (.size a) (.size ^Vector b))))

  Matrix
  (same-shape [a b]
    (and (instance? Matrix b)
         (let [b ^Matrix b]
           (and (= (.numRows a) (.numRows b))
                (= (.numColumns a) (.numColumns b)))))))

(register-implementation imp)

(comment
  
  (m/set-current-implementation :mtj)
  (m/rows (m/matrix :vectorz [[1 2] [3 4]]))
  (m/mmul (m/diagonal-matrix [1 2]) (m/matrix [2 3]))
(let [mat (m/matrix :mtj [1 2 3 4])]
  (first mat)
  )

  )
