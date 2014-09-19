(ns cav.mtj
  "You will find APIs not offered by core.matrix here."
  (:import [java.io Writer Reader]
           [no.uib.cipr.matrix
            DenseVector DenseMatrix]
           [no.uib.cipr.matrix.io
            MatrixVectorReader MatrixVectorWriter MatrixInfo MatrixSize
            VectorInfo VectorSize MatrixInfo$MatrixField VectorInfo$VectorField
            MatrixInfo$MatrixSymmetry]
           [cav.mtj Mat Vec]))

(defn save-matrix
  "Saves a matrix `m` to `writer`. Only dense matrices are supported, structured
   sparse matrices are coerced to dense matrices before saving. This may have
   a memory impact."
  [^Mat m ^Writer writer]
  (let [matrix-info (MatrixInfo. false
                                 MatrixInfo$MatrixField/Real
                                 MatrixInfo$MatrixSymmetry/General)]
    (doto (MatrixVectorWriter. writer)
      (.printMatrixInfo matrix-info)
      (.printMatrixSize (MatrixSize. (.numRows m)
                                     (.numColumns m)
                                     (* (.numRows m) (.numColumns m))))
      (.printArray (.getData (if (instance? DenseMatrix (.mtj m))
                               ^DenseMatrix (.mtj m)
                               (DenseMatrix. (.mtj m) true))))))
  (.flush writer))

(defn save-vector
  "Saves a vector `v` to `writer`. Only dense vectors are supported, vectors
   are coerced to dense vectors before saving. This may have a memory impact."
  [^Vec v ^Writer writer]
  (let [vector-info (VectorInfo. false VectorInfo$VectorField/Real)]
    (doto (MatrixVectorWriter. writer)
      (.printVectorInfo vector-info)
      (.printVectorSize (VectorSize. (.size v) (.size v)) vector-info)
      (.printArray (.getData (if (instance? DenseVector (.mtj v))
                               ^DenseVector (.mtj v)
                               (DenseVector. (.mtj v) true))))))
  (.flush writer))

(defn read-matrix
  "Reads a saved matrix from `reader`."
  [^Reader reader]
  (Mat. (DenseMatrix. (MatrixVectorReader. reader))))

(defn read-vector
  "Reads a saved vector from `reader`."
  [^Reader reader]
  (Vec. (DenseVector. (MatrixVectorReader. reader))))


(comment
  
(require '[clojure.core.matrix :as m])
(require '[cav.mtj.core.matrix])
(require '[clojure.java.io :as io])


  )
