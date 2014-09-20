(ns cav.mtj.core-matrix-test
  (:require [clojure.test :refer :all]
            [clojure.core.matrix :refer [matrix add! equals]]
            [clojure.core.matrix.compliance-tester :refer [compliance-test]]
            [cav.mtj.core.matrix]))

(deftest mtj-compliance-test
  (compliance-test :mtj))

(deftest seq-mutability-test
  (testing "seq elements of MTJ vectors/matrices are also mutable"
    (let [mat (matrix :mtj [[1 2] [3 4]])
          first-row (first mat)]
      (add! first-row 10)
      (is (equals first-row [11 12]))
      (is (equals mat [[11 12] [3 4]])))))
