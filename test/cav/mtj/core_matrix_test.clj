(ns cav.mtj.core-matrix-test
  (:require [clojure.test :refer :all]
            [clojure.core.matrix.compliance-tester :refer [compliance-test]]
            [cav.mtj.core.matrix :refer [imp]]))

(deftest compliance-test
  (compliance-test imp))
