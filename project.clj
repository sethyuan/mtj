(defproject cav/mtj "0.4.1"
  :description "A MTJ clojure wrapper."
  :url "https://github.com/sethyuan/mtj"
  :license {:name "MIT"}
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [com.googlecode.matrix-toolkits-java/mtj "1.0.1"]
                 [com.github.fommil.netlib/all "1.1.2" :extension "pom"]
                 [net.mikera/core.matrix "0.29.1"]]
  :global-vars {*warn-on-reflection* true
                *unchecked-math* true}
  :source-paths ["src/clojure"]
  :java-source-paths ["src/java"]
  :profiles {:aot {:aot :all}
             :dev {:dependencies [[net.mikera/vectorz-clj "0.25.0"]
                                  [clatrix "0.4.0"]]}}
  :jvm-opts ^:replace []
  :plugins [[codox "0.8.10"]]
  :codox {:language :clojure
          :defaults {:doc/format :markdown}
          :src-dir-uri "https://github.com/sethyuan/mtj/blob/master/"
          :src-linenum-anchor-prefix "L"})
