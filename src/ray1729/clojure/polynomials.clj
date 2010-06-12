;; Univariate Polynomials
;; Copyright (c) Ray Miller, 2010. All rights reserved.

(ns
  #^{:author "Ray Miller"
     :doc "Generic arithmetic operations for univariate polynomials - incomplete"}
  ray1729.clojure.polynomials
  (:refer-clojure :exclude (deftype))
  (:use [clojure.contrib.types :only (deftype)]
	[clojure.contrib.generic :only (root-type)])
  (:require [clojure.contrib.generic.arithmetic :as ga]
	    [clojure.contrib.generic.comparison :as gc]))

;
; Polynomials are represented as struct maps with a variable and a
; sorted hash of terms: the key in the hash is the order, and the
; value the corresponding coefficient
;
(defstruct polynomial-struct :variable :terms)

;
; Helper function to handle different input types and construct a sorted-map, 
; eliminating terms with zero coefficient
;
(defn- build-term-list [terms]
  (letfn [(canonicalize-terms
           [terms]
           (cond
             (empty? terms) '()
             (map? (first terms)) (interleave (keys (first terms)) (vals (first terms)))
             (or (vector? (first terms)) (seq? (first terms))) (first terms)
             :else terms))
          (build-term-list
           [accum terms]
           (if (empty? terms) accum
               (let [[order coeff] (first terms)]
                 (if (gc/zero? coeff) 
                   (recur accum (rest terms))
                   (recur (assoc accum order (ga/+ (get accum order 0) coeff)) (rest terms))))))]
    (build-term-list (sorted-map) (partition 2 (canonicalize-terms terms)))))

(deftype ::polynomial polynomial
  (fn [variable & terms] (struct polynomial-struct variable (build-term-list terms)))
  (fn [p] (vals p)))

(derive ::polynomial root-type)

(def variable (accessor polynomial-struct :variable))
(def terms    (accessor polynomial-struct :terms))

(defn degree [p] (key (last (terms p))))

(defn to-fn [p]
  (letfn [(pow [x y] (reduce * (repeat y x)))]
    (fn [u]
      (reduce ga/+ (map (fn [[order coeff]] (ga/* coeff (pow u order))) (terms p))))))

(defn evaluate [p u]
  ((to-fn p) u))

(defn to-str [p]
  (letfn [(plus-or-minus [coeff]
            (if (neg? coeff) " - " " + "))
          (display-coeff [coeff order]
            (cond
              (and (or (= coeff 1) (= coeff -1)) (not (zero? order))) ""
              (neg? coeff) (- coeff)
              :else coeff))
          (x-term [order]
            (cond
              (zero? order) ""
              (= order 1) (variable p)
              :else (str (variable p) "^" order)))
          (leading-term [terms d]
              (let [c (get terms d 0)
                    dc (display-coeff c d)
                    xt (x-term d)]
                (str (if (neg? c) "-" "") dc xt)))]
    (let [d (degree p)]
      (loop [s (leading-term (terms p) d) d (dec d)]
        (if (neg? d) s
          (let [c ((terms p) d)
                pm (plus-or-minus c)
                dc (display-coeff c d)
                xt (x-term d)]
            (recur (if (zero? c) s (str s pm dc xt)) (dec d))))))))

;
; Comparison operations
;

(defmethod gc/zero? ::polynomial
  [p]
  (every? zero? (vals (terms p))))

(defmethod gc/= [::polynomial ::polynomial]
  [p q]
  (and (= (variable p) (variable q)) (= (terms p) (terms q))))

;
; Arithmetic operations
;

(defmethod ga/+ [::polynomial ::polynomial]
  [p q]
  (when (not= (variable p) (variable q))
    (throw (IllegalArgumentException. "addition of polynomials in different variables not supported")))
  (polynomial (variable p) (merge-with ga/+ (terms p) (terms q))))

(defmethod ga/+ [root-type ::polynomial]
  [c p]
  (ga/+ (polynomial (variable p) (merge-with ga/+ {0 c} (terms p)))))

(defmethod ga/+ [::polynomial root-type] 
  [p c] 
  (ga/+ (polynomial (variable p) (merge-with ga/+ (terms p) {0 c}))))

(defmethod ga/- ::polynomial
  [p]
  (polynomial (variable p) (interleave (keys (terms p)) (map ga/- (vals (terms p))))))

(defn- multiply-terms [tp tq]
  (reduce concat (for [[order coeff] tp]
                   (interleave (map #(+ order %) (keys tq)) (map #(ga/* coeff %) (vals tq))))))

(defmethod ga/* [::polynomial ::polynomial]
  [p q]
  (when (not= (variable p) (variable q))
    (throw (IllegalArgumentException. "multiplication of polynomials in different variables not supported")))
  (polynomial (variable p) (multiply-terms (terms p) (terms q))))
  
(defmethod ga/* [root-type ::polynomial]
  [c p]
  (polynomial (variable p) (multiply-terms {0 c} (terms p))))

(defmethod ga/* [::polynomial root-type] 
  [p c]
  (polynomial (variable p) (multiply-terms (terms p) {0 c})))
