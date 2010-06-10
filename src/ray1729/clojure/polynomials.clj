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
; Polynomials are represented as struct maps with a variable
; and hash of terms: the key is the order, and the value the coefficient

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
                 (if (zero? coeff) 
                   (recur accum (rest terms))
                   (recur (assoc accum order (+ (get accum order 0) coeff)) (rest terms))))))]
    (build-term-list (sorted-map) (partition 2 (canonicalize-terms terms)))))

(deftype ::polynomial polynomial
  (fn [variable & terms] (struct polynomial-struct variable (build-term-list terms)))
  (fn [p] (vals p)))

(derive ::polynomial root-type)

(def variable (accessor polynomial-struct :variable))
(def terms    (accessor polynomial-struct :terms))

(defn degree [p] (reduce max (map key (filter #(not (zero? (val %))) (terms p)))))

(defmethod gc/zero? ::polynomial
  [p]
  (every? zero? (vals (terms p))))

(defmethod gc/= [::polynomial ::polynomial]
  [p q]
  (and (= (variable p) (variable q)) (= (terms p) (terms q))))

(defmethod ga/+ [::polynomial ::polynomial]
  [p q]
  (when (not= (variable p) (variable q))
    (throw (IllegalArgumentException. "addition of polynomials in different variables not supported")))
  (polynomial (variable p) (merge-with + (terms p) (terms q))))

(defmethod ga/- ::polynomial
  [p]
  (polynomial (variable p) (interleave (keys (terms p)) (map ga/- (vals (terms p))))))

(defmethod ga/* [::polynomial ::polynomial]
  [p q]
  (when (not= (variable p) (variable q))
    (throw (IllegalArgumentException. "multiplication of polynomials in different variables not supported")))  
  (let [tq (terms q)]
    (reduce ga/+ (for [[order coeff] (terms p)]
                   (polynomial (variable p) 
                    (interleave (map #(+ order %) (keys tq)) (map #(ga/* coeff %) (vals tq))))))))
