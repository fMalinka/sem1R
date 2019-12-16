#ifndef EVALUATE_H
#define EVALUATE_H

#include <iostream>
#include <vector>
#include <string>
#include <list>
#include "ontology.h"


#include <RcppArmadillo.h>
using namespace Rcpp;
#include "error.h"

class Evaluate {
public:
    //constructors
    Evaluate();
    static int evaluateSize(arma::vec *gene, arma::vec *sample);
    static int evaluateXOR(int *gene, int *sample, int nrow, int ncol, arma::mat *data);
    static double evaluateF1score(arma::vec *gene, arma::vec *sample, arma::mat *data);
    static int evaluateXOR(Rcpp::NumericVector *tmpGene, Rcpp::NumericVector *tmpSample, arma::mat *data);  //for significant test
    static int evaluateCorrectPrediction(arma::vec *gene, arma::vec *sample, arma::mat *data);
    static int evaluateCorrectPrediction(Rcpp::NumericVector *tmpGene, Rcpp::NumericVector *tmpSample, arma::mat *data);
    static double evaluateSemantic(arma::vec *gene, int nStart, int nStop, boost::unordered_map<std::string, Ontology *> *ontology);
    static double evaluateSemantic(Rcpp::NumericVector *bitset, boost::unordered_map<std::string, Ontology *> *ontology);   //for significant test

    static void updataDataMatrix(arma::vec *gene, arma::vec *sample, arma::mat *data);
private:
    int x;

};

#endif // EVALUATE_H
