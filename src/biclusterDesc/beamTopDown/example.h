#ifndef EXAMPLE_H
#define EXAMPLE_H

#include "../../ontology.h"
#include "../../node.h"
//#include "mlpack/core.hpp"
//#include "../../nsga2/nsga2.h"
#include "ontologyrules.h"
//#include "conjunction_learning.h"

class Example
{
public:    
    Example(int kbics, arma::mat *armaData, std::vector<Ontology *> *refOntologies);
    Example(arma::mat *testarmaData, boost::unordered_map<std::string, std::vector<std::vector<std::string> > > *rowTestOntologyDesc, boost::unordered_map<std::string, std::vector<std::vector<std::string> > > *colTestOntologyDesc, std::vector<Ontology *> *refOntologies);
    ~Example();

    void findPositiveExamples(double threshold);
    void findNegativeExamples(double threshold);

    //sem1r version!!!!!
    void findPositiveExamples();
    void findNegativeExamples();
    //sem1r version END!!!!!

    Rcpp::DataFrame getPropositionalTable(int bic);
    Rcpp::DataFrame getEnrichTable(int bic, std::vector<Node *> *enrichNodes);
    Rcpp::DataFrame getEnrichTablev2(int bic, std::vector<Node *> *enrichNodes);
    mydataframe getEnrichTablev2bitset(int bic, std::vector<Node *> *enrichNodes);
    mydataframe getEnrichTablev2Sigbitset(int bic, std::vector<Node *> *enrichNodes);
    Rcpp::DataFrame getEnrichTablev2Ancestors(int bic, std::vector<Node *> *enrichNodes, int minth);
    mydataframe getEnrichTablev2Ancestorsbitset(int bic, std::vector<Node *> *enrichNodes, int minth);
    mydataframe getEnrichTablev2AncestorsSigbitset(int bic, std::vector<Node *> *enrichNodes, int minth);
    void getEnrichTablev3THNegative(int bic, std::vector<Node *> *enrichNodes, int minth, Rcpp::DataFrame *enrichPropositionTable);
    Rcpp::DataFrame getEnrichTablev3TH(int bic, std::vector<Node *> *enrichNodes, int minth);
    mydataframe getEnrichTablev3THbitset(int bic, std::vector<Node *> *enrichNodes, int minth);
    void computeFeatureVectors(int bic, std::vector<boost::dynamic_bitset<> > *enrichItems, std::vector<Rcpp::LogicalVector> *featureVector, Rcpp::LogicalVector *featureClass, std::vector<std::string> *featureName);
    void computeFeatureVectorsbitset(int bic, std::vector<boost::dynamic_bitset<> > *enrichItems, std::vector<boost::dynamic_bitset<> > *featureVector, boost::dynamic_bitset<> *featureClass, std::vector<std::string> *featureName);
    std::vector<bottomFeature> initBottomFeatures(std::vector<Node *> *enrichNodes, Rcpp::DataFrame *enrichPropositionTable, int onlyPosFeatures);
    std::vector<bottomFeature> initBottomFeaturesPosNeg(std::vector<Node *> *enrichNodes, std::vector<Node *> *enrichNodesNEG,Rcpp::DataFrame *enrichPropositionTable);
    std::vector<bottomFeature> initBottomFeaturesbitset(std::vector<Node *> *enrichNodes, mydataframe *enrichPropositionTable, int onlyPosFeatures);

    boost::dynamic_bitset<> getEnrichClass(Rcpp::DataFrame *enrichPropositionTable);
    static boost::dynamic_bitset<> completAllChilds(newComplex *feature, std::vector<bottomFeature> *bottomFeatures);
    boost::dynamic_bitset<> completAllChilds(bottomFeature *feature, std::vector<bottomFeature> *bottomFeatures);
    static boost::dynamic_bitset<> completAllChildsLastRule(newComplex *feature, std::vector<bottomFeature> *bottomFeatures);
    static boost::dynamic_bitset<> completAllParents(newComplex *feature, std::vector<bottomFeature> *bottomFeatures);
    boost::dynamic_bitset<> completAllParents(bottomFeature *feature, std::vector<bottomFeature> *bottomFeatures);
    static boost::dynamic_bitset<> getAllSpecifics(newComplex *feature);
    static boost::dynamic_bitset<> getAllGenerals(newComplex *feature);
    static void addCurrNode2CLOSED(newComplex *element2Expand, boost::dynamic_bitset<> *closed);
    static boost::dynamic_bitset<> getRoots(std::vector<bottomFeature> *bottom, boost::dynamic_bitset<> *restriction);
    static boost::dynamic_bitset<> getRoots(std::vector<bottomFeature> *bottom);

    void buildTestingExamples(double positiveTH, double negativeTH);
    boost::dynamic_bitset<> coveredExamplesByRule(newComplexStat *rule);
    boost::dynamic_bitset<> getClassBitMask(){return this->classBitMask;}

private:
    double threshold;
    int kbics;
    arma::mat *armaData;
    //std::vector<paretoSet> *myparetoSet;
    std::vector<Ontology *> *refOntologies;
    arma::Col<size_t> *bicAssign;

    arma::mat sumMatOverAll;
    //std::vector<std::vector<boost::unordered_map<int, bool> > > indexPosExamples;
    std::vector<boost::unordered_map<int, bool> > indexPosExamples;
    //std::vector<Rcpp::LogicalVector> featureVector;
    //Rcpp::LogicalVector featureClass;
    //std::vector<std::string> featureName;

    //std::vector< std::vector< std::vector<boost::dynamic_bitset<>* > > > POSexamples;
    std::vector< std::vector<boost::dynamic_bitset<>* > > POSexamples;
    std::vector< std::vector<boost::dynamic_bitset<>* > > NEGexamples;

    boost::dynamic_bitset<> classBitMask;   //vektor urcujici pozitivni a negativni priklady
    boost::unordered_map<std::string, std::vector<std::vector<std::string> > > *rowTestOntologyDesc;
    boost::unordered_map<std::string, std::vector<std::vector<std::string> > > *colTestOntologyDesc;

    void buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexes, std::vector<Ontology *> *refOntologies);
    void buildTrainingExamples(std::vector< std::vector<boost::dynamic_bitset<>* > > *examples, arma::umat *indexMatrix, std::vector<Ontology *> *refOntologies, std::vector<boost::unordered_map<int, bool> > *indexPosExamples);
    double getHyperGeometricScore(int *successSample, int *successPop, int *failurePop, int *sampleSize);

};

#endif // EXAMPLE_H
