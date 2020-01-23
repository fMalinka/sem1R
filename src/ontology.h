#ifndef ONTOLOGY_H
#define ONTOLOGY_H

#include <iostream>
#include <vector>
#include <string>
#include <list>
#include "hypothesis.h"


//tip pro wapper fcí
//http://r-pkgs.had.co.nz/src.html

#include "error.h"
#include "ontologyParser.h"

#include <RcppArmadillo.h>

#define NDEBUG
#include <boost/dynamic_bitset.hpp>


//#include <Rcpp.h>
//class OntologyParser; //Forward declaration

typedef struct
{
	std::vector<boost::dynamic_bitset<> > data;
	std::vector<std::string> names;
	boost::dynamic_bitset<> predClass;
} mydataframe;

class Ontology {
public:
    //constructors
    Ontology();
    Ontology(std::string name, std::string pathToOntology, Rcpp::List descriptionTerms, int ontoType);
    ~Ontology();
    bool isCorrect() {return this->correct;}
    std::vector<std::vector<std::string> > getDescriptionTerms() {return this->descriptionTerms;}
    OntologyParser* getOntologyParser() {return this->ontology;}
    int isCovered(Node *searched, arma::vec *bitset, int *arrayIndex);
    int setCoverageVector(hypothesis *searched, arma::vec *bitset, int *arrayIndex);
    int getOntoType() {return this->ontoType;}
    std::vector<Node*> findLGNroots();
    std::string getName() {return this->name;}
    void precomputeSemanticPatterns(int vectorsize);
    std::vector<boost::dynamic_bitset<> *> getSemanticPatterns() {return this->SemanticPatterns;}
    boost::dynamic_bitset<>* getSemanticPattern(int ipattern) {return this->SemanticPatterns[ipattern];}
    boost::dynamic_bitset<>* getSemanticPatternTest(int ipattern) {return this->SemanticPatternsTest[ipattern];}
    void printSemanticPattern(boost::dynamic_bitset<> *bitset);
    double getSemanticDistance(int x, int y) { return this->precomputedDistanceMatrix[x][y];}
    std::vector<std::vector<double> > getSemanticSimilarityMatrix() { return this->precomputedDistanceMatrix;}
    std::vector<std::vector<std::string> > convertRList2Vector(Rcpp::List desc);
    void addDescriptionTermsTest(std::vector<std::vector<std::string> > descriptionTermsTest) {this->descriptionTermsTest = descriptionTermsTest;}
    boost::dynamic_bitset<> getTermBitAncestors(Node *searched);    
private:
    std::string name;
    std::string pathToOntology;
    std::vector<std::vector<std::string> > descriptionTerms;
    std::vector<std::vector<std::string> > descriptionTermsTest;
    OntologyParser *ontology;
    bool correct;
    //std::vector<std::vector<std::string> > convertRList2Vector(Rcpp::List desc);
    int findTestAncestors(std::list<Node*> *nodePool, Node *searched);
    int findTestAncestors2(std::list<Node*> *nodePool, hypothesis *searched);
    int findLGNancestors(Node *searched, boost::dynamic_bitset<> *bitset, boost::unordered_map<Node*, int> *mapNode2Int);
    
    int getLevelOfSpecialization(Node *specific, Node *general);
    int ontoType;
    std::vector<boost::dynamic_bitset<> *> SemanticPatterns;
    std::vector<boost::dynamic_bitset<> *> SemanticPatternsTest;
    std::vector<std::vector<double> > precomputedDistanceMatrix;
};

#endif // ONTOLOGY_H
