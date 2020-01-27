#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <ctime>
#include <queue>
#include <functional>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>
#include <limits.h>
#include <math.h>

//#include <RcppArmadilloExtensions/sample.h>
#include "../../ontology.h"
#include "ontologyrules.h"
#include "example.h"

using namespace std;
using namespace boost;

const int BITSET_SIZE = 100000;

enum MODE {TRAIN_MODE, TEST_MODE};

struct conjunction_max
{
    double object;  //score
    std::vector<int> ontoID; //identify each ontology
    Node *node2Specified;
    int intNode2Specified;
    //std::vector<Node*> node; //poiter to node
    //std::vector<boost::dynamic_bitset<> > ontoBitset;
    std::vector<boost::dynamic_bitset<> > bitset;
    std::vector<boost::dynamic_bitset<> > discovered;
    int P;
    int N;
    double bestPathScore;

    bool operator<(const conjunction_max& rhs) const
    {
        return object < rhs.object;
    }
};


struct conjunction_min
{
    double object;
    std::vector<int> ontoID;
    Node *node2Specified;
    int intNode2Specified;
    //std::vector<Node*> node;
    std::vector<boost::dynamic_bitset<> > bitset;
    std::vector<boost::dynamic_bitset<> > discovered;
    int P;
    int N;
    double bestPathScore;

    bool operator<(const conjunction_min& rhs) const
    {
        return object > rhs.object;
    }
};

struct exhaustiveRule
{
    boost::dynamic_bitset<> coveredExample;
    std::vector<boost::dynamic_bitset<> *> candidates;
    double score;
};

const string HELP = "Conjunction positive rule learner - Frantisek Malinka\n\
Parameters:\n\
conjunction_learning --train FILE [--test FILE] [--max_rules UINT] [--max_iteration UINT] [--max_queue UINT]\n\
\t--train FILE - path to the ARFF training dataset\n\
\t--test FILE - path to the ARFF testing dataset\n\
\t--max_rules UINT - upper limit for number of rules. DEFAULT 10. IF 0 then cover all positive examples.\n\
\t--max_iteration UINT - upper limit for number of iterations during the procces of finding a rule. DEFAULT 100.\n\
\t--max_queue UINT - upper limit for number of candidates in a priority queue. DEFAULT 10000. IF 0 then queue is unlimited.\n\
\t--verbose - verbose mode\n\
\t--debug - For debug mode use --debug and --verbose\n\
-\t--generateRFile FILE - generate R file for computing ROC curve and AUC for TEST dataset\n\
";


class BeamTopDown
{
public:
    BeamTopDown(){;}
    BeamTopDown(std::vector<bottomFeature> *bottomFeatures, boost::dynamic_bitset<> *myclass, std::string objective, int verbose);
    BeamTopDown(arma::mat *data, std::vector<Ontology *> *refOntologies, std::vector< std::vector<boost::dynamic_bitset<>* > > *POSexamples, std::vector< std::vector<boost::dynamic_bitset<>* > > *NEGexamples, std::vector<boost::dynamic_bitset<> > *enrichItems);
    int run();
    int runExhaustive(Rcpp::DataFrame *enrichPropositionTable);
    newComplex runExhaustive(int filterTH, int ruleDepth, double signTH);
    newComplex runTopologyVersion(int filterTH, int ruleDepth, double signTH, int minLevel);
    int runTest();
    arma::mat removeCoveredExamples(newComplex actbest, arma::mat *new_armaData, boost::dynamic_bitset<> classBitMask);
    std::string printCoveredExamples(newComplex actbest, arma::mat *new_armaData, boost::dynamic_bitset<> classBitMask, Rcpp::CharacterVector rownames, Rcpp::CharacterVector colnames);
    std::string printOnlyCoveredColExamples(newComplex actbest, arma::mat *new_armaData, boost::dynamic_bitset<> classBitMask, Rcpp::CharacterVector rownames, Rcpp::CharacterVector colnames);
    Rcpp::CharacterVector getOnlyCoveredRowExamples(newComplex actbest, boost::dynamic_bitset<> *classBitMask, arma::mat *new_armaData, Rcpp::CharacterVector rownames);
    Rcpp::CharacterVector getOnlyCoveredRowNegExamples(newComplex actbest, boost::dynamic_bitset<> *classBitMask, arma::mat *new_armaData, Rcpp::CharacterVector rownames);
    Rcpp::CharacterVector getOnlyCoveredColExamples(newComplex actbest, boost::dynamic_bitset<> *classBitMask, arma::mat *new_armaData, Rcpp::CharacterVector colnames);
    Rcpp::CharacterVector getOnlyCoveredColNegExamples(newComplex actbest, boost::dynamic_bitset<> *classBitMask, arma::mat *new_armaData, Rcpp::CharacterVector colnames);
    Rcpp::CharacterVector getOnlyCoveredRowColExamples(newComplex actbest, boost::dynamic_bitset<> *classBitMask, arma::mat *new_armaData, Rcpp::CharacterVector rownames, Rcpp::CharacterVector colnames);
    Rcpp::CharacterVector getOnlyCoveredRowColNegExamples(newComplex actbest,  boost::dynamic_bitset<> *classBitMask, arma::mat *new_armaData, Rcpp::CharacterVector rownames, Rcpp::CharacterVector colnames);

private:
    arma::mat *data;
    std::vector<Ontology *> *refOntologies;
    std::vector< std::vector<boost::dynamic_bitset<>* > > *POSexamples;
    std::vector< std::vector<boost::dynamic_bitset<>* > > *NEGexamples;
    std::vector<boost::dynamic_bitset<> > *enrichItems;

    void initPriorityQueue(std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED);
    double evaluateF1Score(boost::dynamic_bitset<> *bitset, int *ontoID);
    double evaluateF1Score(conjunction_max *item);

    std::vector<bottomFeature> *bottomFeatures;
    boost::dynamic_bitset<> *myclass;
    int verbose;


    bool generateNewConjunctions(std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED);
    bool generateNewConjunctions_COMPLETE(std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED);
    bool isBetter(conjunction_max *bestConjunction, conjunction_max *max_heap_top);
    void checkQueueMaxLimit(std::priority_queue<conjunction_max> *max_heap, std::priority_queue<conjunction_min> *min_heap, int limit);

    void initPriorityQueue2(conjunction_max *toInit, std::priority_queue<conjunction_max> *OPEN_heap, std::priority_queue<conjunction_max> *max_heap, std::vector<boost::dynamic_bitset<> > *CLOSED);


	//pointers to functions
	double (BeamTopDown::*evaluationFunction)(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
	void (BeamTopDown::*evaluationFunctionPotential)(newComplex *complex, boost::dynamic_bitset<> *classVector, double *accscore, double *potencial);

    double evaluateF1Score(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    double evaluateF1Score(newComplex *complex, boost::dynamic_bitset<> *classVector);
    double evaluateF1ScorePotencial(newComplex *complex, boost::dynamic_bitset<> *classVector);
    double evaluateF1ScorePotencial(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    double evaluateEntropyScore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    double evaluateEntropyScore(newComplex *complex, boost::dynamic_bitset<> *classVector);
    double evaluateACCScore(newComplex *complex, boost::dynamic_bitset<> *classVector);
    double evaluateACCScore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    double evaluateAUCcore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    double evaluateLaplaceScore(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    double evaluateLaplaceScore(newComplex *complex, boost::dynamic_bitset<> *classVector);
    double evaluateSignificance(newComplex *complex, boost::dynamic_bitset<> *classVector);
    void evaluateACCScorePotencialBoth(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector, double *accscore, double *potential);
    void evaluateACCScorePotencialBoth(newComplex *complex, boost::dynamic_bitset<> *classVector, double *accscore, double *potencial);
    void evaluateF1ScorePotencialBoth(newComplex *complex, boost::dynamic_bitset<> *classVector, double *f1score, double *potencial);
    void evaluateF1ScorePotencialBoth(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector, double *f1score, double *potential);
    void evaluateAUCScorePotencialBoth(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector, double *score, double *potential);
    void evaluateAUCScorePotencialBoth(newComplex *complex, boost::dynamic_bitset<> *classVector, double *score, double *potencial);
    int computeFP(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    int computeFN(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    int computeTP(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    int computeTN(boost::dynamic_bitset<> *vector, boost::dynamic_bitset<> *classVector);
    void filterComplexes(std::priority_queue<newComplex> *queue, int th);
    std::string createHashKey(std::vector<int> *identicalInt);
    std::vector<bottomFeature> makeComplex(newComplex *toExpand);
    std::string getPrintableFeature(newComplex *toPrint);
    std::string getPrintableCoverage(newComplex *toPrint, boost::dynamic_bitset<> *classVector);
    void getCoverage(newComplex *toPrint, boost::dynamic_bitset<> *classVector, int *npos, int *nneg);
    boost::dynamic_bitset<> getRuleBitset(newComplex *rule);
    void addNonPotencialNodes2BITSET(newComplex *mycand, boost::unordered_map<boost::dynamic_bitset<>, bool> *duplicitiesBITSET);
};

/*
struct statistics
{
	int TP;
	int FP;
	int FN;
	int TN;
	double ACC;
	std::vector<double> confidence;	//confidence for test data
};

std::vector<conjunction_max> train();
void test(std::vector<conjunction_max> *rules);
bool parseCommandLineParametrs(int argc, char* argv[]);
int parseFile(const char* file, vector<boost::dynamic_bitset<> > *featureBitSet, boost::dynamic_bitset<> *classMask, vector<string> *features);
int getdataPositionFomARFF(ifstream *myfile, int *featuresNumber, int *examplesNumber);
bool swap(std::priority_queue<conjunction_max> *max_heap, std::priority_queue<conjunction_min> *min_heap);
bool swap(std::priority_queue<conjunction_min> *min_heap, std::priority_queue<conjunction_max> *max_heap);
void checkQueueMaxLimit(std::priority_queue<conjunction_max> *max_heap, std::priority_queue<conjunction_min> *min_heap);
std::vector<conjunction_max> initPriorityQueue(std::priority_queue<conjunction_max> *max_heap, vector<boost::dynamic_bitset<> > *featureBitSet, vector<string> *features, boost::dynamic_bitset<> *classMask);
int countCover(boost::dynamic_bitset<> *example, boost::dynamic_bitset<> *classMask, int *P, int *N);
bool generateNewConjunctions(std::priority_queue<conjunction_max> *max_heap, vector<boost::dynamic_bitset<> > *featureBitSet, boost::unordered_map<string, int> *CLOSED, boost::dynamic_bitset<> *classMask, std::vector<conjunction_max> *baseTerms);
string getPrintableConjunction(conjunction_max best, std::vector<string> *features);
string generateHashKey(boost::dynamic_bitset<> *best);
void eraseCoveredExamples(conjunction_max best, vector<boost::dynamic_bitset<> > *featureBitSet, boost::dynamic_bitset<> *classMask);
void countPN(boost::dynamic_bitset<> *bitset, int *P, int *N);
bool isBetter(conjunction_max *bestConjunction, conjunction_max *max_heap_top);
void printSettings();
int countTrue(boost::dynamic_bitset<> *bitset);
statistics evaluateDataset(MODE mod, vector<boost::dynamic_bitset<> > *featureBitSet, boost::dynamic_bitset<> *classMask, std::vector<conjunction_max> *rules);
void generateRFile(vector<boost::dynamic_bitset<> > *featureBitSet, boost::dynamic_bitset<> *classMask, std::vector<conjunction_max> *rules, statistics *test);
*/
